#Code for data analysis associated with Woodruff et al. 2024, "The bacteria of a fig microcommunity"


#This resource was a major guide in analyzing these data:
#https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/set-up-and-pre-processing.html

#also this: https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/doc/MicrobiomeWorkshopII.html

#load libraries

library(ape)
library(reshape2)
library(ggplot2)
library(microbiome)
library(phyloseq)
library(microbiomeutilities)
library(RColorBrewer)
library(ggpubr)
library(DT) 
library(data.table)
library(dplyr) 
library(ggforce)
library(cowplot)
library(lemon)
library(patchwork)
library(effsize)

#set working directory
setwd("/Users/gavin/genome/genome/16S_enrivonmental_fig_microbe_Illumina_metabarcoding_2022/phyloseq_attempt_4-23-22/pub_prep/")

#load OTU counts (see line 319 of workflow.sh), taxonomy table (line 333 of workflow.sh), and sample metadata
ps_object <- read_phyloseq(otu.file = "feature-table.csv", 
                    taxonomy.file = "revised_taxonomy.csv", 
                    metadata.file = "fig_microbe_sample_taiwan_2019_metadata_d_4-23-22.csv", 
                    type = 'simple')

#load phylogenetic tree (line 326 of workflow.sh)
treefile <- read.tree("tree.tree")

#add tree to phyloseq object
ps.ng.tax <- merge_phyloseq(ps_object, treefile)

#get number of figs in this study
    #sample metadata
metadat <- read.csv("fig_microbe_sample_taiwan_2019_metadata_d_4-23-22.csv")
unique(metadat$fig_field_label)
length(unique(metadat$fig_field_label))-2
    #38 figs 

#get number of plants in this study
unique(metadat$plant_label)
length(unique(metadat$plant_label))-2
    #12 plants

metadat
as.data.frame(table(metadat$fig_stage))
#  Var1 Freq
#1        11
#2    B   58
#3    C    6
#4    E    6


#get the number of raw reads per sample before OTU clustering and filtering for organellar DNA
    #line 64 workflow.sh
rawreads <- read.table("/Users/gavin/genome/genome/16S_enrivonmental_fig_microbe_Illumina_metabarcoding_2022/read_counts_one_per_sample.tsv",sep='\t', header=TRUE)
rawreads2 <- data.frame(sample_id = rawreads$sample_id, num_paired_end_reads=rawreads$num_paired_end_reads, step= "Raw reads")

#how many paired-end reads total?
sum(rawreads2$num_paired_end_reads)
#[1] 11361994
    #11,361,994 raw paired-end reads

#per-sample summary stats
summary(rawreads2$num_paired_end_reads)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   2983   98301  157812  140272  186237  273672

    #140,272 paired-end reads on average per sample (raw reads)

#get the number of reads associated with OTU clusters
phyloreads <- as.data.frame(sample_sums(ps.ng.tax))
#get data frame with column ID's
phyloreads$sample_id <- rownames(phyloreads)
phyloreads$num_paired_end_reads <- phyloreads[,1]
phyloreads2 <- data.frame(sample_id = phyloreads$sample_id, num_paired_end_reads=phyloreads$num_paired_end_reads, step= "OTU clustering")


#how many paired-end reads total?
sum(phyloreads2$num_paired_end_reads)
#9421644
    ##9,421,644 paired-end reads after OTU clustering

#per-sample summary stats
summary(phyloreads2$num_paired_end_reads)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   2812   83004  129963  116317  157416  235738
    #116,317 paired-end reads on average per sample (after OTU clustering)

#remove mitochondrial and chloroplast reads

nochmi <- subset_taxa(ps.ng.tax, !Family %in% "Mitochondria" & !Order %in% "Chloroplast")

#reads after removing mito and chloro
orgoreads <- as.data.frame(sample_sums(nochmi))
orgoreads$sample_id <- rownames(orgoreads)
orgoreads$num_paired_end_reads <- orgoreads[,1]

orgoreads2 <- data.frame(sample_id = orgoreads$sample_id, num_paired_end_reads=orgoreads$num_paired_end_reads, step= "Remove organelles")


all_reads <- rbind(rawreads2,phyloreads2,orgoreads2)


all_reads$step <- factor(all_reads$step, levels=c("Raw reads","OTU clustering","Remove organelles"))

#merge df's to get per sample fraction of reads
m1 <- merge(rawreads2,phyloreads2,by="sample_id")
m2 <- merge(m1,orgoreads2,by="sample_id")
#get fraction of reads organellar
m2$fra_org <- 1-(m2$num_paired_end_reads/m2$num_paired_end_reads.y)

#total number of reads after OTU clustering
sum(m2$num_paired_end_reads.y)
#[1] 9421644
    #same as before
#total number of reads after removing organelles
sum(m2$num_paired_end_reads)
#[1] 5215832
#total number of organellar reads
sum(m2$num_paired_end_reads.y)-sum(m2$num_paired_end_reads)
#[1] 4205812

#summary stats per sample fraction organellar reads

summary(m2$fra_org)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.09169 0.33299 0.37272 0.64272 0.93930

#summary stats per reads after removing organelles

summary(m2$num_paired_end_reads)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   2794   35889   50897   64393   94201  172997

#make df for supplemental figure
rdata <- all_reads %>% group_by(step) %>% mutate(mean = mean(num_paired_end_reads))
#supplemental figure 1, read histograms
ggplot(rdata, aes(x = num_paired_end_reads)) + geom_histogram(colour="black",binwidth=6000) + facet_rep_wrap(~step, ncol=1) +theme_cowplot() +geom_vline(aes(xintercept = mean, group = step), colour = 'red') + xlab("Number of paired end reads per sample") + ylab("Frequency") +scale_y_continuous(limits=c(0,9), breaks=c(0:9)) + scale_x_continuous(breaks=c(0,50000,100000,150000,200000,250000,300000),labels = scales::comma)


ggsave("Supplemental_Figure_1.png",bg="white",height=7,width=6,units="in")
    #this is supplemental figure 1

#make the supplemental figure to show the distribution of fraction of organellar reads...
#get just the mitochondrial and chloroplast reads

chloro <- subset_taxa(ps.ng.tax, Order %in% "Chloroplast")
mito <- subset_taxa(ps.ng.tax, Family %in% "Mitochondria")
#fraction chloroplast
chlorodf <- as.data.frame(sample_sums(chloro)/sample_sums(ps.ng.tax))
chlorodf$sample.id <- rownames(chlorodf)
names(chlorodf)[names(chlorodf) == 'sample_sums(chloro)/sample_sums(ps.ng.tax)'] <- 'fra_chloro'

#fraction mitochondrion
mitodf <- as.data.frame(sample_sums(mito)/sample_sums(ps.ng.tax))
mitodf$sample.id <- rownames(mitodf)
names(mitodf)[names(mitodf) == 'sample_sums(mito)/sample_sums(ps.ng.tax)'] <- 'fra_mito'

#metadata configured for R (without phyloseq info)
metadat <- read.csv("/Users/gavin/genome/genome/16S_enrivonmental_fig_microbe_Illumina_metabarcoding_2022/phyloseq_attempt_4-23-22/fig_microbe_sample_taiwan_2019_metadata_R_4-23-22.csv")
#merge with the fraction organelle 
mitmeta <- merge(metadat,mitodf,by="sample.id")

chlmitmeta <- merge(mitmeta,chlorodf,by="sample.id")

#melt, rename things
mchlmitmeta <- melt(chlmitmeta,measure.vars = c("fra_mito","fra_chloro"))
levels(mchlmitmeta$variable) <- c("Mitochondria", "Chloroplast")
names(mchlmitmeta)[names(mchlmitmeta) == 'variable'] <- 'mito.chloro'
names(mchlmitmeta)[names(mchlmitmeta) == 'value'] <- 'fraction_reads'

#get summary stats again

chlmitmeta$fra_mito_chloro <- chlmitmeta$fra_mito + chlmitmeta$fra_chloro
summary(chlmitmeta$fra_mito_chloro)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.09169 0.33299 0.37272 0.64272 0.93930

#replace factor level
levels(mchlmitmeta$mito.chloro)[levels(mchlmitmeta$mito.chloro)=="Mitochondria"] <- "Mitochondrion"


#supplemental figure 2
ggplot(mchlmitmeta, aes(fill=mito.chloro, y=fraction_reads, x=reorder(sample.id, -fraction_reads))) + geom_bar(position="stack", stat="identity") +geom_hline(yintercept=0.33299,linetype="dashed",colour="black") + geom_hline(yintercept=0.37272,linetype="dotted",colour="black") + theme_cowplot() + theme(axis.text.x=element_blank()) + scale_fill_brewer(palette="Set1") +ylim(0,1) + xlab("Samples") + ylab("Fraction of reads") + labs(fill="Organelle")


ggsave("Supplemental_Figure_2.png",bg="white",height=6.5,width=7.5,units="in")

#remove control OTUs


#subset the controls
controls <- prune_samples((sample_names(nochmi) %in% c("EXTRNEG","ExtrNeg1","EXTRPOS","GW10","GW11","GW31","GW9","PCRNeg1","PCRNeg2","PCRPos1","PCRPos2")), nochmi)

# Compute prevalence of each feature, store as data.frame (from https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/doc/MicrobiomeWorkshopII.html)

prevdf = apply(X = otu_table(controls),
               MARGIN = ifelse(taxa_are_rows(controls), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(controls),
                    tax_table(controls))


control_OTU <- prevdf[prevdf$Prevalence > 0,]

nrow(control_OTU)
#[1] 122

#there are 122 OTUs in the controls that will ultimately be removed from the fig microbiome analysis.
#
nochmiprevdf = apply(X = otu_table(nochmi),
               MARGIN = ifelse(taxa_are_rows(controls), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
nochmiprevdf = data.frame(Prevalence = nochmiprevdf,
                    TotalAbundance = taxa_sums(nochmi),
                    tax_table(nochmi))

nrow(nochmiprevdf)
#[1] 3443
122/3443
#[1] 0.03543421
#3.5% of all OTU's were seen at least once in the controls.

#these are the OTUs to remove
control_OTU$OTU

#remove control samples

no_controls <- prune_samples(!(sample_names(nochmi) %in% c("EXTRNEG","ExtrNeg1","EXTRPOS","GW10","GW11","GW31","GW9","PCRNeg1","PCRNeg2","PCRPos1","PCRPos2")), nochmi)

#remove OTU's found in controls and controls; remove control samples
ps3 <- subset_taxa(no_controls, !OTU %in% control_OTU$OTU)


#OTU prevalence


prevdf = apply(X = otu_table(ps3),
               MARGIN = ifelse(taxa_are_rows(ps3), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps3),
                    tax_table(ps3))

#get number of samples

nsamples(ps3)
#[1] 70

#get summary stats for prevalence

summary(prevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.000   1.000   1.000   2.976   2.000  65.000
    
    #No OTU is present in all samples

#prevalence for OTUs histogram

ggplot(prevdf, aes(x = Prevalence)) + geom_histogram(colour="black",binwidth=1) + theme_cowplot() +geom_vline(xintercept = 2.976, colour = 'red') + ylab("OTU count") + scale_x_continuous(breaks=seq(0,70,5)) + scale_y_continuous(limits=c(0,2500),breaks=seq(0,2500,500))

ggsave("Supplemental_Figure_3.png",bg="white",height=4.5,width=7.5,units="in")


#try it all with genera

ps3genera <- tax_glom(ps3, "Genus", NArm = TRUE)



ps3generaprevdf = apply(X = otu_table(ps3genera),
               MARGIN = ifelse(taxa_are_rows(ps3genera), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
ps3generaprevdf = data.frame(Prevalence = ps3generaprevdf,
                    TotalAbundance = taxa_sums(ps3genera),
                    tax_table(ps3genera))


#get number of samples

nsamples(ps3genera)
#[1] 70

#get summary stats for prevalence

summary(ps3generaprevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.000   1.000   3.000   8.803   9.000  70.000

#there is are genera with 100% prevalence, what are they?

ps3generaprevdf[ps3generaprevdf$Prevalence == 70,]$Genus
#[1] "Methylobacterium-Methylorubrum"
    #only one

#get number of genera
length(unique(ps3generaprevdf$Genus))
#[1] 363

#306 after removing dubious genera

#plot Supplemental Figure 4
ggplot(ps3generaprevdf, aes(x = Prevalence)) + geom_histogram(colour="black",binwidth=1) + theme_cowplot() +geom_vline(xintercept = 8.803, colour = 'red')  + ylab("Genus count") + scale_x_continuous(breaks=seq(0,70,5))  + scale_y_continuous(limits=c(0,150),breaks=seq(0,150,25))
ggsave("Supplemental_Figure_4.png",bg="white",height=4.5,width=7.5,units="in")

#Is Wolbachia in here?
ps3generaprevdf[ps3generaprevdf$Genus == "Wolbachia",]$Prevalence
#[1] 2

#only found in two samples

#What about other specific genera?
ps3generaprevdf[ps3generaprevdf$Genus == "Burkholderia",]$Prevalence
    #not present


ps3generaprevdf[ps3generaprevdf$Genus == "Ralstonia",]$Prevalence
    #not present


#okay, subset fig suspensions (interior) and fig surface washes (exterior)

interior_samples <- subset_samples(ps3, surface_or_interior == "interior")

surface_samples <- subset_samples(ps3, surface_or_interior == "surface")



#combine counts to family level for all fig samples, fig suspensions, and fig surface washes
allfam <- tax_glom(ps3, "Family", NArm = TRUE)
prevdf = apply(X = otu_table(allfam),
               MARGIN = ifelse(taxa_are_rows(allfam), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(allfam),
                    tax_table(allfam))

write.table(prevdf, "family_prevalence.tsv",quote=FALSE,sep='\t',row.names=FALSE)


#get number of families
length(unique(prevdf$Family))
#[1] 216

#after removing obvious bad groups,
#182 families

intfam <- tax_glom(interior_samples, "Family", NArm = TRUE)
intprevdf = apply(X = otu_table(intfam),
               MARGIN = ifelse(taxa_are_rows(intfam), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
intprevdf = data.frame(Prevalence = intprevdf,
                    TotalAbundance = taxa_sums(intfam),
                    tax_table(intfam))


surffam <- tax_glom(surface_samples, "Family", NArm = TRUE)
surfprevdf = apply(X = otu_table(surffam),
               MARGIN = ifelse(taxa_are_rows(surffam), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
surfprevdf = data.frame(Prevalence = surfprevdf,
                    TotalAbundance = taxa_sums(surffam),
                    tax_table(surffam))

#these are the 80% prevalence levels
all_eighty <- prevdf[prevdf$Prevalence > 55,]

int_eighty <- intprevdf[intprevdf$Prevalence > 30,]

surf_eighty <- surfprevdf[surfprevdf$Prevalence > 25,]

allrbind <- rbind(all_eighty,int_eighty,surf_eighty)

#get just the unique families
allforplot <- subset(prevdf,Family %in% unique(allrbind$Family))
intforplot <- subset(intprevdf,Family %in% unique(allrbind$Family))
surfforplot <- subset(surfprevdf,Family %in% unique(allrbind$Family))

allforplot$Group <- "All fig samples"
intforplot$Group <- "Fig suspensions"
surfforplot$Group <- "Fig surface washes"

allforplot$Percent_prevalence <- ((allforplot$Prevalence/70)*100)
intforplot$Percent_prevalence <- ((intforplot$Prevalence/38)*100)
surfforplot$Percent_prevalence <- ((surfforplot$Prevalence/32)*100)

plot2df <- rbind(allforplot,intforplot,surfforplot)
#transform abundance for plotting
plot2df$LogAbundance <- log(1 + plot2df$TotalAbundance)

plot2df$Family <- as.factor(plot2df$Family)

plot2df$Family <- factor(plot2df$Family, levels=c("Geodermatophilaceae", "Kineosporiaceae", "Microbacteriaceae", "Nocardiaceae", "Chitinophagaceae", "Hymenobacteraceae", "Sphingobacteriaceae", "Spirosomaceae", "Weeksellaceae", "Acetobacteraceae", "Beijerinckiaceae", "Rhizobiaceae", "Rhodobacteraceae", "Sphingomonadaceae", "Xanthobacteraceae", "Comamonadaceae", "Enterobacteriaceae", "Erwiniaceae", "Moraxellaceae", "Xanthomonadaceae"))

plot2df$Class <- as.factor(plot2df$Class)

plot2df$Class <- factor(plot2df$Class, levels=c("Actinobacteria","Bacteroidia","Alphaproteobacteria","Gammaproteobacteria"))



#number of OTU (again)
prevdf = apply(X = otu_table(ps3),
               MARGIN = ifelse(taxa_are_rows(ps3), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps3),
                    tax_table(ps3))
nrow(prevdf)
#[1] 3321
#3,321 OTU in fig samples


#get summary stats for prevalence
summary(prevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.000   1.000   1.000   2.976   2.000  65.000

#what phyla are OTUs associated with?
phy_otu_count_df <- as.data.frame(table(prevdf$Phylum))

phy_otu_count_df$Fra <- phy_otu_count_df$Freq/sum(phy_otu_count_df$Freq)

phy_otu_count_df[order(phy_otu_count_df$Freq, decreasing = TRUE),]

#                           Var1 Freq          Fra
#23               Proteobacteria 1230 0.3703703704
#6                  Bacteroidota  656 0.1975308642
#4              Actinobacteriota  327 0.0984643180
#1                                275 0.0828063836
#22              Planctomycetota  140 0.0421559771
#7              Bdellovibrionota  135 0.0406504065
#26            Verrucomicrobiota   83 0.0249924721
#17                   Firmicutes   67 0.0201746462
#3               Acidobacteriota   63 0.0189701897
#20                  Myxococcota   62 0.0186690756
#2              Abditibacteriota   59 0.0177657332
#11                Cyanobacteria   54 0.0162601626
#9                   Chloroflexi   50 0.0150557061
#5                Armatimonadota   48 0.0144534779
#12                 Deinococcota   25 0.0075278531
#19              Gemmatimonadota   12 0.0036133695
#21              Patescibacteria   12 0.0036133695
#18               Fusobacteriota    5 0.0015055706
#27                        WPS-2    4 0.0012044565
#10                Crenarchaeota    3 0.0009033424
#8              Campilobacterota    2 0.0006022282
#13                 Dependentiae    2 0.0006022282
#14             Desulfobacterota    2 0.0006022282
#24 SAR324_clade(Marine_group_B)    2 0.0006022282
#15              Elusimicrobiota    1 0.0003011141
#16                Euryarchaeota    1 0.0003011141
#25                  Sumerlaeota    1 0.0003011141




#genus supplemental figure
#combine counts to family level for all fig samples, fig suspensions, and fig surface washes
allgen <- tax_glom(ps3, "Genus", NArm = TRUE)
prevdf = apply(X = otu_table(allgen),
               MARGIN = ifelse(taxa_are_rows(allgen), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(allgen),
                    tax_table(allgen))


intgen <- tax_glom(interior_samples, "Genus", NArm = TRUE)
intprevdf = apply(X = otu_table(intgen),
               MARGIN = ifelse(taxa_are_rows(intgen), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
intprevdf = data.frame(Prevalence = intprevdf,
                    TotalAbundance = taxa_sums(intgen),
                    tax_table(intgen))


surfgen <- tax_glom(surface_samples, "Genus", NArm = TRUE)
surfprevdf = apply(X = otu_table(surfgen),
               MARGIN = ifelse(taxa_are_rows(surfgen), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
surfprevdf = data.frame(Prevalence = surfprevdf,
                    TotalAbundance = taxa_sums(surfgen),
                    tax_table(surfgen))

summary(prevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.000   1.000   3.000   8.803   9.000  70.000

#trying to find a prevalence threshold that affords a good balance-- not too few and not too many genera in the figure

all_ten <- prevdf[prevdf$Prevalence > 10,]
nrow(all_ten)
#[1] 97

all_twenty <- prevdf[prevdf$Prevalence > 20,]
nrow(all_twenty)
#[1] 52


all_thirty <- prevdf[prevdf$Prevalence > 30,]
nrow(all_thirty)

nrow(prevdf[prevdf$Prevalence > 35,])
#[1] 27

nrow(prevdf[prevdf$Prevalence > 40,])
#[1] 23


nrow(prevdf[prevdf$Prevalence > 45,])
#[1] 16

#these are the 64% prevalence levels
all_45 <- prevdf[prevdf$Prevalence > 45,]

int_45 <- intprevdf[intprevdf$Prevalence > 24,]

surf_45 <- surfprevdf[surfprevdf$Prevalence > 20,]

allrbind <- rbind(all_45,int_45,surf_45)

#get just the unique families
allforplot <- subset(prevdf,Genus %in% unique(allrbind$Genus))
intforplot <- subset(intprevdf,Genus %in% unique(allrbind$Genus))
surfforplot <- subset(surfprevdf,Genus %in% unique(allrbind$Genus))

allforplot$Group <- "All fig samples"
intforplot$Group <- "Fig suspensions"
surfforplot$Group <- "Fig surface washes"

allforplot$Percent_prevalence <- ((allforplot$Prevalence/70)*100)
intforplot$Percent_prevalence <- ((intforplot$Prevalence/38)*100)
surfforplot$Percent_prevalence <- ((surfforplot$Prevalence/32)*100)

plot2df <- rbind(allforplot,intforplot,surfforplot)
#transform abundance for plotting
plot2df$LogAbundance <- log(1 + plot2df$TotalAbundance)

write.table(plot2df,"genera.tsv",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)


prevdf$Group <- "All fig samples"
intprevdf$Group <- "Fig suspensions"
surfprevdf$Group <- "Fig surface washes"


prevdf$Percent_prevalence <- ((prevdf$Prevalence/70)*100)
intprevdf$Percent_prevalence <- ((intprevdf$Prevalence/38)*100)
surfprevdf$Percent_prevalence <- ((surfprevdf$Prevalence/32)*100)


allgenera <- rbind(prevdf,intprevdf,surfprevdf)


write.table(allgenera[order(allgenera$Percent_prevalence, decreasing = TRUE),], "all_genera.tsv",,sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)

#get factor levels right


plot2df$Genus <- as.factor(plot2df$Genus)

plot2df$Genus <- factor(plot2df$Genus, levels=c("Gordonia", "Kineococcus", "Klenkia", "Nocardioides", "Quadrisphaera", "Chryseobacterium", "Hymenobacter", "Mucilaginibacter", "Spirosoma", "Bdellovibrio", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Aureimonas", "Methylobacterium-Methylorubrum", "Novosphingobium", "Ochrobactrum", "Paracoccus", "Roseomonas", "Sphingomonas", "Acinetobacter", "Massilia", "Raoultella", "Stenotrophomonas", "Xanthomonas"))


plot2df$Class <- as.factor(plot2df$Class)

plot2df$Class <- factor(plot2df$Class, levels=c("Actinobacteria","Bacteroidia","Bdellovibrionia","Alphaproteobacteria","Gammaproteobacteria"))

levels(plot2df$Genus)[levels(plot2df$Genus)=="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Allorhizobium-Neo.-Par.-Rhi."


ggplot(plot2df,aes(y=Percent_prevalence,x=Genus)) + geom_hline(yintercept=64,linetype="dotted") + geom_point(aes(shape=Group, colour=Class),size=2)  + theme_cowplot(font_size=14) + theme(axis.text.x = element_text(angle = 45, hjust=1,face="italic")) + scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) + scale_colour_manual("Class",values=c("yellow3","#7b3294","orange","#bd0026","#3182bd")) + ylab("Prevalence (%)")
#this is supplemental figure 5
ggsave("prevalence_figure_genus.pdf",height=5.6,width=10,units="in")



#upset plot, supplemental figure 6
library(UpSetR)

    #fig families in variables all_eighty , int_eighty , and surf_eighty above
    #C elegans families from Figure 3C of Zhang et al. 2017 https://doi.org/10.3389/fmicb.2017.00485
        #Actinomycetales is not included because it is an order and not a family

upset_list <- list(All_fig_samples=c("Microbacteriaceae", "Kineosporiaceae", "Comamonadaceae", "Enterobacteriaceae", "Moraxellaceae", "Xanthomonadaceae", "Weeksellaceae", "Hymenobacteraceae", "Spirosomaceae", "Rhizobiaceae", "Beijerinckiaceae", "Acetobacteraceae", "Sphingomonadaceae"),Fig_suspensions=c("Microbacteriaceae", "Kineosporiaceae", "Nocardiaceae", "Comamonadaceae", "Enterobacteriaceae", "Erwiniaceae", "Moraxellaceae", "Xanthomonadaceae", "Sphingobacteriaceae", "Chitinophagaceae", "Weeksellaceae", "Hymenobacteraceae", "Spirosomaceae", "Rhizobiaceae", "Rhodobacteraceae", "Beijerinckiaceae", "Xanthobacteraceae", "Acetobacteraceae", "Sphingomonadaceae"),Fig_surface_washes=c("Microbacteriaceae", "Kineosporiaceae", "Geodermatophilaceae", "Comamonadaceae", "Moraxellaceae", "Xanthomonadaceae", "Hymenobacteraceae", "Spirosomaceae", "Rhizobiaceae", "Beijerinckiaceae", "Acetobacteraceae", "Sphingomonadaceae"),C_elegans_natural=c("Xanthomonadaceae", "Pseudomonadaceae", "Moraxellaceae", "Enterobacteriaceae", "Oxalobacteraceae", "Comamonadaceae", "Sphingomonadaceae", "Acetobacteraceae", "Rhodobacteraceae", "Sphingobacteriaceae", "Weeksellaceae", "Flavobacteriaceae", "Microbacteriaceae"))


#pdf(file="upset.pdf", onefile=FALSE) 
upset(fromList(upset_list), order.by = "freq",text.scale = 2)
#dev.off()
#this is supplemental figure 6

x1 <- unlist(upset_list, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
x1
fromList(upset_list)
    #used to get the families defining various intersections


#Ordination with all samples (for supplemental figures 7-8), no chloro/mito

pslog <- transform_sample_counts(nochmi, function(x) log(1 + x))

#PCoA using Bray-Curtis dissimilarity
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")


pslog <- transform_sample_counts(nochmi, function(x) log(1 + x))

#PCoA using Bray-Curtis dissimilarity, including controls
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")

evals <- out.pcoa.log$values[,1]
ord_df <- plot_ordination(pslog, out.pcoa.log, color = "surface_or_interior", justDF = TRUE)
ggplot(ord_df, aes(x=Axis.1,y=Axis.2)) + geom_point(aes(color = surface_or_interior),size=1) + coord_fixed(sqrt(evals[2] / evals[1])) + scale_color_manual(labels = c("Control", "Fig Suspension", "Fig Surface Wash"), values = c("#fde725", "#21918c","#440154")) +labs(colour="Sample Type") + theme_cowplot()  + xlab("Axis 1 (12%)") + ylab("Axis 2 (10.1%)")
    #figure 3b


#PCoA with weighted Unifrac, including controls

out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "surface_or_interior") + labs(col = "surface/interior/\nctrl")+ coord_fixed(sqrt(evals[2] / evals[1])) + scale_color_manual(labels = c("Control", "Fig Suspension", "Fig Surface Wash"), values = c("#fde725", "#21918c","#440154")) + theme_cowplot() + ggtitle("PCoA, weighted Unifrac")
    #this is supplemental figure 7
ggsave("supplemental_figure_7.png",bg="white",height=4.5,width=7.5,units="in")




#PCoA with unweighted Unifrac, including controls

out.uf.log <- ordinate(pslog, method = "PCoA", distance ="unifrac")
evals <- out.uf.log$values$Eigenvalues
plot_ordination(pslog, out.uf.log, color = "surface_or_interior") + labs(col = "surface/interior/\nctrl")+ coord_fixed(sqrt(evals[2] / evals[1])) + scale_color_manual(labels = c("Control", "Fig Suspension", "Fig Surface Wash"), values = c("#fde725", "#21918c","#440154")) + theme_cowplot() + ggtitle("PCoA, unweighted Unifrac")
    #supplemental figure 8
ggsave("supplemental_figure_8.png",bg="white",height=4.5,width=7.5,units="in")




#PCoA using Bray-Curtis dissimilarity, without controls

pslog <- transform_sample_counts(ps3, function(x) log(1 + x))

out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")

evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "surface_or_interior")+ labs(col = "surface/interior")+ coord_fixed(sqrt(evals[2] / evals[1])) + scale_color_manual(labels = c("Fig Suspension", "Fig Surface Wash"), values = c("#21918c","#440154")) + theme_cowplot() + ggtitle("PCoA, Bray-Curtis")
    #supplemental figure 9
ggsave("supplemental_figure_9.png",bg="white",height=4.5,width=7.5,units="in")




#PCoA with weighted Unifrac, without controls

out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "surface_or_interior") + labs(col = "surface/interior/\nctrl")+ coord_fixed(sqrt(evals[2] / evals[1])) + scale_color_manual(labels = c("Fig Suspension", "Fig Surface Wash"), values = c("#21918c","#440154")) + theme_cowplot() + ggtitle("PCoA, weighted Unifrac")
    #this is supplemental figure 9
ggsave("supplemental_figure_10.png",bg="white",height=4.5,width=7.5,units="in")




#PCoA with unweighted Unifrac, without controls

out.uf.log <- ordinate(pslog, method = "PCoA", distance ="unifrac")
evals <- out.uf.log$values$Eigenvalues
plot_ordination(pslog, out.uf.log, color = "surface_or_interior") + labs(col = "surface/interior/\nctrl")+ coord_fixed(sqrt(evals[2] / evals[1])) + scale_color_manual(labels = c("Fig Suspension", "Fig Surface Wash"), values = c("#21918c","#440154")) + theme_cowplot() + ggtitle("PCoA, unweighted Unifrac")
    #supplemental figure 10
ggsave("supplemental_figure_11.png",bg="white",height=4.5,width=7.5,units="in")



#PERMANOVA, fig suspensions v. fig surface washes

#transform
ps_clr <- microbiome::transform(ps3, "clr")
    #centered log-ratio
    #CLR transform applies a pseudocount of
     #min(relative abundance)/2 to exact zero relative abundance entries
     #in OTU table before taking logs.

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
#ADONIS test
surfintperma <- vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$surface_or_interior)

surfintperma
#                                                  Df SumOfSqs      R2      F
#phyloseq::sample_data(ps_clr)$surface_or_interior  1    11501 0.04496 3.2015
#Residual                                          68   244285 0.95504
#Total                                             69   255786 1.00000
#                                                  Pr(>F)
#phyloseq::sample_data(ps_clr)$surface_or_interior  0.001 ***
#Residual
#Total
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#wilcoxon rank-sum tests for finding OTU differentially abundant among surface and interior.
#melt data
mphyseq = psmelt(ps_clr)

mphyseq$OTU <- as.factor(mphyseq$OTU)
#make a bunch of empty variables for the loop
OTU_id <- NULL
effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL
dom  <- NULL
phyl <- NULL
clas <- NULL
orde <- NULL
fami <- NULL
genu <- NULL
spec <- NULL
#do a MWU test for every OTU comparing surf and int; get cohens d effect size; put in a df with stats I want
for (i in levels(mphyseq$OTU)){
 dat <- mphyseq[mphyseq$OTU == i,]
 if (nrow(dat) >10){surface <- dat[dat$surface_or_interior == "surface",]
 interior <- dat[dat$surface_or_interior == "interior",]
 OTU_id <- rbind(OTU_id, i)
 dom  <- rbind(dom, unique(dat$Domain))
 phyl <- rbind(phyl, unique(dat$Phylum))
 clas <- rbind(clas, unique(dat$Class))
 orde <- rbind(orde, unique(dat$Order))
 fami <- rbind(fami, unique(dat$Family))
 genu <- rbind(genu, unique(dat$Genus))
 spec <- rbind(spec, unique(dat$Species))
 effect_size_lower <- rbind(effect_size_lower, cohen.d(interior$Abundance,surface$Abundance)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(interior$Abundance,surface$Abundance)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(interior$Abundance,surface$Abundance)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(interior$Abundance,surface$Abundance)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(interior$Abundance,surface$Abundance)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(interior$Abundance,surface$Abundance)$p.value)}
}

    #effect size is int-surf, so > 0 reveals enriched in interior, <0 reveals enriched in surface

stat_df <- cbind(OTU_id,dom,phyl,clas,orde,fami,genu,spec,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p)
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("OTU_id","Domain","Phylum","Class","Order","Family","Genus","Species","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df <- as.data.frame(stat_df)
#correct for multiple tests
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

write.table(stat_df, "wilcox_tests_OTU_surface_interior.tsv",sep="\t",row.names=FALSE, quote=FALSE)
stat_df <- read.table("wilcox_tests_OTU_surface_interior.tsv",sep='\t',header=TRUE)
#get sig OTU's

sigotu <- subset(stat_df, p.adj < 0.05)
#how many
nrow(sigotu)
#[1] 11


sigotu$OTU <- as.factor(sigotu$OTU_id)

#make a supplemental figure

mphyseq_plot <- mphyseq[mphyseq$OTU %in% sigotu$OTU,]

summary(mphyseq_plot$Abundance)

mphyseq_plot$OTU <- as.factor(mphyseq_plot$OTU)

levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='607aee8b5f8abc3baf12a67738ff12d1'] <- 'Rhodobacteraceae OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='0d44989f98a705edc1bd5c29df6c552d'] <- 'Acinetobacter OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='b37bf03859da8d004dbb9f4f4a220214'] <- 'Xanthobacteraceae OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='4552af3e905206840cf8b67bcaf2693f'] <- 'Raoultella OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='7ce647f7d457226fc3950cd5ec4eefdf'] <- 'Enterobacterales OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='40ffac9cd8b8ef62bba893e2223310f0'] <- 'Comamonadaceae OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='4754344449d3a4bbfa13be1b17979fab'] <- 'Gordonia OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='04ecfad5772d2e09a84a0f5ef460536c'] <- 'Methylorubrum OTU 1'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='3d3f89fadf28973358458a50d436e91e'] <- 'Methylorubrum OTU 2'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='7521edd61e578d0d49cfd01fbf8cc2f5'] <- 'Xanthomonas OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='42c4ad8258e5ae1d878307dd458b5990'] <- 'Paracoccus OTU'

mphyseq_plot$OTU <- factor(mphyseq_plot$OTU, levels=c('Rhodobacteraceae OTU', 'Xanthobacteraceae OTU', 'Acinetobacter OTU', 'Gordonia OTU', 'Enterobacterales OTU', 'Comamonadaceae OTU', 'Raoultella OTU', 'Paracoccus OTU', 'Xanthomonas OTU', 'Methylorubrum OTU 1', 'Methylorubrum OTU 2'))


ggplot(mphyseq_plot, aes(x=OTU,y=Abundance)) + geom_sina(aes(colour=surface_or_interior),scale="width",alpha=0.5) + stat_summary(aes(group=surface_or_interior),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none"))+ scale_color_manual(labels = c("Fig supsensions", "Fig surface\nwashes"),values = c("#D5770C","#760D7D")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic",size=11),legend.title=element_blank()) + scale_y_continuous(breaks=c(-2:13)) + xlab("OTU") + ylab("Transformed abundance")
    #this is supplemental figure 12
ggsave("supplemental_figure_12.pdf",bg="white",height=4.5,width=7.5,units="in",useDingbats=FALSE)

#now do genera



#transform
ps_clr <- microbiome::transform(ps3genera, "clr")

mphyseq = psmelt(ps_clr)

mphyseq$Genus <- as.factor(mphyseq$Genus)
#make a bunch of empty variables for the loop
Genus_id <- NULL
effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL
dom  <- NULL
phyl <- NULL
clas <- NULL
orde <- NULL
fami <- NULL

stat_df <- NULL

#do a MWU test for every OTU comparing surf and int; get cohens d effect size; put in a df with stats I want
for (i in levels(mphyseq$Genus)){
 dat <- mphyseq[mphyseq$Genus == i,]
 if (nrow(dat) >10){surface <- dat[dat$surface_or_interior == "surface",]
 interior <- dat[dat$surface_or_interior == "interior",]
 Genus_id <- rbind(Genus_id, i)
 dom  <- rbind(dom, unique(dat$Domain))
 phyl <- rbind(phyl, unique(dat$Phylum))
 clas <- rbind(clas, unique(dat$Class))
 orde <- rbind(orde, unique(dat$Order))
 fami <- rbind(fami, unique(dat$Family))
 effect_size_lower <- rbind(effect_size_lower, cohen.d(interior$Abundance,surface$Abundance)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(interior$Abundance,surface$Abundance)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(interior$Abundance,surface$Abundance)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(interior$Abundance,surface$Abundance)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(interior$Abundance,surface$Abundance)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(interior$Abundance,surface$Abundance)$p.value)}
}



    #effect size is int-surf, so > 0 reveals enriched in interior, <0 reveals enriched in surface

stat_df <- as.data.frame(cbind(Genus_id,dom,phyl,clas,orde,fami,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p))
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("Genus_id","Domain","Phylum","Class","Order","Family","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)
#correct for multiple tests
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

write.table(stat_df, "wilcox_tests_Genus_surface_interior.tsv",sep="\t",row.names=FALSE, quote=FALSE)

#get sig genera

sigGEN <- subset(stat_df, p.adj < 0.05)
#how many
nrow(sigGEN)
#[1] 16

sigGEN$Genus_id
    #line 412 of taxonomy
    #remove uncultured

sigGEN <- sigGEN[sigGEN$Genus_id != "uncultured",]

#make a supplemental figure

mphyseq_plot <- mphyseq[mphyseq$Genus %in% sigGEN$Genus_id,]

summary(mphyseq_plot$Abundance)

mphyseq_plot$Genus <- as.factor(mphyseq_plot$Genus)

mphyseq_plot$Genus <- factor(mphyseq_plot$Genus, levels=c("Ochrobactrum", "Gordonia", "Raoultella", "Acinetobacter", "Paracoccus", "Xanthomonas", "Kosakonia", "Stenotrophomonas", "uncultured", "Hymenobacter", "Aureimonas", "Kineococcus", "Quadrisphaera", "Nocardioides", "Sphingomonas", "Methylobacterium-Methylorubrum"))


ggplot(mphyseq_plot, aes(x=Genus,y=Abundance)) + geom_sina(aes(colour=surface_or_interior),scale="width",alpha=0.5) + stat_summary(aes(group=surface_or_interior),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none"))+ scale_color_manual(labels = c("Fig supsensions", "Fig surface\nwashes"),values = c("#D5770C","#760D7D")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic"),legend.title=element_blank()) + scale_y_continuous(breaks=c(-2:13)) + xlab("Genus") + ylab("Transformed abundance")

    #this is supplemental figure 13
ggsave("supplemental_figure_13.pdf",bg="white",height=7,width=10,units="in",useDingbats=FALSE)


#allfam
#now do families

#transform
ps_clr <- microbiome::transform(allfam, "clr")

mphyseq = psmelt(ps_clr)

mphyseq$Family <- as.factor(mphyseq$Family)
#make a bunch of empty variables for the loop
Fami <- NULL
effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL
dom  <- NULL
phyl <- NULL
clas <- NULL
orde <- NULL

stat_df <- NULL



#do a MWU test for every OTU comparing surf and int; get cohens d effect size; put in a df with stats I want
for (i in levels(mphyseq$Family)){
 dat <- mphyseq[mphyseq$Family == i,]
 if (nrow(dat) >10){surface <- dat[dat$surface_or_interior == "surface",]
 interior <- dat[dat$surface_or_interior == "interior",]
 Fami <- rbind(Fami, i)
 dom  <- rbind(dom, unique(dat$Domain))
 phyl <- rbind(phyl, unique(dat$Phylum))
 clas <- rbind(clas, unique(dat$Class))
 orde <- rbind(orde, unique(dat$Order))
 effect_size_lower <- rbind(effect_size_lower, cohen.d(interior$Abundance,surface$Abundance)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(interior$Abundance,surface$Abundance)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(interior$Abundance,surface$Abundance)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(interior$Abundance,surface$Abundance)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(interior$Abundance,surface$Abundance)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(interior$Abundance,surface$Abundance)$p.value)}
}



    #effect size is int-surf, so > 0 reveals enriched in interior, <0 reveals enriched in surface

stat_df <- as.data.frame(cbind(Fami,dom,phyl,clas,orde,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p))
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("Family","Domain","Phylum","Class","Order","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)
#correct for multiple tests
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

write.table(stat_df, "wilcox_tests_Family_surface_interior.tsv",sep="\t",row.names=FALSE, quote=FALSE)

#get sig families

sigfam <- subset(stat_df, p.adj < 0.05)
#how many
nrow(sigfam)
#[1] 16

sigfam$Family
    #line 412 of taxonomy
    #remove uncultured

sigfam <- sigfam[sigfam$Family != "uncultured",]

#make a supplemental figure

mphyseq_plot <- mphyseq[mphyseq$Family %in% sigfam$Family,]

summary(mphyseq_plot$Abundance)

mphyseq_plot$Family <- as.factor(mphyseq_plot$Family)

mphyseq_plot$Family <- factor(mphyseq_plot$Family, levels=c("Xanthobacteraceae", "Rhodobacteraceae", "Moraxellaceae", "Nocardiaceae", "Enterobacteriaceae", "Comamonadaceae", "Chitinophagaceae", "Hymenobacteraceae", "Geodermatophilaceae", "Pseudonocardiaceae", "Kineosporiaceae", "Nocardioidaceae", "Sphingomonadaceae", "Microbacteriaceae", "Beijerinckiaceae"))


ggplot(mphyseq_plot, aes(x=Family,y=Abundance)) + geom_sina(aes(colour=surface_or_interior),scale="width",alpha=0.5) + stat_summary(aes(group=surface_or_interior),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none"))+ scale_color_manual(labels = c("Fig supsensions", "Fig surface\nwashes"),values = c("#D5770C","#760D7D")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic"),legend.title=element_blank()) + scale_y_continuous(breaks=c(-2:13)) + xlab("Family") + ylab("Transformed abundance")


    #this is supplemental figure 14
ggsave("supplemental_figure_14.pdf",bg="white",height=7,width=10,units="in",useDingbats=FALSE)


#diversity metrics
#inside/outside

#
sample_sums(ps3)
#   GW1   GW12   GW13   GW14   GW15   GW16   GW17   GW18   GW19    GW2   GW20
# 40824  96209  73532  36394  23641  26205  78326  19477  40911  24375   6116
#  GW21   GW22   GW23   GW24   GW25   GW26   GW27   GW28   GW29    GW3   GW30
# 34740  41511  17912  46192  35361  42355  13290  25341  31509  37133  29383
#  GW32   GW33   GW34   GW35   GW36   GW37   GW38   GW39    GW4   GW40   GW41
# 37181  13825   6578  12252  30405  86099  23987  57922  62652  39163  28859
#  GW42   GW43   GW44   GW45   GW46   GW47   GW48   GW49    GW5   GW50   GW51
# 41427  72559 125935 152913  55156  51954 122305  45481  85128  62933  32815
#  GW52   GW53   GW54   GW55   GW56   GW57   GW58   GW59    GW6   GW60   GW61
# 21361  23253  63182  25158  44176  23357  32624  29790  86585  65532  78247
#  GW62   GW63   GW64   GW65   GW66   GW67   GW68   GW69    GW7   GW70   GW71
# 42601 113888  29863  34646 101479  79046  43929  96536  72833  18020  90736
#  GW72   GW73   GW74    GW8
# 44788  76275  30840  55364

#exclude GW20, way too few reads

ps3_no_small <- prune_samples((!sample_names(ps3) %in% c("GW20","GW34")), ps3)


#Subsample reads
ps_rare <- phyloseq::rarefy_even_depth(ps3_no_small, rngseed = 123, replace = FALSE)


head(phyloseq::sample_sums(ps_rare))
#  GW1  GW12  GW13  GW14  GW15  GW16
#12252 12252 12252 12252 12252 12252

library(picante)

adiv <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps_rare)))), tree = phyloseq::phy_tree(ps_rare))[, 1],
  "surface_or_interior" = phyloseq::sample_data(ps_rare)$surface_or_interior)

adiv$sample_id <- rownames(adiv)
names(adiv)[names(adiv) == "Observed"] <- "Number_of_OTU"

adiv_melt <- melt(adiv, id.vars= c("sample_id","surface_or_interior"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic"))

levels(adiv_melt$variable)[levels(adiv_melt$variable)=="Number_of_OTU"] <- "Number of OTU"

adiv_melt$surface_or_interior <- as.factor(adiv_melt$surface_or_interior)

levels(adiv_melt$surface_or_interior)[levels(adiv_melt$surface_or_interior)=="surface"] <- "Surface washes"
levels(adiv_melt$surface_or_interior)[levels(adiv_melt$surface_or_interior)=="interior"] <- "Suspensions"

#supplemental figure 15
ggplot(adiv_melt, aes(x=surface_or_interior,y=value)) + geom_sina(scale="width") + stat_summary(aes(group=surface_or_interior),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + facet_wrap(~variable,nrow=1,scales="free") +theme_cowplot() + scale_y_continuous(limits = c(0,NA)) + theme(strip.background = element_rect(colour="white", fill="white"))  + xlab("Sample type") +ylab("Diversity metric value")

ggsave("supplemental_figure_15.png",bg="white",height=6,width=9,units="in")


#diversity statistical tests


div.num.otu <- adiv_melt[adiv_melt$variable == "Number of OTU",]
div.shannon <- adiv_melt[adiv_melt$variable == "Shannon",]
div.phylogenetic <- adiv_melt[adiv_melt$variable == "Phylogenetic",]



#num otu

wilcox.test(div.num.otu[div.num.otu$surface_or_interior == "Suspensions",]$value, div.num.otu[div.num.otu$surface_or_interior == "Surface washes",]$value)

#    Wilcoxon rank sum test with continuity correction
#
#data:  div.num.otu[div.num.otu$surface_or_interior == "Suspensions", ]$value and div.num.otu[div.num.otu$surface_or_interior == "Surface washes", ]$value
#W = 662, p-value = 0.2785
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(div.num.otu[div.num.otu$surface_or_interior == "Surface washes",]$value, div.num.otu[div.num.otu$surface_or_interior == "Suspensions",]$value)
#    Wilcoxon rank sum test with continuity correction
#
#data:  div.num.otu[div.num.otu$surface_or_interior == "Surface washes", ]$value and div.num.otu[div.num.otu$surface_or_interior == "Suspensions", ]$value
#W = 485, p-value = 0.2785
#alternative hypothesis: true location shift is not equal to 0



#shannon
wilcox.test(div.shannon[div.shannon$surface_or_interior == "Suspensions",]$value, div.shannon[div.shannon$surface_or_interior == "Surface washes",]$value)
#    Wilcoxon rank sum exact test
#
#data:  div.shannon[div.shannon$surface_or_interior == "Suspensions", ]$value and div.shannon[div.shannon$surface_or_interior == "Surface washes", ]$value
#W = 482, p-value = 0.2644
#alternative hypothesis: true location shift is not equal to 0


wilcox.test(div.shannon[div.shannon$surface_or_interior == "Surface washes",]$value, div.shannon[div.shannon$surface_or_interior == "Suspensions",]$value)

#    Wilcoxon rank sum exact test
#
#data:  div.shannon[div.shannon$surface_or_interior == "Surface washes", ]$value and div.shannon[div.shannon$surface_or_interior == "Suspensions", ]$value
#W = 665, p-value = 0.2644
#alternative hypothesis: true location shift is not equal to 0



#phylogenetic
wilcox.test(div.phylogenetic[div.phylogenetic$surface_or_interior == "Surface washes",]$value, div.phylogenetic[div.phylogenetic$surface_or_interior == "Suspensions",]$value)

#    Wilcoxon rank sum exact test
#
#data:  div.phylogenetic[div.phylogenetic$surface_or_interior == "Surface washes", ]$value and div.phylogenetic[div.phylogenetic$surface_or_interior == "Suspensions", ]$value
#W = 440, p-value = 0.1017
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(div.phylogenetic[div.phylogenetic$surface_or_interior == "Suspensions",]$value, div.phylogenetic[div.phylogenetic$surface_or_interior == "Surface washes",]$value)

#    Wilcoxon rank sum exact test
#
#data:  div.phylogenetic[div.phylogenetic$surface_or_interior == "Suspensions", ]$value and div.phylogenetic[div.phylogenetic$surface_or_interior == "Surface washes", ]$value
#W = 707, p-value = 0.1017
#alternative hypothesis: true location shift is not equal to 0


#okay, now, nematode occupancy

#just get interior samples
ps4 <- subset_samples(ps3, surface_or_interior == "interior")

#transform
ps_clr <- microbiome::transform(ps4, "clr")


#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 

#ADONIS test PERMANOVA
vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$worms_present)
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999
#
#vegan::adonis2(formula = clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$worms_present)
#                                            Df SumOfSqs      R2     F Pr(>F)
#phyloseq::sample_data(ps_clr)$worms_present  1     4708 0.03498 1.305   0.04 *
#Residual                                    36   129879 0.96502
#Total                                       37   134587 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#OTU


ps_clr <- microbiome::transform(ps4, "clr")


mphyseq = psmelt(ps_clr)

library(effsize)


mphyseq$OTU <- as.factor(mphyseq$OTU)

OTU_id <- NULL
effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL
dom  <- NULL
phyl <- NULL
clas <- NULL
orde <- NULL
fami <- NULL
genu <- NULL
spec <- NULL

for (i in levels(mphyseq$OTU)){
 dat <- mphyseq[mphyseq$OTU == i,]
 if (nrow(dat) >10){worms_yes <- dat[dat$worms_present == "yes",]
 worms_no <- dat[dat$worms_present == "no",]
 OTU_id <- rbind(OTU_id, i)
 dom  <- rbind(dom, unique(dat$Domain))
 phyl <- rbind(phyl, unique(dat$Phylum))
 clas <- rbind(clas, unique(dat$Class))
 orde <- rbind(orde, unique(dat$Order))
 fami <- rbind(fami, unique(dat$Family))
 genu <- rbind(genu, unique(dat$Genus))
 spec <- rbind(spec, unique(dat$Species))
 effect_size_lower <- rbind(effect_size_lower, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(worms_yes$Abundance,worms_no$Abundance)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(worms_yes$Abundance,worms_no$Abundance)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(worms_yes$Abundance,worms_no$Abundance)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(worms_yes$Abundance,worms_no$Abundance)$p.value)}
}

stat_df <- cbind(OTU_id,dom,phyl,clas,orde,fami,genu,spec,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p)
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("OTU_id","Domain","Phylum","Class","Order","Family","Genus","Species","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df <- as.data.frame(stat_df)
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")


stat_df[stat_df$p.adj < 0.05,]

stat_df[stat_df$wilcox_p < 0.05,]

nrow(stat_df[stat_df$wilcox_p < 0.05,])

write.table(stat_df, "wilcox_tests_OTU_worm_presence.tsv",sep="\t",row.names=FALSE, quote=FALSE)
    #This is supplemental table 9


#ZERO sig OTU's after correcting for multiple testing
#let's get top ten by effect size
stat_df$asb_effsize <- abs(stat_df$effect_size)


stat_dfES <- stat_df[order(stat_df$asb_effsize, decreasing = TRUE),]
stat_dfES_top_ten <- head(stat_dfES,10)

mphyseq_plot <- mphyseq[mphyseq$OTU %in% stat_dfES_top_ten$OTU_id,]


#mphyseq_plot <- mphyseq[mphyseq$OTU %in% c("51d810a3fed597ec8151df6d57268336","34ed5dee580cb7e4333f93fa7bbb2354","40ffac9cd8b8ef62bba893e2223310f0","708bb3f490daf2f0465a4095a5c010e4","777c50185200d35116740a777e462069","8b414fe98a85dff09990951214b966b5","b37bf03859da8d004dbb9f4f4a220214","dcba105f35d8ebc9e22269c7491ad3a7","ecedae16c5190f358d2b40c47a943fa0","fc7b257064027e1bc2b6fbfb4584677d","4552af3e905206840cf8b67bcaf2693f"),]

summary(mphyseq_plot$Abundance)

stat_dfES_top_ten[order(stat_dfES_top_ten$effect_size, decreasing = TRUE),]

stat_dfES_top_ten[order(stat_dfES_top_ten$effect_size, decreasing = TRUE),]$OTU_id

mphyseq_plot$OTU <- as.factor(mphyseq_plot$OTU)

levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='51d810a3fed597ec8151df6d57268336'] <- 'Ochrobactrum OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='34ed5dee580cb7e4333f93fa7bbb2354'] <- 'Kosakonia OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='40ffac9cd8b8ef62bba893e2223310f0'] <- 'Comamonadaceae OTU 1'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='708bb3f490daf2f0465a4095a5c010e4'] <- 'Unassigned OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='777c50185200d35116740a777e462069'] <- 'Acinetobacter OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='8b414fe98a85dff09990951214b966b5'] <- 'Aureimonas OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='b37bf03859da8d004dbb9f4f4a220214'] <- 'Xanthobacteraceae OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='dcba105f35d8ebc9e22269c7491ad3a7'] <- 'Stenotrophomonas OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='ecedae16c5190f358d2b40c47a943fa0'] <- 'Chitinophagaceae OTU'
levels(mphyseq_plot$OTU)[levels(mphyseq_plot$OTU)=='fc7b257064027e1bc2b6fbfb4584677d'] <- 'Comamonadaceae OTU 2'

mphyseq_plot$OTU <- factor(mphyseq_plot$OTU, levels=c('Ochrobactrum OTU', 'Kosakonia OTU', 'Stenotrophomonas OTU', 'Comamonadaceae OTU 1', 'Xanthobacteraceae OTU', 'Unassigned OTU', 'Comamonadaceae OTU 2', 'Aureimonas OTU', 'Chitinophagaceae OTU', 'Acinetobacter OTU'))



ggplot(mphyseq_plot, aes(x=OTU,y=Abundance)) + geom_sina(aes(colour=worms_present),scale="width") + stat_summary(aes(group=worms_present),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none")) + scale_color_manual(values = c("#fabb01", "#0081c7")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic",size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15)) + scale_y_continuous(breaks=c(-2:13)) + xlab("OTU") + ylab("Transformed abundance") +labs(colour="Worms\npresent?")

ggsave("supplemental_figure_16.pdf", height=5,width=7.5, units="in",bg="white",useDingbats=FALSE)
    #that is supplemental figure 16



#genera, supplemental table 10


ps_clr <- transform_sample_counts(ps4, function(x){x / sum(x)})


ps_clr_genus <- tax_glom(ps_clr, "Genus", NArm = TRUE)


mphyseq = psmelt(ps_clr_genus)

library(effsize)


mphyseq$OTU <- as.factor(mphyseq$OTU)

mphyseq$Genus <- as.factor(mphyseq$Genus)

genus <- NULL
effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL


for (i in levels(mphyseq$Genus)){
 dat <- mphyseq[mphyseq$Genus == i,]
 if (nrow(dat) >10){worms_yes <- dat[dat$worms_present == "yes",]
 worms_no <- dat[dat$worms_present == "no",]
 genus <- rbind(genus, i)
 effect_size_lower <- rbind(effect_size_lower, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(worms_yes$Abundance,worms_no$Abundance)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(worms_yes$Abundance,worms_no$Abundance)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(worms_yes$Abundance,worms_no$Abundance,exact = FALSE)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(worms_yes$Abundance,worms_no$Abundance,exact = FALSE)$p.value)}
}

stat_df <- cbind(genus,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p)
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("genus","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df <- as.data.frame(stat_df)
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)

stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

write.table(stat_df, "wilcox_tests_Genus_worm_presence.tsv",sep="\t",row.names=FALSE, quote=FALSE)

#what is up with the NAs?
mphyseq[mphyseq$Genus == "Actinomyces",]
mphyseq[mphyseq$Genus == "Anaerocolumna",]$Abundance
mphyseq[mphyseq$Genus == "Anoxybacillus",]$Abundance
mphyseq[mphyseq$Genus == "Campylobacter",]$Abundance
mphyseq[mphyseq$Genus == "Candidatus_Protochlamydia",]$Abundance
mphyseq[mphyseq$Genus == "Cetobacterium",]$Abundance
mphyseq[mphyseq$Genus == "Chalicogloea_CCALA_975",]$Abundance
mphyseq[mphyseq$Genus == "Clostridium_sensu_stricto_12",]$Abundance
mphyseq[mphyseq$Genus == "Cohnella",]$Abundance
#okay, where all figs have zero abundance. these must be exterior-specific taxa?

#family, supplemental table 11

ps_clr_family <- tax_glom(ps_clr, "Family", NArm = TRUE)


mphyseq = psmelt(ps_clr_family)

library(effsize)


mphyseq$OTU <- as.factor(mphyseq$OTU)

mphyseq$Family <- as.factor(mphyseq$Family)

family <- NULL
effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL


for (i in levels(mphyseq$Family)){
 dat <- mphyseq[mphyseq$Family == i,]
 if (nrow(dat) >10){worms_yes <- dat[dat$worms_present == "yes",]
 worms_no <- dat[dat$worms_present == "no",]
 family <- rbind(family, i)
 effect_size_lower <- rbind(effect_size_lower, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(worms_yes$Abundance,worms_no$Abundance)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(worms_yes$Abundance,worms_no$Abundance)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(worms_yes$Abundance,worms_no$Abundance,exact = FALSE)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(worms_yes$Abundance,worms_no$Abundance,exact = FALSE)$p.value)}
}

stat_df <- cbind(family,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p)
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("family","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df <- as.data.frame(stat_df)
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

stat_df[stat_df$p.adj < 0.1,]
write.table(stat_df, "wilcox_tests_Family_worm_presence.tsv",sep="\t",row.names=FALSE, quote=FALSE)

#####
#order, supplemental table 12

ps_clr_order <- tax_glom(ps_clr, "Order", NArm = TRUE)


mphyseq = psmelt(ps_clr_order)

library(effsize)


mphyseq$OTU <- as.factor(mphyseq$OTU)

mphyseq$Order <- as.factor(mphyseq$Order)

order <- NULL
effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL


for (i in levels(mphyseq$Order)){
 dat <- mphyseq[mphyseq$Order == i,]
 if (nrow(dat) >10){worms_yes <- dat[dat$worms_present == "yes",]
 worms_no <- dat[dat$worms_present == "no",]
 order <- rbind(order, i)
 effect_size_lower <- rbind(effect_size_lower, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(worms_yes$Abundance,worms_no$Abundance)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(worms_yes$Abundance,worms_no$Abundance)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(worms_yes$Abundance,worms_no$Abundance,exact = FALSE)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(worms_yes$Abundance,worms_no$Abundance,exact = FALSE)$p.value)}
}

stat_df <- cbind(order,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p)
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("order","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df <- as.data.frame(stat_df)
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

stat_df[stat_df$p.adj < 0.1,]


write.table(stat_df, "wilcox_tests_Order_worm_presence.tsv",sep="\t",row.names=FALSE, quote=FALSE)



#class, supplemental table 13



ps_clr_class <- tax_glom(ps_clr, "Class", NArm = TRUE)


mphyseq = psmelt(ps_clr_class)

library(effsize)


mphyseq$OTU <- as.factor(mphyseq$OTU)

mphyseq$Class <- as.factor(mphyseq$Class)

class <- NULL
effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL


for (i in levels(mphyseq$Class)){
 dat <- mphyseq[mphyseq$Class == i,]
 if (nrow(dat) >10){worms_yes <- dat[dat$worms_present == "yes",]
 worms_no <- dat[dat$worms_present == "no",]
 class <- rbind(class, i)
 effect_size_lower <- rbind(effect_size_lower, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(worms_yes$Abundance,worms_no$Abundance)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(worms_yes$Abundance,worms_no$Abundance)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(worms_yes$Abundance,worms_no$Abundance,exact = FALSE)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(worms_yes$Abundance,worms_no$Abundance,exact = FALSE)$p.value)}
}

stat_df <- cbind(class,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p)
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("class","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df <- as.data.frame(stat_df)
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

stat_df[stat_df$p.adj < 0.1,]


write.table(stat_df, "wilcox_tests_Class_worm_presence.tsv",sep="\t",row.names=FALSE, quote=FALSE)


#phylum, supplemental table 14

ps_clr_phylum <- tax_glom(ps_clr, "Phylum", NArm = TRUE)


mphyseq = psmelt(ps_clr_phylum)

library(effsize)


mphyseq$OTU <- as.factor(mphyseq$OTU)

mphyseq$Phylum <- as.factor(mphyseq$Phylum)

phylum_id <- NULL
effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL
dom  <- NULL
phyl <- NULL
clas <- NULL
orde <- NULL
fami <- NULL
genu <- NULL
spec <- NULL

for (i in levels(mphyseq$Phylum)){
 dat <- mphyseq[mphyseq$Phylum == i,]
 if (nrow(dat) >10){worms_yes <- dat[dat$worms_present == "yes",]
 worms_no <- dat[dat$worms_present == "no",]
 phylum_id <- rbind(phylum_id, i)
 dom  <- rbind(dom, unique(dat$Domain))
 phyl <- rbind(phyl, unique(dat$Phylum))
 clas <- rbind(clas, unique(dat$Class))
 orde <- rbind(orde, unique(dat$Order))
 fami <- rbind(fami, unique(dat$Family))
 genu <- rbind(genu, unique(dat$Genus))
 spec <- rbind(spec, unique(dat$Species))
 effect_size_lower <- rbind(effect_size_lower, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(worms_yes$Abundance,worms_no$Abundance)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(worms_yes$Abundance,worms_no$Abundance)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(worms_yes$Abundance,worms_no$Abundance)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(worms_yes$Abundance,worms_no$Abundance,exact = FALSE)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(worms_yes$Abundance,worms_no$Abundance,exact = FALSE)$p.value)}
}

stat_df <- cbind(phylum_id,dom,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p)
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("phylum_id","Domain","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df <- as.data.frame(stat_df)
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

stat_df[stat_df$p.adj < 0.1,]
#ZERO sig Phyla after adjustment

write.table(stat_df, "wilcox_tests_Phylum_worm_presence.tsv",sep="\t",row.names=FALSE, quote=FALSE)


#ordinations with interior (fig suspension) samples only, colored by nematode occupancy


pslog <- transform_sample_counts(ps4, function(x) log(1 + x))

#PCoA using Bray-Curtis dissimilarity
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")

evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "worms_present") + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot() + scale_color_manual( values = c("#fabb01", "#0081c7"))  + ggtitle("PCoA, Bray-Curtis")
    #this is figure 4b



#PCoA with weighted Unifrac , supplemental figure 17

out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "worms_present") + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot() + scale_color_manual( values = c("#fabb01", "#0081c7")) + theme_cowplot() + ggtitle("PCoA, weighted Unifrac")


ggsave("supplemental_figure_17.pdf", height=4.5,width=7.5, units="in",bg="white",useDingbats=FALSE)



#PCoA with unweighted Unifrac, supplemental figure 18

out.uf.log <- ordinate(pslog, method = "PCoA", distance ="unifrac")
evals <- out.uf.log$values$Eigenvalues
plot_ordination(pslog, out.uf.log, color = "worms_present") + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot() + scale_color_manual( values = c("#fabb01", "#0081c7")) + theme_cowplot() + ggtitle("PCoA, unweighted Unifrac")
ggsave("supplemental_figure_18.pdf", height=4.5,width=7.5, units="in",bg="white",useDingbats=FALSE)


#CCA , supplemental figure 19

out.cca.log <- ordinate(pslog, method = "CCA")

plot_ordination(pslog, out.cca.log , color = "worms_present") + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot() + scale_color_manual( values = c("#fabb01", "#0081c7"))  + theme_cowplot() + ggtitle("Canonical Correspondence Analysis")
ggsave("supplemental_figure_19.pdf", height=4.5,width=7.5, units="in",bg="white",useDingbats=FALSE)




#alpha diversity depending on nematode occupancy





#diversity metrics
#worms present/absent

sample_sums(ps4)
#   GW1   GW13   GW14   GW15   GW17   GW19    GW2   GW21   GW22   GW23   GW25
# 40824  73532  36394  23641  78326  40911  24375  34740  41511  17912  35361
#  GW27   GW29    GW3   GW32   GW34   GW36   GW38    GW4   GW42   GW43   GW45
# 13290  31509  37133  37181   6578  30405  23987  62652  41427  72559 152913
#  GW47   GW48   GW50   GW51   GW52   GW55   GW58   GW61   GW63   GW64   GW65
# 51954 122305  62933  32815  21361  25158  32624  78247 113888  29863  34646
#  GW66   GW68   GW70   GW73   GW74
#101479  43929  18020  76275  30840

ps3_no_small <- prune_samples((!sample_names(ps4) %in% c("GW20","GW34")), ps4)


ps_rare <- phyloseq::rarefy_even_depth(ps3_no_small, rngseed = 123, replace = FALSE)

#ordination is the arrangement or ‘ordering’ of species and/or sample units along gradients.

head(phyloseq::sample_sums(ps_rare))
#13290 13290 13290 13290 13290 13290

adiv <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps_rare)))), tree = phyloseq::phy_tree(ps_rare))[, 1],
  "worms_present" = phyloseq::sample_data(ps_rare)$worms_present)

adiv$sample_id <- rownames(adiv)
names(adiv)[names(adiv) == "Observed"] <- "Number_of_OTU"

adiv_melt <- melt(adiv, id.vars= c("sample_id","worms_present"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic"))

levels(adiv_melt$variable)[levels(adiv_melt$variable)=="Number_of_OTU"] <- "Number of OTU"

ggplot(adiv_melt, aes(x=worms_present,y=value)) + geom_sina(scale="width") + stat_summary(aes(group=worms_present),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + facet_wrap(~variable,nrow=1,scales="free") +theme_cowplot() + scale_y_continuous(limits = c(0,NA)) + xlab("Worms Present?") +ylab("Diversity measure") + theme(strip.background = element_rect(colour="white", fill="white"))

ggsave("worms_present_diversity_sina_2024.pdf", height=4,width=7, units="in",bg="white",useDingbats=FALSE)
    #this is figure 5

#diversity statistical tests

div.num.otu <- adiv_melt[adiv_melt$variable == "Number of OTU",]
div.shannon <- adiv_melt[adiv_melt$variable == "Shannon",]
div.phylogenetic <- adiv_melt[adiv_melt$variable == "Phylogenetic",]




#num otu

wilcox.test(div.num.otu[div.num.otu$worms_present == "yes",]$value, div.num.otu[div.num.otu$worms_present == "no",]$value)

#  Wilcoxon rank sum test with continuity correction
#
#data:  div.num.otu[div.num.otu$worms_present == "yes", ]$value and div.num.otu[div.num.otu$worms_present == "no", ]$value
#W = 141.5, p-value = 0.3781
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(div.shannon[div.shannon$worms_present == "yes",]$value, div.shannon[div.shannon$worms_present == "no",]$value)

#> wilcox.test(div.shannon[div.shannon$worms_present == "yes",]$value, div.shannon[div.shannon$worms_present == "no",]$value)
#
#  Wilcoxon rank sum exact test
#
#data:  div.shannon[div.shannon$worms_present == "yes", ]$value and div.shannon[div.shannon$worms_present == "no", ]$value
#W = 148, p-value = 0.4987
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(div.phylogenetic[div.phylogenetic$worms_present == "yes",]$value, div.phylogenetic[div.phylogenetic$worms_present == "no",]$value)

#  Wilcoxon rank sum exact test
#
#data:  div.phylogenetic[div.phylogenetic$worms_present == "yes", ]$value and div.phylogenetic[div.phylogenetic$worms_present == "no", ]$value
#W = 146, p-value = 0.4612
#alternative hypothesis: true location shift is not equal to 0




###foundress number



adiv <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps_rare)))), tree = phyloseq::phy_tree(ps_rare))[, 1],
  "foundress_number" = phyloseq::sample_data(ps_rare)$foundress_number)

adiv$sample_id <- rownames(adiv)
names(adiv)[names(adiv) == "Observed"] <- "Number_of_OTU"

adiv_melt <- melt(adiv, id.vars= c("sample_id","foundress_number"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic"))

levels(adiv_melt$variable)[levels(adiv_melt$variable)=="Number_of_OTU"] <- "Number of OTU"

str(adiv_melt)

adiv_melt$foundress_number <- as.numeric(adiv_melt$foundress_number)

ggplot(adiv_melt, aes(x=foundress_number,y=value)) + geom_point() + facet_wrap(~variable,nrow=1,scales="free") +theme_cowplot() + stat_smooth() + scale_y_continuous(limits = c(0,NA)) + xlab("Foundress number") +ylab("Diversity measure") + theme(strip.background = element_rect(colour="white", fill="white")) + scale_y_continuous(oob=scales::rescale_none)
    #this is figure 6A


#linear models, foundress number



adiv$foundress_number <- as.numeric(adiv$foundress_number)

summary(lm(Number_of_OTU ~ foundress_number, data=adiv))

#Call:
#lm(formula = Number_of_OTU ~ foundress_number, data = adiv)
#
#Residuals:
#    Min      1Q  Median      3Q     Max
#-141.32  -84.32  -34.75   61.16  456.16
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)       194.837     25.682   7.586  6.8e-09 ***
#foundress_number  -11.521      4.801  -2.400   0.0219 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 126 on 35 degrees of freedom
#Multiple R-squared:  0.1413,    Adjusted R-squared:  0.1168
#F-statistic:  5.76 on 1 and 35 DF,  p-value: 0.02185

summary(lm(Shannon ~ foundress_number, data=adiv))

#Call:
#lm(formula = Shannon ~ foundress_number, data = adiv)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-2.00545 -0.65114  0.01302  0.44890  1.69981
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)       3.12836    0.19776  15.819   <2e-16 ***
#foundress_number -0.06030    0.03697  -1.631    0.112
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.9703 on 35 degrees of freedom
#Multiple R-squared:  0.07065,   Adjusted R-squared:  0.04409
#F-statistic: 2.661 on 1 and 35 DF,  p-value: 0.1118

summary(lm(Phylogenetic ~ foundress_number, data=adiv))

#Call:
#lm(formula = Phylogenetic ~ foundress_number, data = adiv)
#
#Residuals:
#    Min      1Q  Median      3Q     Max
#-6.0236 -3.1945  0.0568  2.1393 15.0269
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)        9.7292     0.8969  10.847 9.68e-13 ***
#foundress_number  -0.4294     0.1677  -2.561   0.0149 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 4.401 on 35 degrees of freedom
#Multiple R-squared:  0.1578,    Adjusted R-squared:  0.1338
#F-statistic:  6.56 on 1 and 35 DF,  p-value: 0.0149

#looking at evenness

eve_df <- evenness(ps3_no_small)
eve_df$sample_id <- rownames(eve_df)

eve_adiv <- merge(adiv,eve_df,by="sample_id")


eve_adiv$foundress_number <- as.numeric(eve_adiv$foundress_number)


ggplot(eve_adiv, aes(x=foundress_number,y=Phylogenetic)) + geom_point() + stat_smooth() +theme_cowplot() 

evefounddf <- data.frame(sample_id=eve_adiv$sample_id,foundress_number=eve_adiv$foundress_number,camargo=eve_adiv$camargo,pielou=eve_adiv$pielou,simpson=eve_adiv$simpson,evar=eve_adiv$evar,bulla=eve_adiv$bulla)


eve_melt <- melt(evefounddf, id.vars= c("sample_id","foundress_number"),measure.vars = c("camargo", "pielou","simpson","evar","bulla"))


ggplot(eve_melt, aes(x=foundress_number,y=value)) + geom_point() + facet_wrap(~variable,scales="free") +theme_cowplot() + stat_smooth() + scale_y_continuous(limits = c(0,NA)) + xlab("Foundress number") +ylab("Evenness measure") + theme(strip.background = element_rect(colour="white", fill="white")) + scale_y_continuous(oob=scales::rescale_none)

ggsave("evenness.pdf",useDingbats=FALSE)
    #this is a supplemental figure
summary(lm(camargo ~ foundress_number, data=evefounddf))
#Call:
#lm(formula = camargo ~ foundress_number, data = evefounddf)
#
#Residuals:
#      Min        1Q    Median        3Q       Max
#-0.059369 -0.002893  0.003430  0.009225  0.012440
#
#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)
#(Intercept)      0.9875599  0.0029474 335.066   <2e-16 ***
#foundress_number 0.0010477  0.0005509   1.902   0.0655 .
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.01446 on 35 degrees of freedom
#Multiple R-squared:  0.09366,   Adjusted R-squared:  0.06776
#F-statistic: 3.617 on 1 and 35 DF,  p-value: 0.06546
summary(lm(pielou ~ foundress_number, data=evefounddf))
#Call:
#lm(formula = pielou ~ foundress_number, data = evefounddf)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-0.32177 -0.06073  0.02026  0.11530  0.23308
#
#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.606137   0.029850  20.306   <2e-16 ***
#foundress_number -0.002618   0.005580  -0.469    0.642
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1465 on 35 degrees of freedom
#Multiple R-squared:  0.006252,  Adjusted R-squared:  -0.02214
#F-statistic: 0.2202 on 1 and 35 DF,  p-value: 0.6418
summary(lm(simpson ~ foundress_number, data=evefounddf))
#Call:
#lm(formula = simpson ~ foundress_number, data = evefounddf)
#
#Residuals:
#      Min        1Q    Median        3Q       Max
#-0.076257 -0.033178 -0.006922  0.019618  0.137131
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)      0.081488   0.010406   7.831 3.34e-09 ***
#foundress_number 0.001473   0.001945   0.757    0.454
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.05106 on 35 degrees of freedom
#Multiple R-squared:  0.01611,   Adjusted R-squared:  -0.012
#F-statistic: 0.5732 on 1 and 35 DF,  p-value: 0.4541
summary(lm(evar ~ foundress_number, data=evefounddf))
#Call:
#lm(formula = evar ~ foundress_number, data = evefounddf)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-0.09521 -0.04041  0.01113  0.03235  0.10013
#
#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.227211   0.010645  21.345   <2e-16 ***
#foundress_number -0.005269   0.001990  -2.648    0.012 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.05223 on 35 degrees of freedom
#Multiple R-squared:  0.1669,    Adjusted R-squared:  0.1431
#F-statistic: 7.014 on 1 and 35 DF,  p-value: 0.01205
summary(lm(bulla ~ foundress_number, data=evefounddf))
#Call:
#lm(formula = bulla ~ foundress_number, data = evefounddf)
#
#Residuals:
#      Min        1Q    Median        3Q       Max
#-0.195793 -0.066910 -0.001223  0.070614  0.170240
#
#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.298914   0.019169  15.594   <2e-16 ***
#foundress_number -0.003221   0.003583  -0.899    0.375
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.09405 on 35 degrees of freedom
#Multiple R-squared:  0.02257,   Adjusted R-squared:  -0.005359
#F-statistic: 0.8081 on 1 and 35 DF,  p-value: 0.3748





ps4 <- subset_samples(ps3, surface_or_interior == "interior")


prevdf = apply(X = otu_table(ps4),
               MARGIN = ifelse(taxa_are_rows(ps4), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps4),
                    tax_table(ps4))

nrow(prevdf)
#[1] 3321


pslog <- transform_sample_counts(ps4, function(x) log(1 + x))

#PCoA using Bray-Curtis dissimilarity
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")


evals <- out.pcoa.log$values[,1]

ord_df <-plot_ordination(pslog, out.pcoa.log, color = "foundress_number", shape="worms_present", justDF = TRUE)
ord_df$foundress_number <- as.numeric(ord_df$foundress_number)

ggplot(ord_df, aes(x=Axis.1,y=Axis.2)) + geom_point(aes(size = foundress_number, shape=worms_present, colour=foundress_number),alpha=0.6) + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot()  + xlab("Axis 1 (19.4%)") + ylab("Axis 2 (8.1%)")  + scale_size_continuous(limits=c(0,20),breaks=c(0,5,10,15,20)) + scale_color_continuous(limits=c(0,20),breaks=c(0,5,10,15,20)) + guides(color=guide_legend(title="Foundress\nnumber"), size = guide_legend(title="Foundress\nnumber"), shape=guide_legend(title="Worms\npresent?")) + theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),legend.title = element_text(size = 14))
    #this is figure 6b

#PCoA with weighted Unifrac

out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
evals <- out.wuf.log$values$Eigenvalues

plot_ordination(pslog, out.wuf.log, shape = "worms_present") + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot() + scale_color_manual( values = c("#fabb01", "#0081c7")) + theme_cowplot() + ggtitle("PCoA, weighted Unifrac")

ord_df <-plot_ordination(pslog, out.wuf.log, color = "foundress_number", shape="worms_present", justDF = TRUE)


ord_df$foundress_number <- as.numeric(ord_df$foundress_number)


ggplot(ord_df, aes(x=Axis.1,y=Axis.2)) + geom_point(aes(size = foundress_number, shape=worms_present, colour=foundress_number),alpha=0.6) + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot()  + xlab("Axis 1 (30.2%)") + ylab("Axis 2 (13.4%)")  + scale_size_continuous(limits=c(0,20),breaks=c(0,5,10,15,20)) + scale_color_continuous(limits=c(0,20),breaks=c(0,5,10,15,20)) + guides(color=guide_legend(title="Foundress\nnumber"), size = guide_legend(title="Foundress\nnumber"), shape=guide_legend(title="Worms\npresent?")) + theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),legend.title = element_text(size = 14)) + ggtitle("PCoA, weighted Unifrac")

ggsave("supplemental_figure_20.pdf", height=4.5,width=7.5, units="in",bg="white",useDingbats=FALSE)
    #supplemental figure 20


#PCoA with unweighted Unifrac , supplemental figure 21

out.uf.log <- ordinate(pslog, method = "PCoA", distance ="unifrac")
evals <- out.uf.log$values$Eigenvalues
plot_ordination(pslog, out.uf.log, color = "worms_present") + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot() + scale_color_manual( values = c("#fabb01", "#0081c7")) + theme_cowplot() + ggtitle("PCoA, unweighted Unifrac")


ord_df <-plot_ordination(pslog, out.uf.log, color = "foundress_number", shape="worms_present", justDF = TRUE)


ord_df$foundress_number <- as.numeric(ord_df$foundress_number)


ggplot(ord_df, aes(x=Axis.1,y=Axis.2)) + geom_point(aes(size = foundress_number, shape=worms_present, colour=foundress_number),alpha=0.6) + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot()  + xlab("Axis 1 (15.3%)") + ylab("Axis 2 (8.1%)")  + scale_size_continuous(limits=c(0,20),breaks=c(0,5,10,15,20)) + scale_color_continuous(limits=c(0,20),breaks=c(0,5,10,15,20)) + guides(color=guide_legend(title="Foundress\nnumber"), size = guide_legend(title="Foundress\nnumber"), shape=guide_legend(title="Worms\npresent?")) + theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),legend.title = element_text(size = 14)) + ggtitle("PCoA, weighted Unifrac")



ggsave("supplemental_figure_21.pdf", height=4.5,width=7.5, units="in",bg="white",useDingbats=FALSE)


#CCA  , supplemental figure 22

out.cca.log <- ordinate(pslog, method = "CCA")

plot_ordination(pslog, out.cca.log , color = "worms_present") + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot() + scale_color_manual( values = c("#fabb01", "#0081c7"))  + theme_cowplot() + ggtitle("Canonical Correspondence Analysis")


ord_df <-plot_ordination(pslog, out.cca.log, color = "foundress_number", shape="worms_present", justDF = TRUE)


ord_df$foundress_number <- as.numeric(ord_df$foundress_number)


ggplot(ord_df, aes(x=CA1,y=CA2)) + geom_point(aes(size = foundress_number, shape=worms_present, colour=foundress_number),alpha=0.6) + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot()  + xlab("Axis 1 (6.6%)") + ylab("Axis 2 (5.3%)")  + scale_size_continuous(limits=c(0,20),breaks=c(0,5,10,15,20)) + scale_color_continuous(limits=c(0,20),breaks=c(0,5,10,15,20)) + guides(color=guide_legend(title="Foundress\nnumber"), size = guide_legend(title="Foundress\nnumber"), shape=guide_legend(title="Worms\npresent?")) + theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),legend.title = element_text(size = 14)) + ggtitle("Canonical Correspondence Analysis")


ggsave("supplemental_figure_22.pdf", height=4.5,width=7.5, units="in",bg="white",useDingbats=FALSE)


#write.table(table(tax_table(ps3)[, "Family"]), "family_table.tsv",sep='\t')



#####rarefaction and diversity

##per-sample summary stats
#summary(rawreads2$num_paired_end_reads)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##   2983   98301  157812  140272  186237  273672


#Subsample reads
ps_250 <- vegan::rarefy(ps3, 250)



otu.rare = otu_table(ps3)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)

# we will use vegan rarecurve 
otu.rarecurve = vegan::rarecurve(otu.rare, step = 500, label = T)

str(otu.rarecurve)


otu.rarecurve = vegan::rarecurve(otu.rare, step = 500, label = T,tidy = TRUE)


a <- ggplot(otu.rarecurve, aes(x=Sample,y=Species,group=Site)) + geom_line() + ylab("Number of OTU") + xlab("Number of sampled reads")


b <- ggplot(otu.rarecurve, aes(x=Sample,y=Species,group=Site)) + geom_line() +geom_vline(xintercept=12252,colour="red",linetype="dotted") +xlim(0,25000) + ylab("Number of OTU") + xlab("Number of sampled reads")

a
ggsave("supplemental_figure_rarefaction_all.png",bg="white",height=7,width=6,units="in")
    #Supplemental Figure 23
b
ggsave("supplemental_figure_rarefaction_25000.png",bg="white",height=7,width=6,units="in")
    #Supplemental Figure 24


#maybe that's it
