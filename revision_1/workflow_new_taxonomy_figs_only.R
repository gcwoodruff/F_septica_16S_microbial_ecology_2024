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
setwd("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/")

#load OTU counts (see line 319 of workflow.sh), taxonomy table (line 333 of workflow.sh), and sample metadata
ps_object <- read_phyloseq(otu.file = "/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/07_qiime_export_dada2/exported-table/feature-table.csv", 
                    taxonomy.file = "/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/05_qiime_feature-classifier_classify-consensus-vsearch/export/taxonomy_revised_fig_samples_only_5-2025.tsv", 
                    metadata.file = "/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/fig_microbe_sample_taiwan_2019_metadata_e_5-7-25.csv", 
                    type = 'simple')

#load phylogenetic tree (line 326 of workflow.sh || line 290 of workflow.sh following revisions)
treefile <- read.tree("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/06_qiime_phylogeny_align-to-tree-mafft-fasttree/export/rooted-tree.qza/tree.tree")


ps.ng.tax <- merge_phyloseq(ps_object, treefile)

#get number of figs in this study
    #sample metadata
metadat <- read.csv("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/fig_microbe_sample_taiwan_2019_metadata_e_5-7-25_2.csv")

metadat$plant_label <- as.factor(metadat$plant_label)

surf <- metadat[metadat$surface_or_interior == "surface",]
inte <-metadat[metadat$surface_or_interior == "interior",]

as.data.frame(table(surf$plant_label))
#   Var1 Freq
#1     9    1
#2    10    1
#3    11    1
#4    13    4
#5    14    3
#6    15    5
#7    16    1
#8    17    3
#9    18    1
#10   19    3
#11   20    3
#12   21    6

as.data.frame(table(inte$plant_label))
#   Var1 Freq
#1     9    1
#2    10    1
#3    11    1
#4    13    6
#5    14    3
#6    15    7
#7    16    1
#8    17    5
#9    18    1
#10   19    3
#11   20    3
#12   21    6


as.data.frame(table(metadat$fig_field_label))
#2  10-1    2
#3  11-1    2
#4  13-1    2
#5  13-2    2
#6  13-3    2
#7  13-4    1
#8  13-5    1
#9  13-6    2
#10 14-1    2
#11 14-2    2
#12 14-3    2
#13 15-1    2
#14 15-2    2
#15 15-3    2
#16 15-4    1
#17 15-5    1
#18 15-6    2
#19 15-7    2
#20 16-1    2
#21 17-1    1
#22 17-2    2
#23 17-3    2
#24 17-4    1
#25 17-5    2
#26 18-1    2
#27 19-1    2
#28 19-2    2
#29 19-3    2
#30 20-1    2
#31 20-2    2
#32 20-3    2
#33 21-1    2
#34 21-2    2
#35 21-3    2
#36 21-4    2
#37 21-5    2
#38 21-6    2
#39  9-1    2

as.data.frame(table(inte$worms_present))

unique(metadat$fig_field_label)
length(unique(metadat$fig_field_label))-2
#38 figs

#get number of plants in this study

unique(metadat$plant_label)
length(unique(metadat$plant_label))-2
    #12 plants

metadat
as.data.frame(table(metadat$fig_stage))
#2           B   58
#3           C    6
#5           E    6



#foundress number

#   Var1 Freq
#1     0   11
#2     1    5
#3     2    8
#4     3    4
#5     4    4
#6     6    2
#7    11    1
#8    12    1
#9    15    1
#10   19    1


as.data.frame(table(inte$foundress_number))


rawreads <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/read_counts_one_per_sample.tsv",sep='\t', header=TRUE)
rawreads2 <- data.frame(sample_id = rawreads$sample_id, num_paired_end_reads=rawreads$num_paired_end_reads, step= "Raw reads")


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
phyloreads2 <- data.frame(sample_id = phyloreads$sample_id, num_paired_end_reads=phyloreads$num_paired_end_reads, step= "ASV clustering")





#how many paired-end reads total?
sum(phyloreads2$num_paired_end_reads)
#9459003
    ##9,459,003 paired-end reads after ASV clustering

#per-sample summary stats
summary(phyloreads2$num_paired_end_reads)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   2800   82032  127450  116778  153961  226560
	#116,778 paired-end reads on average per sample (after ASV clustering)

#remove mitochondrial and chloroplast reads

nochmi <- subset_taxa(ps.ng.tax, !Family %in% "Mitochondria" & !Order %in% "Chloroplast")

#reads after removing mito and chloro
orgoreads <- as.data.frame(sample_sums(nochmi))
orgoreads$sample_id <- rownames(orgoreads)
orgoreads$num_paired_end_reads <- orgoreads[,1]

orgoreads2 <- data.frame(sample_id = orgoreads$sample_id, num_paired_end_reads=orgoreads$num_paired_end_reads, step= "Remove organelles")


all_reads <- rbind(rawreads2,phyloreads2,orgoreads2)


all_reads$step <- factor(all_reads$step, levels=c("Raw reads","ASV clustering","Remove organelles"))

#merge df's to get per sample fraction of reads
m1 <- merge(rawreads2,phyloreads2,by="sample_id")
m2 <- merge(m1,orgoreads2,by="sample_id")
#get fraction of reads organellar
m2$fra_org <- 1-(m2$num_paired_end_reads/m2$num_paired_end_reads.y)

#total number of reads after ASV clustering
sum(m2$num_paired_end_reads.y)

#[1] 9459003



#total number of reads after removing organelles
sum(m2$num_paired_end_reads)
#[1] 5515057
#total number of organellar reads
sum(m2$num_paired_end_reads.y)-sum(m2$num_paired_end_reads)
#[1] 3943946

#summary stats per sample fraction organellar reads

summary(m2$fra_org)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.07677 0.31664 0.35908 0.59831 0.93693

#summary stats per reads after removing organelles

summary(m2$num_paired_end_reads)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   2800   37088   52301   68087   99902  194311

#make df for supplemental figure
rdata <- all_reads %>% group_by(step) %>% mutate(mean = mean(num_paired_end_reads))
#supplemental figure 1, read histograms
ggplot(rdata, aes(x = num_paired_end_reads)) + geom_histogram(colour="black",binwidth=6000) + facet_rep_wrap(~step, ncol=1) +theme_cowplot() +geom_vline(aes(xintercept = mean, group = step), colour = 'red') + xlab("Number of paired end reads per sample") + ylab("Frequency") +scale_y_continuous(limits=c(0,9), breaks=c(0:9)) + scale_x_continuous(breaks=c(0,50000,100000,150000,200000,250000,300000),labels = scales::comma)


ggsave("Supplemental_Figure_histograms_organellar_reads_revisions_5-13-25_FORMER_SUPP_FIG_1.pdf",bg="white",height=7,width=6,units="in",useDingbats=FALSE)


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
metadat <- read.csv("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/fig_microbe_sample_taiwan_2019_metadata_R_4-23-22.csv")
#merge with the fraction organelle 
mitmeta <- merge(metadat,mitodf,by="sample.id")

chlmitmeta <- merge(mitmeta,chlorodf,by="sample.id")

#melt, rename things
mchlmitmeta <- reshape2::melt(chlmitmeta,measure.vars = c("fra_mito","fra_chloro"))
levels(mchlmitmeta$variable) <- c("Mitochondria", "Chloroplast")
names(mchlmitmeta)[names(mchlmitmeta) == 'variable'] <- 'mito.chloro'
names(mchlmitmeta)[names(mchlmitmeta) == 'value'] <- 'fraction_reads'
summary(chlmitmeta$fra_mito_chloro)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.07677 0.31664 0.35908 0.59831 0.93693



#replace factor level
levels(mchlmitmeta$mito.chloro)[levels(mchlmitmeta$mito.chloro)=="Mitochondria"] <- "Mitochondrion"


#supplemental figure 2
ggplot(mchlmitmeta, aes(fill=mito.chloro, y=fraction_reads, x=reorder(sample.id, -fraction_reads))) + geom_bar(position="stack", stat="identity") +geom_hline(yintercept=0.31664,linetype="dashed",colour="black") + geom_hline(yintercept=0.35908,linetype="dotted",colour="black") + theme_cowplot() + theme(axis.text.x=element_blank()) + scale_fill_brewer(palette="Set1") +ylim(0,1) + xlab("Samples") + ylab("Fraction of reads") + labs(fill="Organelle")


ggsave("Supplemental_Figure_fraction_organellar_reads_revisions_5-13-25_FORMER_Supplemental_Figure_2.pdf",bg="white",height=6.5,width=7.5,units="in",useDingbats=FALSE)



#remove control ASVs


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
#[1] 125


#there are 125 ASVs in the controls that will ultimately be removed from the fig microbiome analysis.
#
nochmiprevdf = apply(X = otu_table(nochmi),
               MARGIN = ifelse(taxa_are_rows(controls), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
nochmiprevdf = data.frame(Prevalence = nochmiprevdf,
                    TotalAbundance = taxa_sums(nochmi),
                    tax_table(nochmi))

nrow(nochmiprevdf)
#[1] 2970
	#we lost some ASV's
125/2970
#[1] 0.04208754
#4.2% of all ASVs were seen at least once in the controls.

#these are the OTUs to remove
control_OTU$OTU

#remove control samples

no_controls <- prune_samples(!(sample_names(nochmi) %in% c("EXTRNEG","ExtrNeg1","EXTRPOS","GW10","GW11","GW31","GW9","PCRNeg1","PCRNeg2","PCRPos1","PCRPos2")), nochmi)



#remove ASV's found in controls and controls; remove control samples
ps3 <- subset_taxa(no_controls, !OTU %in% control_OTU$OTU)



#save.image(file = "figs_only_workspace_8-11-2025.RData")
#load("figs_only_workspace_8-11-2025.RData")

#write.table(as.data.frame(otu_table(ps3)),"figs_only_no_control_otu_table_8-2025.tsv",quote=FALSE,sep="\t")


#write.table(as.data.frame(tax_table(ps3)),"figs_only_no_control_tax_table_8-2025.tsv",quote=FALSE,sep="\t")




#ASV prevalence




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


#get number of ASV
nrow(prevdf)
#[1] 2845


summary(prevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.000   1.000   1.000   3.111   2.000  67.000
    
    #No ASV is present in all samples

#prevalence for ASVs histogram

ggplot(prevdf, aes(x = Prevalence)) + geom_histogram(colour="black",binwidth=1) + theme_cowplot() +geom_vline(xintercept = 3.111, colour = 'red') + ylab("ASV count") + scale_x_continuous(breaks=seq(0,70,5)) + scale_y_continuous(limits=c(0,2000),breaks=seq(0,2000,250))

ggsave("REVISIONS_ASV_prevalence_histogram_FORMER_Supplemental_Figure_3.pdf",bg="white",height=4.5,width=7.5,units="in",useDingbats=FALSE)

####return to this!
#########
####
####
#########
####
####
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
#  1.000   1.000   3.000   8.655   9.250  68.000

#No genera with 100% prevalence!
#What are the 68/70 prevalence genera?

ps3generaprevdf[ps3generaprevdf$Prevalence == 68,]$Genus
#[1] "Methylobacterium-Methylorubrum"
#[2] "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
#[3] "Sphingomonas"


#get number of genera
length(unique(ps3generaprevdf$Genus))
#[1] 373

#check for bad genera (i.e., "uncultured" or something)
#321 after removing dubious genera

ggplot(ps3generaprevdf, aes(x = Prevalence)) + geom_histogram(colour="black",binwidth=1) + theme_cowplot() +geom_vline(xintercept = 8.655, colour = 'red')  + ylab("Genus count") + scale_x_continuous(breaks=seq(0,70,5))  + scale_y_continuous(limits=c(0,150),breaks=seq(0,150,25))
ggsave("REVISIONS_genus_prevalence_histogram_FORMER+Supplemental_Figure_4.png",bg="white",height=4.5,width=7.5,units="in")


#Is Wolbachia in here?
ps3generaprevdf[ps3generaprevdf$Genus == "Wolbachia",]$Prevalence
#[1] 2

#only found in two samples

#What about other specific genera?
ps3generaprevdf[ps3generaprevdf$Genus == "Burkholderia",]$Prevalence
    #not present


ps3generaprevdf[ps3generaprevdf$Genus == "Ralstonia",]$Prevalence
    #not present



#combine counts to family level for all fig samples, fig suspensions, and fig surface washes
allfam <- tax_glom(ps3, "Family", NArm = TRUE)
prevdf = apply(X = otu_table(allfam),
               MARGIN = ifelse(taxa_are_rows(allfam), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(allfam),
                    tax_table(allfam))

write.table(prevdf, "figs_only_family_prevalence.tsv",quote=FALSE,sep='\t',row.names=FALSE)


#get number of families
length(unique(prevdf$Family))
#[1] 221

#after removing dubious families,
#185 families




interior_samples <- subset_samples(ps3, surface_or_interior == "interior")

surface_samples <- subset_samples(ps3, surface_or_interior == "surface")






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

cat(unique(plot2df$Family), sep = "\n")
#Spirosomaceae
#Hymenobacteraceae
#Chitinophagaceae
#Sphingobacteriaceae
#Weeksellaceae
#Xanthomonadaceae
#Oxalobacteraceae
#Comamonadaceae
#Enterobacteriaceae
#Erwiniaceae
#Geodermatophilaceae
#Kineosporiaceae
#Microbacteriaceae
#Nocardiaceae
#Beijerinckiaceae
#Acetobacteraceae
#Xanthobacteraceae
#Rhizobiaceae
#Rhodobacteraceae
#Sphingomonadaceae


#transform abundance for plotting
plot2df$LogAbundance <- log(1 + plot2df$TotalAbundance)

plot2df$Family <- as.factor(plot2df$Family)

length(unique(plot2df$Family))

unique(plot2df$Class)

length(unique(all_eighty$Family))
[1] 12


length(unique(int_eighty$Family))

length(unique(surf_eighty$Family))

plot2df[plot2df$Family == "Microbacteriaceae",]

#getting the fraction of ASVs classified as various phyla....


prevdf = apply(X = otu_table(ps3),
               MARGIN = ifelse(taxa_are_rows(ps3), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps3),
                    tax_table(ps3))

phy_df <- as.data.frame(table(prevdf$Phylum))

phy_df$perc_ASV <- (phy_df$Freq/2845)*100
phy_df
#                           Var1 Freq    perc_ASV
#1                                256  8.99824253
#2              Abditibacteriota   41  1.44112478
#3               Acidobacteriota   55  1.93321617
#4              Actinobacteriota  271  9.52548330
#5                Armatimonadota   50  1.75746924
#6                  Bacteroidota  521 18.31282953
#7              Bdellovibrionota  124  4.35852373
#8              Campilobacterota    2  0.07029877
#9                   Chloroflexi   52  1.82776801
#10                Crenarchaeota    3  0.10544815
#11                Cyanobacteria   45  1.58172232
#12                 Deinococcota   24  0.84358524
#13                 Dependentiae    2  0.07029877
#14             Desulfobacterota    2  0.07029877
#15              Elusimicrobiota    1  0.03514938
#16                Euryarchaeota    1  0.03514938
#17                   Firmicutes   73  2.56590510
#18               Fusobacteriota    5  0.17574692
#19              Gemmatimonadota   15  0.52724077
#20             Latescibacterota    1  0.03514938
#21                  Myxococcota   59  2.07381371
#22                 Nitrospirota    1  0.03514938
#23              Patescibacteria   13  0.45694200
#24              Planctomycetota  137  4.81546573
#25               Proteobacteria 1011 35.53602812
#26 SAR324_clade(Marine_group_B)    2  0.07029877
#27                  Sumerlaeota    3  0.10544815
#28            Verrucomicrobiota   71  2.49560633
#29                        WPS-2    4  0.14059754



##okay, ordinations for figs + controls ; Aitchison distances.......


ps_clr <- microbiome::transform(nochmi, "clr")

#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + geom_bar(stat="identity", fill = "blue") + labs(x = "\nAxis", y = "Proportion of Variance\n")

ggsave("REVISIONS_scree_plot_fig_samples_and_controls_ordination_Aitchison.png",bg="white",height=4.5,width=7.5,units="in")



#Examine eigenvalues and % prop. variance explained
head(ord_clr$CA$eig)                                                  

#      PC1       PC2       PC3       PC4       PC5       PC6
#333.39521 217.63127 193.65895 148.44787 102.99301  97.03873


sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))     

#       PC1        PC2        PC3        PC4        PC5
#0.10526087 0.06871141 0.06114278 0.04686855 0.03251737

clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(nochmi, ord_clr, color="control_interior_surface") + 
  geom_point(size = 1) +
  coord_fixed(clr2 / clr1)




ord_df <-plot_ordination(nochmi, ord_clr, color = "control_interior_surface", justDF = TRUE)


ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=control_interior_surface),size=1) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (10.5%)") + ylab("Axis 2 (6.9%)")  

ord_df$control_interior_surface <- as.factor(ord_df$control_interior_surface)

levels(ord_df$control_interior_surface)


ord_df$control_interior_surface <- factor(ord_df$control_interior_surface, levels=c("surface", "interior", "fig_negative", "extraction_negative", "pcr_negative","extraction_positive","pcr_positive"))


levels(ord_df$control_interior_surface)[match("surface",levels(ord_df$control_interior_surface))] <- "Fig surface wash"
levels(ord_df$control_interior_surface)[match("interior",levels(ord_df$control_interior_surface))] <- "Fig suspension"

levels(ord_df$control_interior_surface)[match("fig_negative",levels(ord_df$control_interior_surface))] <- "M9 buffer neg. ctrl."
levels(ord_df$control_interior_surface)[match("extraction_negative",levels(ord_df$control_interior_surface))] <- "DNA ext. neg. ctrl."
levels(ord_df$control_interior_surface)[match("pcr_negative",levels(ord_df$control_interior_surface))] <- "PCR neg. ctrl."
levels(ord_df$control_interior_surface)[match("extraction_positive",levels(ord_df$control_interior_surface))] <- "DNA ext. pos. ctrl."
levels(ord_df$control_interior_surface)[match("pcr_positive",levels(ord_df$control_interior_surface))] <- "PCR pos. ctrl."



levels(ord_df$control_interior_surface)[match("M9 Buffer neg. ctrl.",levels(ord_df$control_interior_surface))] <- "M9 buffer neg. ctrl."

summary(ord_df$PC1)
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
#-3.9093764 -1.5399061  0.0009561  0.0000000  1.7827019  6.4236353
summary(ord_df$PC2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-5.4697 -1.5858  0.4444  0.0000  1.0958  7.5676


ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=control_interior_surface),size=1) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (10.5%)") + ylab("Axis 2 (6.9%)")  + scale_colour_manual(values = c("tomato1","royalblue3","gray22","gray32","gray42","gray52","gray62")) +labs(colour="Sample Type") + scale_x_continuous(limits=c(-4,7),breaks=c(seq(-4,7,by=1))) + scale_y_continuous(limits=c(-6,8),breaks=c(seq(-6,8,by=1)))

ggsave("REVISIONS_fig_samples_controls_ordination_Aitchison.pdf",bg="white",height=4.5,width=7.5,units="in",useDingbats=FALSE)


#okay, remove controls



ps_clr <- microbiome::transform(ps3, "clr")

#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + geom_bar(stat="identity", fill = "blue") + labs(x = "\nAxis", y = "Proportion of Variance\n")

ggsave("REVISIONS_scree_plot_fig_samples_only_ordination_Aitchison.png",bg="white",height=4.5,width=7.5,units="in")



#Examine eigenvalues and % prop. variance explained
head(ord_clr$CA$eig)                                                  

#     PC1      PC2      PC3      PC4      PC5      PC6
#281.7568 225.7669 165.3672 160.8611 116.7727 105.2526


sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))     

#       PC1        PC2        PC3        PC4        PC5
#0.08689936 0.06963096 0.05100251 0.04961274 0.03601500

clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ps3, ord_clr, color="surface_or_interior") + 
  geom_point(size = 1) +
  coord_fixed(clr2 / clr1)




ord_df <-plot_ordination(ps3, ord_clr, color = "surface_or_interior", justDF = TRUE)


ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=surface_or_interior),size=1) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (8.7%)") + ylab("Axis 2 (7.0%)")  

ord_df$surface_or_interior <- as.factor(ord_df$surface_or_interior)

levels(ord_df$surface_or_interior)


ord_df$surface_or_interior <- factor(ord_df$surface_or_interior, levels=c("surface", "interior"))


levels(ord_df$surface_or_interior)[match("surface",levels(ord_df$surface_or_interior))] <- "Fig surface wash"
levels(ord_df$surface_or_interior)[match("interior",levels(ord_df$surface_or_interior))] <- "Fig suspension"


summary(ord_df$PC1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-3.0364 -1.7292 -0.9847  0.0000  0.6383  9.4621
summary(ord_df$PC2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-5.1997 -2.2181 -0.3562  0.0000  1.9946  4.6843


ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=surface_or_interior),size=1) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (8.7%)") + ylab("Axis 2 (7.0%)")  + scale_colour_manual(values = c("tomato1","royalblue3")) +labs(colour="Sample Type") + scale_x_continuous(limits=c(-4,10),breaks=c(seq(-4,10,by=1))) + scale_y_continuous(limits=c(-6,5),breaks=c(seq(-6,5,by=1)))

ggsave("REVISIONS_fig_samples_only_ordination_Aitchison.pdf",bg="white",height=4.5,width=7.5,units="in",useDingbats=FALSE)




a <- ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=surface_or_interior),size=2) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (8.7%)") + ylab("Axis 2 (7.0%)")  + scale_colour_manual(values = c("tomato1","royalblue3")) +labs(colour="Sample Type") + scale_x_continuous(limits=c(-4,10),breaks=c(seq(-4,10,by=1))) + scale_y_continuous(limits=c(-6,5),breaks=c(seq(-6,5,by=1)))


library(rcartocolor)

ord_df$plant_label <- as.factor(ord_df$plant_label)

clrs <- carto_pal(12, "Prism")

ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(shape=surface_or_interior,colour=plant_label),size=2) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (8.7%)") + ylab("Axis 2 (7.0%)")  + scale_color_manual(values=clrs) +labs(colour="Plant ID",shape="Sample Type") + scale_x_continuous(limits=c(-4,10),breaks=c(seq(-4,10,by=1))) + scale_y_continuous(limits=c(-6,5),breaks=c(seq(-6,5,by=1)))



ggsave("REVISIONS_fig_samples_only_ordination_Aitchison_surf_int_plant_labels.pdf",bg="white",height=4.5,width=7.5,units="in",useDingbats=FALSE)


b <- ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(shape=surface_or_interior,colour=plant_label),size=2) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (8.7%)") + ylab("Axis 2 (7.0%)")  + scale_color_manual(values=clrs) +labs(colour="Plant ID",shape="Sample Type") + scale_x_continuous(limits=c(-4,10),breaks=c(seq(-4,10,by=1))) + scale_y_continuous(limits=c(-6,5),breaks=c(seq(-6,5,by=1)))

#interior samples only, worm occupancy

ps4 <- subset_samples(ps3, surface_or_interior == "interior")



ps_clr <- microbiome::transform(ps4, "clr")

#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + geom_bar(stat="identity", fill = "blue") + labs(x = "\nAxis", y = "Proportion of Variance\n")

ggsave("REVISIONS_scree_plot_interior_fig_samples_only_ordination_Aitchison.png",bg="white",height=4.5,width=7.5,units="in")



#Examine eigenvalues and % prop. variance explained
head(ord_clr$CA$eig)                                                  

#     PC1      PC2      PC3      PC4      PC5      PC6
#420.5848 265.8617 211.3597 199.3661 156.7677 146.7852


sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))     

#       PC1        PC2        PC3        PC4        PC5
#0.13363104 0.08447139 0.06715464 0.06334395 0.04980929

clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ps4, ord_clr, color="worms_present") + 
  geom_point(size = 1) +
  coord_fixed(clr2 / clr1)




ord_df <-plot_ordination(ps4, ord_clr, color = "worms_present", justDF = TRUE)


ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=worms_present),size=1) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (13.4%)") + ylab("Axis 2 (8.4%)")  

ord_df$worms_present <- as.factor(ord_df$worms_present)

levels(ord_df$worms_present)


ord_df$worms_present <- factor(ord_df$worms_present, levels=c("yes", "no"))


levels(ord_df$worms_present)[match("yes",levels(ord_df$worms_present))] <- "Caenorhabditis animals present"
levels(ord_df$worms_present)[match("no",levels(ord_df$worms_present))] <- "Caenorhabditis animals absent"


summary(ord_df$PC1)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-10.3179  -0.4944   1.2901   0.0000   1.8063   2.4026
summary(ord_df$PC2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-8.1482 -1.6061  0.5889  0.0000  1.6623  9.3209


ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=worms_present),size=1) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (13.4%)") + ylab("Axis 2 (8.4%)")  + scale_color_manual(values = c("#fabb01", "#0081c7")) +labs(colour="Nematode occupancy") + scale_x_continuous(limits=c(-11,3),breaks=c(seq(-11,3,by=1))) + scale_y_continuous(limits=c(-9,10),breaks=c(seq(-9,10,by=1)))

ggsave("REVISIONS_fig_interior_samples_only_ordination_Aitchison.pdf",bg="white",height=4.5,width=7.5,units="in",useDingbats=FALSE)



c <- ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=worms_present),size=2) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (13.4%)") + ylab("Axis 2 (8.4%)")  + scale_color_manual(values = c("#fabb01", "#0081c7")) +labs(colour="Nematode occupancy") + scale_x_continuous(limits=c(-11,3),breaks=c(seq(-11,3,by=1))) + scale_y_continuous(limits=c(-9,10),breaks=c(seq(-9,10,by=1)))



ord_df$plant_label <- as.factor(ord_df$plant_label)
ord_df$foundress_number <- as.integer(ord_df$foundress_number)

ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=plant_label,shape=worms_present,size=foundress_number),alpha=0.85) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (13.4%)") + ylab("Axis 2 (8.4%)")  + scale_color_manual(values=clrs) +labs(colour="Plant ID",size="Foundress\nnumber",shape="Nematode\noccupancy") + scale_x_continuous(limits=c(-11,3),breaks=c(seq(-11,3,by=1))) + scale_y_continuous(limits=c(-9,10),breaks=c(seq(-9,10,by=1)))


ggsave("REVISIONS_fig_interior_samples_only_ordination_Aitchison_plant_id_foundress_number_worm_occupancy.pdf",bg="white",height=4.5,width=7.5,units="in",useDingbats=FALSE)



d <- ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=plant_label,shape=worms_present,size=foundress_number),alpha=0.85) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (13.4%)") + ylab("Axis 2 (8.4%)")  + scale_color_manual(values=clrs) +labs(colour="Plant ID",size="Foundress\nnumber",shape="Nematode\noccupancy") + scale_x_continuous(limits=c(-11,3),breaks=c(seq(-11,3,by=1))) + scale_y_continuous(limits=c(-9,10),breaks=c(seq(-9,10,by=1)))


library(patchwork)

(a+c)/((b+d))

ggsave("REVISIONS_provisional_figure_4_5-13-25.pdf",bg="white",height=10,width=15,units="in",useDingbats=FALSE)



#permanova
###all fig samples, fig samples only; interior v. exterior , plant_label )

sample_data(ps3)$plant_label <- as.factor(sample_data(ps3)$plant_label)

#transform
ps_clr <- microbiome::transform(ps3, "clr")
    #centered log-ratio
    #CLR transform applies a pseudocount of
     #min(relative abundance)/2 to exact zero relative abundance entries
     #in OTU table before taking logs.

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 


#ADONIS test all factors
allperma <- vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$surface_or_interior + phyloseq::sample_data(ps_clr)$plant_label + phyloseq::sample_data(ps_clr)$foundress_number + phyloseq::sample_data(ps_clr)$worms_present)
allperma
#         Df SumOfSqs      R2      F Pr(>F)
#Model    22   105830 0.47304 1.9178  0.001 ***
#Residual 47   117891 0.52696
#Total    69   223721 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


allperma <- vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$surface_or_interior + phyloseq::sample_data(ps_clr)$plant_label + phyloseq::sample_data(ps_clr)$foundress_number + phyloseq::sample_data(ps_clr)$worms_present,by="terms")
allperma
#                                                  Df SumOfSqs      R2      F
#phyloseq::sample_data(ps_clr)$surface_or_interior  1    11261 0.05033 4.4893
#phyloseq::sample_data(ps_clr)$plant_label         11    71227 0.31838 2.5815
#phyloseq::sample_data(ps_clr)$foundress_number     9    20246 0.09049 0.8968
#phyloseq::sample_data(ps_clr)$worms_present        1     3096 0.01384 1.2343
#Residual                                          47   117891 0.52696
#Total                                             69   223721 1.00000
#                                                  Pr(>F)
#phyloseq::sample_data(ps_clr)$surface_or_interior  0.001 ***
#phyloseq::sample_data(ps_clr)$plant_label          0.001 ***
#phyloseq::sample_data(ps_clr)$foundress_number     0.848
#phyloseq::sample_data(ps_clr)$worms_present        0.115
#Residual
#Total
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




allperma <- vegan::adonis2(clr_dist_matrix ~  phyloseq::sample_data(ps_clr)$plant_label + phyloseq::sample_data(ps_clr)$surface_or_interior + phyloseq::sample_data(ps_clr)$foundress_number + phyloseq::sample_data(ps_clr)$worms_present,by="terms")
allperma
#                                                  Df SumOfSqs      R2      F
#phyloseq::sample_data(ps_clr)$plant_label         11    71769 0.32080 2.6011
#phyloseq::sample_data(ps_clr)$surface_or_interior  1    10719 0.04791 4.2734
#phyloseq::sample_data(ps_clr)$foundress_number     9    20246 0.09049 0.8968
#phyloseq::sample_data(ps_clr)$worms_present        1     3096 0.01384 1.2343
#Residual                                          47   117891 0.52696
#Total                                             69   223721 1.00000
#                                                  Pr(>F)
#phyloseq::sample_data(ps_clr)$plant_label          0.001 ***
#phyloseq::sample_data(ps_clr)$surface_or_interior  0.001 ***
#phyloseq::sample_data(ps_clr)$foundress_number     0.855
#phyloseq::sample_data(ps_clr)$worms_present        0.125
#Residual
#Total
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


ps4 <- subset_samples(ps3, surface_or_interior == "interior")
sample_data(ps4)$plant_label <- as.factor(sample_data(ps4)$plant_label)




#transform
ps_clr <- microbiome::transform(ps4, "clr")

    #centered log-ratio
    #CLR transform applies a pseudocount of
     #min(relative abundance)/2 to exact zero relative abundance entries
     #in OTU table before taking logs.

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
#ADONIS test
intallperm <- vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$plant_label + phyloseq::sample_data(ps_clr)$foundress_number + phyloseq::sample_data(ps_clr)$worms_present)
intallperm
#         Df SumOfSqs      R2      F Pr(>F)
#Model    21    68922 0.59185 1.1048  0.228
#Residual 16    47530 0.40815
#Total    37   116452 1.00000



intallperm <- vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$plant_label + phyloseq::sample_data(ps_clr)$foundress_number + phyloseq::sample_data(ps_clr)$worms_present,by="terms")
intallperm

#                                               Df SumOfSqs      R2      F
#phyloseq::sample_data(ps_clr)$plant_label      11    50920 0.43726 1.5583
#phyloseq::sample_data(ps_clr)$foundress_number  9    15985 0.13727 0.5979
#phyloseq::sample_data(ps_clr)$worms_present     1     2016 0.01731 0.6787
#Residual                                       16    47530 0.40815
#Total                                          37   116452 1.00000
#                                               Pr(>F)
#phyloseq::sample_data(ps_clr)$plant_label       0.013 *
#phyloseq::sample_data(ps_clr)$foundress_number  0.999
#phyloseq::sample_data(ps_clr)$worms_present     0.898
#Residual
#Total
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$worms_present,by="terms")
#                                            Df SumOfSqs      R2      F Pr(>F)
#phyloseq::sample_data(ps_clr)$worms_present  1     4016 0.03448 1.2858  0.051 .
#Residual                                    36   112436 0.96552
#Total                                       37   116452 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#no plant 17


noplant17 <- subset_samples(ps3, plant_label != "17")



ps_clr <- microbiome::transform(noplant17, "clr")
    #centered log-ratio
    #CLR transform applies a pseudocount of
     #min(relative abundance)/2 to exact zero relative abundance entries
     #in OTU table before taking logs.

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 


#ADONIS test all factors
allperma <- vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$surface_or_interior + phyloseq::sample_data(ps_clr)$plant_label + phyloseq::sample_data(ps_clr)$foundress_number + phyloseq::sample_data(ps_clr)$worms_present)
#allperma
#         Df SumOfSqs      R2      F Pr(>F)
#Model    21    88093 0.51472 2.0203  0.001 ***
#Residual 40    83056 0.48528
#Total    61   171149 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


allperma <- vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$surface_or_interior + phyloseq::sample_data(ps_clr)$plant_label + phyloseq::sample_data(ps_clr)$foundress_number + phyloseq::sample_data(ps_clr)$worms_present,by="terms")
#allperma
#                                                  Df SumOfSqs      R2      F
#phyloseq::sample_data(ps_clr)$surface_or_interior  1    10781 0.06299 5.1923
#phyloseq::sample_data(ps_clr)$plant_label         10    56677 0.33116 2.7296
#phyloseq::sample_data(ps_clr)$foundress_number     9    17649 0.10312 0.9444
#phyloseq::sample_data(ps_clr)$worms_present        1     2985 0.01744 1.4377
#Residual                                          40    83056 0.48528
#Total                                             61   171149 1.00000
#                                                  Pr(>F)
#phyloseq::sample_data(ps_clr)$surface_or_interior  0.001 ***
#phyloseq::sample_data(ps_clr)$plant_label          0.001 ***
#phyloseq::sample_data(ps_clr)$foundress_number     0.736
#phyloseq::sample_data(ps_clr)$worms_present        0.023 *
#Residual
#Total
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#alpha diversity among individual plants


adiv <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps3, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps3, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps3)))), tree = phyloseq::phy_tree(ps3))[, 1],
  "plant_label" = phyloseq::sample_data(ps3)$plant_label)





adiv$sample_id <- rownames(adiv)
names(adiv)[names(adiv) == "Observed"] <- "Number_of_OTU"

adiv_melt <- reshape2::melt(adiv, id.vars= c("sample_id","plant_label"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic"))

levels(adiv_melt$variable)[levels(adiv_melt$variable)=="Number_of_OTU"] <- "Number of ASV"

plant_label_plot_df <- adiv_melt

ggplot(plant_label_plot_df, aes(x=plant_label,y=value)) + geom_sina(scale="width",size=0.5) + stat_summary(aes(group=plant_label),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + facet_wrap(~variable,nrow=1,scales="free") +theme_cowplot(font_size=9) + background_grid(major = c("x"),minor = c("none"))+ scale_y_continuous(limits = c(0,NA)) + xlab("Individual Plant ID") +ylab("Diversity measure") + theme(strip.background = element_rect(colour="white", fill="white"))

ggsave("REVISIONS_supplemental_figure_plant_id_alpha_diversity_9-4-2025.pdf",bg="white",height=4,width=8,units="in",useDingbats=FALSE)

summary(aov(Number_of_OTU ~ plant_label,data=adiv))
#            Df Sum Sq Mean Sq F value   Pr(>F)
#plant_label 11 498007   45273   9.822 9.69e-10 ***
#Residuals   58 267345    4609
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


summary(aov(Shannon ~ plant_label,data=adiv))
#            Df Sum Sq Mean Sq F value   Pr(>F)
#plant_label 11  31.10  2.8270   4.804 3.22e-05 ***
#Residuals   58  34.13  0.5885
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(aov(Phylogenetic ~ plant_label,data=adiv))
#            Df Sum Sq Mean Sq F value   Pr(>F)
#plant_label 11  737.5   67.04   10.97 1.36e-10 ***
#Residuals   58  354.6    6.11
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




#alpha diversity, surface v interior, alright...






adiv <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps3, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps3, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps3)))), tree = phyloseq::phy_tree(ps3))[, 1],
  "surface_or_interior" = phyloseq::sample_data(ps3)$surface_or_interior)

adiv$sample_id <- rownames(adiv)
names(adiv)[names(adiv) == "Observed"] <- "Number_of_OTU"

adiv_melt <- reshape2::melt(adiv, id.vars= c("sample_id","surface_or_interior"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic"))

levels(adiv_melt$variable)[levels(adiv_melt$variable)=="Number_of_OTU"] <- "Number of ASV"

adiv_melt$surface_or_interior <- as.factor(adiv_melt$surface_or_interior)

levels(adiv_melt$surface_or_interior)[levels(adiv_melt$surface_or_interior)=="surface"] <- "Surface washes"
levels(adiv_melt$surface_or_interior)[levels(adiv_melt$surface_or_interior)=="interior"] <- "Suspensions"

#supplemental figure 15
ggplot(adiv_melt, aes(x=surface_or_interior,y=value)) + geom_sina(scale="width") + stat_summary(aes(group=surface_or_interior),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + facet_wrap(~variable,nrow=1,scales="free") +theme_cowplot() + scale_y_continuous(limits = c(0,NA)) + theme(strip.background = element_rect(colour="white", fill="white"))  + xlab("Sample type") +ylab("Diversity metric value")

ggsave("REVISIONS_supplemental_figure_16_alpha_diversity_fig_surface_wash_v_fig_suspension_no_rarefaction.png",bg="white",height=6,width=9,units="in")


#diversity statistical tests


div.num.otu <- adiv_melt[adiv_melt$variable == "Number of ASV",]
div.shannon <- adiv_melt[adiv_melt$variable == "Shannon",]
div.phylogenetic <- adiv_melt[adiv_melt$variable == "Phylogenetic",]



#num otu

wilcox.test(div.num.otu[div.num.otu$surface_or_interior == "Suspensions",]$value, div.num.otu[div.num.otu$surface_or_interior == "Surface washes",]$value)

#    Wilcoxon rank sum test with continuity correction
#
#data:  div.num.otu[div.num.otu$surface_or_interior == "Suspensions", ]$value and div.num.otu[div.num.otu$surface_or_interior == "Surface washes", ]$value
#W = 698, p-value = 0.2913
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(div.num.otu[div.num.otu$surface_or_interior == "Surface washes",]$value, div.num.otu[div.num.otu$surface_or_interior == "Suspensions",]$value)
#    Wilcoxon rank sum test with continuity correction
#
#data:  div.num.otu[div.num.otu$surface_or_interior == "Surface washes", ]$value and div.num.otu[div.num.otu$surface_or_interior == "Suspensions", ]$value
#W = 518, p-value = 0.2913
#alternative hypothesis: true location shift is not equal to 0



#shannon
wilcox.test(div.shannon[div.shannon$surface_or_interior == "Suspensions",]$value, div.shannon[div.shannon$surface_or_interior == "Surface washes",]$value)
#    Wilcoxon rank sum exact test
#
#data:  div.shannon[div.shannon$surface_or_interior == "Suspensions", ]$value and div.shannon[div.shannon$surface_or_interior == "Surface washes", ]$value
#W = 574, p-value = 0.6947
#alternative hypothesis: true location shift is not equal to 0


wilcox.test(div.shannon[div.shannon$surface_or_interior == "Surface washes",]$value, div.shannon[div.shannon$surface_or_interior == "Suspensions",]$value)

#    Wilcoxon rank sum exact test
#
#data:  div.shannon[div.shannon$surface_or_interior == "Surface washes", ]$value and div.shannon[div.shannon$surface_or_interior == "Suspensions", ]$value
#W = 642, p-value = 0.6947
#alternative hypothesis: true location shift is not equal to 0



#phylogenetic
wilcox.test(div.phylogenetic[div.phylogenetic$surface_or_interior == "Surface washes",]$value, div.phylogenetic[div.phylogenetic$surface_or_interior == "Suspensions",]$value)

#    Wilcoxon rank sum exact test
#
#data:  div.phylogenetic[div.phylogenetic$surface_or_interior == "Surface washes", ]$value and div.phylogenetic[div.phylogenetic$surface_or_interior == "Suspensions", ]$value
#W = 456, p-value = 0.07396
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(div.phylogenetic[div.phylogenetic$surface_or_interior == "Suspensions",]$value, div.phylogenetic[div.phylogenetic$surface_or_interior == "Surface washes",]$value)

#    Wilcoxon rank sum exact test
#
#data:  div.phylogenetic[div.phylogenetic$surface_or_interior == "Suspensions", ]$value and div.phylogenetic[div.phylogenetic$surface_or_interior == "Surface washes", ]$value
#W = 760, p-value = 0.07396
#alternative hypothesis: true location shift is not equal to 0





#worm occupancy

adiv <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps4, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps4, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps4)))), tree = phyloseq::phy_tree(ps4))[, 1],
  "worms_present" = phyloseq::sample_data(ps4)$worms_present)

adiv$sample_id <- rownames(adiv)
names(adiv)[names(adiv) == "Observed"] <- "Number_of_OTU"

adiv_melt <- reshape2::melt(adiv, id.vars= c("sample_id","worms_present"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic"))

levels(adiv_melt$variable)[levels(adiv_melt$variable)=="Number_of_OTU"] <- "Number of ASV"

ggplot(adiv_melt, aes(x=worms_present,y=value)) + geom_sina(scale="width") + stat_summary(aes(group=worms_present),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + facet_wrap(~variable,nrow=1,scales="free") +theme_cowplot() + scale_y_continuous(limits = c(0,NA)) + xlab("Worms Present?") +ylab("Diversity measure") + theme(strip.background = element_rect(colour="white", fill="white"))

ggsave("worms_present_diversity_sina_2024.pdf", height=4,width=7, units="in",bg="white",useDingbats=FALSE)
    #this is figure 5

#diversity statistical tests

div.num.otu <- adiv_melt[adiv_melt$variable == "Number of ASV",]
div.shannon <- adiv_melt[adiv_melt$variable == "Shannon",]
div.phylogenetic <- adiv_melt[adiv_melt$variable == "Phylogenetic",]




#num otu

wilcox.test(div.num.otu[div.num.otu$worms_present == "yes",]$value, div.num.otu[div.num.otu$worms_present == "no",]$value)

#    Wilcoxon rank sum test with continuity correction
#
#data:  div.num.otu[div.num.otu$worms_present == "yes", ]$value and div.num.otu[div.num.otu$worms_present == "no", ]$value
#W = 161.5, p-value = 0.589
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(div.shannon[div.shannon$worms_present == "yes",]$value, div.shannon[div.shannon$worms_present == "no",]$value)

#    Wilcoxon rank sum exact test
#
#data:  div.shannon[div.shannon$worms_present == "yes", ]$value and div.shannon[div.shannon$worms_present == "no", ]$value
#W = 157, p-value = 0.5063
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(div.phylogenetic[div.phylogenetic$worms_present == "yes",]$value, div.phylogenetic[div.phylogenetic$worms_present == "no",]$value)

#    Wilcoxon rank sum exact test
#
#data:  div.phylogenetic[div.phylogenetic$worms_present == "yes", ]$value and div.phylogenetic[div.phylogenetic$worms_present == "no", ]$value
#W = 165, p-value = 0.665
#alternative hypothesis: true location shift is not equal to 0




adiv <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps4, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps4, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps4)))), tree = phyloseq::phy_tree(ps4))[, 1],
  "foundress_number" = phyloseq::sample_data(ps4)$foundress_number,
  "plant_label" = phyloseq::sample_data(ps4)$plant_label)

num.otu <- aggregate(Observed ~ plant_label, FUN=mean,data=adiv)
shann <- aggregate(Shannon ~ plant_label, FUN=mean,data=adiv)
phylog <- aggregate(Phylogenetic ~ plant_label, FUN=mean,data=adiv)

adiv$foundress_number <- as.numeric(adiv$foundress_number)

found_num <- aggregate(foundress_number ~ plant_label, FUN=mean,data=adiv)

merg.otu <- merge(num.otu,found_num,by="plant_label")
merg.shann <- merge(shann,found_num,by="plant_label")
merg.phylog <- merge(phylog,found_num,by="plant_label")

a <- ggplot(merg.otu,aes(x=foundress_number,y=Observed)) + geom_point() +theme_cowplot() +xlab("Foundress Number") + ylab("Number of ASV") + ylim(0,400) + ggtitle("a")


b <- ggplot(merg.shann,aes(x=foundress_number,y=Shannon)) + geom_point()  +theme_cowplot() +xlab("Foundress Number") + ylab("Shannon diversity") + ylim(0,5)+ ggtitle("b")



c <- ggplot(merg.phylog,aes(x=foundress_number,y=Phylogenetic)) + geom_point()  +theme_cowplot() +xlab("Foundress Number") + ylab("Phylogenetic diversity") + ylim(0,17)+ ggtitle("c")


summary(lm(Observed ~ foundress_number,data=merg.otu))
#Call:
#lm(formula = Observed ~ foundress_number, data = merg.otu)
#
#Residuals:
#    Min      1Q  Median      3Q     Max
#-129.48  -51.34  -30.54   40.22  205.53
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)        205.54      44.01   4.670 0.000881 ***
#foundress_number   -19.84      13.29  -1.493 0.166395
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 102.5 on 10 degrees of freedom
#Multiple R-squared:  0.1822,    Adjusted R-squared:  0.1004
#F-statistic: 2.228 on 1 and 10 DF,  p-value: 0.1664


summary(lm(Shannon ~ foundress_number,data=merg.shann))
#Call:
#lm(formula = Shannon ~ foundress_number, data = merg.shann)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-1.56372 -0.28593  0.06182  0.55079  1.08518
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)        3.2091     0.3491   9.194 3.41e-06 ***
#foundress_number  -0.0736     0.1054  -0.698    0.501
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.813 on 10 degrees of freedom
#Multiple R-squared:  0.0465,    Adjusted R-squared:  -0.04885
#F-statistic: 0.4877 on 1 and 10 DF,  p-value: 0.5009


summary(lm(Phylogenetic ~ foundress_number,data=merg.phylog))
#Call:
#lm(formula = Phylogenetic ~ foundress_number, data = merg.phylog)
#
#Residuals:
#    Min      1Q  Median      3Q     Max
#-4.9049 -1.5174 -0.0192  0.9100  6.9335
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)       10.3350     1.4362   7.196 2.94e-05 ***
#foundress_number  -0.8802     0.4337  -2.030   0.0698 .
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3.345 on 10 degrees of freedom
#Multiple R-squared:  0.2918,    Adjusted R-squared:  0.221
#F-statistic:  4.12 on 1 and 10 DF,  p-value: 0.06982

a+b+c

ggsave("REVISIONS_alpha_diversity_foundress_number_averaged_by_plant.pdf", height=4,width=7, units="in",bg="white",useDingbats=FALSE)



adiv$sample_id <- rownames(adiv)
names(adiv)[names(adiv) == "Observed"] <- "Number_of_OTU"

adiv_melt <- reshape2::melt(adiv, id.vars= c("sample_id","plant_label"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic","foundress_number"))

levels(adiv_melt$variable)[levels(adiv_melt$variable)=="Number_of_OTU"] <- "Number of ASV"


adiv_melt_fn <- adiv_melt[adiv_melt$value == "foundress_number"]




ggplot(adiv, aes(x=plant_label,y=foundress_number)) + geom_sina(scale="width",size=0.5) + stat_summary(aes(group=plant_label),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) +theme_cowplot(font_size=9) + background_grid(major = c("x"),minor = c("none"))+ scale_y_continuous(limits = c(0,NA)) + xlab("Individual Plant ID") +ylab("Foundress number") + theme(strip.background = element_rect(colour="white", fill="white")) +ylim(-1,20)


ggsave("REVISIONS_foundress_number_by_plant.pdf", height=4,width=7, units="in",bg="white",useDingbats=FALSE)

summary(aov(foundress_number ~ plant_label,data=adiv))

#> summary(aov(foundress_number ~ plant_label,data=adiv))
#            Df Sum Sq Mean Sq F value Pr(>F)
#plant_label 11  314.7   28.61   1.963 0.0769 .
#Residuals   26  378.9   14.57
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#co-occurence network analysis stuff...
#https://netcomi.de/articles/netcomi

library(NetCoMi)
library(phyloseq)

#okay, let's do it at the genus level.
#including exterior and interior to begin




####
ps3genera <- tax_glom(ps3, "Genus", NArm = TRUE)



ps3generaprevdf = apply(X = otu_table(ps3genera),
               MARGIN = ifelse(taxa_are_rows(ps3genera), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
ps3generaprevdf = data.frame(Prevalence = ps3generaprevdf,
                    TotalAbundance = taxa_sums(ps3genera),
                    tax_table(ps3genera))

ps3genera_renamed <- renameTaxa(ps3genera, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")



taxtabr <- tax_table(ps3genera_renamed)

net_spring <- netConstruct(ps3genera_renamed,
                           taxRank = "Genus",
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 50),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           measure = "aitchison",
                           measurePar = list(nlambda=10, 
                                             rep.num=10,
                                             Rmethod = "approx"),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 2,
                           seed = 666)

head(net_spring$edgelist1)

   v1   v2      diss      adja
1 GW1 GW12 0.5292302 0.4707698
2 GW1 GW13 0.4639369 0.5360631
3 GW1 GW14 0.2072251 0.7927749
4 GW1 GW15 0.2794186 0.7205814
5 GW1 GW16 0.7276904 0.2723096
6 GW1 GW17 0.7849122 0.2150878


str(net_spring)



plotHeat(net_spring$assoMat1, textUpp = "none", textLow = "none")


####trying this again........... https://chiliubio.github.io/microeco_tutorial/meconetcomp-package.html

library(microeco)
library(meconetcomp)
# use pipe operator in magrittr package
library(magrittr)
library(igraph)
library(ggplot2)
theme_set(theme_bw())
# load soil amplicon sequencing dataset
data(soil_amp)

#okay, need my data in the microtable format...
ps3genera@otu_table

taxtab <- as.data.frame(tax_table(ps3genera))

otutab <- as.data.frame(ps3genera@otu_table)

infotab <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/fig_microbe_sample_taiwan_2019_metadata_e_5-7-25_4.tsv", header = TRUE, sep = "\t", row.names = 1,
               as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
               check.names = FALSE)

infotab$SampleID <- rownames(infotab)

infotab <- subset(infotab, select=c(22,1:21))

taxtab %<>% tidy_taxonomy


mt <- microtable$new(sample_table = infotab, otu_table = otutab, tax_table = taxtab, phylo_tree = treefile)



# use clone to get a deep copy of soil_amp (R6 object)
tmp <- clone(mt)

# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
    #the pairwise correlation coefficients are in tmp$res_cor_p as a list! (see manual entry for trans_network())


# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
# put the network into the list

mt_network <- list()


mt_network$all <- tmp




tmp <- clone(mt)
# change sample_table directly
tmp$sample_table %<>% subset(surface_or_interior == "surface")
# trim all files in the object
tmp$tidy_dataset()
# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
# put the network into the list
mt_network$surface <- tmp
# select samples of "TW" group
tmp <- clone(mt)
tmp$sample_table %<>% subset(surface_or_interior == "interior")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
mt_network$interior <- tmp


#11.1 Network modularity for all networks

#The function cal_module in meconetcomp package is designed to partition modules for all the networks in the list.

mt_network %<>% cal_module(undirected_method = "cluster_fast_greedy")



#11.2 Network topological attributes for all networks

#we extracted all the res_network_attr tables in the networks and merged them into one final table by using cal_network_attr function in meconetcomp package.

tmp <- cal_network_attr(mt_network)
# tmp is a data.frame object

write.table(tmp,"meconetcomp_network_attributes_figs_only.tsv",quote=FALSE,sep="\t")

#11.3 Node and edge properties extraction for all networks
#
#The get_node_table and get_edge_table functions of meconetcomp package can be used to directly extract node and edge properties for all the networks. The return table is stored in each network object.

mt_network %<>% get_node_table(node_roles = TRUE) %>% get_edge_table

#11.4 Compare nodes across networks
#
#The nodes in all the networks can be converted to a new microtable object by using the node_comp function of meconetcomp package. Then, it is easy to analyse the nodes overlap with trans_venn class.
# obtain the node distributions by searching the res_node_table in the object
tmp <- node_comp(mt_network, property = "name")
# obtain nodes intersection
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g1 <- tmp1$plot_venn(fill_color = FALSE)
ggsave("meconetcomp_network_attributes_figs_only_node_overlap.pdf", g1, width = 7, height = 6)
# calculate jaccard distance to reflect the overall differences of networks
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard


#               all   surface  interior
#all      0.0000000 0.4067797 0.3709677
#surface  0.4067797 0.0000000 0.4782609
#interior 0.3709677 0.4782609 0.0000000


#11.5 Compare edges across networks
#
#The pipeline of studying edges overlap is similar with the above operations of nodes comparison. The edge_comp function of meconetcomp package is used to convert edges distribution to a new microtable object.

# get the edge distributions across networks
tmp <- edge_comp(mt_network)
# obtain edges intersection
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g1 <- tmp1$plot_venn(fill_color = FALSE)
ggsave("meconetcomp_network_attributes_figs_only_edge_overlap.pdf", g1, width = 7, height = 6)
# calculate jaccard distance
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard

#               all   surface  interior
#all      0.0000000 0.7153285 0.7991361
#surface  0.7153285 0.0000000 0.9313929
#interior 0.7991361 0.9313929 0.0000000


11.6 Extract overlapped edges of networks to a new network

Then we extracted the subset of edges according to the intersections of edges across networks, which can be accomplished with the subset_network function in meconetcomp package.

# first obtain edges distribution and intersection
tmp <- edge_comp(mt_network)
tmp1 <- trans_venn$new(tmp)
# convert intersection result to a microtable object
tmp2 <- tmp1$trans_comm()
# extract the intersection of all the three networks
# please use colnames(tmp2$otu_table) to find the required name
Intersec_all <- subset_network(mt_network, venn = tmp2, name = "all")
# Intersec_all is a trans_network object
# for example, save Intersec_all as gexf format
Intersec_all$save_network("meconetcomp_network_attributes_figs_only_Intersec_all.gexf")




#11.7 Compare phylogenetic distances of paired nodes in edges
#
#The edge_node_distance class (R6 class) in meconetcomp package is designed to compare the distribution of distance values of paired nodes in all the edges across networks. Here, we indicated the phylogenetic distance distributions and performed the differential test among networks. The input parameter dis_matrix can be any symmetric matrix with both the column names and row names (i.e. feature names). So it is also feasible to compare other properties of features, such as Levin’s niche overlap.

# filter useless features to speed up the calculation
node_names <- unique(unlist(lapply(mt_network, function(x){colnames(x$data_abund)})))
mt_amp <- microeco::clone(mt)
mt_amp$otu_table <- mt_amp$otu_table[node_names, ]
mt_amp$tidy_dataset()
# obtain phylogenetic distance matrix
phylogenetic_distance_mt <- as.matrix(cophenetic(mt_amp$phylo_tree))
# use both the positive and negative labels
tmp <- edge_node_distance$new(network_list = mt_network, dis_matrix = phylogenetic_distance_mt, label = c("+", "-"))
tmp$cal_diff(method = "anova")
# visualization
g1 <- tmp$plot(add = "none", add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")
g1
ggsave("meconetcomp_network_attributes_figs_only_phylo_distance.pdf", g1, width = 7, height = 6)

# show different modules with at least 10 nodes and positive edges
tmp <- edge_node_distance$new(network_list = mt_network, dis_matrix = phylogenetic_distance_mt, 
    label = "+", with_module = TRUE, module_thres = 10)
tmp$cal_diff(method = "anova")
g1 <- tmp$plot(add = "none", add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")
ggsave("meconetcomp_network_attributes_figs_only_phylo_distance_modules.pdf", g1, width = 8, height = 6)

#11.8 Compare node sources of edges across networks
#
#To know which taxa constitute the nodes in edges is important in understanding species co-occurrence patterns and answering ecological questions. In this part, as an instance, we used edge_tax_comp function of meconetcomp package to get the sums of node sources (at Phylum level) in the positive edges. In other words, how many linked nodes of positive edges come from different phyla or the same phyla. Then, to make the results comparable, the ratio was calculated with the positive edge number as denominator.

mt_network_edgetax <- edge_tax_comp(mt_network, taxrank = "Phylum", label = "+", rel = TRUE)
# filter the features with small number
mt_network_edgetax <- mt_network_edgetax[apply(mt_network_edgetax, 1, mean) > 0.01, ]
# visualization
g1 <- pheatmap::pheatmap(mt_network_edgetax, display_numbers = TRUE)
ggsave("meconetcomp_network_attributes_figs_only_edge_tax_comp.pdf", g1, width = 7, height = 7)


mt_network_edgetax <- edge_tax_comp(mt_network, taxrank = "Genus", label = "+", rel = TRUE)
# filter the features with small number
mt_network_edgetax <- mt_network_edgetax[apply(mt_network_edgetax, 1, mean) > 0.01, ]
# visualization
g1 <- pheatmap::pheatmap(mt_network_edgetax, display_numbers = TRUE)
ggsave("meconetcomp_network_attributes_figs_only_edge_tax_comp_Genus.pdf", g1, width = 7, height = 7)
#######
#trying something else..........


str(ps3genera)

head(ps3genera@otu_table)
head(ps3genera@tax_table)

ps3genera_otu_table <- as.data.frame(ps3genera@otu_table)
ps3genera_tax_table <- as.data.frame(ps3genera@tax_table)

ps3genera_otu_table$asvID <- rownames(ps3genera_otu_table)
ps3genera_tax_table$asvID <- rownames(ps3genera_tax_table)


ps3genera_df_merge <- merge(ps3genera_otu_table,ps3genera_tax_table)

#get rid of genera with prevalence of one.....

genera_prev_one <- ps3generaprevdf[ps3generaprevdf$Prevalence < 2,]$Genus

ps3genera_df_merge_no_prev_one <- ps3genera_df_merge[!ps3genera_df_merge$Genus %in% genera_prev_one, ]


intsamps <- metadat[metadat$surface_or_interior == "interior",]$sample.id

surfsamps <- metadat[metadat$surface_or_interior == "surface",]$sample.id

ps3genera_df_merge_no_prev_one_SURF <- ps3genera_df_merge_no_prev_one[,names(ps3genera_df_merge_no_prev_one)[(names(ps3genera_df_merge_no_prev_one) %in% surfsamps)]]

ps3genera_df_merge_no_prev_one_INT <- ps3genera_df_merge_no_prev_one[,names(ps3genera_df_merge_no_prev_one)[(names(ps3genera_df_merge_no_prev_one) %in% intsamps)]]

ps3genera_df_merge_no_prev_one_SURFb <- cbind(ps3genera_df_merge_no_prev_one$asvID,ps3genera_df_merge_no_prev_one_SURF,ps3genera_df_merge_no_prev_one[,c(72:79)])


ps3genera_df_merge_no_prev_one_INTb <- cbind(ps3genera_df_merge_no_prev_one$asvID,ps3genera_df_merge_no_prev_one_INT,ps3genera_df_merge_no_prev_one[,c(72:79)])

#ps3genera_df_merge_no_prev_one_all_print <- ps3genera_df_merge_no_prev_one[,c(1,77,2:71)]
#ps3genera_df_merge_no_prev_one_SURFb_print <- ps3genera_df_merge_no_prev_one_SURFb[,c(1,39,2:33)]
#ps3genera_df_merge_no_prev_one_INTb_print <- ps3genera_df_merge_no_prev_one_INTb[,c(1,45,2:39)]

#names(ps3genera_df_merge_no_prev_one_SURFb_print)[names(ps3genera_df_merge_no_prev_one_SURFb_print) == 'ps3genera_df_merge_no_prev_one$asvID'] <- 'asvID'
#names(ps3genera_df_merge_no_prev_one_INTb_print)[names(ps3genera_df_merge_no_prev_one_INTb_print) == 'ps3genera_df_merge_no_prev_one$asvID'] <- 'asvID'


ps3genera_df_merge_no_prev_one_all_print <- ps3genera_df_merge_no_prev_one[,c(77,2:71)]
ps3genera_df_merge_no_prev_one_SURFb_print <- ps3genera_df_merge_no_prev_one_SURFb[,c(39,2:33)]
ps3genera_df_merge_no_prev_one_INTb_print <- ps3genera_df_merge_no_prev_one_INTb[,c(45,2:39)]


write.table(ps3genera_df_merge_no_prev_one_all_print,"/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/sparcc_tables/ps3genera_df_merge_no_prev_one_all_print.tsv",quote=FALSE,sep="\t",row.names=FALSE)

write.table(ps3genera_df_merge_no_prev_one_SURFb_print,"/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/sparcc_tables/ps3genera_df_merge_no_prev_one_SURFb_print.tsv",quote=FALSE,sep="\t",row.names=FALSE)

write.table(ps3genera_df_merge_no_prev_one_all_print,"/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/sparcc_tables/ps3genera_df_merge_no_prev_one_INTb_print.tsv",quote=FALSE,sep="\t",row.names=FALSE)





####
#okay, trying the tutorial data set

data("amgut1.filt") # ASV count matrix
data("amgut2.filt.phy") # phyloseq objext





amgut_genus <- tax_glom(amgut2.filt.phy, taxrank = "Rank6")

# Rename taxonomic table and make Rank6 (genus) unique
amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Rank6")


net_spring <- netConstruct(amgut_genus_renamed,
                           taxRank = "Rank6",
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 50),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           measure = "spring",
                           measurePar = list(nlambda=10, 
                                             rep.num=10,
                                             Rmethod = "approx"),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 2,
                           seed = 123456)
#okay... it works



net_spring <- netConstruct(ps3genera_renamed,
                           taxRank = "Genus",
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 50),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           measure = "spring",
                           measurePar = list(nlambda=10, 
                                             rep.num=10,
                                             Rmethod = "approx"),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 2,
                           seed = 666)



head(net_spring$edgelist1)


str(net_spring)
#1                                       Enterobacter     Enterococcus
#2 Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium     Ochrobactrum
#3                                   Chryseobacterium      Sphingobium
#4                                   Chryseobacterium Stenotrophomonas
#5                                        Sphingobium Aurantisolimonas
#6                     Methylobacterium-Methylorubrum     Sphingomonas
#        asso      diss      adja
#1 0.01335519 0.7023691 0.2976309
#2 0.10304288 0.6696854 0.3303146
#3 0.04458882 0.6911625 0.3088375
#4 0.06497269 0.6837497 0.3162503
#5 0.03679723 0.6939751 0.3060249
#6 0.23898259 0.6168539 0.3831461
    #okay, aitchinson option just doesn't work

plotHeat(net_spring$assoMat1, textUpp = "none", textLow = "none",type="upper")



props_spring <- netAnalyze(net_spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, 
                           normDeg = FALSE)


summary(props_spring, numbNodes = 5L)

#Component sizes
#```````````````
#size: 45 1
#   #:  1 5
#______________________________
#Global network properties
#`````````````````````````
#Largest connected component (LCC):
#
#Relative LCC size         0.90000
#Clustering coefficient    0.20400
#Modularity                0.56543
#Positive edge percentage 98.43750
#Edge density              0.06465
#Natural connectivity      0.02794
#Vertex connectivity       1.00000
#Edge connectivity         1.00000
#Average dissimilarity*    0.97940
#Average path length**     3.63355
#
#Whole network:
#
#Number of components      6.00000
#Clustering coefficient    0.20400
#Modularity                0.56543
#Positive edge percentage 98.43750
#Edge density              0.05224
#Natural connectivity      0.02467
#-----
#*: Dissimilarity = 1 - edge weight
#**: Path length = Units with average dissimilarity
#
#______________________________
#Clusters
#- In the whole network
#- Algorithm: cluster_fast_greedy
#````````````````````````````````
#
#name: 0  1 2 3 4 5 6
#   #: 5 12 9 5 8 7 4
#
#______________________________
#Hubs
#- In alphabetical/numerical order
#- Based on empirical quantiles of centralities
#```````````````````````````````````````````````
# Abditibacterium
# Larkinella
# Spirosoma
#
#______________________________
#Centrality measures
#- In decreasing order
#- Centrality of disconnected components is zero
#````````````````````````````````````````````````
#Degree (unnormalized):
#
#Larkinella       10
#Hymenobacter      6
#Aurantisolimonas  6
#Kineococcus       5
#Spirosoma         5
#
#Betweenness centrality (normalized):
#
#Aurantisolimonas 0.51586
#Sphingobium      0.40592
#Chryseobacterium 0.38372
#Stenotrophomonas 0.35941
#Taibaiella       0.34144
#
#Closeness centrality (normalized):
#
#Larkinella       0.64708
#Aurantisolimonas 0.57777
#Devosia          0.56653
#Armatimonadales  0.54835
#Hymenobacter     0.54582
#
#Eigenvector centrality (normalized):
#
#Larkinella      1.00000
#Abditibacterium 0.61074
#Spirosoma       0.60732
#Armatimonadales 0.55740
#Hymenobacter    0.53305



p <- plot(props_spring, 
          labelScale = FALSE,
          nodeColor = "cluster", 
          nodeSize = "eigenvector",
          title1 = "Network on genus level with SPRING associations", 
          showTitle = TRUE,
          cexTitle = 2.3,
          cexLabels = 1.5)

legend(0.7, 1.1, cex = 2.2, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)




net_pears <- netConstruct(ps3genera,  
                          measure = "pearson",
                          normMethod = "clr",
                          zeroMethod = "multRepl",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)

props_pears <- netAnalyze(net_pears, 
                          clustMethod = "cluster_fast_greedy")


plot(props_pears, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     title1 = "Network on Genus level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)



net_pears <- netConstruct(ps3,  
                          measure = "pearson",
                          normMethod = "clr",
                          zeroMethod = "multRepl",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)

props_pears <- netAnalyze(net_pears, 
                          clustMethod = "cluster_fast_greedy")


plot(props_pears, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     title1 = "Network on ASV level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)












##okay try community assembly, neutral/null models with iCAMP


install.packages("iCAMP")
library(iCAMP)

"/Users/gavin/Downloads/iCAMP1-master/Examples/SimpleOTU"


asvtab <- as.data.frame(otu_table(ps3))

tasvtab <- t(asvtab)

asvtab$SpeciesID <- rownames(asvtab)

asvtab <- subset(asvtab, select=c(71,1:70))

comm=tasvtab

tree=read.tree(file = "/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/06_qiime_phylogeny_align-to-tree-mafft-fasttree/export/rooted-tree.qza/tree.tree")
clas=read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/05_qiime_feature-classifier_classify-consensus-vsearch/export/taxonomy_revised_fig_samples_only_5-2025_tab_delim_icamp_2.csv", header = TRUE, sep = "\t", row.names = 1,as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",check.names = FALSE)
treat=read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/fig_microbe_sample_taiwan_2019_metadata_e_5-7-25_4.tsv", header = TRUE, sep = "\t", row.names = 1,
               as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
               check.names = FALSE)

env=read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/fig_microbe_sample_taiwan_2019_metadata_e_5-7-25_4.tsv", header = TRUE, sep = "\t", row.names = 1,
               as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
               check.names = FALSE)


library(iCAMP)
library(ape)

# 4 # match sample IDs in OTU table and treatment information table
sampid.check=match.name(rn.list=list(comm=comm,treat=treat))
# sampid.check=match.name(rn.list=list(comm=comm,treat=treat)) # if you do not have env.file
# for the example data, the output should be "All match very well".
# for your data files, if you have not matched their IDs, the unmatched samples will be removed.
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE] # if some unmatched samples were removed, some OTUs may become ghosts, then you may use this line to remove them if necessary.
#env=sampid.check$env # skip this if you do not have env.file

spid.check=match.name(cn.list=list(comm=comm),rn.list=list(clas=clas),tree.list=list(tree=tree))
# for the example data, the output should be "All match very well".
# for your data files, if you have not matched the IDs before, the unmatched OTUs will be removed.
comm=spid.check$comm
clas=spid.check$clas
tree=spid.check$tree





prefix="Test"  # prefix of the output file names. usually use a project ID.
rand.time=100  # randomization time, 1000 is usually enough. For example test, you may set as 100 or less to save time.
nworker=4 # nworker is thread number for parallel computing, which depends on the CPU core number of your computer.
memory.G=50 # to set the memory size as you need (but should be less than the available space in your hard disk), so that calculation of large tree will not be limited by physical memory. unit is Gb.

save.wd="/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/icamp_out_figs_only/"

# 6 # calculate pairwise phylogenetic distance matrix.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large phylogenetic distance matrix. 
setwd(save.wd)
if(!file.exists("pd.desc")) 
{
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
  # output files:
  # path.rda: a R object to list all the nodes and  edge lengthes from root to every tip. saved in R data format. an intermediate output when claculating phylogenetic distance matrix.
  # pd.bin: BIN file (backingfile) generated by function big.matrix in R package bigmemory. This is the big matrix storing pairwise phylogenetic distance values. By using this bigmemory format file, we will not need memory but hard disk when calling big matrix for calculation.
  # pd.desc: the DESC file (descriptorfile) to hold the backingfile (pd.bin) description.
  # pd.taxon.name.csv: comma delimited csv file storing the IDs of tree tips (OTUs), serving as the row/column names of the big phylogenetic distance matrix.
}else{
  # if you already calculated the phylogenetic distance matrix in a previous run
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
}

env_num_only <- data.frame(latitude=env$latitude,longitude=env$longitude,day_of_aug_2019_picked=env$day_of_aug_2019_picked,day_of_aug_2019_dissected=env$day_of_aug_2019_dissected,foundress_number=env$foundress_number,row.names = rownames(env))

# 7 # assess niche preference difference between species
# env is required for this step.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large niche difference matrix. 
setwd(save.wd)
niche.dif=iCAMP::dniche(env = env_num_only,comm = comm,method = "niche.value",
                        nworker = nworker,out.dist=FALSE,bigmemo=TRUE,
                        nd.wd=save.wd)


# 8 # within-bin phylogenetic signal assessment.
# For real data, you may try several different settings of binning, and choose the one leading to the best within-bin phylogenetic signal.
# env is required for this step.
# 8.1 # try phylogenetic binning using current setttings.
ds = 0.2 # setting can be changed to explore the best choice
bin.size.limit = 5 # setting can be changed to explore the best choice. # here set as 5 just for the small example dataset. For real data, usually try 12 to 48.

# The tree for taxa.binphy.big must be a rooted tree.
if(!ape::is.rooted(tree))
{
  tree.rt=iCAMP::midpoint.root.big(tree = tree, pd.desc = pd.big$pd.file,
                                   pd.spname = pd.big$tip.label,pd.wd = pd.big$pd.wd,
                                   nworker = nworker)
  tree=tree.rt$tree
}
phylobin=taxa.binphy.big(tree = tree, pd.desc = pd.big$pd.file,pd.spname = pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, ds = ds, bin.size.limit = bin.size.limit,
                         nworker = nworker)


# 8.2 # test within-bin phylogenetic signal.
sp.bin=phylobin$sp.bin[,3,drop=FALSE]
sp.ra=colMeans(comm/rowSums(comm))
abcut=3 # you may remove some species, if they are too rare to perform reliable correlation test.
commc=comm[,colSums(comm)>=abcut,drop=FALSE]
dim(commc)
spname.use=colnames(commc)
binps=iCAMP::ps.bin(sp.bin = sp.bin,sp.ra = sp.ra,spname.use = spname.use,
                    pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                    nd.list = niche.dif$nd,nd.spname = niche.dif$names,ndbig.wd = niche.dif$nd.wd,
                    cor.method = "pearson",r.cut = 0.1, p.cut = 0.05, min.spn = 5)
if(file.exists(paste0(prefix,".PhyloSignalSummary.csv"))){appendy=TRUE;col.namesy=FALSE}else{appendy=FALSE;col.namesy=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binps$Index),file = paste0(prefix,".PhyloSignalSummary.csv"),
            append = appendy, quote=FALSE, sep=",", row.names = FALSE,col.names = col.namesy)
if(file.exists(paste0(prefix,".PhyloSignalDetail.csv"))){appendy2=TRUE;col.namesy2=FALSE}else{appendy2=FALSE;col.namesy2=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binID=rownames(binps$detail),binps$detail),file = paste0(prefix,".PhyloSignalDetail.csv"),
            append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)
# since this example small data is randomly generated, the correlation should be very weak.
# usually, you are looking for a binning setting lead to higher RAsig.abj (relative abundance of bins with significant phylogenetic signal) and relative high meanR (mean correlation coefficient across bins).
# see help document of the function "ps.bin" for the meaning of output.





# 9 # iCAMP analysis
# 9.1 # without omitting small bins.
# commonly use # set sig.index as Confidence instead of SES.RC (betaNRI/NTI + RCbray)
bin.size.limit = 5 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
sig.index="Confidence" # see other options in help document of icamp.big.
icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                       prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
# there are quite a few parameters in this function, please check the help document of "icamp.big".
# output files:
# Test.iCAMP.detail.rda: the object "icres" saved in R data format. it is a list object. The first element bNRIiRCa is the result of relative importance of each assembly process in each pairwise comparison (each turnover). The second element "detail" including binning information (named taxabin), phylogenetic and taxonomic metrics results in each bin (named like bNRIi, RCa, etc.), relative abundance of each bin (bin.weight), relative importance of each process in each turnover between communities (processes), input settings (setting), and input community data matrix (comm). See help document of the function icamp.big for more details.

detail.null=TRUE
bin.size.limit = 5 
sig.index="SES.RC" # this is traditional way, with assumption that null values of phylogenetic metrics follow normal distribution. 
prefixb="TestB"

icres2=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                       prefix = prefixb, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = detail.null, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)





# 9.2.2 # normality test
nntest=iCAMP::null.norm(icamp.output=icres2, p.norm.cut=0.05, detail.out=FALSE)
# output shows non-normal distribution ratio in each bin, i.e. the proportion of turnovers which have null values significantly deviated from normal distribution.
# if some ratio values are very high, may need to change to use "Confidence" as sig.index.


# 9.2.3 # change sig.index to "Confidence".
icres3=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "Confidence", detail.save = TRUE, detail.null = FALSE, conf.cut = 0.975)
head(icres3$CbMPDiCBraya)


# 9.2.4 # change sig.index to "RC" for both phylogenetic and taxonomic metrics.
icres4=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "RC", detail.save = TRUE, detail.null = FALSE, rc.cut = 0.95)
head(icres4$RCbMPDiRCbraya)

# 9.2.5 # the function can also change the significance threshold.
icres5=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "SES.RC", detail.save = TRUE, detail.null = FALSE, ses.cut = 1.64, rc.cut = 0.9)
head(icres5$bNRIiRCbraya)



# 9.5 # input community matrix as relative abundances (values < 1) rather than counts
comra=comm/rowSums(comm)
prefixra=paste0(prefix,"RA")
bin.size.limit = 5 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
icres6=iCAMP::icamp.big(comm=comra,tree=tree,pd.desc=pd.big$pd.file, pd.spname=pd.big$tip.label, pd.wd=pd.big$pd.wd,
                        rand=rand.time,prefix=prefixra,ds=0.2,pd.cut=NA,sp.check=TRUE,
                        phylo.rand.scale="within.bin",taxa.rand.scale="across.all",
                        phylo.metric="bMPD",sig.index="Confidence",
                        bin.size.limit=bin.size.limit,nworker=nworker,memory.G=memory.G,
                        rtree.save=FALSE,detail.save=TRUE,qp.save=FALSE,detail.null=FALSE,
                        ignore.zero=TRUE,output.wd=save.wd,correct.special=TRUE,unit.sum=rowSums(comra),
                        special.method="depend",ses.cut = 1.96,rc.cut = 0.95,conf.cut=0.975,
                        omit.option="no",meta.ab=NULL, taxo.metric="bray", transform.method=NULL,
                        logbase=2, dirichlet=TRUE)






#evenness


ps4 <- subset_samples(ps3, surface_or_interior == "interior")

ps_no_small <- prune_samples((!sample_names(ps4) %in% c("GW20","GW34")), ps4)





eve_df <- evenness(ps_no_small)
eve_df$sample_id <- rownames(eve_df)


adiv <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps_no_small, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps_no_small, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps_no_small)))), tree = phyloseq::phy_tree(ps_no_small))[, 1],
  "foundress_number" = phyloseq::sample_data(ps_no_small)$foundress_number)


adiv$sample_id <- rownames(adiv)

eve_adiv <- merge(adiv,eve_df,by="sample_id")


eve_adiv$foundress_number <- as.numeric(eve_adiv$foundress_number)



evefounddf <- data.frame(sample_id=eve_adiv$sample_id,foundress_number=eve_adiv$foundress_number,camargo=eve_adiv$camargo,pielou=eve_adiv$pielou,simpson=eve_adiv$simpson,evar=eve_adiv$evar,bulla=eve_adiv$bulla)


eve_melt <- reshape2::melt(evefounddf, id.vars= c("sample_id","foundress_number"),measure.vars = c("camargo", "pielou","simpson","evar","bulla"))



ggplot(eve_melt, aes(x=foundress_number,y=value)) + geom_point() + facet_wrap(~variable,scales="free") +theme_cowplot() + stat_smooth() + scale_y_continuous(limits = c(0,NA)) + xlab("Foundress number") +ylab("Evenness measure") + theme(strip.background = element_rect(colour="white", fill="white")) + scale_y_continuous(oob=scales::rescale_none)





asymmodcamargo <- nls(camargo ~ SSasymp(foundress_number, Asym, R0, lrc), data = evefounddf)

summary(asymmodcamargo)

#Parameters:
#     Estimate Std. Error t value Pr(>|t|)
#Asym 0.994654   0.005349  185.94   <2e-16 ***
#R0   0.978998   0.006230  157.15   <2e-16 ***
#lrc  0.470230   2.241405    0.21    0.835
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.02067 on 34 degrees of freedom
#
#Number of iterations to convergence: 3
#Achieved convergence tolerance: 1.652e-06

pred_x <- seq(0, 20, length.out = 37)

#pred_y <- predict(asymmodNumber_of_OTU, pred_x)
    #not sure what happened here
    #Asym+(R0-Asym)*exp(-exp(lrc)*input)
pred_y2 <- (0.994654+(0.978998-0.994654)*exp(-exp(-0.470230)*pred_x))
    #this is the function R says SSasymp is trying to fit, I just plugged in the fitted parameters

predicted_data <- data.frame(x = pred_x,y=pred_y2)
observed_data <- data.frame(x=evefounddf$foundress_number, y=evefounddf$camargo)

se_upper <- ((0.994654+0.005349)+((0.978998+0.006230)-(0.994654+0.005349))*exp(-exp((-0.470230))*pred_x))

se_lower <- ((0.994654-0.005349)+((0.978998-0.006230)-(0.994654-0.005349))*exp(-exp((-0.470230))*pred_x))

se_high_df <- data.frame(x = pred_x,y=se_upper)
se_low_df <- data.frame(x = pred_x,y=se_lower)

ggplot(predicted_data, aes(x=x,y=y)) + geom_line(colour="blue") + geom_line(data=se_high_df,colour="blue",linetype="dotted") + geom_line(data=se_low_df,colour="blue",linetype="dotted") + geom_point(data = observed_data) +theme_cowplot() + xlab("Foundress number") + ylab("Camargo evenness") + ylim(0.85,1.01) + ggtitle("a")

a <- ggplot(predicted_data, aes(x=x,y=y)) + geom_line(colour="blue") + geom_line(data=se_high_df,colour="blue",linetype="dotted") + geom_line(data=se_low_df,colour="blue",linetype="dotted") + geom_point(data = observed_data) +theme_cowplot() + xlab("Foundress number") + ylab("Camargo evenness") + ylim(0.85,1.01) + ggtitle("a")





asymmodpielou <- nls(pielou ~ SSasymp(foundress_number, Asym, R0, lrc), data = evefounddf)

#Error in nls(y ~ cbind(1 - exp(-exp(lrc) * x), exp(-exp(lrc) * x)), data = xy,  :
#  singular gradient
#alright

ggplot(evefounddf, aes(x=foundress_number,y=pielou)) + geom_point() +theme_cowplot() + xlab("Foundress number") + ylab("Pielou evenness") + ggtitle("b")

b <- ggplot(evefounddf, aes(x=foundress_number,y=pielou)) + geom_point() +theme_cowplot() + xlab("Foundress number") + ylab("Shannon-Weaver evenness") + ggtitle("b")













asymmodsimpson <- nls(simpson ~ SSasymp(foundress_number, Asym, R0, lrc), data = evefounddf)

summary(asymmodsimpson)

#Formula: simpson ~ SSasymp(foundress_number, Asym, R0, lrc)
#
#Parameters:
#     Estimate Std. Error t value Pr(>|t|)
#Asym  0.13069    0.03368   3.880 0.000456 ***
#R0    0.09543    0.01829   5.219 8.94e-06 ***
#lrc  -1.14903    2.23420  -0.514 0.610374
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.06315 on 34 degrees of freedom
#
#Number of iterations to convergence: 24
#Achieved convergence tolerance: 7.247e-06

pred_x <- seq(0, 20, length.out = 37)

#pred_y <- predict(asymmodNumber_of_OTU, pred_x)
    #not sure what happened here
    #Asym+(R0-Asym)*exp(-exp(lrc)*input)
pred_y2 <- (0.13069+(0.09543-0.13069)*exp(-exp(-1.14903)*pred_x))
    #this is the function R says SSasymp is trying to fit, I just plugged in the fitted parameters

predicted_data <- data.frame(x = pred_x,y=pred_y2)
observed_data <- data.frame(x=evefounddf$foundress_number, y=evefounddf$simpson)

se_upper <- ((0.13069+0.03368)+((0.09543+0.01829)-(0.13069+0.03368))*exp(-exp((-1.14903))*pred_x))

se_lower <- ((0.13069-0.03368)+((0.09543-0.01829)-(0.13069-0.03368))*exp(-exp((-1.14903))*pred_x))

se_high_df <- data.frame(x = pred_x,y=se_upper)
se_low_df <- data.frame(x = pred_x,y=se_lower)

ggplot(predicted_data, aes(x=x,y=y)) + geom_line(colour="blue") + geom_line(data=se_high_df,colour="blue",linetype="dotted") + geom_line(data=se_low_df,colour="blue",linetype="dotted") + geom_point(data = observed_data) +theme_cowplot() + xlab("Foundress number") + ylab("Simpson evenness") + ylim(0,0.3) + ggtitle("c")

c <- ggplot(predicted_data, aes(x=x,y=y)) + geom_line(colour="blue") + geom_line(data=se_high_df,colour="blue",linetype="dotted") + geom_line(data=se_low_df,colour="blue",linetype="dotted") + geom_point(data = observed_data) +theme_cowplot() + xlab("Foundress number") + ylab("Simpson evenness") + ylim(0,0.3) + ggtitle("c")






asymmodevar <- nls(evar ~ SSasymp(foundress_number, Asym, R0, lrc), data = evefounddf)

summary(asymmodevar)

#> summary(asymmodevar)
#
#Formula: evar ~ SSasymp(foundress_number, Asym, R0, lrc)
#
#Parameters:
#     Estimate Std. Error t value Pr(>|t|)
#Asym  0.18513    0.02664   6.949 5.19e-08 ***
#R0    0.27700    0.01680  16.491  < 2e-16 ***
#lrc  -0.86526    0.73763  -1.173    0.249
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.0572 on 34 degrees of freedom
#
#Number of iterations to convergence: 8
#Achieved convergence tolerance: 3.751e-06

pred_x <- seq(0, 20, length.out = 37)

#pred_y <- predict(asymmodNumber_of_OTU, pred_x)
    #not sure what happened here
    #Asym+(R0-Asym)*exp(-exp(lrc)*input)
pred_y2 <- (0.18513+(0.27700-0.18513)*exp(-exp(-0.86526)*pred_x))
    #this is the function R says SSasymp is trying to fit, I just plugged in the fitted parameters

predicted_data <- data.frame(x = pred_x,y=pred_y2)
observed_data <- data.frame(x=evefounddf$foundress_number, y=evefounddf$evar)

se_upper <- ((0.18513+0.02664)+((0.27700+0.01680)-(0.18513+0.02664))*exp(-exp((-0.86526))*pred_x))

se_lower <- ((0.18513-0.02664)+((0.27700-0.01680)-(0.18513-0.02664))*exp(-exp((-0.86526))*pred_x))

se_high_df <- data.frame(x = pred_x,y=se_upper)
se_low_df <- data.frame(x = pred_x,y=se_lower)

ggplot(predicted_data, aes(x=x,y=y)) + geom_line(colour="blue") + geom_line(data=se_high_df,colour="blue",linetype="dotted") + geom_line(data=se_low_df,colour="blue",linetype="dotted") + geom_point(data = observed_data) +theme_cowplot() + xlab("Foundress number") + ylab("Evar index")  + ggtitle("d")

d <- ggplot(predicted_data, aes(x=x,y=y)) + geom_line(colour="blue") + geom_line(data=se_high_df,colour="blue",linetype="dotted") + geom_line(data=se_low_df,colour="blue",linetype="dotted") + geom_point(data = observed_data) +theme_cowplot() + xlab("Foundress number") + ylab("Evar index")  + ggtitle("d")




asymmodbulla <- nls(bulla ~ SSasymp(foundress_number, Asym, R0, lrc), data = evefounddf)

summary(asymmodbulla)
#Error in nls(y ~ cbind(1 - exp(-exp(lrc) * x), exp(-exp(lrc) * x)), data = xy,  :
#  singular gradient

ggplot(evefounddf, aes(x=foundress_number,y=bulla)) + geom_point() +theme_cowplot() + xlab("Foundress number") + ylab("Bulla index") + ggtitle("e")

e <- ggplot(evefounddf, aes(x=foundress_number,y=bulla)) + geom_point() +theme_cowplot() + xlab("Foundress number") + ylab("Bulla index") + ggtitle("e")

library(patchwork)

(a+b+c)/(d+e)


ggsave("Supplemental_Figure_evenness_fig_interiors_8-18-2025.pdf",bg="white",height=7,width=10,units="in",useDingbats=FALSE)



