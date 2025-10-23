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
setwd("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/")

#load ASV counts (see line 319 of workflow.sh), taxonomy table (line 333 of workflow.sh), and sample metadata
ps_object <- read_phyloseq(otu.file = "feature-table.csv", 
                    taxonomy.file = "taxonomy_revised_2.csv", 
                    metadata.file = "inopinata_and_elegans_metadata_3.csv", 
                    type = 'simple')

#load phylogenetic tree (line 326 of workflow.sh)
treefile <- read.tree("tree.tree")

#add tree to phyloseq object
ps.ng.tax <- merge_phyloseq(ps_object, treefile)

#samples i want to keep

samplestoexclude <- scan("samples_to_exclude.txt",sep = "\n",what="character")

#remove elegans samples that ARE NOT substrates



ps.ng.tax1 <- prune_samples(!(sample_names(ps.ng.tax) %in% samplestoexclude), ps.ng.tax)


#remove mitochondrial and chloroplast reads

nochmi <- subset_taxa(ps.ng.tax1, !Family %in% "Mitochondria" & !Order %in% "Chloroplast")

##reads after removing mito and chloro
#orgoreads <- as.data.frame(sample_sums(nochmi))
#orgoreads$sample_id <- rownames(orgoreads)
#orgoreads$num_paired_end_reads <- orgoreads[,1]
#
#orgoreads2 <- data.frame(sample_id = orgoreads$sample_id, num_paired_end_reads=orgoreads$num_paired_end_reads, step= "Remove organelles")
#

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
#[1] 150

#there are 150 ASVs in the controls that will ultimately be removed from the fig microbiome analysis.
#
nochmiprevdf = apply(X = otu_table(nochmi),
               MARGIN = ifelse(taxa_are_rows(controls), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
nochmiprevdf = data.frame(Prevalence = nochmiprevdf,
                    TotalAbundance = taxa_sums(nochmi),
                    tax_table(nochmi))

nrow(nochmiprevdf)
#[1] 58457
150/58457
#[1] 0.002565989

#these are the OTUs to remove
control_OTU$OTU

#remove control samples

no_controls <- prune_samples(!(sample_names(nochmi) %in% c("EXTRNEG","ExtrNeg1","EXTRPOS","GW10","GW11","GW31","GW9","PCRNeg1","PCRNeg2","PCRPos1","PCRPos2")), nochmi)

#remove OTU's found in controls and controls; remove control samples
ps3 <- subset_taxa(no_controls, !OTU %in% control_OTU$OTU)

#save.image(file = "elegans_included_workspace_8-11-2025.RData")
#load("elegans_included_workspace_8-11-2025.RData")

as.data.frame(table(sam_data(ps.ng.tax1)$category))

#                          Var1 Freq
#1            dirksen_substrate   64
#2             samuel_substrate   60
#3 woodruff_extraction_negative    2
#4 woodruff_extraction_positive    1
#5        woodruff_fig_interior   38
#6         woodruff_fig_surface   32
#7         woodruff_M9_negative    4
#8        woodruff_pcr_negative    2
#9        woodruff_pcr_positive    2

#write.table(as.data.frame(otu_table(ps3)),"elegans_no_control_otu_table_8-25-2025.tsv",quote=FALSE,sep="\t")


#write.table(as.data.frame(tax_table(ps3)),"elegans_no_control_tax_table_8-25-2025.tsv",quote=FALSE,sep="\t")

#ASV prevalence


prevdf = apply(X = otu_table(ps3),
               MARGIN = ifelse(taxa_are_rows(ps3), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps3),
                    tax_table(ps3))

#get number of samples

nrow(prevdf)
#[1] 58457
length(unique(prevdf$OTU))
#[1] 58457



nsamples(ps3)
#[1] 194

#get summary stats for prevalence

summary(prevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   1.000   1.000   2.247   1.000 157.000
    
    #No OTU is present in all samples

#prevalence for OTUs histogram

#ggplot(prevdf, aes(x = Prevalence)) + geom_histogram(colour="black",binwidth=1) + theme_cowplot() +geom_vline(xintercept = 2.247, colour = 'red') + ylab("OTU count") 

#ggsave("Supplemental_Figure_3_but_with_elegans_samples_7-30-25.png",bg="white",height=4.5,width=7.5,units="in")


#just the elegans samples


elegans_substrates <- subset_samples(ps3, nematode_species == "C. elegans")



#ASV prevalence


prevdf2 = apply(X = otu_table(elegans_substrates),
               MARGIN = ifelse(taxa_are_rows(elegans_substrates), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf2 = data.frame(Prevalence = prevdf2,
                    TotalAbundance = taxa_sums(elegans_substrates),
                    tax_table(elegans_substrates))

prevdf2nozero <- prevdf2[prevdf2$Prevalence > 0,]

nrow(prevdf2nozero)
#[1] 49467
length(unique(prevdf2nozero$OTU))
#[1] 49467
    #okay! 



nsamples(elegans_substrates)
#[1] 124

summary(prevdf2nozero$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.000   1.000   1.000   2.423   2.000  94.000

prevdf2nozero[prevdf2nozero$Prevalence > 89,]

inopinata_substrates <- subset_samples(ps3, nematode_species == "C. inopinata")

prevdf3 = apply(X = otu_table(inopinata_substrates),
               MARGIN = ifelse(taxa_are_rows(inopinata_substrates), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf3 = data.frame(Prevalence = prevdf3,
                    TotalAbundance = taxa_sums(inopinata_substrates),
                    tax_table(inopinata_substrates))


prevdf3nozero <- prevdf3[prevdf3$Prevalence > 0,]

nrow(prevdf3nozero)
#[1] 3662
length(unique(prevdf3nozero$OTU))
#[1] 3662

49467/3662
#13.50819 -- elegans substrates have 13x as many ASVs as inopinata substrates?

elegans_substrates_genera <- tax_glom(elegans_substrates, "Genus", NArm = TRUE)

elegans_substrates_generaprevdf = apply(X = otu_table(elegans_substrates_genera),
               MARGIN = ifelse(taxa_are_rows(elegans_substrates_genera), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
elegans_substrates_generaprevdf = data.frame(Prevalence = elegans_substrates_generaprevdf,
                    TotalAbundance = taxa_sums(elegans_substrates_genera),
                    tax_table(elegans_substrates_genera))


summary(elegans_substrates_generaprevdf$Prevalence)
#> summary(elegans_substrates_generaprevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.00    2.00    7.00   18.47   28.75  115.00

elegans_substrates_generaprevdf[elegans_substrates_generaprevdf$Prevalence > 111,]



elegans_substrates_family <- tax_glom(elegans_substrates, "Family", NArm = TRUE)

elegans_substrates_familyprevdf = apply(X = otu_table(elegans_substrates_family),
               MARGIN = ifelse(taxa_are_rows(elegans_substrates_family), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
elegans_substrates_familyprevdf = data.frame(Prevalence = elegans_substrates_familyprevdf,
                    TotalAbundance = taxa_sums(elegans_substrates_family),
                    tax_table(elegans_substrates_family))

summary(elegans_substrates_familyprevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.00    3.00   15.00   27.11   47.00  123.00

elegans_substrates_familyprevdf[elegans_substrates_familyprevdf$Prevalence > 99,]$Family



#by paper...

dirksen_subs <- subset_samples(ps3, category == "dirksen_substrate")

samuel_subs <- subset_samples(ps3, category == "samuel_substrate")

dirksen_subsprevdf = apply(X = otu_table(dirksen_subs),
               MARGIN = ifelse(taxa_are_rows(dirksen_subs), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
dirksen_subsprevdf = data.frame(Prevalence = dirksen_subsprevdf,
                    TotalAbundance = taxa_sums(dirksen_subs),
                    tax_table(dirksen_subs))




samuel_subsprevdf = apply(X = otu_table(samuel_subs),
               MARGIN = ifelse(taxa_are_rows(samuel_subs), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
samuel_subsprevdf = data.frame(Prevalence = samuel_subsprevdf,
                    TotalAbundance = taxa_sums(samuel_subs),
                    tax_table(samuel_subs))

#dirksen ASVs



dirksen_subsprevdfnozero <- dirksen_subsprevdf[dirksen_subsprevdf$Prevalence > 0,]

nrow(dirksen_subsprevdfnozero)
#[1] 47365
length(unique(dirksen_subsprevdfnozero$OTU))
#[1] 47365

summary(dirksen_subsprevdfnozero$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.000   1.000   1.000   2.277   2.000  59.000

nsamples(dirksen_subs)
#[1] 64

samuel_subsprevdfnozero <- samuel_subsprevdf[samuel_subsprevdf$Prevalence > 0,]

nrow(samuel_subsprevdfnozero)
#[1] 4746
length(unique(samuel_subsprevdfnozero$OTU))
#[1] 4746


summary(samuel_subsprevdfnozero$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.000   1.000   1.000   2.524   2.000  56.000

nsamples(samuel_subs)
#[1] 60

#just plot one ordination first

#no chloro/mito

pslog <- transform_sample_counts(nochmi, function(x) log(1 + x))

#PCoA using Bray-Curtis dissimilarity, including controls
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")

evals <- out.pcoa.log$values[,1]
ord_df <- plot_ordination(pslog, out.pcoa.log, color = "category", justDF = TRUE)

evals[1]/sum(evals)
#[1] 0.1395958
evals[2]/sum(evals)
#[1] 0.08870025
metdattt <- read.csv("inopinata_and_elegans_metadata_2.csv",header=TRUE)
ord_df$sample.id <- rownames(ord_df)
ord_dfm <- merge(ord_df,metdattt)

nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

ggplot(ord_dfm, aes(x=Axis.1,y=Axis.2)) + geom_point(aes(color = category),size=1) + scale_colour_manual(values = mycolors) + coord_fixed(sqrt(evals[2] / evals[1])) + theme_cowplot() + xlab("Axis 1 (14.0%)") + ylab("Axis 2 (8.9%)")
ggsave("PCoA_bray_curtis_with_elegans_substrates_7-30-25.png",bg="white",height=4.5,width=7.5,units="in")



    #original  figure 3b

#+ scale_color_manual(labels = c("Control", "Fig Suspension", "Fig Surface Wash"), values = c("#fde725", "#21918c","#440154")) +labs(colour="Sample Type") 

#PCoA with weighted Unifrac, including controls

#out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
#evals <- out.wuf.log$values$Eigenvalues

#evals[1]/sum(evals)
#evals[2]/sum(evals)

#plot_ordination(pslog, out.wuf.log, color = category) + coord_fixed(sqrt(evals[2] / evals[1])) + scale_color_manual(values = mycolors) + theme_cowplot() + xlab("Axis 1 (14.0%)") + ylab("Axis 2 (8.9%)") + ggtitle("PCoA, weighted Unifrac")
    #this is supplemental figure 7
#ggsave("PCoA_weighted_unifrac_with_elegans_substrates_7-30-25.png",bg="white",height=4.5,width=7.5,units="in")



as.data.frame(table(sam_data(ps3)$category))
#                   Var1 Freq
#1     dirksen_substrate   64
#2      samuel_substrate   60
#3 woodruff_fig_interior   38
#4  woodruff_fig_surface   32

sam_data(ps3)$SampleID

###Aitchison distances (https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/)

#okay, want to add some metadata here--

sample_data(ps3)

metadata2d <- read.csv("inopinata_and_elegans_metadata_2d.csv",header=TRUE,row.names=1)

as.data.frame(table(metadata2d$substrate_type))
#              Var1 Freq
#1          Compost   62
#2 fig surface wash   32
#3   fig suspension   38
#4            Fruit   52
#5             Stem    6
#6           Vector    4


dirk <- metadata2d[metadata2d$study == "Dirksen",]
samu <- metadata2d[metadata2d$study == "Samuel",]


as.data.frame(table(dirk$substrate_type))



as.data.frame(table(samu$substrate_type))

sample_data(ps3)$substrate_type <- metadata2d$substrate_type
sample_data(ps3)$year <- metadata2d$year
sample_data(ps3)$year2 <- as.factor(metadata2d$year)
sample_data(ps3)$study <- metadata2d$study
ps_clr <- microbiome::transform(ps3, "clr")

#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + geom_bar(stat="identity", fill = "blue") + labs(x = "\nAxis", y = "Proportion of Variance\n")

#Examine eigenvalues and % prop. variance explained
head(ord_clr$CA$eig)                                                  

#      PC1       PC2       PC3       PC4       PC5       PC6
#1516.7109  784.8807  503.9165  342.1571  333.2347  265.3946


sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))     

#       PC1        PC2        PC3        PC4        PC5
#0.12346202 0.06389020 0.04101938 0.02785199 0.02712569
#

clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ps3, ord_clr, color="category") + 
  geom_point(size = 1) +
  coord_fixed(clr2 / clr1)

phyloseq::plot_ordination(ps3, ord_clr, color="substrate_type") + 
  geom_point(size = 1) +
  coord_fixed(clr2 / clr1)

phyloseq::plot_ordination(ps3, ord_clr, color="year2") + 
  geom_point(size = 1) +
  coord_fixed(clr2 / clr1)


ord_df <-plot_ordination(ps3, ord_clr, color = "study", justDF = TRUE)


summary(ord_df$PC1)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-7.9198  0.1717  1.1871  0.0000  1.5939  2.3046


summary(ord_df$PC2)

#summary(ord_df$PC2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-5.9248 -2.0216 -0.8189  0.0000  1.7753  7.5160


a <- ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=study),size=1) + coord_fixed(clr2 / clr1) + theme_cowplot()  + xlab("Axis 1 (12.3%)") + ylab("Axis 2 (6.4%)") + scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a")) + scale_x_continuous(limits=c(-8,3),breaks=seq(-8,3,by=1)) + scale_y_continuous(limits=c(-6,8),breaks=seq(-6,8,by=1))



ord_df2 <-plot_ordination(ps3, ord_clr, color = "substrate_type", justDF = TRUE)


b <- ggplot(ord_df2, aes(x=PC1,y=PC2)) + geom_point(aes(colour=substrate_type),size=1) + coord_fixed(clr2 / clr1) + theme_cowplot() + scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))  + xlab("Axis 1 (12.3%)") + ylab("Axis 2 (6.4%)") + scale_x_continuous(limits=c(-8,3),breaks=seq(-8,3,by=1)) + scale_y_continuous(limits=c(-6,8),breaks=seq(-6,8,by=1))


ord_df3 <-plot_ordination(ps3, ord_clr, color = "year2", justDF = TRUE)

c <- ggplot(ord_df3, aes(x=PC1,y=PC2)) + geom_point(aes(colour=year2),size=1) + coord_fixed(clr2 / clr1) + theme_cowplot() + scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"))  + xlab("Axis 1 (12.3%)") + ylab("Axis 2 (6.4%)") + scale_x_continuous(limits=c(-8,3),breaks=seq(-8,3,by=1)) + scale_y_continuous(limits=c(-6,8),breaks=seq(-6,8,by=1))

d <- ggplot(ord_df2, aes(x=PC1,y=PC2)) + geom_point(aes(colour=substrate_type),size=1) + theme_cowplot() + scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))  + scale_x_continuous(limits=c(-1,3),breaks=seq(-1,3,by=1)) + scale_y_continuous(limits=c(-6,8),breaks=seq(-6,8,by=1))

((a+b)/(c+d))


ggsave("figure_3_7-30-25.pdf",bg="white",height=10,width=15,units="in",useDingbats=FALSE)
    #to be figure 3

############
############
############
#here, on 7/30/25 jumped to ~line 781 for figure 2
############
############
############

summary(ord_df$PC1)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-1.7545 -1.0493 -0.9154  0.0000 -0.5047  8.6618

ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=category),size=0.8) + theme_cowplot()  + xlab("Axis 1 (10.8%)") + ylab("Axis 2 (5.8%)") +xlim(-2,3)  

ggsave("ordination_aitchinson_2.png",bg="white")


ggplot(ord_df, aes(x=PC1,y=PC2)) + geom_point(aes(colour=category),size=0.8) + theme_cowplot()  + xlab("Axis 1 (10.8%)") + ylab("Axis 2 (5.8%)") +xlim(-2,1)  


ggsave("ordination_aitchinson_3.png",bg="white")



+ coord_fixed(clr2 / clr1)



+ scale_size_continuous(limits=c(0,20),breaks=c(0,5,10,15,20)) + scale_color_continuous(limits=c(0,20),breaks=c(0,5,10,15,20)) + guides(color=guide_legend(title="Foundress\nnumber"), size = guide_legend(title="Foundress\nnumber"), shape=guide_legend(title="Worms\npresent?")) + theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),legend.title = element_text(size = 14))


###Aitchison distances (https://microbiome.github.io/course_2021_radboud/beta-diversity.html)

# List of packages that we need from cran and bioc 
cran_pkg <- c("BiocManager", "bookdown", "dplyr", "ecodist", "ggplot2", 
              "gridExtra", "kableExtra", "knitr", "scales", "vegan")
bioc_pkg <- c("ANCOMBC", "ape", "DESeq2",  "DirichletMultinomial", "mia", "miaViz")

# Gets those packages that are already installed
#cran_pkg_already_installed <- cran_pkg[ cran_pkg %in% installed.packages() ]
#bioc_pkg_already_installed <- bioc_pkg[ bioc_pkg %in% installed.packages() ]
#
## Gets those packages that need to be installed
#cran_pkg_to_be_installed <- setdiff(cran_pkg, cran_pkg_already_installed)
#bioc_pkg_to_be_installed <- setdiff(bioc_pkg, bioc_pkg_already_installed)
## If there are packages that need to be installed, installs them from CRAN
#if( length(cran_pkg_to_be_installed) ) {
#   install.packages(cran_pkg_to_be_installed)
#}
## If there are packages that need to be installed, installs them from Bioconductor
#if( length(bioc_pkg_to_be_installed) ) {
#   BiocManager::install(bioc_pkg_to_be_installed, ask = F)
#}

# Reorders bioc packages, so that mia and miaViz are first
bioc_pkg <- c(bioc_pkg[ bioc_pkg %in% c("mia", "miaViz") ], 
              bioc_pkg[ !bioc_pkg %in% c("mia", "miaViz") ] ) 

# Loading all packages into session. Returns true if package was successfully loaded.
loaded <- sapply(c(bioc_pkg, cran_pkg), require, character.only = TRUE)
as.data.frame(loaded)


cfpps3 <- convertFromPhyloseq(ps3)

# Does clr transformation. Pseudocount is added, because data contains zeros. 
tse <- transformCounts(ps3, method = "clr", pseudocount = 1)

# Gets clr table
clr_assay <- assays(tse)$clr

# Transposes it to get taxa to columns
clr_assay <- t(clr_assay)

# Calculates Euclidean distances between samples. Because taxa is in columns,
# it is used to compare different samples.
euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")

# Does principal coordinate analysis
euclidean_pcoa <- ecodist::pco(euclidean_dist)

# Creates a data frame from principal coordinates
euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1], 
                                pcoa2 = euclidean_pcoa$vectors[,2])









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
#[1] 310

#get summary stats for prevalence

summary(ps3generaprevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.00    2.25   10.00   25.77   37.00  285.00

#285/310
#[1] 0.9193548


#there is a genus with 91% prevalence, what is it?

ps3generaprevdf[ps3generaprevdf$Prevalence == 285,]$Genus
#[1] "Sphingomonas"
    #only one

#get number of genera
length(unique(ps3generaprevdf$Genus))
#[1] 1473

#306 after removing dubious genera

#plot Supplemental Figure 4
ggplot(ps3generaprevdf, aes(x = Prevalence)) + geom_histogram(colour="black",binwidth=1) + theme_cowplot() +geom_vline(xintercept = 8.803, colour = 'red')  + ylab("Genus count") 


ggsave("Supplemental_Figure_4_with_elegans.png",bg="white",height=4.5,width=7.5,units="in")

#Is Wolbachia in here?
ps3generaprevdf[ps3generaprevdf$Genus == "Wolbachia",]$Prevalence
#[1] 16



#number of OTU (again)
prevdf = apply(X = otu_table(ps3),
               MARGIN = ifelse(taxa_are_rows(ps3), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps3),
                    tax_table(ps3))
nrow(prevdf)
#[1] 58307
#ASVs across elegans and inopinata samples


#get summary stats for prevalence
summary(prevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.000   1.000   1.000   2.597   2.000 219.000

#what phyla are OTUs associated with?
phy_otu_count_df <- as.data.frame(table(prevdf$Phylum))

phy_otu_count_df$Fra <- phy_otu_count_df$Freq/sum(phy_otu_count_df$Freq)

phy_otu_count_df[order(phy_otu_count_df$Freq, decreasing = TRUE),]

#                            Var1  Freq          Fra
#51                Proteobacteria 16770 0.2876155522
#7                   Bacteroidota  7957 0.1364673195
#50               Planctomycetota  7166 0.1229011954
#60             Verrucomicrobiota  3661 0.0627883445
#4               Actinobacteriota  2832 0.0485704975
#3                Acidobacteriota  2810 0.0481931844
#28                    Firmicutes  2762 0.0473699556
#1                                 2467 0.0423105288
#12                   Chloroflexi  2456 0.0421218722
#43                   Myxococcota  1940 0.0332721629
#9               Bdellovibrionota  1700 0.0291560190
#49               Patescibacteria  1242 0.0213010445
#31               Gemmatimonadota   978 0.0167732862
#21                  Dependentiae   500 0.0085752997
#17                 Cyanobacteria   421 0.0072204024
#22              Desulfobacterota   398 0.0068259386
#6                 Armatimonadota   371 0.0063628724
#55                   Sumerlaeota   224 0.0038417343
#53  SAR324_clade(Marine_group_B)   167 0.0028641501
#37              Latescibacterota   152 0.0026068911
#23               Elusimicrobiota   143 0.0024525357
#2               Abditibacteriota   142 0.0024353851
#36               Hydrogenedentes   100 0.0017150599
#54                 Spirochaetota    95 0.0016293069
#27                Fibrobacterota    90 0.0015435539
#20                  Deinococcota    83 0.0014234998
#45                         NB1-j    82 0.0014063492
#44                 Nanoarchaeota    52 0.0008918312
#61                         WPS-2    52 0.0008918312
#8                  Basidiomycota    50 0.0008575300
#48                  Nitrospirota    49 0.0008403794
#11              Campilobacterota    48 0.0008232288
#63                           WS2    42 0.0007203252
#16                 Crenarchaeota    39 0.0006688734
#41             Methylomirabilota    37 0.0006345722
#40                        MBNT15    31 0.0005316686
#26                       FCPU426    27 0.0004630662
#29                Fusobacteriota    24 0.0004116144
#34                 Halobacterota    17 0.0002915602
#65                  Zixibacteria    16 0.0002744096
#24             Entotheonellaeota    15 0.0002572590
#25                 Euryarchaeota    10 0.0001715060
#38              Margulisbacteria    10 0.0001715060
#18                  Dadabacteria     9 0.0001543554
#52                       RCP2-54     9 0.0001543554
#46                     Nematozoa     8 0.0001372048
#56                       Sva0485     8 0.0001372048
#19              Deferribacterota     6 0.0001029036
#33              Halanaerobiaeota     6 0.0001029036
#15                Cloacimonadota     4 0.0000686024
#57                  Synergistota     4 0.0000686024
#58              Thermoplasmatota     3 0.0000514518
#59                  Thermotogota     3 0.0000514518
#5                    Apicomplexa     2 0.0000343012
#13                   Chlorophyta     2 0.0000343012
#14                    Ciliophora     2 0.0000343012
#30                         GAL15     2 0.0000343012
#32                          GN01     2 0.0000343012
#64                           WS4     2 0.0000343012
#10                 Caldisericota     1 0.0000171506
#35                 Heterolobosea     1 0.0000171506
#39 Marinimicrobia_(SAR406_clade)     1 0.0000171506
#42                      Mollusca     1 0.0000171506
#47                  Nitrospinota     1 0.0000171506
#62                           WS1     1 0.0000171506
#66                 Zoopagomycota     1 0.0000171506




#genus supplemental figure


allgen <- tax_glom(ps3, "Genus", NArm = TRUE)
prevdf = apply(X = otu_table(allgen),
               MARGIN = ifelse(taxa_are_rows(allgen), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(allgen),
                    tax_table(allgen))



inopinata_substrates <- subset_samples(ps3, nematode_species == "C. inopinata")

elegans_substrates <- subset_samples(ps3, nematode_species == "C. elegans")




inopinatagen <- tax_glom(inopinata_substrates, "Genus", NArm = TRUE)
inopinataprevdf = apply(X = otu_table(inopinatagen),
               MARGIN = ifelse(taxa_are_rows(inopinatagen), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
inopinataprevdf = data.frame(Prevalence = inopinataprevdf,
                    TotalAbundance = taxa_sums(inopinatagen),
                    tax_table(inopinatagen))


elegansgen <- tax_glom(elegans_substrates, "Genus", NArm = TRUE)
elegansprevdf = apply(X = otu_table(elegansgen),
               MARGIN = ifelse(taxa_are_rows(elegansgen), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
elegansprevdf = data.frame(Prevalence = elegansprevdf,
                    TotalAbundance = taxa_sums(elegansgen),
                    tax_table(elegansgen))



interior_samples <- subset_samples(ps3, category == "woodruff_fig_interior")

surface_samples <- subset_samples(ps3, category == "woodruff_fig_surface")

figintgen <- tax_glom(interior_samples, "Genus", NArm = TRUE)
figintprevdf = apply(X = otu_table(figintgen),
               MARGIN = ifelse(taxa_are_rows(figintgen), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
figintprevdf = data.frame(Prevalence = figintprevdf,
                    TotalAbundance = taxa_sums(figintgen),
                    tax_table(figintgen))


figextgen <- tax_glom(surface_samples, "Genus", NArm = TRUE)
figextprevdf = apply(X = otu_table(figextgen),
               MARGIN = ifelse(taxa_are_rows(figextgen), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
figextprevdf = data.frame(Prevalence = figextprevdf,
                    TotalAbundance = taxa_sums(figextgen),
                    tax_table(figextgen))



summary(prevdf$Prevalence)
#summary(prevdf$Prevalence)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.00    2.25   10.00   25.77   37.00  285.00

#80%?
0.8*310
#[1] 248

eighty <- prevdf[prevdf$Prevalence > 247,]
nrow(eighty)
#[1] 4
eighty
#                                 Prevalence TotalAbundance   Domain
#74b7dc5aa537fcaf4d0046dc7715b9bf        253         971116 Bacteria
#399a96716238bd9fb94415d1f390461b        259         538186 Bacteria
#a538e0c72f940bfe5957d13756d7b83b        285        1021915 Bacteria
#925e56c65cae6ca645de5b98f0dd5028        260         985516 Bacteria
#                                         Phylum               Class
#74b7dc5aa537fcaf4d0046dc7715b9bf Proteobacteria Gammaproteobacteria
#399a96716238bd9fb94415d1f390461b   Bacteroidota         Bacteroidia
#a538e0c72f940bfe5957d13756d7b83b Proteobacteria Alphaproteobacteria
#925e56c65cae6ca645de5b98f0dd5028 Proteobacteria Gammaproteobacteria
#                                            Order            Family
#74b7dc5aa537fcaf4d0046dc7715b9bf  Xanthomonadales  Xanthomonadaceae
#399a96716238bd9fb94415d1f390461b Flavobacteriales     Weeksellaceae
#a538e0c72f940bfe5957d13756d7b83b Sphingomonadales Sphingomonadaceae
#925e56c65cae6ca645de5b98f0dd5028  Pseudomonadales  Pseudomonadaceae
#                                            Genus Species  OTU
#74b7dc5aa537fcaf4d0046dc7715b9bf Stenotrophomonas    <NA> <NA>
#399a96716238bd9fb94415d1f390461b Chryseobacterium    <NA> <NA>
#a538e0c72f940bfe5957d13756d7b83b     Sphingomonas    <NA> <NA>
#925e56c65cae6ca645de5b98f0dd5028      Pseudomonas    <NA> <NA>



#70%?
0.7*310
#[1] 217

seventy <- prevdf[prevdf$Prevalence > 216,]
nrow(seventy)
#[1] 8
seventy
#                                 Prevalence TotalAbundance   Domain
#74b7dc5aa537fcaf4d0046dc7715b9bf        253         971116 Bacteria
#f1f2ca4dce7d38d3b9803546b3632575        218         362162 Bacteria
#4461d2a06acfa9eb83d74e43087f3725        220         429591 Bacteria
#7769cac4c6bfbe7e3cddd33636c22519        224         218388 Bacteria
#399a96716238bd9fb94415d1f390461b        259         538186 Bacteria
#a538e0c72f940bfe5957d13756d7b83b        285        1021915 Bacteria
#f5a8edcca6d3ed5107cb077e7b84c16b        245         484645 Bacteria
#925e56c65cae6ca645de5b98f0dd5028        260         985516 Bacteria
#                                         Phylum               Class
#74b7dc5aa537fcaf4d0046dc7715b9bf Proteobacteria Gammaproteobacteria
#f1f2ca4dce7d38d3b9803546b3632575   Bacteroidota         Bacteroidia
#4461d2a06acfa9eb83d74e43087f3725   Bacteroidota         Bacteroidia
#7769cac4c6bfbe7e3cddd33636c22519 Proteobacteria Alphaproteobacteria
#399a96716238bd9fb94415d1f390461b   Bacteroidota         Bacteroidia
#a538e0c72f940bfe5957d13756d7b83b Proteobacteria Alphaproteobacteria
#f5a8edcca6d3ed5107cb077e7b84c16b Proteobacteria Alphaproteobacteria
#925e56c65cae6ca645de5b98f0dd5028 Proteobacteria Gammaproteobacteria
#                                              Order              Family
#74b7dc5aa537fcaf4d0046dc7715b9bf    Xanthomonadales    Xanthomonadaceae
#f1f2ca4dce7d38d3b9803546b3632575 Sphingobacteriales Sphingobacteriaceae
#4461d2a06acfa9eb83d74e43087f3725   Flavobacteriales   Flavobacteriaceae
#7769cac4c6bfbe7e3cddd33636c22519        Rhizobiales    Beijerinckiaceae
#399a96716238bd9fb94415d1f390461b   Flavobacteriales       Weeksellaceae
#a538e0c72f940bfe5957d13756d7b83b   Sphingomonadales   Sphingomonadaceae
#f5a8edcca6d3ed5107cb077e7b84c16b        Rhizobiales        Rhizobiaceae
#925e56c65cae6ca645de5b98f0dd5028    Pseudomonadales    Pseudomonadaceae
#                                                                              Genus
#74b7dc5aa537fcaf4d0046dc7715b9bf                                   Stenotrophomonas
#f1f2ca4dce7d38d3b9803546b3632575                                   Sphingobacterium
#4461d2a06acfa9eb83d74e43087f3725                                     Flavobacterium
#7769cac4c6bfbe7e3cddd33636c22519                     Methylobacterium-Methylorubrum
#399a96716238bd9fb94415d1f390461b                                   Chryseobacterium
#a538e0c72f940bfe5957d13756d7b83b                                       Sphingomonas
#f5a8edcca6d3ed5107cb077e7b84c16b Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium
#925e56c65cae6ca645de5b98f0dd5028                                        Pseudomonas
#                                 Species  OTU
#74b7dc5aa537fcaf4d0046dc7715b9bf    <NA> <NA>
#f1f2ca4dce7d38d3b9803546b3632575    <NA> <NA>
#4461d2a06acfa9eb83d74e43087f3725    <NA> <NA>
#7769cac4c6bfbe7e3cddd33636c22519    <NA> <NA>
#399a96716238bd9fb94415d1f390461b    <NA> <NA>
#a538e0c72f940bfe5957d13756d7b83b    <NA> <NA>
#f5a8edcca6d3ed5107cb077e7b84c16b    <NA> <NA>
#925e56c65cae6ca645de5b98f0dd5028    <NA> <NA>


metdat2 <- read.csv("inopinata_and_elegans_metadata_2.csv",header=TRUE)

as.data.frame(table(metdat2$nematode_species))
#          Var1 Freq
#1   C. elegans  241
#2 C. inopinata   70
#3      neither   11


as.data.frame(table(metdat2$category))
#                          Var1 Freq
#1            dirksen_substrate  181
#2             samuel_substrate   60
#3 woodruff_extraction_negative    2
#4 woodruff_extraction_positive    1
#5        woodruff_fig_interior   38
#6         woodruff_fig_surface   32
#7         woodruff_M9_negative    4
#8        woodruff_pcr_negative    2
#9        woodruff_pcr_positive    2


seventy_all <- prevdf[prevdf$Prevalence > 216,]

0.7*241
#[1] 168.7
seventy_elegans <- elegansprevdf[elegansprevdf$Prevalence > 168,]
seventy_fig_all <- inopinataprevdf[inopinataprevdf$Prevalence > 48,]
seventy_fig_interior <- figintprevdf[figintprevdf$Prevalence > 26,]
seventy_fig_exterior <- figextprevdf[figextprevdf$Prevalence > 22,]

seventy_elegans$Genus
#[1] "Rhodococcus"
#[2] "Stenotrophomonas"
#[3] "Sphingobacterium"
#[4] "Flavobacterium"
#[5] "Chryseobacterium"
#[6] "Sphingomonas"
#[7] "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
#[8] "Pseudomonas"

seventy_fig_all$Genus
# [1] "Klenkia"
# [2] "Kineococcus"
# [3] "Quadrisphaera"
# [4] "Nocardioides"
# [5] "Massilia"
# [6] "Xanthomonas"
# [7] "Spirosoma"
# [8] "Deinococcus"
# [9] "Mucilaginibacter"
#[10] "Aureimonas"
#[11] "Methylobacterium-Methylorubrum"
#[12] "Hymenobacter"
#[13] "Chryseobacterium"
#[14] "Roseomonas"
#[15] "Sphingomonas"
#[16] "Novosphingobium"
#[17] "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
#[18] "Bdellovibrio"

seventy_fig_interior$Genus
# [1] "Klenkia"
# [2] "Gordonia"
# [3] "Kineococcus"
# [4] "Quadrisphaera"
# [5] "Nocardioides"
# [6] "Massilia"
# [7] "Xanthomonas"
# [8] "Spirosoma"
# [9] "Deinococcus"
#[10] "Mucilaginibacter"
#[11] "Aureimonas"
#[12] "Ochrobactrum"
#[13] "Methylobacterium-Methylorubrum"
#[14] "Hymenobacter"
#[15] "Chryseobacterium"
#[16] "Paracoccus"
#[17] "Roseomonas"
#[18] "Sphingomonas"
#[19] "Novosphingobium"
#[20] "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
#[21] "Azorhizobium"
#[22] "Bdellovibrio"

seventy_fig_exterior$Genus
# [1] "Klenkia"
# [2] "Kineococcus"
# [3] "Quadrisphaera"
# [4] "Nocardioides"
# [5] "Massilia"
# [6] "Spirosoma"
# [7] "Aureimonas"
# [8] "Methylobacterium-Methylorubrum"
# [9] "Hymenobacter"
#[10] "Sphingomonas"
#[11] "Novosphingobium"
#[12] "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
#[13] "Bdellovibrio"

seventy_elegans_genus <- seventy_elegans
seventy_fig_all_genus <- seventy_fig_all
seventy_fig_interior_genus <- seventy_fig_interior
seventy_fig_exterior_genus <- seventy_fig_exterior



#family supplemental figure
#starting again, 7-30-25

allfam <- tax_glom(ps3, "Family", NArm = TRUE)
prevdffam = apply(X = otu_table(allfam),
               MARGIN = ifelse(taxa_are_rows(allfam), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
prevdffam = data.frame(Prevalence = prevdffam,
                    TotalAbundance = taxa_sums(allfam),
                    tax_table(allfam))

inopinatafam <- tax_glom(inopinata_substrates, "Family", NArm = TRUE)
inopinatafamdf = apply(X = otu_table(inopinatafam),
               MARGIN = ifelse(taxa_are_rows(inopinatafam), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
inopinatafamdf = data.frame(Prevalence = inopinatafamdf,
                    TotalAbundance = taxa_sums(inopinatafam),
                    tax_table(inopinatafam))


elegansfam <- tax_glom(elegans_substrates, "Family", NArm = TRUE)
elegansfamdf = apply(X = otu_table(elegansfam),
               MARGIN = ifelse(taxa_are_rows(elegansfam), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
elegansfamdf = data.frame(Prevalence = elegansfamdf,
                    TotalAbundance = taxa_sums(elegansfam),
                    tax_table(elegansfam))

figintfam <- tax_glom(interior_samples, "Family", NArm = TRUE)
figintfamdf = apply(X = otu_table(figintfam),
               MARGIN = ifelse(taxa_are_rows(figintfam), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
figintfamdf = data.frame(Prevalence = figintfamdf,
                    TotalAbundance = taxa_sums(figintfam),
                    tax_table(figintfam))


figextfam <- tax_glom(surface_samples, "Family", NArm = TRUE)
figextfamdf = apply(X = otu_table(figextfam),
               MARGIN = ifelse(taxa_are_rows(figextfam), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
figextfamdf = data.frame(Prevalence = figextfamdf,
                    TotalAbundance = taxa_sums(figextfam),
                    tax_table(figextfam))

#7-30-25: skipping to 80% prevalence (line 1190)

seventy_all_fam <- prevdffam[prevdffam$Prevalence > 216,]
seventy_elegans_fam <- elegansfamdf[elegansfamdf$Prevalence > 168,]
seventy_fig_all_fam <- inopinatafamdf[inopinatafamdf$Prevalence > 48,]
seventy_fig_interior_fam <- figintfamdf[figintfamdf$Prevalence > 26,]
seventy_fig_exterior_fam <- figextfamdf[figextfamdf$Prevalence > 22,]

cat(x, sep = "\n")

cat(seventy_all_fam$Family, sep = "\n")
#Nocardiaceae
#Microbacteriaceae
#Comamonadaceae
#Oxalobacteraceae
#Xanthomonadaceae
#Sphingobacteriaceae
#Flavobacteriaceae
#Acetobacteraceae
#Rhodobacteraceae
#Beijerinckiaceae
#Spirosomaceae
#Weeksellaceae
#Sphingomonadaceae
#Rhizobiaceae
#Pseudomonadaceae
#Erwiniaceae
#Enterobacteriaceae
cat(seventy_elegans_fam$Family, sep = "\n")
#Nocardiaceae
#Microbacteriaceae
#Alcaligenaceae
#Comamonadaceae
#Oxalobacteraceae
#Xanthomonadaceae
#Sphingobacteriaceae
#Flavobacteriaceae
#Acetobacteraceae
#Rhodobacteraceae
#Rhizobiaceae
#Caulobacteraceae
#Beijerinckiaceae
#Spirosomaceae
#Weeksellaceae
#Sphingomonadaceae
#Pseudomonadaceae
#Erwiniaceae
#Enterobacteriaceae
cat(seventy_fig_all_fam$Family, sep = "\n")
#Geodermatophilaceae
#Nocardiaceae
#Pseudonocardiaceae
#Microbacteriaceae
#Kineosporiaceae
#Nocardioidaceae
#Comamonadaceae
#Oxalobacteraceae
#Xanthomonadaceae
#Spirosomaceae
#Deinococcaceae
#Sphingobacteriaceae
#Chitinophagaceae
#Beijerinckiaceae
#Hymenobacteraceae
#Weeksellaceae
#Rhodobacteraceae
#Acetobacteraceae
#Sphingomonadaceae
#Rhizobiaceae
#Enterobacteriaceae
#Erwiniaceae
#Bdellovibrionaceae
cat(seventy_fig_interior_fam$Family, sep = "\n")
#Geodermatophilaceae
#Nocardiaceae
#Pseudonocardiaceae
#Microbacteriaceae
#Kineosporiaceae
#Nocardioidaceae
#Comamonadaceae
#Oxalobacteraceae
#Xanthomonadaceae
#Spirosomaceae
#Deinococcaceae
#Chitinophagaceae
#Sphingobacteriaceae
#Beijerinckiaceae
#Hymenobacteraceae
#Weeksellaceae
#Rhodobacteraceae
#Acetobacteraceae
#Sphingomonadaceae
#Rhizobiaceae
#Xanthobacteraceae
#Enterobacteriaceae
#Erwiniaceae
#Bdellovibrionaceae
cat(seventy_fig_exterior_fam$Family, sep = "\n")
#Geodermatophilaceae
#Pseudonocardiaceae
#Microbacteriaceae
#Kineosporiaceae
#Nocardioidaceae
#Comamonadaceae
#Oxalobacteraceae
#Xanthomonadaceae
#Spirosomaceae
#Beijerinckiaceae
#Hymenobacteraceae
#Weeksellaceae
#Acetobacteraceae
#Sphingomonadaceae
#Rhizobiaceae
#Enterobacteriaceae
#Erwiniaceae
#Bdellovibrionaceae

library(UpSetR)

    #fig families in variables all_eighty , int_eighty , and surf_eighty above
    #C elegans families from Figure 3C of Zhang et al. 2017 https://doi.org/10.3389/fmicb.2017.00485
        #Actinomycetales is not included because it is an order and not a family

upset_list <- list(all_samples=c("Nocardiaceae","Microbacteriaceae","Comamonadaceae","Oxalobacteraceae","Xanthomonadaceae","Sphingobacteriaceae","Flavobacteriaceae","Acetobacteraceae","Rhodobacteraceae","Beijerinckiaceae","Spirosomaceae","Weeksellaceae","Sphingomonadaceae","Rhizobiaceae","Pseudomonadaceae","Erwiniaceae","Enterobacteriaceae"),elegans_samples=c("Nocardiaceae","Microbacteriaceae","Alcaligenaceae","Comamonadaceae","Oxalobacteraceae","Xanthomonadaceae","Sphingobacteriaceae","Flavobacteriaceae","Acetobacteraceae","Rhodobacteraceae","Rhizobiaceae","Caulobacteraceae","Beijerinckiaceae","Spirosomaceae","Weeksellaceae","Sphingomonadaceae","Pseudomonadaceae","Erwiniaceae","Enterobacteriaceae"),fig_homogenates=c("Geodermatophilaceae","Nocardiaceae","Pseudonocardiaceae","Microbacteriaceae","Kineosporiaceae","Nocardioidaceae","Comamonadaceae","Oxalobacteraceae","Xanthomonadaceae","Spirosomaceae","Deinococcaceae","Chitinophagaceae","Sphingobacteriaceae","Beijerinckiaceae","Hymenobacteraceae","Weeksellaceae","Rhodobacteraceae","Acetobacteraceae","Sphingomonadaceae","Rhizobiaceae","Xanthobacteraceae","Enterobacteriaceae","Erwiniaceae","Bdellovibrionaceae"),fig_surface_washes=c("Geodermatophilaceae","Pseudonocardiaceae","Microbacteriaceae","Kineosporiaceae","Nocardioidaceae","Comamonadaceae","Oxalobacteraceae","Xanthomonadaceae","Spirosomaceae","Beijerinckiaceae","Hymenobacteraceae","Weeksellaceae","Acetobacteraceae","Sphingomonadaceae","Rhizobiaceae","Enterobacteriaceae","Erwiniaceae","Bdellovibrionaceae"))

upset_list <- list(elegans_samples=c("Nocardiaceae","Microbacteriaceae","Alcaligenaceae","Comamonadaceae","Oxalobacteraceae","Xanthomonadaceae","Sphingobacteriaceae","Flavobacteriaceae","Acetobacteraceae","Rhodobacteraceae","Rhizobiaceae","Caulobacteraceae","Beijerinckiaceae","Spirosomaceae","Weeksellaceae","Sphingomonadaceae","Pseudomonadaceae","Erwiniaceae","Enterobacteriaceae"),fig_homogenates=c("Geodermatophilaceae","Nocardiaceae","Pseudonocardiaceae","Microbacteriaceae","Kineosporiaceae","Nocardioidaceae","Comamonadaceae","Oxalobacteraceae","Xanthomonadaceae","Spirosomaceae","Deinococcaceae","Chitinophagaceae","Sphingobacteriaceae","Beijerinckiaceae","Hymenobacteraceae","Weeksellaceae","Rhodobacteraceae","Acetobacteraceae","Sphingomonadaceae","Rhizobiaceae","Xanthobacteraceae","Enterobacteriaceae","Erwiniaceae","Bdellovibrionaceae"),fig_surface_washes=c("Geodermatophilaceae","Pseudonocardiaceae","Microbacteriaceae","Kineosporiaceae","Nocardioidaceae","Comamonadaceae","Oxalobacteraceae","Xanthomonadaceae","Spirosomaceae","Beijerinckiaceae","Hymenobacteraceae","Weeksellaceae","Acetobacteraceae","Sphingomonadaceae","Rhizobiaceae","Enterobacteriaceae","Erwiniaceae","Bdellovibrionaceae"))


#pdf(file="upset.pdf", onefile=FALSE) 
upset(fromList(upset_list), order.by = "freq",text.scale = 2)
#dev.off()
#this is supplemental figure 6

x1 <- unlist(upset_list, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
x1
fromList(upset_list)
    #used to get the families defining various intersections

#all: 310
#elegans: 241
#inopinata: 70
#fig interior: 38
#fig surface: 32

#7-30-25

sample_data(ps3)


bigtab <- as.data.frame(table(sample_data(ps3)$nematode_species))
smalltab <-as.data.frame(table(sample_data(ps3)$category))

bigtab$eighty <- bigtab$Freq*0.8
smalltab$eighty <- smalltab$Freq*0.8

bigtab
#          Var1 Freq eighty
#1   C. elegans  124   99.2
#2 C. inopinata   70   56.0

smalltab
#                   Var1 Freq eighty
#1     dirksen_substrate   64   51.2
#2      samuel_substrate   60   48.0
#3 woodruff_fig_interior   38   30.4
#4  woodruff_fig_surface   32   25.6

nrow(sample_data(ps3))
#[1] 194
nrow(sample_data(ps3))*0.8
#[1] 155.2
eighty_all_fam <- prevdffam[prevdffam$Prevalence > 155,]
eighty_elegans_fam <- elegansfamdf[elegansfamdf$Prevalence > 99,]
eighty_fig_all_fam <- inopinatafamdf[inopinatafamdf$Prevalence > 55,]
eighty_fig_interior_fam <- figintfamdf[figintfamdf$Prevalence > 30,]
eighty_fig_exterior_fam <- figextfamdf[figextfamdf$Prevalence > 25,]


cat(eighty_all_fam$Family, sep = "\n")
#Nocardiaceae
#Microbacteriaceae
#Nocardioidaceae
#Comamonadaceae
#Oxalobacteraceae
#Xanthomonadaceae
#Spirosomaceae
#Sphingobacteriaceae
#Acetobacteraceae
#Beijerinckiaceae
#Weeksellaceae
#Sphingomonadaceae
#Rhizobiaceae
cat(eighty_elegans_fam$Family, sep = "\n")
#Nocardiaceae
#Microbacteriaceae
#Nocardioidaceae
#Comamonadaceae
#Oxalobacteraceae
#Rhodanobacteraceae
#Xanthomonadaceae
#Sphingobacteriaceae
#Chitinophagaceae
#Flavobacteriaceae
#Acetobacteraceae
#Rhodobacteraceae
#Caulobacteraceae
#Beijerinckiaceae
#Spirosomaceae
#Weeksellaceae
#Sphingomonadaceae
#Rhizobiaceae
#Pseudomonadaceae
cat(eighty_fig_all_fam$Family, sep = "\n")
#Geodermatophilaceae
#Pseudonocardiaceae
#Microbacteriaceae
#Kineosporiaceae
#Nocardioidaceae
#Comamonadaceae
#Oxalobacteraceae
#Xanthomonadaceae
#Spirosomaceae
#Beijerinckiaceae
#Hymenobacteraceae
#Weeksellaceae
#Acetobacteraceae
#Sphingomonadaceae
#Rhizobiaceae
#Enterobacteriaceae
#Bdellovibrionaceae
cat(eighty_fig_interior_fam$Family, sep = "\n")
#Nocardiaceae
#Pseudonocardiaceae
#Microbacteriaceae
#Kineosporiaceae
#Nocardioidaceae
#Comamonadaceae
#Oxalobacteraceae
#Xanthomonadaceae
#Spirosomaceae
#Chitinophagaceae
#Sphingobacteriaceae
#Beijerinckiaceae
#Hymenobacteraceae
#Weeksellaceae
#Rhodobacteraceae
#Acetobacteraceae
#Sphingomonadaceae
#Rhizobiaceae
#Xanthobacteraceae
#Enterobacteriaceae
#Bdellovibrionaceae
cat(eighty_fig_exterior_fam$Family, sep = "\n")
#Geodermatophilaceae
#Pseudonocardiaceae
#Microbacteriaceae
#Kineosporiaceae
#Nocardioidaceae
#Comamonadaceae
#Xanthomonadaceae
#Spirosomaceae
#Beijerinckiaceae
#Hymenobacteraceae
#Acetobacteraceae
#Sphingomonadaceae
#Rhizobiaceae

cat(unique(c(eighty_fig_all_fam$Family,eighty_fig_interior_fam$Family,eighty_fig_exterior_fam$Family)),sep = "\n")
Geodermatophilaceae
Pseudonocardiaceae
Microbacteriaceae
Kineosporiaceae
Nocardioidaceae
Comamonadaceae
Oxalobacteraceae
Xanthomonadaceae
Spirosomaceae
Beijerinckiaceae
Hymenobacteraceae
Weeksellaceae
Acetobacteraceae
Sphingomonadaceae
Rhizobiaceae
Enterobacteriaceae
Bdellovibrionaceae
Nocardiaceae
Chitinophagaceae
Sphingobacteriaceae
Rhodobacteraceae
Xanthobacteraceae

upset_list <- list(elegans_samples=eighty_elegans_fam$Family, fig_homogenates=eighty_fig_interior_fam$Family, fig_surface_washes=eighty_fig_exterior_fam$Family)


upset(fromList(upset_list), order.by = "freq",text.scale = 2)

allrbind <- rbind(eighty_all_fam,eighty_elegans_fam,eighty_fig_all_fam,eighty_fig_interior_fam,eighty_fig_exterior_fam)

#get just the unique families
allforplot <- subset(prevdffam,Family %in% unique(allrbind$Family))
elgforplot <- subset(elegansfamdf,Family %in% unique(allrbind$Family))
allfigfforplot <- subset(inopinatafamdf,Family %in% unique(allrbind$Family))
inffigforplot <- subset(figintfamdf,Family %in% unique(allrbind$Family))
extfigforplot <- subset(figextfamdf,Family %in% unique(allrbind$Family))

allforplot$Group <- "All samples"
elgforplot$Group <- "C. elegans substrates"
allfigfforplot$Group <- "All fig samples"
inffigforplot$Group <- "Fig suspensions"
extfigforplot$Group <- "Fig surface washes"

#all: 194
#elegans: 124
#inopinata: 70
#fig interior: 38
#fig surface: 32

allforplot$Percent_prevalence <- ((allforplot$Prevalence/194)*100)
elgforplot$Percent_prevalence <- ((elgforplot$Prevalence/124)*100)
allfigfforplot$Percent_prevalence <- ((allfigfforplot$Prevalence/70)*100)
inffigforplot$Percent_prevalence <- ((inffigforplot$Prevalence/38)*100)
extfigforplot$Percent_prevalence <- ((extfigforplot$Prevalence/32)*100)


plot2df <- rbind(allforplot,elgforplot,allfigfforplot,inffigforplot,extfigforplot)
#transform abundance for plotting
plot2df$LogAbundance <- log(1 + plot2df$TotalAbundance)

write.table(plot2df,"plot2df_family_table_7-30-25.csv",sep=",")

plot2df$Family <- as.factor(plot2df$Family)
plot2df$Family <- factor(plot2df$Family, levels=c("Geodermatophilaceae", "Kineosporiaceae", "Microbacteriaceae", "Nocardiaceae", "Nocardioidaceae", "Pseudonocardiaceae", "Chitinophagaceae", "Flavobacteriaceae", "Hymenobacteraceae", "Sphingobacteriaceae", "Spirosomaceae", "Weeksellaceae", "Bdellovibrionaceae", "Acetobacteraceae", "Beijerinckiaceae", "Caulobacteraceae", "Rhizobiaceae", "Rhodobacteraceae", "Sphingomonadaceae", "Xanthobacteraceae", "Comamonadaceae", "Enterobacteriaceae", "Oxalobacteraceae", "Pseudomonadaceae", "Rhodanobacteraceae", "Xanthomonadaceae"))


unique(plot2df$Class)
#[1] "Actinobacteria"      "Gammaproteobacteria" "Bacteroidia"
#[4] "Alphaproteobacteria" "Bdellovibrionia"
plot2df$Class <- as.factor(plot2df$Class)

plot2df$Class <- factor(plot2df$Class, levels=c("Actinobacteria","Bacteroidia","Bdellovibrionia","Alphaproteobacteria","Gammaproteobacteria"))


unique(plot2df$Group)

plot2df$Group <- as.factor(plot2df$Group)

plot2df$Group <- factor(plot2df$Group, levels=c("All samples","C. elegans substrates","All fig samples","Fig suspensions","Fig surface washes"))


ggplot(plot2df,aes(y=Percent_prevalence,x=Family)) + geom_hline(yintercept=80,linetype="dotted") + geom_point(aes(shape=Group, colour=Class),size=2)  +theme_half_open(font_size =14) + background_grid(major = c("x"),minor = c("none")) + theme(axis.text.x = element_text(angle = 45, hjust=1,face="italic")) + scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) + scale_colour_manual("Class",values=c("yellow3","#7b3294","#bd0026","darkblue","#3182bd")) + scale_shape_manual(values=c(16,15,1,0,2)) + ylab("Prevalence (%)")

ggsave("elegans_figs_family_prevalence_7-30-25.pdf",height=5.6,width=10,units="in")


ggplot(plot2df,aes(x=Percent_prevalence,y=Family)) + geom_vline(xintercept=80,linetype="dotted") + geom_point(aes(shape=Group, colour=Class),size=2)  +theme_half_open(font_size =14) + background_grid(major = c("y"),minor = c("none")) + theme(axis.text.y = element_text(face="italic")) + scale_x_continuous(limits=c(0,100),breaks=seq(0,100,10)) + scale_colour_manual("Class",values=c("yellow3","#7b3294","#bd0026","darkblue","#3182bd")) + scale_shape_manual(values=c(16,15,1,0,2)) + xlab("Prevalence (%)")


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


####doing this again, 7-31-25


#diversity metrics
#inside/outside

#
sample_sums(ps3)
#ERR1307313 ERR1307319 ERR1307320 ERR1307322 ERR1307326 ERR1307327 ERR1307328
#     42658      60281      43003      46656      45222      52070      47486
#ERR1307329 ERR1307330 ERR1307332 ERR1307333 ERR1307334 ERR1307335 ERR1307336
#     30895      51752      47225      46941      56422      54697      47859
#ERR1307372 ERR1307373 ERR1307380 ERR1307385 ERR1307386 ERR1307388 ERR1307391
#     70544      56433      47026      69972      47830      61780      47530
#ERR1307392 ERR1307397 ERR1307402 ERR1307408 ERR1307409 ERR1307412 ERR1307415
#     52579      52204      50500      73433      43145      43465      38865
#ERR1307416 ERR1307417 ERR1307418 ERR1307419 ERR1307420 ERR1307421 ERR1307422
#     49410      57331      49934      49047      47983      42538      77465
#ERR1307423 ERR1307424 ERR1307425 ERR1307426 ERR1307427 ERR1307428 ERR1307429
#     44034      44914      38414      34675      38423      47905      48279
#ERR1307430 ERR1307431 ERR1307432 ERR1307433 ERR1307434 ERR1307435 ERR1307436
#     49333      36437      65577      44527      43700      36928      44908
#ERR1307437 ERR1307438 ERR1307439 ERR1307440 ERR1307441 ERR1307442 ERR1307443
#     39657      45975      42168      37561      68224      36130      76522
#ERR1307444 ERR1307445 ERR1307446 ERR1307447 ERR1307448 ERR1307449 ERR1307450
#     55084      49679      52601      49146      59190      55406      64219
#ERR1307451        GW1       GW12       GW13       GW14       GW15       GW16
#     54670      34490      71239      71576      42444      25420      31165
#      GW17       GW18       GW19        GW2       GW20       GW21       GW22
#     31642      20413      30935      24564       8400      15950      26015
#      GW23       GW24       GW25       GW26       GW27       GW28       GW29
#     13826      73470      33737      45089      11128      24718      19799
#       GW3       GW30       GW32       GW33       GW34       GW35       GW36
#     46735      28019      32638      17479       7279      14895      41980
#      GW37       GW38       GW39        GW4       GW40       GW41       GW42
#     78718      22088      66635      61244      41816      32834      36479
#      GW43       GW44       GW45       GW46       GW47       GW48       GW49
#     60357      80578      74688      59266      48792      63037      45625
#       GW5       GW50       GW51       GW52       GW53       GW54       GW55
#     69756      57413      26515      16282      14788      75070      19555
#      GW56       GW57       GW58       GW59        GW6       GW60       GW61
#     48311      24824      17794      32183      62902      77843      69398
#      GW62       GW63       GW64       GW65       GW66       GW67       GW68
#     38682      79767      27747       9890      62487      98966      34484
#      GW69        GW7       GW70       GW71       GW72       GW73       GW74
#    109059      67349      16067     107523      45322      62264      23066
#       GW8 SRR5094338 SRR5094346 SRR5094349 SRR5094351 SRR5094352 SRR5094356
#     55762      64349      39415     132661       5918       2005      77688
#SRR5094361 SRR5094362 SRR5094366 SRR5094370 SRR5094374 SRR5094376 SRR5094384
#     48124      41692      80994      26786     218567      50680      59898
#SRR5094388 SRR5094389 SRR5094400 SRR5094401 SRR5094404 SRR5094405 SRR5094408
#      1206      29729     115678     104610       9056      60795     126422
#SRR5094415 SRR5094419 SRR5094420 SRR5094421 SRR5094424 SRR5094428 SRR5094432
#     61593      77326      90942        779      43885      76081      30417
#SRR5094436 SRR5094437 SRR5094438 SRR5094440 SRR5094443 SRR5094445 SRR5094449
#    230702       2219       5784     114991          0     180066      67656
#SRR5094460 SRR5094469 SRR5094471 SRR5094472 SRR5094475 SRR5094478 SRR5094479
#      9793       1844      82142        849      43936      54278     230663
#SRR5094480 SRR5094485 SRR5094489 SRR5094515 SRR5094517 SRR5094525 SRR5094528
#     52010     102072       9145      65966       1210      41572      50000
#SRR5094531 SRR5094532 SRR5094534 SRR5094536 SRR5094551 SRR5094552 SRR5094554
#     36175      29433      51529     156484      57629      20957         22
#SRR5094565 SRR5094566 SRR5094569 SRR5094574 SRR5094576
#     36997      52994      95369       9548      12296

#exclude some with way too few reads

rcou <- sample_sums(ps3)

names(which(rcou < 9000))
# [1] "GW20"       "GW34"       "SRR5094351" "SRR5094352" "SRR5094388"
# [6] "SRR5094421" "SRR5094437" "SRR5094438" "SRR5094443" "SRR5094469"
#[11] "SRR5094472" "SRR5094517" "SRR5094554"
ps3_no_small <- prune_samples((!sample_names(ps3) %in% c("GW20","GW34","SRR5094351","SRR5094352","SRR5094388","SRR5094421","SRR5094437","SRR5094438","SRR5094443","SRR5094469","SRR5094472","SRR5094517","SRR5094554")), ps3)

rcou2 <- sample_sums(ps3_no_small)

summary(rcou2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   9056   36175   47983   54259   63037  230702

#Subsample reads
#ps_rare <- phyloseq::rarefy_even_depth(ps3_no_small, rngseed = 123, replace = FALSE)

set.seed(666)
rseeds <- runif(1000, min=0, max=1000000)
rseeds <- round(rseeds)

#head(phyloseq::sample_sums(ps_rare))
#  GW1  GW12  GW13  GW14  GW15  GW16
#12252 12252 12252 12252 12252 12252

library(picante)

big_div_df <- NULL

for (i in rseeds){
 ps_rare <- phyloseq::rarefy_even_depth(ps3_no_small, rngseed = i, replace = FALSE)
 adiv <- data.frame(
   "Number_of_OTU" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
   "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
   "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps_rare)))), tree = phyloseq::phy_tree(ps_rare))[, 1],
   "category" = phyloseq::sample_data(ps_rare)$category,
   "rngseed" = i)
 big_div_df <- rbind(big_div_df,adiv)
}

#write.table(big_div_df,"diversity_metrics_rarefy_1000_times_4-2025.tsv",quote=FALSE,sep='\t')

write.table(big_div_df,"diversity_metrics_rarefy_1000_times_7-2025.tsv",quote=FALSE,sep='\t')


#big_div_df <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/diversity_metrics_rarefy_1000_times_4-2025.tsv",sep='\t',header=TRUE)

big_div_df <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/diversity_metrics_rarefy_1000_times_7-2025.tsv",sep='\t',header=TRUE)


big_div_df$sample_id <- rownames(big_div_df)

shannon_means <- aggregate(Shannon ~ sample_id, FUN=mean, data=big_div_df)

Phylogenetic_means <- aggregate(Phylogenetic ~ sample_id, FUN=mean, data=big_div_df)

Number_of_OTU_means <- aggregate(Observed ~ sample_id, FUN=mean, data=big_div_df)

readcount <- as.data.frame(rcou2) 

readcount$sample_id <- rownames(readcount)

metdat2$sample_id <- metdat2$sample.id


readmerge <- merge(readcount,metdat2)

names(readmerge)[names(readmerge) == 'rcou2'] <- 'read_count'

divmerge <- merge(readmerge,shannon_means)

divmerge2 <- merge(divmerge, Phylogenetic_means)

divmerge3 <- merge(divmerge2, Number_of_OTU_means)
names(divmerge3)[names(divmerge3) == "Observed"] <- "Number_of_OTU"


write.table(divmerge3,"alpha_diversity_1000_rarify_means_8-2025.tsv",sep='\t',quote=FALSE,row.names=FALSE)
divmerge3 <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/alpha_diversity_1000_rarify_means_8-2025.tsv",sep='\t',header=TRUE)

ps3_no_small






library(picante)

adiv_actual <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps3_no_small, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps3_no_small, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps3_no_small)))), tree = phyloseq::phy_tree(ps3_no_small))[, 1],
  "category" = phyloseq::sample_data(ps3_no_small)$category,
  "nematode_species" = phyloseq::sample_data(ps3_no_small)$nematode_species)


#okay, I did an investigation. Singletons re prevalence are included in this data. Two (and greater)-read ASVs are included, as well. However, Qiime2 removes singletons wrt ASV by READ (i.e., ASVs that are represented only by a single read)! I kinda agree this is appropriate. What do we do with ASVs where there is just one read? Two-read (three read, four read, etc.) ASVs DO exist in this data. So low abundance is merely cut off at two reads. I am fine ignoring this warning.

adiv_actual$sample_id <- rownames(adiv_actual)
names(adiv_actual)[names(adiv_actual) == "Observed"] <- "Number_of_ASV"

#write.table(adiv_actual,"alpha_diversity_no_rarefy_8-2025.tsv",sep="\t",quote=FALSE,row.names=FALSE)
adiv_actual <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/alpha_diversity_no_rarefy_8-2025.tsv",sep='\t',header=TRUE)


adiv_actual_melt <- reshape2::melt(adiv_actual, id.vars= c("sample_id","category","nematode_species"),measure.vars = c("Number_of_ASV", "Shannon","Phylogenetic"))

levels(adiv_actual_melt$variable)[levels(adiv_actual_melt$variable)=="Number_of_ASV"] <- "Number of ASV"

adiv_actual_melt$category <- as.factor(adiv_actual_melt$category)

levels(adiv_actual_melt$category)[levels(adiv_actual_melt$category)=="dirksen_substrate"] <- "Dirksen substrate"
levels(adiv_actual_melt$category)[levels(adiv_actual_melt$category)=="samuel_substrate"] <- "Samuel substrate"
levels(adiv_actual_melt$category)[levels(adiv_actual_melt$category)=="woodruff_fig_surface"] <- "Fig surface wash"
levels(adiv_actual_melt$category)[levels(adiv_actual_melt$category)=="woodruff_fig_interior"] <- "Fig suspension"

ggplot(adiv_actual_melt, aes(x=category,y=value)) + geom_sina(scale="width",aes(colour=nematode_species),size=0.85,alpha=0.75) + stat_summary(aes(group=category),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + facet_wrap(~variable,ncol=1,scales="free") +theme_cowplot() + scale_y_continuous(limits = c(0,NA)) + scale_color_brewer(palette="Set1") + theme(strip.background = element_rect(colour="white", fill="white"))  + xlab("Sample type") +ylab("Diversity metric value") + labs(colour="Nematode species") + theme(legend.text = element_text(face ="italic"))


ggsave("no_rarefy_alpha_diversity_sina_8-2025.pdf",unit="in",height=7,width=10)


adiv_actual_merge <- merge(readcount,adiv_actual)

names(adiv_actual_merge)[names(adiv_actual_merge) == 'rcou2'] <- 'read_count'



adiv_actual_melt2 <- reshape2::melt(adiv_actual_merge, id.vars= c("sample_id","category","nematode_species","read_count"),measure.vars = c("Number_of_ASV", "Shannon","Phylogenetic"))

levels(adiv_actual_melt2$variable)[levels(adiv_actual_melt2$variable)=="Number_of_ASV"] <- "Number of ASV"

adiv_actual_melt2$category <- as.factor(adiv_actual_melt2$category)

levels(adiv_actual_melt2$category)[levels(adiv_actual_melt2$category)=="dirksen_substrate"] <- "Dirksen substrate"
levels(adiv_actual_melt2$category)[levels(adiv_actual_melt2$category)=="samuel_substrate"] <- "Samuel substrate"
levels(adiv_actual_melt2$category)[levels(adiv_actual_melt2$category)=="woodruff_fig_surface"] <- "Fig surface wash"
levels(adiv_actual_melt2$category)[levels(adiv_actual_melt2$category)=="woodruff_fig_interior"] <- "Fig suspension"



ggplot(adiv_actual_melt2, aes(x=read_count,y=value)) + geom_point(aes(colour=category),size=0.85,alpha=0.75) + facet_wrap(~variable,ncol=1,scales="free") +theme_cowplot() + scale_color_brewer(palette="Set1") + theme(strip.background = element_rect(colour="white", fill="white"))  + xlab("Read count") +ylab("Diversity metric value") + labs(colour="Sample type") + scale_x_continuous(labels = scales::label_comma())


ggsave("no_rarefy_alpha_diversity_read_count_scatterplot_8-25-25.pdf",unit="in",height=7,width=10)










divmerge3_melt <- reshape2::melt(divmerge3, id.vars= c("sample_id","category","nematode_species"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic"))


levels(divmerge3_melt$variable)[levels(divmerge3_melt$variable)=="Number_of_OTU"] <- "Number of ASV"

divmerge3_melt$category <- as.factor(divmerge3_melt$category)

levels(divmerge3_melt$category)[levels(divmerge3_melt$category)=="dirksen_substrate"] <- "Dirksen substrate"
levels(divmerge3_melt$category)[levels(divmerge3_melt$category)=="samuel_substrate"] <- "Samuel substrate"
levels(divmerge3_melt$category)[levels(divmerge3_melt$category)=="woodruff_fig_surface"] <- "Fig surface wash"
levels(divmerge3_melt$category)[levels(divmerge3_melt$category)=="woodruff_fig_interior"] <- "Fig suspension"


ggplot(divmerge3_melt, aes(x=category,y=value)) + geom_sina(scale="width",aes(colour=nematode_species),size=0.85,alpha=0.75) + stat_summary(aes(group=category),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + facet_wrap(~variable,ncol=1,scales="free") +theme_cowplot() + scale_y_continuous(limits = c(0,NA)) + scale_color_brewer(palette="Set1") + theme(strip.background = element_rect(colour="white", fill="white"))  + xlab("Sample type") +ylab("Diversity metric value") + labs(colour="Nematode species") + theme(legend.text = element_text(face ="italic"))

ggsave("rarefy_1000_alpha_diversity_sina_8-25-25.pdf",unit="in",height=7,width=10)





#####stats comparing alpha diversity


#divmerge3 <- read.table("alpha_diversity_1000_rarify_means.tsv",header=TRUE,sep="\t")

#adiv_actual <- read.table("alpha_diversity_no_rarefy.tsv",sep="\t",header=TRUE)

rcdf <- data.frame(sample_id=divmerge3$sample_id,read_count=divmerge3$read_count)

adiv_actual2 <- merge(adiv_actual,rcdf)


#okay, no rarefaction

library(tidyverse)
library(ggpubr)
library(rstatix)

#substrate categories

adiv_actual2 %>% kruskal_test(Number_of_ASV ~ category)
#  .y.               n statistic    df        p method
#* <chr>         <int>     <dbl> <int>    <dbl> <chr>
#1 Number_of_ASV   181      124.     3 8.41e-27 Kruskal-Wallis

as.data.frame(adiv_actual2 %>% dunn_test(Number_of_ASV ~ category, p.adjust.method = "BH") )

#            .y.                group1                group2 n1 n2  statistic
#1 Number_of_ASV     dirksen_substrate      samuel_substrate 64 49 -7.9959222
#2 Number_of_ASV     dirksen_substrate woodruff_fig_interior 64 37 -8.2274272
#3 Number_of_ASV     dirksen_substrate  woodruff_fig_surface 64 31 -9.1465964
#4 Number_of_ASV      samuel_substrate woodruff_fig_interior 49 37 -0.8326145
#5 Number_of_ASV      samuel_substrate  woodruff_fig_surface 49 31 -2.1075238
#6 Number_of_ASV woodruff_fig_interior  woodruff_fig_surface 37 31 -1.2416272
#             p        p.adj p.adj.signif
#1 1.286076e-15 2.572152e-15         ****
#2 1.912773e-16 5.738319e-16         ****
#3 5.875556e-20 3.525334e-19         ****
#4 4.050622e-01 4.050622e-01           ns
#5 3.507220e-02 5.260830e-02           ns
#6 2.143741e-01 2.572490e-01           ns


adiv_actual2 %>% kruskal_test(Shannon ~ category)
#  .y.         n statistic    df        p method
#* <chr>   <int>     <dbl> <int>    <dbl> <chr>
#1 Shannon   181      116.     3 6.96e-25 Kruskal-Wallis

as.data.frame(adiv_actual2 %>% dunn_test(Shannon ~ category, p.adjust.method = "BH"))
#      .y.                group1                group2 n1 n2  statistic
#1 Shannon     dirksen_substrate      samuel_substrate 64 49 -8.4653508
#2 Shannon     dirksen_substrate woodruff_fig_interior 64 37 -8.4735215
#3 Shannon     dirksen_substrate  woodruff_fig_surface 64 31 -7.6206465
#4 Shannon      samuel_substrate woodruff_fig_interior 49 37 -0.6568331
#5 Shannon      samuel_substrate  woodruff_fig_surface 49 31 -0.2642281
#6 Shannon woodruff_fig_interior  woodruff_fig_surface 37 31  0.3384915
#             p        p.adj p.adj.signif
#1 2.553829e-17 7.661487e-17         ****
#2 2.380839e-17 7.661487e-17         ****
#3 2.524082e-14 5.048165e-14         ****
#4 5.112883e-01 7.669324e-01           ns
#5 7.916041e-01 7.916041e-01           ns
#6 7.349928e-01 7.916041e-01           ns

adiv_actual2 %>% kruskal_test(Phylogenetic ~ category)
## A tibble: 1 × 6
#  .y.              n statistic    df        p method
#* <chr>        <int>     <dbl> <int>    <dbl> <chr>
#1 Phylogenetic   181      123.     3 1.55e-26 Kruskal-Wallis

as.data.frame(adiv_actual2 %>% dunn_test(Phylogenetic ~ category, p.adjust.method = "BH"))
#           .y.                group1                group2 n1 n2 statistic
#1 Phylogenetic     dirksen_substrate      samuel_substrate 64 49 -8.763494
#2 Phylogenetic     dirksen_substrate woodruff_fig_interior 64 37 -7.318911
#3 Phylogenetic     dirksen_substrate  woodruff_fig_surface 64 31 -9.045605
#4 Phylogenetic      samuel_substrate woodruff_fig_interior 49 37  0.697870
#5 Phylogenetic      samuel_substrate  woodruff_fig_surface 49 31 -1.376330
#6 Phylogenetic woodruff_fig_interior  woodruff_fig_surface 37 31 -1.921467
#             p        p.adj p.adj.signif
#1 1.892903e-18 5.678710e-18         ****
#2 2.499922e-13 4.999844e-13         ****
#3 1.488385e-19 8.930311e-19         ****
#4 4.852585e-01 4.852585e-01           ns
#5 1.687196e-01 2.024635e-01           ns
#6 5.467281e-02 8.200922e-02           ns




adiv_actual2 %>% kruskal_test(Number_of_ASV ~ nematode_species)
#  .y.               n statistic    df        p method
#* <chr>         <int>     <dbl> <int>    <dbl> <chr>
#1 Number_of_ASV   181      59.0     1 1.58e-14 Kruskal-Wallis

as.data.frame(adiv_actual2 %>% dunn_test(Number_of_ASV ~ nematode_species, p.adjust.method = "BH") )
#            .y.     group1       group2  n1 n2 statistic            p
#1 Number_of_ASV C. elegans C. inopinata 113 68  -7.68066 1.582708e-14
#         p.adj p.adj.signif
#1 1.582708e-14         ****

wilcox.test(adiv_actual2[adiv_actual2$nematode_species == "C. elegans",]$Number_of_ASV,adiv_actual2[adiv_actual2$nematode_species == "C. inopinata",]$Number_of_ASV)

#    Wilcoxon rank sum test with continuity correction
#
#data:  adiv_actual2[adiv_actual2$nematode_species == "C. elegans", ]$Number_of_ASV and adiv_actual2[adiv_actual2$nematode_species == "C. inopinata", ]$Number_of_ASV
#W = 6464, p-value = 1.601e-14
#alternative hypothesis: true location shift is not equal to 0


adiv_actual2 %>% kruskal_test(Shannon ~ nematode_species)
## A tibble: 1 × 6
#  .y.         n statistic    df        p method
#* <chr>   <int>     <dbl> <int>    <dbl> <chr>
#1 Shannon   181      43.8     1 3.66e-11 Kruskal-Wallis

as.data.frame(adiv_actual2 %>% dunn_test(Shannon ~ nematode_species, p.adjust.method = "BH") )
#      .y.     group1       group2  n1 n2 statistic            p        p.adj
#1 Shannon C. elegans C. inopinata 113 68 -6.617249 3.659446e-11 3.659446e-11
#  p.adj.signif
#1         ****


wilcox.test(adiv_actual2[adiv_actual2$nematode_species == "C. elegans",]$Shannon,adiv_actual2[adiv_actual2$nematode_species == "C. inopinata",]$Shannon)

#    Wilcoxon rank sum test with continuity correction
#
#data:  adiv_actual2[adiv_actual2$nematode_species == "C. elegans", ]$Shannon and adiv_actual2[adiv_actual2$nematode_species == "C. inopinata", ]$Shannon
#W = 6101, p-value = 3.696e-11
#alternative hypothesis: true location shift is not equal to 0


adiv_actual2 %>% kruskal_test(Phylogenetic ~ nematode_species)
#  .y.              n statistic    df     p method
#* <chr>        <int>     <dbl> <int> <dbl> <chr>
#1 Phylogenetic   293      1.25     1 0.264 Kruskal-Wallis
as.data.frame(adiv_actual2 %>% dunn_test(Phylogenetic ~ nematode_species, p.adjust.method = "BH") )

#           .y.     group1       group2  n1 n2 statistic         p     p.adj
#1 Phylogenetic C. elegans C. inopinata 225 68 -1.117191 0.2639128 0.2639128
#  p.adj.signif
#1           ns


wilcox.test(adiv_actual2[adiv_actual2$nematode_species == "C. elegans",]$Phylogenetic,adiv_actual2[adiv_actual2$nematode_species == "C. inopinata",]$Phylogenetic)

#    Wilcoxon rank sum test with continuity correction
#
#data:  adiv_actual2[adiv_actual2$nematode_species == "C. elegans", ]$Phylogenetic and adiv_actual2[adiv_actual2$nematode_species == "C. inopinata", ]$Phylogenetic
#W = 6074, p-value = 6.289e-11
#alternative hypothesis: true location shift is not equal to 0




aggregate(Number_of_ASV ~ nematode_species, FUN=summary, data=adiv_actual)
#  nematode_species Number_of_ASV.Min. Number_of_ASV.1st Qu.
#1       C. elegans            26.0000              198.0000
#2     C. inopinata            28.0000               77.5000
#  Number_of_ASV.Median Number_of_ASV.Mean Number_of_ASV.3rd Qu.
#1             805.0000          1057.1239             1874.0000
#2             119.5000           163.2353              197.7500
#  Number_of_ASV.Max.
#1          3214.0000
#2           843.0000

aggregate(Shannon ~ nematode_species, FUN=summary, data=adiv_actual)
#  nematode_species Shannon.Min. Shannon.1st Qu. Shannon.Median Shannon.Mean
#1       C. elegans     0.820531        3.546060       5.464499     5.062682
#2     C. inopinata     1.063852        2.459432       3.352614     3.179015
#  Shannon.3rd Qu. Shannon.Max.
#1        6.615591     7.409771
#2        3.902805     5.726773

aggregate(Phylogenetic ~ nematode_species, FUN=summary, data=adiv_actual)
#  nematode_species Phylogenetic.Min. Phylogenetic.1st Qu. Phylogenetic.Median
#1       C. elegans          9.507724            25.225062           60.003690
#2     C. inopinata         13.312480            17.601942           21.767647
#  Phylogenetic.Mean Phylogenetic.3rd Qu. Phylogenetic.Max.
#1         73.610620           123.417195        177.044051
#2         24.246995            26.966993         61.778419



cohen.d(adiv_actual2[adiv_actual2$nematode_species == "C. elegans",]$Number_of_ASV,adiv_actual2[adiv_actual2$nematode_species == "C. inopinata",]$Number_of_ASV)

#Cohen's d
#
#d estimate: 1.265355 (large)
#95 percent confidence interval:
#    lower     upper
#0.9352851 1.5954251

cohen.d(adiv_actual2[adiv_actual2$nematode_species == "C. elegans",]$Shannon,adiv_actual2[adiv_actual2$nematode_species == "C. inopinata",]$Shannon)
#Cohen's d
#
#d estimate: 1.250734 (large)
#95 percent confidence interval:
#    lower     upper
#0.9212637 1.5802036


cohen.d(adiv_actual2[adiv_actual2$nematode_species == "C. elegans",]$Phylogenetic,adiv_actual2[adiv_actual2$nematode_species == "C. inopinata",]$Phylogenetic)

#Cohen's d
#
#d estimate: 1.213316 (large)
#95 percent confidence interval:
#    lower     upper
#0.8853544 1.5412772



as.data.frame(adiv_actual2 %>% cohens_d(Number_of_ASV ~ category))

#            .y.                group1                group2   effsize n1 n2
#1 Number_of_ASV     dirksen_substrate      samuel_substrate 2.9559689 64 49
#2 Number_of_ASV     dirksen_substrate woodruff_fig_interior 3.0524958 64 37
#3 Number_of_ASV     dirksen_substrate  woodruff_fig_surface 3.3222900 64 31
#4 Number_of_ASV      samuel_substrate woodruff_fig_interior 0.1871137 49 37
#5 Number_of_ASV      samuel_substrate  woodruff_fig_surface 0.7953111 49 31
#6 Number_of_ASV woodruff_fig_interior  woodruff_fig_surface 0.6107700 37 31
#   magnitude
#1      large
#2      large
#3      large
#4 negligible
#5   moderate
#6   moderate

as.data.frame(adiv_actual2 %>% cohens_d(Shannon ~ category))

#      .y.                group1                group2     effsize n1 n2
#1 Shannon     dirksen_substrate      samuel_substrate  3.13834169 64 49
#2 Shannon     dirksen_substrate woodruff_fig_interior  3.38985070 64 37
#3 Shannon     dirksen_substrate  woodruff_fig_surface  3.47821321 64 31
#4 Shannon      samuel_substrate woodruff_fig_interior  0.21251011 49 37
#5 Shannon      samuel_substrate  woodruff_fig_surface  0.19332826 49 31
#6 Shannon woodruff_fig_interior  woodruff_fig_surface -0.02486423 37 31
#   magnitude
#1      large
#2      large
#3      large
#4      small
#5 negligible
#6 negligible


as.data.frame(adiv_actual2 %>% cohens_d(Phylogenetic ~ category))


#           .y.                group1                group2    effsize n1 n2
#1 Phylogenetic     dirksen_substrate      samuel_substrate  3.0324134 64 49
#2 Phylogenetic     dirksen_substrate woodruff_fig_interior  3.0180820 64 37
#3 Phylogenetic     dirksen_substrate  woodruff_fig_surface  3.3713329 64 31
#4 Phylogenetic      samuel_substrate woodruff_fig_interior -0.1140361 49 37
#5 Phylogenetic      samuel_substrate  woodruff_fig_surface  0.5811127 49 31
#6 Phylogenetic woodruff_fig_interior  woodruff_fig_surface  0.8038818 37 31
#   magnitude
#1      large
#2      large
#3      large
#4 negligible
#5   moderate
#6      large








ggplot(adiv_actual_melt, aes(x=nematode_species,y=value)) + geom_sina(scale="width",aes(colour=nematode_species),size=0.85,alpha=0.75) + stat_summary(aes(group=nematode_species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black") + facet_wrap(~variable,ncol=1,scales="free") +theme_cowplot() + scale_y_continuous(limits = c(0,NA)) + scale_color_brewer(palette="Set1") + theme(strip.background = element_rect(colour="white", fill="white"))  + xlab("Nematode species") +ylab("Diversity metric value") + labs(colour="Nematode species") + theme(legend.text = element_text(face ="italic"),axis.text.x = element_text(face = "italic"))


ggsave("no_rarefy_alpha_diversity_sina_nematode_species_8-26-2025.pdf",unit="in",height=7,width=10)




#read counts

summary(lm(Number_of_ASV~read_count,data=adiv_actual2))

#Call:
#lm(formula = Number_of_ASV ~ read_count, data = adiv_actual2)
#
#Residuals:
#   Min     1Q Median     3Q    Max
#-607.8 -456.6 -338.2  104.6 2728.6
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)  6.556e+02  7.336e+01   8.937   <2e-16 ***
#read_count  -2.224e-03  9.401e-04  -2.366   0.0187 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 732.9 on 291 degrees of freedom
#Multiple R-squared:  0.01887,   Adjusted R-squared:  0.0155
#F-statistic: 5.596 on 1 and 291 DF,  p-value: 0.01866

#-0.002224 ASV/read! a negative correlation! 

summary(lm(Shannon~read_count,data=adiv_actual2))
#Call:
#lm(formula = Shannon ~ read_count, data = adiv_actual2)
#
#Residuals:
#    Min      1Q  Median      3Q     Max
#-3.6087 -1.4858 -0.1765  1.5048  4.0651
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)  4.312e+00  1.873e-01  23.022  < 2e-16 ***
#read_count  -1.264e-05  2.400e-06  -5.266 2.71e-07 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.871 on 291 degrees of freedom
#Multiple R-squared:  0.08702,   Adjusted R-squared:  0.08388
#F-statistic: 27.73 on 1 and 291 DF,  p-value: 2.71e-07

#-.00001264e ASV/read! a negative correlation! 

summary(lm(Phylogenetic~read_count,data=adiv_actual2))
#Call:
#lm(formula = Phylogenetic ~ read_count, data = adiv_actual2)
#
#Residuals:
#    Min      1Q  Median      3Q     Max
#-41.018 -26.456 -18.565   7.003 136.103
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)  5.193e+01  4.210e+00  12.337  < 2e-16 ***
#read_count  -1.436e-04  5.394e-05  -2.663  0.00818 **
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 42.06 on 291 degrees of freedom
#Multiple R-squared:  0.02378,   Adjusted R-squared:  0.02043
#F-statistic:  7.09 on 1 and 291 DF,  p-value: 0.008184

# -.0001436 ASV/read! a negative correlation! 

#this seems to me to be driven by these high read count, moderate diversity samples in the dirksen and samuel data sets. basically, no real impact of read count on diversity, which is surprising to me. i think we can use the un-rarefied data in discussions of alpha diversity.

#with rarefaction


divmerge3 %>% kruskal_test(Number_of_OTU ~ category)
#  .y.               n statistic    df     p method
#* <chr>         <int>     <dbl> <int> <dbl> <chr>
#1 Number_of_OTU   293      5.45     3 0.141 Kruskal-Wallis

as.data.frame(divmerge3 %>% dunn_test(Number_of_OTU ~ category, p.adjust.method = "BH") )
#            .y.                group1                group2  n1 n2   statistic
#1 Number_of_OTU     dirksen_substrate      samuel_substrate 176 49 -0.08096966
#2 Number_of_OTU     dirksen_substrate woodruff_fig_interior 176 37 -0.79176059
#3 Number_of_OTU     dirksen_substrate  woodruff_fig_surface 176 31 -2.25755562
#4 Number_of_OTU      samuel_substrate woodruff_fig_interior  49 37 -0.59742091
#5 Number_of_OTU      samuel_substrate  woodruff_fig_surface  49 31 -1.85912032
#6 Number_of_OTU woodruff_fig_interior  woodruff_fig_surface  37 31 -1.21787972
#           p     p.adj p.adj.signif
#1 0.93546608 0.9354661           ns
#2 0.42850029 0.6427504           ns
#3 0.02397338 0.1438403           ns
#4 0.55022640 0.6602717           ns
#5 0.06301009 0.1890303           ns
#6 0.22326968 0.4465394           ns


divmerge3 %>% kruskal_test(Shannon ~ category)
## A tibble: 1 × 6
#  .y.         n statistic    df     p method
#* <chr>   <int>     <dbl> <int> <dbl> <chr>
#1 Shannon   293     0.399     3 0.941 Kruskal-Wallis

as.data.frame(divmerge3 %>% dunn_test(Shannon ~ category, p.adjust.method = "BH"))
#      .y.                group1                group2  n1 n2   statistic
#1 Shannon     dirksen_substrate      samuel_substrate 176 49  0.26850956
#2 Shannon     dirksen_substrate woodruff_fig_interior 176 37 -0.50635030
#3 Shannon     dirksen_substrate  woodruff_fig_surface 176 31 -0.07065688
#4 Shannon      samuel_substrate woodruff_fig_interior  49 37 -0.61960280
#5 Shannon      samuel_substrate  woodruff_fig_surface  49 31 -0.24895671
#6 Shannon woodruff_fig_interior  woodruff_fig_surface  37 31  0.31958299
#          p     p.adj p.adj.signif
#1 0.7883071 0.9436708           ns
#2 0.6126108 0.9436708           ns
#3 0.9436708 0.9436708           ns
#4 0.5355193 0.9436708           ns
#5 0.8033943 0.9436708           ns
#6 0.7492845 0.9436708           ns

divmerge3 %>% kruskal_test(Phylogenetic ~ category)
#  .y.              n statistic    df     p method
#* <chr>        <int>     <dbl> <int> <dbl> <chr>
#1 Phylogenetic   293      4.70     3 0.195 Kruskal-Wallis

as.data.frame(divmerge3 %>% dunn_test(Phylogenetic ~ category, p.adjust.method = "BH"))

#           .y.                group1                group2  n1 n2  statistic
#1 Phylogenetic     dirksen_substrate      samuel_substrate 176 49 -1.1582260
#2 Phylogenetic     dirksen_substrate woodruff_fig_interior 176 37  0.3597346
#3 Phylogenetic     dirksen_substrate  woodruff_fig_surface 176 31 -1.8221967
#4 Phylogenetic      samuel_substrate woodruff_fig_interior  49 37  1.1576939
#5 Phylogenetic      samuel_substrate  woodruff_fig_surface  49 31 -0.7313971
#6 Phylogenetic woodruff_fig_interior  woodruff_fig_surface  37 31 -1.7249114
#           p     p.adj p.adj.signif
#1 0.24677181 0.3704835           ns
#2 0.71904562 0.7190456           ns
#3 0.06842514 0.2536303           ns
#4 0.24698899 0.3704835           ns
#5 0.46453664 0.5574440           ns
#6 0.08454345 0.2536303           ns












####compare genera, elegans and inopinata substrates
    #relative abundances first (this sucks just do clr transformed)
    #microbiome::transform(ps4, "clr")


ps_genus <- tax_glom(ps3, "Genus", NArm = TRUE)


ps_clr_genus <- microbiome::transform(ps_genus, "clr")



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
 if (nrow(dat) >10){elegans <- dat[dat$nematode_species == "C. elegans",]
 inopinata <- dat[dat$nematode_species == "C. inopinata",]
 genus <- rbind(genus, i)
 effect_size_lower <- rbind(effect_size_lower, cohen.d(na.omit(elegans$Abundance),na.omit(inopinata$Abundance))$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(na.omit(elegans$Abundance),na.omit(inopinata$Abundance))$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(na.omit(elegans$Abundance),na.omit(inopinata$Abundance))$conf.int[2]) 
 magn <- rbind(magn, cohen.d(na.omit(elegans$Abundance),na.omit(inopinata$Abundance))$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(elegans$Abundance,inopinata$Abundance,exact = FALSE)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(elegans$Abundance,inopinata$Abundance,exact = FALSE)$p.value)}
}

#bleh <- mphyseq[mphyseq$Genus == "Pantoea",]
#blehelegans <- bleh[bleh$nematode_species == "C. elegans",]
#blehinopinata <- bleh[bleh$nematode_species == "C. inopinata",]
#cohen.d(na.omit(blehelegans$Abundance),na.omit(blehinopinata$Abundance))
#blehelegans$Abundance


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

write.table(stat_df, "wilcox_tests_Genus_worm_presence_eleg_inop.tsv",sep="\t",row.names=FALSE, quote=FALSE)

sigdf <- stat_df[stat_df$p.adj < 0.05,]
write.table(sigdf, "wilcox_tests_Genus_worm_presence_eleg_inop_padj_less_0.05.tsv",sep="\t",row.names=FALSE, quote=FALSE)




mphyseqbigeff <- mphyseq[mphyseq$Genus %in% c("Quadrisphaera", "Kineococcus", "Klenkia", "Aureimonas", "Spirosoma", "Methylobacterium-Methylorubrum", "Deinococcus", "Hymenobacter", "Azorhizobium", "Actinomycetospora","Duganella", "Luteolibacter", "Pseudomonas", "Cellvibrio", "Stenotrophomonas", "Rickettsiella", "Achromobacter", "Flavobacterium", "Microbacterium", "Sphingobacterium"), ]


mphyseqbigeff$Genus <- as.factor(mphyseqbigeff$Genus)

mphyseqbigeff$Genus <- factor(mphyseqbigeff$Genus, levels=c("Quadrisphaera", "Kineococcus", "Klenkia", "Aureimonas", "Spirosoma", "Methylobacterium-Methylorubrum", "Deinococcus", "Hymenobacter", "Azorhizobium", "Actinomycetospora","Duganella", "Luteolibacter", "Pseudomonas", "Cellvibrio", "Stenotrophomonas", "Rickettsiella", "Achromobacter", "Flavobacterium", "Microbacterium", "Sphingobacterium"))



ggplot(mphyseqbigeff, aes(x=Genus,y=Abundance)) + geom_sina(aes(shape=category,colour=nematode_species),scale="width",alpha=0.5) + stat_summary(aes(group=category,colour=nematode_species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1,position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic"),legend.title=element_blank()) + scale_colour_brewer(palette="Set1") + xlab("Genus") + ylab("Transformed abundance")



ggsave("supplemental_figure_genus_elegans_inopinata.png",bg="white",height=8,width=15,units="in")






#######

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


setwd("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay")

#load ASV counts (see line 319 of workflow.sh), taxonomy table (line 333 of workflow.sh), and sample metadata
ps_object <- read_phyloseq(otu.file = "feature-table.csv", 
                    taxonomy.file = "taxonomy_revised_2.csv", 
                    metadata.file = "inopinata_and_elegans_metadata_2.csv", 
                    type = 'simple')

#load phylogenetic tree (line 326 of workflow.sh)
treefile <- read.tree("tree.tree")

#add tree to phyloseq object
ps.ng.tax <- merge_phyloseq(ps_object, treefile)

#remove mitochondrial and chloroplast reads

nochmi <- subset_taxa(ps.ng.tax, !Family %in% "Mitochondria" & !Order %in% "Chloroplast")

##reads after removing mito and chloro
#orgoreads <- as.data.frame(sample_sums(nochmi))
#orgoreads$sample_id <- rownames(orgoreads)
#orgoreads$num_paired_end_reads <- orgoreads[,1]
#
#orgoreads2 <- data.frame(sample_id = orgoreads$sample_id, num_paired_end_reads=orgoreads$num_paired_end_reads, step= "Remove organelles")
#

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
#[1] 150

#there are 150 ASVs in the controls that will ultimately be removed from the fig microbiome analysis.
#
nochmiprevdf = apply(X = otu_table(nochmi),
               MARGIN = ifelse(taxa_are_rows(controls), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
nochmiprevdf = data.frame(Prevalence = nochmiprevdf,
                    TotalAbundance = taxa_sums(nochmi),
                    tax_table(nochmi))

nrow(nochmiprevdf)
#[1] 58457
150/58457
#[1] 0.002565989

#these are the OTUs to remove
control_OTU$OTU

#remove control samples

no_controls <- prune_samples(!(sample_names(nochmi) %in% c("EXTRNEG","ExtrNeg1","EXTRPOS","GW10","GW11","GW31","GW9","PCRNeg1","PCRNeg2","PCRPos1","PCRPos2")), nochmi)

#remove OTU's found in controls and controls; remove control samples
ps3 <- subset_taxa(no_controls, !OTU %in% control_OTU$OTU)



#transform
ps_clr <- microbiome::transform(ps3, "clr")


#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 

#ADONIS test PERMANOVA
vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$category)

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999
#
#vegan::adonis2(formula = clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$category)
#          Df SumOfSqs      R2      F Pr(>F)
#Model      3   280318 0.11823 8.4918  0.001 ***
#Residual 190  2090656 0.88177
#Total    193  2370974 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

library(pairwiseAdonis)

pairwise.adonis(clr_dist_matrix, phyloseq::sample_data(ps_clr)$category)


#                                          pairs Df SumsOfSqs   F.Model
#1    dirksen_substrate vs woodruff_fig_interior  1 150852.51  8.867092
#2     dirksen_substrate vs woodruff_fig_surface  1 125736.14  7.127452
#3         dirksen_substrate vs samuel_substrate  1 133845.99  8.917579
#4 woodruff_fig_interior vs woodruff_fig_surface  1  11209.09  2.936922
#5     woodruff_fig_interior vs samuel_substrate  1  57919.05 12.859245
#6      woodruff_fig_surface vs samuel_substrate  1  42026.77  9.713602
#          R2 p.value p.adjusted sig
#1 0.08144878   0.001      0.006   *
#2 0.07047989   0.001      0.006   *
#3 0.06811598   0.001      0.006   *
#4 0.04140188   0.001      0.006   *
#5 0.11812727   0.001      0.006   *
#6 0.09741501   0.001      0.006   *



vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$nematode_species)

#vegan::adonis2(formula = clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$nematode_species)
#          Df SumOfSqs      R2      F Pr(>F)
#Model      1   135263 0.05705 11.616  0.001 ***
#Residual 192  2235711 0.94295
#Total    193  2370974 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

library(pairwiseAdonis)

pairwise.adonis(clr_dist_matrix, phyloseq::sample_data(ps_clr)$nematode_species)


vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$study + phyloseq::sample_data(ps_clr)$nematode_species + phyloseq::sample_data(ps_clr)$substrate_type)





vegan::adonis2(formula = clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$nematode_species + phyloseq::sample_data(ps_clr)$study + phyloseq::sample_data(ps_clr)$substrate_type)
#          Df SumOfSqs      R2      F Pr(>F)
#Model      6   326409 0.13767 4.9757  0.001 ***
#Residual 187  2044565 0.86233
#Total    193  2370974 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

vegan::adonis2(formula = clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$nematode_species + phyloseq::sample_data(ps_clr)$study + phyloseq::sample_data(ps_clr)$substrate_type,by="terms")
#                                                Df SumOfSqs      R2       F
#phyloseq::sample_data(ps_clr)$nematode_species   1   135263 0.05705 12.3714
#phyloseq::sample_data(ps_clr)$study              1   133846 0.05645 12.2418
#phyloseq::sample_data(ps_clr)$substrate_type     4    57300 0.02417  1.3102
#Residual                                       187  2044565 0.86233
#Total                                          193  2370974 1.00000
#                                               Pr(>F)
#phyloseq::sample_data(ps_clr)$nematode_species  0.001 ***
#phyloseq::sample_data(ps_clr)$study             0.001 ***
#phyloseq::sample_data(ps_clr)$substrate_type    0.100 .
#Residual
#Total
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


vegan::adonis2(formula = clr_dist_matrix ~  phyloseq::sample_data(ps_clr)$study + phyloseq::sample_data(ps_clr)$substrate_type + phyloseq::sample_data(ps_clr)$nematode_species,by="terms")




vegan::adonis2(formula = clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$study)
#          Df SumOfSqs     R2      F Pr(>F)
#Model      2   269109 0.1135 12.227  0.001 ***
#Residual 191  2101865 0.8865
#Total    193  2370974 1.0000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

vegan::adonis2(formula = clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$nematode_species)
#           Df SumOfSqs      R2      F Pr(>F)
#Model      1   135263 0.05705 11.616  0.001 ***
#Residual 192  2235711 0.94295
#Total    193  2370974 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

vegan::adonis2(formula = clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$substrate_type)

#          Df SumOfSqs      R2      F Pr(>F)
#Model      5   306750 0.12938 5.5875  0.001 ***
#Residual 188  2064224 0.87062
#Total    193  2370974 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pairwise.adonis(clr_dist_matrix, phyloseq::sample_data(ps_clr)$substrate_type)
> pairwise.adonis(clr_dist_matrix, phyloseq::sample_data(ps_clr)$substrate_type)
#                                pairs Df  SumsOfSqs    F.Model         R2
#1                    Compost vs Fruit  1 130132.592  8.3026874 0.06901498
#2           Compost vs fig suspension  1 151800.860  8.9329727 0.08353806
#3         Compost vs fig surface wash  1 126822.806  7.1918641 0.07250458
#4                   Compost vs Vector  1  22135.520  0.9273965 0.01428359
#5                     Compost vs Stem  1  30676.274  1.3063518 0.01940904
#6             Fruit vs fig suspension  1  57847.730 12.9657074 0.12841694
#7           Fruit vs fig surface wash  1  42309.819  9.9233069 0.10795202
#8                     Fruit vs Vector  1   9425.098  1.9970313 0.03566316
#9                       Fruit vs Stem  1  14981.870  3.0275872 0.05129106
#10 fig suspension vs fig surface wash  1  11209.088  2.9369222 0.04140188
#11           fig suspension vs Vector  1  14448.686  3.5077380 0.08062331
#12             fig suspension vs Stem  1  24413.779  5.4826948 0.11546722
#13         fig surface wash vs Vector  1  13054.764  3.6452308 0.09683115
#14           fig surface wash vs Stem  1  22042.401  5.5097396 0.13273366
#15                     Vector vs Stem  1   9564.972  1.5534943 0.16261006
#   p.value p.adjusted sig
#1    0.001      0.015   .
#2    0.001      0.015   .
#3    0.001      0.015   .
#4    0.552      1.000
#5    0.097      1.000
#6    0.001      0.015   .
#7    0.001      0.015   .
#8    0.044      0.660
#9    0.001      0.015   .
#10   0.001      0.015   .
#11   0.001      0.015   .
#12   0.001      0.015   .
#13   0.001      0.015   .
#14   0.001      0.015   .
#15   0.069      1.000



adiv_melt <- melt(adiv_actual, id.vars= c("sample_id","surface_or_interior"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic"))

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
