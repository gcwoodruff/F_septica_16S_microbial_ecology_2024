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

#remove mitochondrial and chloroplast reads

nochmi <- subset_taxa(ps.ng.tax, !Family %in% "Mitochondria" & !Order %in% "Chloroplast")

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


no_controls <- prune_samples(!(sample_names(nochmi) %in% c("EXTRNEG","ExtrNeg1","EXTRPOS","GW10","GW11","GW31","GW9","PCRNeg1","PCRNeg2","PCRPos1","PCRPos2")), nochmi)

#remove ASV's found in controls and controls; remove control samples
ps3 <- subset_taxa(no_controls, !OTU %in% control_OTU$OTU)


#tired of generating this phyloseq object
#saveRDS(ps3, "ps3_figs_only.rds")

ps3 <- readRDS("ps3_figs_only.rds")

#genera only for the moment
ps3genera <- tax_glom(ps3, "Genus", NArm = TRUE)


#https://chiliubio.github.io/microeco_tutorial/meconetcomp-package.html
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
#ps3genera@otu_table

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

# use clone to get a deep copy of soil_amp (R6 object)
tmp <- clone(mt)

# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
    #the pairwise correlation coefficients are in tmp$res_cor_p as a list! (see manual entry for trans_network())


# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5,COR_p_adjust = "fdr")
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_all_genus.gexf")
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
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_surface_genus.gexf")

# put the network into the list
mt_network$surface <- tmp
# select samples of "TW" group
tmp <- clone(mt)
tmp$sample_table %<>% subset(surface_or_interior == "interior")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_interior_genus.gexf")

mt_network$interior <- tmp



str(mt_network$all$res_cor_p)

all_cor_df <- as.data.frame(mt_network$all$res_cor_p$cor)

write.table(all_cor_df,"all_cor_df_genus.tsv",sep="\t")

all_p_df <- as.data.frame(mt_network$all$res_cor_p$p)

all_cor_df$OTU1 <- rownames(all_cor_df)
all_p_df$OTU1 <- rownames(all_p_df)

library(reshape2)

all_cor_df_melt <- reshape2::melt(all_cor_df)
all_p_df_melt <- reshape2::melt(all_p_df)

names(all_cor_df_melt)[names(all_cor_df_melt) == 'variable'] <- 'OTU2'
names(all_p_df_melt)[names(all_p_df_melt) == 'variable'] <- 'OTU2'

names(all_cor_df_melt)[names(all_cor_df_melt) == 'value'] <- 'cor'
names(all_p_df_melt)[names(all_p_df_melt) == 'value'] <- 'p'

all_cor_df_melt$OTU1_OTU2 <- paste(all_cor_df_melt$OTU1,all_cor_df_melt$OTU2)

all_p_df_melt$OTU1_OTU2 <- paste(all_p_df_melt$OTU1,all_p_df_melt$OTU2)

all_cor_df_melt_ii <- data.frame(OTU1_OTU2=all_cor_df_melt$OTU1_OTU2,cor=all_cor_df_melt$cor)
all_p_df_melt_ii <- data.frame(OTU1_OTU2=all_p_df_melt$OTU1_OTU2,p=all_p_df_melt$p)

p_cor_df_melt_merge <- merge(all_cor_df_melt_ii,all_p_df_melt_ii)


library(tidyr)

#split into two
p_cor_df_melt_merge_sep <- separate(p_cor_df_melt_merge, OTU1_OTU2, into = c("OTU1","OTU2"), sep = " (?=[^ ]+$)")

head(p_cor_df_melt_merge_sep)

p_cor_df_melt_merge_sep_sort <- p_cor_df_melt_merge_sep[order(p_cor_df_melt_merge_sep$cor, decreasing = TRUE),]  



library(ggplot2)
library(cowplot)
ggplot(p_cor_df_melt_merge_sep_sort, aes(x=OTU1, y=OTU2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal()

taxonomy_df = read.csv("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/05_qiime_feature-classifier_classify-consensus-vsearch/export/taxonomy_revised_fig_samples_only_5-2025.tsv", header=TRUE)


#okay, replace ASV id with Genus id
taxonomy_df_OTU_cor_OTU1 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU1),]

taxonomy_df_OTU_cor_OTU2 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU2),]


p_cor_df_melt_merge_sep_sort$Genus_1 <- "a"

p_cor_df_melt_merge_sep_sort$Genus_2 <- "b"
p_cor_df_melt_merge_sep_sort$OTU1 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU1)
p_cor_df_melt_merge_sep_sort$OTU2 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU2)

p_cor_df_melt_merge_sep_sort2 <- p_cor_df_melt_merge_sep_sort

head(unique(p_cor_df_melt_merge_sep_sort2$OTU1))

the_otu <- "0b7339d9f641ca23818474b57ca212d8"

the_otu_tax <- taxonomy_df_OTU_cor_OTU1[taxonomy_df_OTU_cor_OTU1$ASV == the_otu,]

the_otu_tax$Genus

p_cor_df_melt_merge_sep_sort2$OTU1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu]


p_cor_df_melt_merge_sep_sort2$Genus_1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu]


for(i in unique(p_cor_df_melt_merge_sep_sort2$OTU1)){
	the_otu <- i
	the_otu_tax <- taxonomy_df_OTU_cor_OTU1[taxonomy_df_OTU_cor_OTU1$ASV == the_otu,]
	p_cor_df_melt_merge_sep_sort2$Genus_1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu] <- the_otu_tax$Genus
	p_cor_df_melt_merge_sep_sort2$Genus_2[p_cor_df_melt_merge_sep_sort2$OTU2 == the_otu] <- the_otu_tax$Genus

}


ggplot(p_cor_df_melt_merge_sep_sort2, aes(x=Genus_1, y=Genus_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

#can i even try to do a upper triangle heat map...  https://albert-rapp.de/posts/ggplot2-tips/24_correlation_heat_map/24_correlation_heat_map.html

p_cor_df_melt_merge_sep_sort2$Genus_1 <- as.factor(p_cor_df_melt_merge_sep_sort2$Genus_1)
p_cor_df_melt_merge_sep_sort2$Genus_2 <- as.factor(p_cor_df_melt_merge_sep_sort2$Genus_2)


p_cor_df_melt_merge_sep_sort3 <- p_cor_df_melt_merge_sep_sort2 |> 
  mutate(
    lvl_a = as.numeric(Genus_1),
    lvl_b = as.numeric(Genus_2),
    cor = if_else(lvl_a < lvl_b, cor, NA)
  ) 



ggplot(p_cor_df_melt_merge_sep_sort3, aes(x=Genus_1, y=Genus_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation",na.value = 'white') + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

p_cor_df_melt_merge_sep_sort4 <- p_cor_df_melt_merge_sep_sort3[p_cor_df_melt_merge_sep_sort3$Genus_1 != "env.OPS_17",]

p_cor_df_melt_merge_sep_sort5 <- p_cor_df_melt_merge_sep_sort4[p_cor_df_melt_merge_sep_sort4$Genus_2 != "env.OPS_17",]


levels(p_cor_df_melt_merge_sep_sort5$Genus_1)[match("Chroococcidiopsis_SAG_2023",levels(p_cor_df_melt_merge_sep_sort5$Genus_1))] <- "Chroococcidiopsis"
levels(p_cor_df_melt_merge_sep_sort5$Genus_2)[match("Chroococcidiopsis_SAG_2023",levels(p_cor_df_melt_merge_sep_sort5$Genus_2))] <- "Chroococcidiopsis"


p_cor_df_melt_merge_sep_sort8 <- p_cor_df_melt_merge_sep_sort5[p_cor_df_melt_merge_sep_sort5$Genus_1 != "uncultured",]
p_cor_df_melt_merge_sep_sort9 <- p_cor_df_melt_merge_sep_sort8[p_cor_df_melt_merge_sep_sort8$Genus_2 != "uncultured",]


ggplot(p_cor_df_melt_merge_sep_sort9, aes(x=Genus_1, y=Genus_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation",na.value = 'white') + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

ggsave("all_genus_heatmap.pdf",useDingbats=FALSE,bg="white")
#omit redundant rows (ie, Genus 2:Genus 1 pair is the same as Genus 1: Genus 2 pair)
p_cor_df_melt_merge_sep_sort6 <- na.omit(p_cor_df_melt_merge_sep_sort5)
#get only significant correlations
p_cor_df_melt_merge_sep_sort7 <- p_cor_df_melt_merge_sep_sort6[p_cor_df_melt_merge_sep_sort6$p < 0.05,]


write.table(p_cor_df_melt_merge_sep_sort6,"all_cor_df_genus_melt.tsv",sep="\t")

nrow(p_cor_df_melt_merge_sep_sort6)
#[1] 2079

nrow(p_cor_df_melt_merge_sep_sort7)
#820 significant pairwise genus associations

nrow(p_cor_df_melt_merge_sep_sort7[p_cor_df_melt_merge_sep_sort7$cor > 0,])
#[1] 751


nrow(p_cor_df_melt_merge_sep_sort7[p_cor_df_melt_merge_sep_sort7$cor < 0,])
#69

head(p_cor_df_melt_merge_sep_sort6)

#get top ten
top10corgenus <- head(p_cor_df_melt_merge_sep_sort7,10)
#get bottom ten
bottom10corgenus <- tail(p_cor_df_melt_merge_sep_sort7,10)
#bind
corpairs <- rbind(top10corgenus,bottom10corgenus)

#get pair for plotting
corpairs$genpair <- paste(corpairs$Genus_1,corpairs$Genus_2,sep="/")



ggplot(corpairs, aes(x=cor,y=reorder(genpair, cor))) + geom_col(fill="lightblue") + theme_cowplot() + background_grid() + scale_x_continuous(limits=c(-1,1),breaks=c(seq(-1,1,by=0.1))) +xlab("Spearman correlation coefficient") + ylab("Genus pair")

ggsave("genus_correlation_top_bottom_ten.pdf",useDingbats=FALSE,bg="white")



#okay, just trying to gut check this by looking at some scatterplots

str(ps3genera)


ps3generaprevdf = apply(X = otu_table(ps3genera),
               MARGIN = ifelse(taxa_are_rows(ps3genera), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this df
ps3generaprevdf = data.frame(Prevalence = ps3generaprevdf,
                    TotalAbundance = taxa_sums(ps3genera),
                    tax_table(ps3genera))
	#okay, not what i want
head(otutab)
	#this looks like what i want but i need genera here!

otutab$OTU <- as.factor(rownames(otutab))
otutab$Genus <- "a"
taxonomy_df[taxonomy_df$ASV == "bc3676c17839094c4fa8e33268905268",]

otutab2 <- otutab

#add genera....
for(i in unique(otutab2$OTU)){
	the_otu <- i
	the_otu_tax <- taxonomy_df[taxonomy_df$ASV == the_otu,]
	otutab2$Genus[otutab2$OTU == the_otu] <- the_otu_tax$Genus
}

otutab2_melt <- reshape2::melt(otutab2,id.vars="Genus")



#wait, let's transform first?

ps_clr <- microbiome::transform(ps3genera, "clr")

otutab3 <- as.data.frame(ps_clr@otu_table)

otutab3$OTU <- as.factor(rownames(otutab))
otutab3$Genus <- "a"
taxonomy_df[taxonomy_df$ASV == "bc3676c17839094c4fa8e33268905268",]

otutab4 <- otutab3

#add genera....
for(i in unique(otutab4$OTU)){
	the_otu <- i
	the_otu_tax <- taxonomy_df[taxonomy_df$ASV == the_otu,]
	otutab4$Genus[otutab4$OTU == the_otu] <- the_otu_tax$Genus
}

as.data.frame(table(otutab4$Genus))
	#okay, get rid of "uncultured"

otutab5 <- otutab4[otutab4$Genus != "uncultured",]
otutab5 <- subset(otutab5, select = -c(OTU))
otutab5_melt <- reshape2::melt(otutab5,id.vars="Genus")

MethylMeth <- otutab5_melt[otutab5_melt$Genus=="Methylobacterium-Methylorubrum",]
Sphingo <- otutab5_melt[otutab5_melt$Genus=="Sphingomonas",]

nrow(MethylMeth)
nrow(Sphingo)

MethylMeth <- subset(MethylMeth, select = -c(Genus))
Sphingo <- subset(Sphingo, select = -c(Genus))

MeSpmerge <- merge(MethylMeth,Sphingo,by="variable")

names(MeSpmerge)[names(MeSpmerge) == 'value.x'] <- 'Methylobacterium_Methylorubrum'
names(MeSpmerge)[names(MeSpmerge) == 'value.y'] <- 'Sphingomonas'


ggplot(MeSpmerge, aes(x=Methylobacterium_Methylorubrum,y=Sphingomonas)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Methylobacterium/Methylorubrum transformed abundance") + ylab("Sphingomonas transformed abundance") +scale_x_continuous(limits=c(-1,10),breaks=c(seq(-1,10,1))) +scale_y_continuous(limits=c(-1,10),breaks=c(seq(-1,10,1)))


ggsave("Methylobacterium_Methylorubrum_Sphingomonas_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")


#one that isn't correlated
gyah <- p_cor_df_melt_merge_sep_sort5[p_cor_df_melt_merge_sep_sort5$Genus_1 == "Methylobacterium-Methylorubrum",]
gyah2 <- gyah[gyah$Genus_2 == "Serratia",]


Serra <- otutab5_melt[otutab5_melt$Genus=="Serratia",]
Serra <- subset(Serra, select = -c(Genus))

MeSamerge <- merge(MethylMeth,Serra,by="variable")


names(MeSamerge)[names(MeSamerge) == 'value.x'] <- 'Methylobacterium_Methylorubrum'
names(MeSamerge)[names(MeSamerge) == 'value.y'] <- 'Serratia'


ggplot(MeSamerge, aes(x=Methylobacterium_Methylorubrum,y=Serratia)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Methylobacterium/Methylorubrum transformed abundance") + ylab("Serratia transformed abundance") +scale_x_continuous(limits=c(-3,10),breaks=c(seq(-3,10,1))) +scale_y_continuous(limits=c(-3,10),breaks=c(seq(-3,10,1)))


ggsave("Methylobacterium_Methylorubrum_Serratia_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")

Ochrobactrum/Quadrisphaera

nrow(MethylMeth)
nrow(Serra)





Ochro <- otutab5_melt[otutab5_melt$Genus=="Ochrobactrum",]
Quadr <- otutab5_melt[otutab5_melt$Genus=="Quadrisphaera",]

nrow(Ochro)
nrow(Quadr)

Ochro <- subset(Ochro, select = -c(Genus))
Quadr <- subset(Quadr, select = -c(Genus))

OchQuamerge <- merge(Ochro,Quadr,by="variable")

names(OchQuamerge)[names(OchQuamerge) == 'value.x'] <- 'Ochrobactrum'
names(OchQuamerge)[names(OchQuamerge) == 'value.y'] <- 'Quadrisphaera'


ggplot(OchQuamerge, aes(x=Ochrobactrum,y=Quadrisphaera)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Ochrobactrum transformed abundance") + ylab("Quadrisphaera transformed abundance") +scale_x_continuous(limits=c(-3,10),breaks=c(seq(-3,10,1))) +scale_y_continuous(limits=c(-3,10),breaks=c(seq(-3,10,1)))



ggsave("Ochrobactrum_Quadrisphaera_transformed_abundance_scatterplot.pdf",useDingbats=FALSE,bg="white")



MethOchmerge <- merge(MethylMeth,Ochro,by="variable")


names(MethOchmerge)[names(MethOchmerge) == 'value.x'] <- 'Methylobacterium_Methylorubrum'
names(MethOchmerge)[names(MethOchmerge) == 'value.y'] <- 'Ochrobactrum'


ggplot(MethOchmerge, aes(x=Methylobacterium_Methylorubrum,y=Ochrobactrum)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Methylobacterium/Methylorubrum transformed abundance") + ylab("Ochrobactrum transformed abundance") +scale_x_continuous(limits=c(-3,10),breaks=c(seq(-3,10,1))) +scale_y_continuous(limits=c(-3,10),breaks=c(seq(-3,10,1)))


ggsave("Methylobacterium_Methylorubrum_Ochrobactrum_transformed_abundance_scatterplot.pdf",useDingbats=FALSE,bg="white")

#okay-- OTU/ASV level



taxtab <- as.data.frame(tax_table(ps3))

otutab <- as.data.frame(ps3@otu_table)

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

# use clone to get a deep copy of soil_amp (R6 object)
tmp <- clone(mt)

# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
    #the pairwise correlation coefficients are in tmp$res_cor_p as a list! (see manual entry for trans_network())


# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_all_ASV.gexf")
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
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_surface_ASV.gexf")

# put the network into the list
mt_network$surface <- tmp
# select samples of "TW" group
tmp <- clone(mt)
tmp$sample_table %<>% subset(surface_or_interior == "interior")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_interior_ASV.gexf")

mt_network$interior <- tmp



str(mt_network$all$res_cor_p)

all_cor_df <- as.data.frame(mt_network$all$res_cor_p$cor)

write.table(all_cor_df,"all_cor_df_ASV.tsv",sep="\t")

all_p_df <- as.data.frame(mt_network$all$res_cor_p$p)

all_cor_df$OTU1 <- rownames(all_cor_df)
all_p_df$OTU1 <- rownames(all_p_df)

library(reshape2)

all_cor_df_melt <- reshape2::melt(all_cor_df)
all_p_df_melt <- reshape2::melt(all_p_df)

names(all_cor_df_melt)[names(all_cor_df_melt) == 'variable'] <- 'OTU2'
names(all_p_df_melt)[names(all_p_df_melt) == 'variable'] <- 'OTU2'

names(all_cor_df_melt)[names(all_cor_df_melt) == 'value'] <- 'cor'
names(all_p_df_melt)[names(all_p_df_melt) == 'value'] <- 'p'

all_cor_df_melt$OTU1_OTU2 <- paste(all_cor_df_melt$OTU1,all_cor_df_melt$OTU2)

all_p_df_melt$OTU1_OTU2 <- paste(all_p_df_melt$OTU1,all_p_df_melt$OTU2)

all_cor_df_melt_ii <- data.frame(OTU1_OTU2=all_cor_df_melt$OTU1_OTU2,cor=all_cor_df_melt$cor)
all_p_df_melt_ii <- data.frame(OTU1_OTU2=all_p_df_melt$OTU1_OTU2,p=all_p_df_melt$p)

p_cor_df_melt_merge <- merge(all_cor_df_melt_ii,all_p_df_melt_ii)


library(tidyr)

#split into two
p_cor_df_melt_merge_sep <- separate(p_cor_df_melt_merge, OTU1_OTU2, into = c("OTU1","OTU2"), sep = " (?=[^ ]+$)")

head(p_cor_df_melt_merge_sep)

p_cor_df_melt_merge_sep_sort <- p_cor_df_melt_merge_sep[order(p_cor_df_melt_merge_sep$cor, decreasing = TRUE),]  



library(ggplot2)
library(cowplot)
ggplot(p_cor_df_melt_merge_sep_sort, aes(x=OTU1, y=OTU2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal()

taxonomy_df = read.csv("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/05_qiime_feature-classifier_classify-consensus-vsearch/export/taxonomy_revised_fig_samples_only_5-2025.tsv", header=TRUE)


#okay, replace ASV id with "Genus OTU N" id
taxonomy_df_OTU_cor_OTU1 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU1),]

taxonomy_df_OTU_cor_OTU2 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU2),]

taxonomy_df_OTU_cor_OTU1 <- data.frame(OTU =taxonomy_df_OTU_cor_OTU1$ASV,Genus=taxonomy_df_OTU_cor_OTU1$Genus)


taxonomy_df_OTU_cor_OTU1$OTU <- as.factor(taxonomy_df_OTU_cor_OTU1$OTU)

taxonomy_df_OTU_cor_OTU1$Genus <- as.factor(taxonomy_df_OTU_cor_OTU1$Genus)

levels(taxonomy_df_OTU_cor_OTU1$Genus)[match("",levels(taxonomy_df_OTU_cor_OTU1$Genus))] <- "Unknown"


levels(taxonomy_df_OTU_cor_OTU1$Genus)[match("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",levels(taxonomy_df_OTU_cor_OTU1$Genus))] <- "Allo.-Neo.-Para.-Rhizobium"


taxonomy_df_OTU_cor_OTU1_sort <- taxonomy_df_OTU_cor_OTU1[order(taxonomy_df_OTU_cor_OTU1$Genus),]

bigdf <- NULL
testsubset<- NULL



for(i in levels(taxonomy_df_OTU_cor_OTU1_sort$Genus)){
	testsubset <- taxonomy_df_OTU_cor_OTU1_sort[taxonomy_df_OTU_cor_OTU1_sort$Genus == i,]
	testsubset$otu_to_paste <- "ASV"
	testsubset$running_count <- seq(1:nrow(testsubset))
	testsubset$LABEL <- paste(testsubset$Genus,testsubset$otu_to_paste,testsubset$running_count)
	bigdf <- rbind(bigdf,testsubset)
}

#
#taxonomy_df_OTU_cor_OTU2 <- data.frame(OTU =taxonomy_df_OTU_cor_OTU2$ASV,Genus=taxonomy_df_OTU_cor_OTU2$Genus)
#
#taxonomy_df_OTU_cor_OTU2$OTU <- as.factor(taxonomy_df_OTU_cor_OTU2$OTU)
#
#taxonomy_df_OTU_cor_OTU2$Genus <- as.factor(taxonomy_df_OTU_cor_OTU2$Genus)
#
#levels(taxonomy_df_OTU_cor_OTU2$Genus)[match("",levels(taxonomy_df_OTU_cor_OTU2$Genus))] <- "Unknown"
#
#levels(taxonomy_df_OTU_cor_OTU2$Genus)[match("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",levels(taxonomy_df_OTU_cor_OTU2$Genus))] <- "Allo.-Neo.-Para.-Rhizobium"
#
#taxonomy_df_OTU_cor_OTU2_sort <- taxonomy_df_OTU_cor_OTU2[order(taxonomy_df_OTU_cor_OTU2$Genus),]
#
#bigdf2 <- NULL
#testsubset<- NULL
#
#
#for(i in levels(taxonomy_df_OTU_cor_OTU2_sort$Genus)){
	#testsubset <- taxonomy_df_OTU_cor_OTU2_sort[taxonomy_df_OTU_cor_OTU2_sort$Genus == i,]
	#testsubset$otu_to_paste <- "ASV"
	#testsubset$running_count <- seq(1:nrow(testsubset))
	#testsubset$LABEL <- paste(testsubset$Genus,testsubset$otu_to_paste,testsubset$running_count)
	#bigdf2 <- rbind(bigdf2,testsubset)
#}
#
#



p_cor_df_melt_merge_sep_sort$Label_1 <- "a"

p_cor_df_melt_merge_sep_sort$Label_2 <- "b"
p_cor_df_melt_merge_sep_sort$OTU1 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU1)
p_cor_df_melt_merge_sep_sort$OTU2 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU2)

p_cor_df_melt_merge_sep_sort2 <- p_cor_df_melt_merge_sep_sort

head(unique(p_cor_df_melt_merge_sep_sort2$OTU1))

the_otu <- "01ef1f83c86823200c8c0e7597cd6a3a"

the_otu_tax <- bigdf[bigdf$OTU == the_otu,]

the_otu_tax$Genus

p_cor_df_melt_merge_sep_sort2$OTU1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu]


p_cor_df_melt_merge_sep_sort2$Genus_1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu]


for(i in unique(p_cor_df_melt_merge_sep_sort2$OTU1)){
	the_otu <- i
	the_otu_tax <- bigdf[bigdf$OTU == the_otu,]
	p_cor_df_melt_merge_sep_sort2$Label_1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu] <- the_otu_tax$LABEL
	p_cor_df_melt_merge_sep_sort2$Label_2[p_cor_df_melt_merge_sep_sort2$OTU2 == the_otu] <- the_otu_tax$LABEL

}


ggplot(p_cor_df_melt_merge_sep_sort2, aes(x=Label_1, y=Label_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + xlab("ASV 1") + ylab("ASV 2")

#can i even try to do a upper triangle heat map...  https://albert-rapp.de/posts/ggplot2-tips/24_correlation_heat_map/24_correlation_heat_map.html

p_cor_df_melt_merge_sep_sort2$Label_1 <- as.factor(p_cor_df_melt_merge_sep_sort2$Label_1)
p_cor_df_melt_merge_sep_sort2$Label_2 <- as.factor(p_cor_df_melt_merge_sep_sort2$Label_2)


p_cor_df_melt_merge_sep_sort3 <- p_cor_df_melt_merge_sep_sort2 |> 
  mutate(
    lvl_a = as.numeric(Label_1),
    lvl_b = as.numeric(Label_2),
    cor = if_else(lvl_a < lvl_b, cor, NA)
  ) 



ggplot(p_cor_df_melt_merge_sep_sort3, aes(x=Label_1, y=Label_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation",na.value = 'white') + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))  + xlab("ASV 1") + ylab("ASV 2")


ggsave("all_ASV_heatmap.pdf",useDingbats=FALSE,bg="white",units="in",height=16,width=16)
#omit redundant rows (ie, thing 2:thing 1 pair is the same as thing 1: thing 2 pair)
p_cor_df_melt_merge_sep_sort6 <- na.omit(p_cor_df_melt_merge_sep_sort3)
#get only significant correlations
p_cor_df_melt_merge_sep_sort7 <- p_cor_df_melt_merge_sep_sort6[p_cor_df_melt_merge_sep_sort6$p < 0.05,]

#get top ten
top10corASV <- head(p_cor_df_melt_merge_sep_sort7,10)
#get bottom ten
bottom10corASV <- tail(p_cor_df_melt_merge_sep_sort7,10)
#bind
corpairs <- rbind(top10corASV,bottom10corASV)

#get pair for plotting
corpairs$genpair <- paste(corpairs$Label_1,corpairs$Label_2,sep="/")



ggplot(corpairs, aes(x=cor,y=reorder(genpair, cor))) + geom_col(fill="lightblue") + theme_cowplot() + background_grid() + scale_x_continuous(limits=c(-1,1),breaks=c(seq(-1,1,by=0.1))) +xlab("Spearman correlation coefficient") + ylab("ASV pair")

ggsave("ASV_correlation_top_bottom_ten.pdf",useDingbats=FALSE,bg="white")



#okay, just trying to gut check this by looking at some scatterplots

	#okay, not what i want
head(otutab)

#wait, let's transform first?

ps_clr <- microbiome::transform(ps3, "clr")

otutab3 <- as.data.frame(ps_clr@otu_table)



otutab3$OTU <- as.factor(rownames(otutab))
otutab3$Label <- "a"
otutab4 <- otutab3

#add labels....
for(i in unique(otutab4$OTU)){
	the_otu <- i
	the_otu_tax <- bigdf[bigdf$ASV == the_otu,]
	otutab4$Label[otutab4$OTU == the_otu] <- the_otu_tax$LABEL
}

as.data.frame(table(otutab4$Genus))
	#okay, get rid of "uncultured"

otutab5 <- otutab4[otutab4$Genus != "uncultured",]
otutab5 <- subset(otutab5, select = -c(OTU))
otutab5_melt <- reshape2::melt(otutab5,id.vars="Genus")

MethylMeth <- otutab5_melt[otutab5_melt$Genus=="Methylobacterium-Methylorubrum",]
Sphingo <- otutab5_melt[otutab5_melt$Genus=="Sphingomonas",]

nrow(MethylMeth)
nrow(Sphingo)

MethylMeth <- subset(MethylMeth, select = -c(Genus))
Sphingo <- subset(Sphingo, select = -c(Genus))

MeSpmerge <- merge(MethylMeth,Sphingo,by="variable")

names(MeSpmerge)[names(MeSpmerge) == 'value.x'] <- 'Methylobacterium_Methylorubrum'
names(MeSpmerge)[names(MeSpmerge) == 'value.y'] <- 'Sphingomonas'


ggplot(MeSpmerge, aes(x=Methylobacterium_Methylorubrum,y=Sphingomonas)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Methylobacterium/Methylorubrum transformed abundance") + ylab("Sphingomonas transformed abundance") +scale_x_continuous(limits=c(-1,10),breaks=c(seq(-1,10,1))) +scale_y_continuous(limits=c(-1,10),breaks=c(seq(-1,10,1)))


ggsave("Methylobacterium_Methylorubrum_Sphingomonas_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")


#one that isn't correlated
gyah <- p_cor_df_melt_merge_sep_sort5[p_cor_df_melt_merge_sep_sort5$Genus_1 == "Methylobacterium-Methylorubrum",]
gyah2 <- gyah[gyah$Genus_2 == "Serratia",]


Serra <- otutab5_melt[otutab5_melt$Genus=="Serratia",]
Serra <- subset(Serra, select = -c(Genus))

MeSamerge <- merge(MethylMeth,Serra,by="variable")


names(MeSamerge)[names(MeSamerge) == 'value.x'] <- 'Methylobacterium_Methylorubrum'
names(MeSamerge)[names(MeSamerge) == 'value.y'] <- 'Serratia'


ggplot(MeSamerge, aes(x=Methylobacterium_Methylorubrum,y=Serratia)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Methylobacterium/Methylorubrum transformed abundance") + ylab("Serratia transformed abundance") +scale_x_continuous(limits=c(-3,10),breaks=c(seq(-3,10,1))) +scale_y_continuous(limits=c(-3,10),breaks=c(seq(-3,10,1)))


ggsave("Methylobacterium_Methylorubrum_Serratia_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")

Ochrobactrum/Quadrisphaera

nrow(MethylMeth)
nrow(Serra)

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
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_surface.gexf")

# put the network into the list
mt_network$surface <- tmp
# select samples of "TW" group
tmp <- clone(mt)
tmp$sample_table %<>% subset(surface_or_interior == "interior")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_interior.gexf")

mt_network$interior <- tmp



str(mt_network$all$res_cor_p)

all_cor_df <- as.data.frame(mt_network$all$res_cor_p$cor)

write.table(all_cor_df,"all_cor_df.tsv",sep="\t")

all_p_df <- as.data.frame(mt_network$all$res_cor_p$p)

all_cor_df$OTU1 <- rownames(all_cor_df)
all_p_df$OTU1 <- rownames(all_p_df)

library(reshape2)

all_cor_df_melt <- reshape2::melt(all_cor_df)
all_p_df_melt <- reshape2::melt(all_p_df)

names(all_cor_df_melt)[names(all_cor_df_melt) == 'variable'] <- 'OTU2'
names(all_p_df_melt)[names(all_p_df_melt) == 'variable'] <- 'OTU2'

names(all_cor_df_melt)[names(all_cor_df_melt) == 'value'] <- 'cor'
names(all_p_df_melt)[names(all_p_df_melt) == 'value'] <- 'p'

all_cor_df_melt$OTU1_OTU2 <- paste(all_cor_df_melt$OTU1,all_cor_df_melt$OTU2)

all_p_df_melt$OTU1_OTU2 <- paste(all_p_df_melt$OTU1,all_p_df_melt$OTU2)

all_cor_df_melt_ii <- data.frame(OTU1_OTU2=all_cor_df_melt$OTU1_OTU2,cor=all_cor_df_melt$cor)
all_p_df_melt_ii <- data.frame(OTU1_OTU2=all_p_df_melt$OTU1_OTU2,p=all_p_df_melt$p)

p_cor_df_melt_merge <- merge(all_cor_df_melt_ii,all_p_df_melt_ii)


library(tidyr)

#split into two
p_cor_df_melt_merge_sep <- separate(p_cor_df_melt_merge, OTU1_OTU2, into = c("OTU1","OTU2"), sep = " (?=[^ ]+$)")

head(p_cor_df_melt_merge_sep)

p_cor_df_melt_merge_sep_sort <- p_cor_df_melt_merge_sep[order(p_cor_df_melt_merge_sep$cor, decreasing = TRUE),]  



library(ggplot2)
library(cowplot)
ggplot(p_cor_df_melt_merge_sep_sort, aes(x=OTU1, y=OTU2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal()

taxonomy_df = read.csv("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/05_qiime_feature-classifier_classify-consensus-vsearch/export/taxonomy_revised_fig_samples_only_5-2025.tsv", header=TRUE)


#okay, replace ASV id with Genus id
taxonomy_df_OTU_cor_OTU1 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU1),]

taxonomy_df_OTU_cor_OTU2 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU2),]


p_cor_df_melt_merge_sep_sort$Genus_1 <- "a"

p_cor_df_melt_merge_sep_sort$Genus_2 <- "b"
p_cor_df_melt_merge_sep_sort$OTU1 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU1)
p_cor_df_melt_merge_sep_sort$OTU2 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU2)

p_cor_df_melt_merge_sep_sort2 <- p_cor_df_melt_merge_sep_sort

head(unique(p_cor_df_melt_merge_sep_sort2$OTU1))

the_otu <- "0b7339d9f641ca23818474b57ca212d8"

the_otu_tax <- taxonomy_df_OTU_cor_OTU1[taxonomy_df_OTU_cor_OTU1$ASV == the_otu,]

the_otu_tax$Genus

p_cor_df_melt_merge_sep_sort2$OTU1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu]


p_cor_df_melt_merge_sep_sort2$Genus_1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu]


for(i in unique(p_cor_df_melt_merge_sep_sort2$OTU1)){
	the_otu <- i
	the_otu_tax <- taxonomy_df_OTU_cor_OTU1[taxonomy_df_OTU_cor_OTU1$ASV == the_otu,]
	p_cor_df_melt_merge_sep_sort2$Genus_1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu] <- the_otu_tax$Genus
	p_cor_df_melt_merge_sep_sort2$Genus_2[p_cor_df_melt_merge_sep_sort2$OTU2 == the_otu] <- the_otu_tax$Genus

}


ggplot(p_cor_df_melt_merge_sep_sort2, aes(x=Genus_1, y=Genus_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

#can i even try to do a upper triangle heat map...  https://albert-rapp.de/posts/ggplot2-tips/24_correlation_heat_map/24_correlation_heat_map.html

p_cor_df_melt_merge_sep_sort2$Genus_1 <- as.factor(p_cor_df_melt_merge_sep_sort2$Genus_1)
p_cor_df_melt_merge_sep_sort2$Genus_2 <- as.factor(p_cor_df_melt_merge_sep_sort2$Genus_2)


p_cor_df_melt_merge_sep_sort3 <- p_cor_df_melt_merge_sep_sort2 |> 
  mutate(
    lvl_a = as.numeric(Genus_1),
    lvl_b = as.numeric(Genus_2),
    cor = if_else(lvl_a < lvl_b, cor, NA)
  ) 



ggplot(p_cor_df_melt_merge_sep_sort3, aes(x=Genus_1, y=Genus_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation",na.value = 'white') + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

p_cor_df_melt_merge_sep_sort4 <- p_cor_df_melt_merge_sep_sort3[p_cor_df_melt_merge_sep_sort3$Genus_1 != "env.OPS_17",]

p_cor_df_melt_merge_sep_sort5 <- p_cor_df_melt_merge_sep_sort4[p_cor_df_melt_merge_sep_sort4$Genus_2 != "env.OPS_17",]


levels(p_cor_df_melt_merge_sep_sort5$Genus_1)[match("Chroococcidiopsis_SAG_2023",levels(p_cor_df_melt_merge_sep_sort5$Genus_1))] <- "Chroococcidiopsis"
levels(p_cor_df_melt_merge_sep_sort5$Genus_2)[match("Chroococcidiopsis_SAG_2023",levels(p_cor_df_melt_merge_sep_sort5$Genus_2))] <- "Chroococcidiopsis"


p_cor_df_melt_merge_sep_sort8 <- p_cor_df_melt_merge_sep_sort5[p_cor_df_melt_merge_sep_sort5$Genus_1 != "uncultured",]
p_cor_df_melt_merge_sep_sort9 <- p_cor_df_melt_merge_sep_sort8[p_cor_df_melt_merge_sep_sort8$Genus_2 != "uncultured",]


ggplot(p_cor_df_melt_merge_sep_sort9, aes(x=Genus_1, y=Genus_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation",na.value = 'white') + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

ggsave("all_genus_heatmap.pdf",useDingbats=FALSE,bg="white")
#omit redundant rows (ie, Genus 2:Genus 1 pair is the same as Genus 1: Genus 2 pair)
p_cor_df_melt_merge_sep_sort6 <- na.omit(p_cor_df_melt_merge_sep_sort5)
#get only significant correlations
p_cor_df_melt_merge_sep_sort7 <- p_cor_df_melt_merge_sep_sort6[p_cor_df_melt_merge_sep_sort6$p < 0.05,]

#get top ten
top10corgenus <- head(p_cor_df_melt_merge_sep_sort7,10)
#get bottom ten
bottom10corgenus <- tail(p_cor_df_melt_merge_sep_sort7,10)
#bind
corpairs <- rbind(top10corgenus,bottom10corgenus)

#get pair for plotting
corpairs$genpair <- paste(corpairs$Genus_1,corpairs$Genus_2,sep="/")



ggplot(corpairs, aes(x=cor,y=reorder(genpair, cor))) + geom_col(fill="lightblue") + theme_cowplot() + background_grid() + scale_x_continuous(limits=c(-1,1),breaks=c(seq(-1,1,by=0.1))) +xlab("Spearman correlation coefficient") + ylab("Genus pair")

ggsave("genus_correlation_top_bottom_ten.pdf",useDingbats=FALSE,bg="white")



#okay, just trying to gut check this by looking at some scatterplots

head(otutab)

#wait, let's transform first?

ps_clr <- microbiome::transform(ps3, "clr")

otutab3 <- as.data.frame(ps_clr@otu_table)

otutab3$OTU <- as.factor(rownames(otutab))
otutab3$LABEL <- "a"
otutab4 <- otutab3[otutab3$OTU %in% unique(p_cor_df_melt_merge_sep_sort$OTU1),]

#add labels
for(i in unique(otutab4$OTU)){
	the_otu <- i
	the_otu_tax <- bigdf[bigdf$OTU == the_otu,]
	otutab4$LABEL[otutab4$OTU == the_otu] <- the_otu_tax$LABEL
}

otutab5 <- subset(otutab4, select = -c(OTU))


otutab5_melt <- reshape2::melt(otutab5,id.vars="LABEL")

#Kineococcus ASV 1/Quadrisphaera ASV 1

Kineo1 <- otutab5_melt[otutab5_melt$LABEL=="Kineococcus ASV 1",]
Quadri1 <- otutab5_melt[otutab5_melt$LABEL=="Quadrisphaera ASV 1",]

nrow(Kineo1)
nrow(Quadri1)

Kineo1 <- subset(Kineo1, select = -c(LABEL))
Quadri1 <- subset(Quadri1, select = -c(LABEL))

Kineo1Quadri1merge <- merge(Kineo1,Quadri1,by="variable")

names(Kineo1Quadri1merge)[names(Kineo1Quadri1merge) == 'value.x'] <- 'Kineococcus_ASV_1'
names(Kineo1Quadri1merge)[names(Kineo1Quadri1merge) == 'value.y'] <- 'Quadrisphaera_ASV_1'


ggplot(Kineo1Quadri1merge, aes(x=Kineococcus_ASV_1,y=Quadrisphaera_ASV_1)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Kineococcus ASV 1 transformed abundance") + ylab("Quadrisphaera ASV 1 transformed abundance") +scale_x_continuous(limits=c(-3,10),breaks=c(seq(-3,10,1))) +scale_y_continuous(limits=c(-3,10),breaks=c(seq(-3,10,1)))


ggsave("Kineococcus_ASV_1_Quadrisphaera_ASV_1_transformed_abundance_scatterplot.pdf",useDingbats=FALSE,bg="white")

#Allo.−Neo.−Para.−Rhizobium ASV 5/Methylobacterium−Methylorubrum ASV 2

Allo5 <- otutab5_melt[otutab5_melt$LABEL=="Allo.-Neo.-Para.-Rhizobium ASV 5",]
Methyl2 <- otutab5_melt[otutab5_melt$LABEL=="Methylobacterium-Methylorubrum ASV 2",]

nrow(Allo5)
nrow(Methyl2)

Allo5 <- subset(Allo5, select = -c(LABEL))
Methyl2 <- subset(Methyl2, select = -c(LABEL))

Allo5Methyl2merge <- merge(Allo5,Methyl2,by="variable")

names(Allo5Methyl2merge)[names(Allo5Methyl2merge) == 'value.x'] <- "Allo_Neo_Para_Rhizobium_ASV_5"
names(Allo5Methyl2merge)[names(Allo5Methyl2merge) == 'value.y'] <- 'Methylobacterium_Methylorubrum_ASV_2'


ggplot(Allo5Methyl2merge, aes(x=Allo_Neo_Para_Rhizobium_ASV_5,y=Methylobacterium_Methylorubrum_ASV_2)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Allo.-Neo.-Para.-Rhizobium ASV 5 transformed abundance") + ylab("Methylobacterium-Methylorubrum ASV 2 transformed abundance") +scale_x_continuous(limits=c(-3,11),breaks=c(seq(-3,11,1))) +scale_y_continuous(limits=c(-3,11),breaks=c(seq(-3,11,1)))



ggsave("Allo_Neo_Para_Rhizobium_ASV_5_Methylobacterium_Methylorubrum_ASV_2_transformed_abundance_scatterplot.pdf",useDingbats=FALSE,bg="white")



#Family level


ps3family <- tax_glom(ps3, "Family", NArm = TRUE)


taxtab <- as.data.frame(tax_table(ps3family))

otutab <- as.data.frame(ps3family@otu_table)

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

# use clone to get a deep copy of soil_amp (R6 object)
tmp <- clone(mt)

# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
    #the pairwise correlation coefficients are in tmp$res_cor_p as a list! (see manual entry for trans_network())


# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_all_family.gexf")
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
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_surface_family.gexf")

# put the network into the list
mt_network$surface <- tmp
# select samples of "TW" group
tmp <- clone(mt)
tmp$sample_table %<>% subset(surface_or_interior == "interior")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_interior_family.gexf")

mt_network$interior <- tmp



str(mt_network$all$res_cor_p)

all_cor_df <- as.data.frame(mt_network$all$res_cor_p$cor)

write.table(all_cor_df,"all_cor_df_family.tsv",sep="\t")

all_p_df <- as.data.frame(mt_network$all$res_cor_p$p)

all_cor_df$OTU1 <- rownames(all_cor_df)
all_p_df$OTU1 <- rownames(all_p_df)

library(reshape2)

all_cor_df_melt <- reshape2::melt(all_cor_df)
all_p_df_melt <- reshape2::melt(all_p_df)

names(all_cor_df_melt)[names(all_cor_df_melt) == 'variable'] <- 'OTU2'
names(all_p_df_melt)[names(all_p_df_melt) == 'variable'] <- 'OTU2'

names(all_cor_df_melt)[names(all_cor_df_melt) == 'value'] <- 'cor'
names(all_p_df_melt)[names(all_p_df_melt) == 'value'] <- 'p'

all_cor_df_melt$OTU1_OTU2 <- paste(all_cor_df_melt$OTU1,all_cor_df_melt$OTU2)

all_p_df_melt$OTU1_OTU2 <- paste(all_p_df_melt$OTU1,all_p_df_melt$OTU2)

all_cor_df_melt_ii <- data.frame(OTU1_OTU2=all_cor_df_melt$OTU1_OTU2,cor=all_cor_df_melt$cor)
all_p_df_melt_ii <- data.frame(OTU1_OTU2=all_p_df_melt$OTU1_OTU2,p=all_p_df_melt$p)

p_cor_df_melt_merge <- merge(all_cor_df_melt_ii,all_p_df_melt_ii)


library(tidyr)

#split into two
p_cor_df_melt_merge_sep <- separate(p_cor_df_melt_merge, OTU1_OTU2, into = c("OTU1","OTU2"), sep = " (?=[^ ]+$)")

head(p_cor_df_melt_merge_sep)

p_cor_df_melt_merge_sep_sort <- p_cor_df_melt_merge_sep[order(p_cor_df_melt_merge_sep$cor, decreasing = TRUE),]  



library(ggplot2)
library(cowplot)
ggplot(p_cor_df_melt_merge_sep_sort, aes(x=OTU1, y=OTU2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal()

taxonomy_df = read.csv("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/05_qiime_feature-classifier_classify-consensus-vsearch/export/taxonomy_revised_fig_samples_only_5-2025.tsv", header=TRUE)


#okay, replace ASV id with Genus id
taxonomy_df_OTU_cor_OTU1 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU1),]

taxonomy_df_OTU_cor_OTU2 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU2),]


p_cor_df_melt_merge_sep_sort$Family_1 <- "a"

p_cor_df_melt_merge_sep_sort$Family_2 <- "b"
p_cor_df_melt_merge_sep_sort$OTU1 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU1)
p_cor_df_melt_merge_sep_sort$OTU2 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU2)

p_cor_df_melt_merge_sep_sort2 <- p_cor_df_melt_merge_sep_sort


for(i in unique(p_cor_df_melt_merge_sep_sort2$OTU1)){
	the_otu <- i
	the_otu_tax <- taxonomy_df_OTU_cor_OTU1[taxonomy_df_OTU_cor_OTU1$ASV == the_otu,]
	p_cor_df_melt_merge_sep_sort2$Family_1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu] <- the_otu_tax$Family
	p_cor_df_melt_merge_sep_sort2$Family_2[p_cor_df_melt_merge_sep_sort2$OTU2 == the_otu] <- the_otu_tax$Family

}


ggplot(p_cor_df_melt_merge_sep_sort2, aes(x=Family_1, y=Family_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

#can i even try to do a upper triangle heat map...  https://albert-rapp.de/posts/ggplot2-tips/24_correlation_heat_map/24_correlation_heat_map.html

p_cor_df_melt_merge_sep_sort2$Family_1 <- as.factor(p_cor_df_melt_merge_sep_sort2$Family_1)
p_cor_df_melt_merge_sep_sort2$Family_2 <- as.factor(p_cor_df_melt_merge_sep_sort2$Family_2)


p_cor_df_melt_merge_sep_sort3 <- p_cor_df_melt_merge_sep_sort2 |> 
  mutate(
    lvl_a = as.numeric(Family_1),
    lvl_b = as.numeric(Family_2),
    cor = if_else(lvl_a < lvl_b, cor, NA)
  ) 



ggplot(p_cor_df_melt_merge_sep_sort3, aes(x=Family_1, y=Family_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation",na.value = 'white') + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

p_cor_df_melt_merge_sep_sort4 <- p_cor_df_melt_merge_sep_sort3[p_cor_df_melt_merge_sep_sort3$Family_1 != "env.OPS_17",]

p_cor_df_melt_merge_sep_sort5 <- p_cor_df_melt_merge_sep_sort4[p_cor_df_melt_merge_sep_sort4$Family_2 != "env.OPS_17",]


ggplot(p_cor_df_melt_merge_sep_sort5, aes(x=Family_1, y=Family_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation",na.value = 'white') + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

ggsave("all_family_heatmap.pdf",useDingbats=FALSE,bg="white")
#omit redundant rows (ie, Genus 2:Genus 1 pair is the same as Genus 1: Genus 2 pair)
p_cor_df_melt_merge_sep_sort6 <- na.omit(p_cor_df_melt_merge_sep_sort5)
#get only significant correlations
p_cor_df_melt_merge_sep_sort7 <- p_cor_df_melt_merge_sep_sort6[p_cor_df_melt_merge_sep_sort6$p < 0.05,]

#get top ten
top10corgenus <- head(p_cor_df_melt_merge_sep_sort7,10)
#get bottom ten
bottom10corgenus <- tail(p_cor_df_melt_merge_sep_sort7,10)
#bind
corpairs <- rbind(top10corgenus,bottom10corgenus)

#get pair for plotting
corpairs$genpair <- paste(corpairs$Family_1,corpairs$Family_2,sep="/")



ggplot(corpairs, aes(x=cor,y=reorder(genpair, cor))) + geom_col(fill="lightblue") + theme_cowplot() + background_grid() + scale_x_continuous(limits=c(-1,1),breaks=c(seq(-1,1,by=0.1))) +xlab("Spearman correlation coefficient") + ylab("Family pair")

ggsave("family_correlation_top_bottom_ten.pdf",useDingbats=FALSE,bg="white")



#okay, just trying to gut check this by looking at some scatterplots


#wait, let's transform first?

ps_clr <- microbiome::transform(ps3family, "clr")

otutab3 <- as.data.frame(ps_clr@otu_table)

otutab3$OTU <- as.factor(rownames(otutab))
otutab3$Family <- "a"

otutab4 <- otutab3

#add families....
for(i in unique(otutab4$OTU)){
	the_otu <- i
	the_otu_tax <- taxonomy_df[taxonomy_df$ASV == the_otu,]
	otutab4$Family[otutab4$OTU == the_otu] <- the_otu_tax$Family
}

otutab5 <- subset(otutab4, select = -c(OTU))
otutab5_melt <- reshape2::melt(otutab5,id.vars="Family")

#top Beijerinckiaceae/Kineosporiaceae
#bottom Enterobacteriaceae/Kineosporiaceae



Beij <- otutab5_melt[otutab5_melt$Family=="Beijerinckiaceae",]
Kineo <- otutab5_melt[otutab5_melt$Family=="Kineosporiaceae",]

nrow(Beij)
nrow(Kineo)

Beij <- subset(Beij, select = -c(Family))
Kineo <- subset(Kineo, select = -c(Family))

BeijKineo <- merge(Beij,Kineo,by="variable")

names(BeijKineo)[names(BeijKineo) == 'value.x'] <- 'Beijerinckiaceae'
names(BeijKineo)[names(BeijKineo) == 'value.y'] <- 'Kineosporiaceae'


ggplot(BeijKineo, aes(x=Beijerinckiaceae,y=Kineosporiaceae)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Beijerinckiaceae transformed abundance") + ylab("Kineosporiaceae transformed abundance") +scale_x_continuous(limits=c(-2,10),breaks=c(seq(-2,10,1))) +scale_y_continuous(limits=c(-2,10),breaks=c(seq(-2,10,1)))


ggsave("Beijerinckiaceae_Kineosporiaceae_FAMILY_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")




Entero <- otutab5_melt[otutab5_melt$Family=="Enterobacteriaceae",]

nrow(Entero)
nrow(Kineo)

Entero <- subset(Entero, select = -c(Family))

EnteroKineo <- merge(Entero,Kineo,by="variable")

names(EnteroKineo)[names(EnteroKineo) == 'value.x'] <- 'Enterobacteriaceae'
names(EnteroKineo)[names(EnteroKineo) == 'value.y'] <- 'Kineosporiaceae'


ggplot(EnteroKineo, aes(x=Enterobacteriaceae,y=Kineosporiaceae)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Enterobacteriaceae transformed abundance") + ylab("Kineosporiaceae transformed abundance") +scale_x_continuous(limits=c(-2,12),breaks=c(seq(-2,12,1))) +scale_y_continuous(limits=c(-2,12),breaks=c(seq(-2,12,1)))


ggsave("Enterobacteriaceae_Kineosporiaceae_FAMILY_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")




#order level



ps3order <- tax_glom(ps3, "Order", NArm = TRUE)


taxtab <- as.data.frame(tax_table(ps3order))

otutab <- as.data.frame(ps3order@otu_table)

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

# use clone to get a deep copy of soil_amp (R6 object)
tmp <- clone(mt)

# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
    #the pairwise correlation coefficients are in tmp$res_cor_p as a list! (see manual entry for trans_network())


# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_all_order.gexf")
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
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_surface_order.gexf")

# put the network into the list
mt_network$surface <- tmp
# select samples of "TW" group
tmp <- clone(mt)
tmp$sample_table %<>% subset(surface_or_interior == "interior")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_interior_order.gexf")

mt_network$interior <- tmp



str(mt_network$all$res_cor_p)

all_cor_df <- as.data.frame(mt_network$all$res_cor_p$cor)

write.table(all_cor_df,"all_cor_df_order.tsv",sep="\t")

all_p_df <- as.data.frame(mt_network$all$res_cor_p$p)

all_cor_df$OTU1 <- rownames(all_cor_df)
all_p_df$OTU1 <- rownames(all_p_df)

library(reshape2)

all_cor_df_melt <- reshape2::melt(all_cor_df)
all_p_df_melt <- reshape2::melt(all_p_df)

names(all_cor_df_melt)[names(all_cor_df_melt) == 'variable'] <- 'OTU2'
names(all_p_df_melt)[names(all_p_df_melt) == 'variable'] <- 'OTU2'

names(all_cor_df_melt)[names(all_cor_df_melt) == 'value'] <- 'cor'
names(all_p_df_melt)[names(all_p_df_melt) == 'value'] <- 'p'

all_cor_df_melt$OTU1_OTU2 <- paste(all_cor_df_melt$OTU1,all_cor_df_melt$OTU2)

all_p_df_melt$OTU1_OTU2 <- paste(all_p_df_melt$OTU1,all_p_df_melt$OTU2)

all_cor_df_melt_ii <- data.frame(OTU1_OTU2=all_cor_df_melt$OTU1_OTU2,cor=all_cor_df_melt$cor)
all_p_df_melt_ii <- data.frame(OTU1_OTU2=all_p_df_melt$OTU1_OTU2,p=all_p_df_melt$p)

p_cor_df_melt_merge <- merge(all_cor_df_melt_ii,all_p_df_melt_ii)


library(tidyr)

#split into two
p_cor_df_melt_merge_sep <- separate(p_cor_df_melt_merge, OTU1_OTU2, into = c("OTU1","OTU2"), sep = " (?=[^ ]+$)")

head(p_cor_df_melt_merge_sep)

p_cor_df_melt_merge_sep_sort <- p_cor_df_melt_merge_sep[order(p_cor_df_melt_merge_sep$cor, decreasing = TRUE),]  



library(ggplot2)
library(cowplot)
ggplot(p_cor_df_melt_merge_sep_sort, aes(x=OTU1, y=OTU2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal()

taxonomy_df = read.csv("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/05_qiime_feature-classifier_classify-consensus-vsearch/export/taxonomy_revised_fig_samples_only_5-2025.tsv", header=TRUE)


#okay, replace ASV id with Order id
taxonomy_df_OTU_cor_OTU1 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU1),]

taxonomy_df_OTU_cor_OTU2 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU2),]


p_cor_df_melt_merge_sep_sort$Order_1 <- "a"

p_cor_df_melt_merge_sep_sort$Order_2 <- "b"
p_cor_df_melt_merge_sep_sort$OTU1 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU1)
p_cor_df_melt_merge_sep_sort$OTU2 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU2)

p_cor_df_melt_merge_sep_sort2 <- p_cor_df_melt_merge_sep_sort


for(i in unique(p_cor_df_melt_merge_sep_sort2$OTU1)){
	the_otu <- i
	the_otu_tax <- taxonomy_df_OTU_cor_OTU1[taxonomy_df_OTU_cor_OTU1$ASV == the_otu,]
	p_cor_df_melt_merge_sep_sort2$Order_1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu] <- the_otu_tax$Order
	p_cor_df_melt_merge_sep_sort2$Order_2[p_cor_df_melt_merge_sep_sort2$OTU2 == the_otu] <- the_otu_tax$Order

}


ggplot(p_cor_df_melt_merge_sep_sort2, aes(x=Order_1, y=Order_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

#can i even try to do a upper triangle heat map...  https://albert-rapp.de/posts/ggplot2-tips/24_correlation_heat_map/24_correlation_heat_map.html

p_cor_df_melt_merge_sep_sort2$Order_1 <- as.factor(p_cor_df_melt_merge_sep_sort2$Order_1)
p_cor_df_melt_merge_sep_sort2$Order_2 <- as.factor(p_cor_df_melt_merge_sep_sort2$Order_2)


p_cor_df_melt_merge_sep_sort3 <- p_cor_df_melt_merge_sep_sort2 |> 
  mutate(
    lvl_a = as.numeric(Order_1),
    lvl_b = as.numeric(Order_2),
    cor = if_else(lvl_a < lvl_b, cor, NA)
  ) 



ggplot(p_cor_df_melt_merge_sep_sort3, aes(x=Order_1, y=Order_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation",na.value = 'white') + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))


ggsave("all_order_heatmap.pdf",useDingbats=FALSE,bg="white")
#omit redundant rows (ie, Genus 2:Genus 1 pair is the same as Genus 1: Genus 2 pair)
p_cor_df_melt_merge_sep_sort6 <- na.omit(p_cor_df_melt_merge_sep_sort3)
#get only significant correlations
p_cor_df_melt_merge_sep_sort7 <- p_cor_df_melt_merge_sep_sort6[p_cor_df_melt_merge_sep_sort6$p < 0.05,]

#get top ten
top10corgenus <- head(p_cor_df_melt_merge_sep_sort7,10)
#get bottom ten
bottom10corgenus <- tail(p_cor_df_melt_merge_sep_sort7,10)
#bind
corpairs <- rbind(top10corgenus,bottom10corgenus)

#get pair for plotting
corpairs$genpair <- paste(corpairs$Order_1,corpairs$Order_2,sep="/")



ggplot(corpairs, aes(x=cor,y=reorder(genpair, cor))) + geom_col(fill="lightblue") + theme_cowplot() + background_grid() + scale_x_continuous(limits=c(-1,1),breaks=c(seq(-1,1,by=0.1))) +xlab("Spearman correlation coefficient") + ylab("Order pair")

ggsave("order_correlation_top_bottom_ten.pdf",useDingbats=FALSE,bg="white")



#okay, just trying to gut check this by looking at some scatterplots


#wait, let's transform first?

ps_clr <- microbiome::transform(ps3order, "clr")

otutab3 <- as.data.frame(ps_clr@otu_table)

otutab3$OTU <- as.factor(rownames(otutab))
otutab3$Order <- "a"

otutab4 <- otutab3

#add families....
for(i in unique(otutab4$OTU)){
	the_otu <- i
	the_otu_tax <- taxonomy_df[taxonomy_df$ASV == the_otu,]
	otutab4$Order[otutab4$OTU == the_otu] <- the_otu_tax$Order
}

otutab5 <- subset(otutab4, select = -c(OTU))
otutab5_melt <- reshape2::melt(otutab5,id.vars="Order")

#top Frankiales/Kineosporiales
#bottom Micrococcales/Rhodobacterales



Frank <- otutab5_melt[otutab5_melt$Order=="Frankiales",]
Kineo <- otutab5_melt[otutab5_melt$Order=="Kineosporiales",]

nrow(Frank)
nrow(Kineo)

Frank <- subset(Frank, select = -c(Order))
Kineo <- subset(Kineo, select = -c(Order))

FrankKineo <- merge(Frank,Kineo,by="variable")

names(FrankKineo)[names(FrankKineo) == 'value.x'] <- 'Frankiales'
names(FrankKineo)[names(FrankKineo) == 'value.y'] <- 'Kineosporiales'


ggplot(FrankKineo, aes(x=Frankiales,y=Kineosporiales)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Frankiales transformed abundance") + ylab("Kineosporiales transformed abundance") +scale_x_continuous(limits=c(-2,10),breaks=c(seq(-2,10,1))) +scale_y_continuous(limits=c(-2,10),breaks=c(seq(-2,10,1)))


ggsave("Frankiales_Kineosporiales_ORDER_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")



Micro <- otutab5_melt[otutab5_melt$Order=="Micrococcales",]
Rhodo <- otutab5_melt[otutab5_melt$Order=="Rhodobacterales",]

nrow(Micro)
nrow(Rhodo)

Micro <- subset(Micro, select = -c(Order))
Rhodo <- subset(Rhodo, select = -c(Order))

MicroRhodo <- merge(Micro,Rhodo,by="variable")

names(MicroRhodo)[names(MicroRhodo) == 'value.x'] <- 'Micrococcales'
names(MicroRhodo)[names(MicroRhodo) == 'value.y'] <- 'Rhodobacterales'


ggplot(MicroRhodo, aes(x=Micrococcales,y=Rhodobacterales)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Micrococcales transformed abundance") + ylab("Rhodobacterales transformed abundance") +scale_x_continuous(limits=c(-2,10),breaks=c(seq(-2,10,1))) +scale_y_continuous(limits=c(-2,10),breaks=c(seq(-2,10,1)))


ggsave("Micrococcales_Rhodobacterales_ORDER_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")


#class level



ps3class <- tax_glom(ps3, "Class", NArm = TRUE)


taxtab <- as.data.frame(tax_table(ps3class))

otutab <- as.data.frame(ps3class@otu_table)

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

# use clone to get a deep copy of soil_amp (R6 object)
tmp <- clone(mt)

# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
    #the pairwise correlation coefficients are in tmp$res_cor_p as a list! (see manual entry for trans_network())


# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_all_class.gexf")
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
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_surface_class.gexf")

# put the network into the list
mt_network$surface <- tmp
# select samples of "TW" group
tmp <- clone(mt)
tmp$sample_table %<>% subset(surface_or_interior == "interior")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_interior_class.gexf")

mt_network$interior <- tmp



str(mt_network$all$res_cor_p)

all_cor_df <- as.data.frame(mt_network$all$res_cor_p$cor)

write.table(all_cor_df,"all_cor_df_class.tsv",sep="\t")

all_p_df <- as.data.frame(mt_network$all$res_cor_p$p)

all_cor_df$OTU1 <- rownames(all_cor_df)
all_p_df$OTU1 <- rownames(all_p_df)

library(reshape2)

all_cor_df_melt <- reshape2::melt(all_cor_df)
all_p_df_melt <- reshape2::melt(all_p_df)

names(all_cor_df_melt)[names(all_cor_df_melt) == 'variable'] <- 'OTU2'
names(all_p_df_melt)[names(all_p_df_melt) == 'variable'] <- 'OTU2'

names(all_cor_df_melt)[names(all_cor_df_melt) == 'value'] <- 'cor'
names(all_p_df_melt)[names(all_p_df_melt) == 'value'] <- 'p'

all_cor_df_melt$OTU1_OTU2 <- paste(all_cor_df_melt$OTU1,all_cor_df_melt$OTU2)

all_p_df_melt$OTU1_OTU2 <- paste(all_p_df_melt$OTU1,all_p_df_melt$OTU2)

all_cor_df_melt_ii <- data.frame(OTU1_OTU2=all_cor_df_melt$OTU1_OTU2,cor=all_cor_df_melt$cor)
all_p_df_melt_ii <- data.frame(OTU1_OTU2=all_p_df_melt$OTU1_OTU2,p=all_p_df_melt$p)

p_cor_df_melt_merge <- merge(all_cor_df_melt_ii,all_p_df_melt_ii)


library(tidyr)

#split into two
p_cor_df_melt_merge_sep <- separate(p_cor_df_melt_merge, OTU1_OTU2, into = c("OTU1","OTU2"), sep = " (?=[^ ]+$)")

head(p_cor_df_melt_merge_sep)

p_cor_df_melt_merge_sep_sort <- p_cor_df_melt_merge_sep[order(p_cor_df_melt_merge_sep$cor, decreasing = TRUE),]  



library(ggplot2)
library(cowplot)
ggplot(p_cor_df_melt_merge_sep_sort, aes(x=OTU1, y=OTU2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal()

taxonomy_df = read.csv("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/05_qiime_feature-classifier_classify-consensus-vsearch/export/taxonomy_revised_fig_samples_only_5-2025.tsv", header=TRUE)


#okay, replace ASV id with Order id
taxonomy_df_OTU_cor_OTU1 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU1),]

taxonomy_df_OTU_cor_OTU2 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU2),]


p_cor_df_melt_merge_sep_sort$Class_1 <- "a"

p_cor_df_melt_merge_sep_sort$Class_2 <- "b"
p_cor_df_melt_merge_sep_sort$OTU1 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU1)
p_cor_df_melt_merge_sep_sort$OTU2 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU2)

p_cor_df_melt_merge_sep_sort2 <- p_cor_df_melt_merge_sep_sort


for(i in unique(p_cor_df_melt_merge_sep_sort2$OTU1)){
  the_otu <- i
  the_otu_tax <- taxonomy_df_OTU_cor_OTU1[taxonomy_df_OTU_cor_OTU1$ASV == the_otu,]
  p_cor_df_melt_merge_sep_sort2$Class_1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu] <- the_otu_tax$Class
  p_cor_df_melt_merge_sep_sort2$Class_2[p_cor_df_melt_merge_sep_sort2$OTU2 == the_otu] <- the_otu_tax$Class

}


ggplot(p_cor_df_melt_merge_sep_sort2, aes(x=Class_1, y=Class_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

#can i even try to do a upper triangle heat map...  https://albert-rapp.de/posts/ggplot2-tips/24_correlation_heat_map/24_correlation_heat_map.html

p_cor_df_melt_merge_sep_sort2$Class_1 <- as.factor(p_cor_df_melt_merge_sep_sort2$Class_1)
p_cor_df_melt_merge_sep_sort2$Class_2 <- as.factor(p_cor_df_melt_merge_sep_sort2$Class_2)


p_cor_df_melt_merge_sep_sort3 <- p_cor_df_melt_merge_sep_sort2 |> 
  mutate(
    lvl_a = as.numeric(Class_1),
    lvl_b = as.numeric(Class_2),
    cor = if_else(lvl_a < lvl_b, cor, NA)
  ) 



ggplot(p_cor_df_melt_merge_sep_sort3, aes(x=Class_1, y=Class_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation",na.value = 'white') + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))


ggsave("all_class_heatmap.pdf",useDingbats=FALSE,bg="white")
#omit redundant rows (ie, Genus 2:Genus 1 pair is the same as Genus 1: Genus 2 pair)
p_cor_df_melt_merge_sep_sort6 <- na.omit(p_cor_df_melt_merge_sep_sort3)
#get only significant correlations
p_cor_df_melt_merge_sep_sort7 <- p_cor_df_melt_merge_sep_sort6[p_cor_df_melt_merge_sep_sort6$p < 0.05,]

#get top ten
top10corgenus <- head(p_cor_df_melt_merge_sep_sort7,10)
#get bottom ten
bottom10corgenus <- tail(p_cor_df_melt_merge_sep_sort7,10)
#bind
corpairs <- rbind(top10corgenus,bottom10corgenus)

#get pair for plotting
corpairs$genpair <- paste(corpairs$Class_1,corpairs$Class_2,sep="/")



ggplot(corpairs, aes(x=cor,y=reorder(genpair, cor))) + geom_col(fill="lightblue") + theme_cowplot() + background_grid() + scale_x_continuous(limits=c(-1,1),breaks=c(seq(-1,1,by=0.1))) +xlab("Spearman correlation coefficient") + ylab("Class pair")

ggsave("class_correlation_top_bottom_ten.pdf",useDingbats=FALSE,bg="white")



#okay, just trying to gut check this by looking at some scatterplots


#wait, let's transform first?

ps_clr <- microbiome::transform(ps3class, "clr")

otutab3 <- as.data.frame(ps_clr@otu_table)

otutab3$OTU <- as.factor(rownames(otutab))
otutab3$Class <- "a"

otutab4 <- otutab3

#add families....
for(i in unique(otutab4$OTU)){
  the_otu <- i
  the_otu_tax <- taxonomy_df[taxonomy_df$ASV == the_otu,]
  otutab4$Class[otutab4$OTU == the_otu] <- the_otu_tax$Class
}

otutab5 <- subset(otutab4, select = -c(OTU))
otutab5_melt <- reshape2::melt(otutab5,id.vars="Class")

#top Actinobacteria/Alphaproteobacteria
#bottom Deinococci/Gammaproteobacteria



Frank <- otutab5_melt[otutab5_melt$Class=="Actinobacteria",]
Kineo <- otutab5_melt[otutab5_melt$Class=="Alphaproteobacteria",]

nrow(Frank)
nrow(Kineo)

Frank <- subset(Frank, select = -c(Class))
Kineo <- subset(Kineo, select = -c(Class))

FrankKineo <- merge(Frank,Kineo,by="variable")

names(FrankKineo)[names(FrankKineo) == 'value.x'] <- 'Actinobacteria'
names(FrankKineo)[names(FrankKineo) == 'value.y'] <- 'Alphaproteobacteria'


ggplot(FrankKineo, aes(x=Actinobacteria,y=Alphaproteobacteria)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Actinobacteria transformed abundance") + ylab("Alphaproteobacteria transformed abundance") +scale_x_continuous(limits=c(-2,11),breaks=c(seq(-2,11,1))) +scale_y_continuous(limits=c(-2,11),breaks=c(seq(-2,11,1)))


ggsave("Actinobacteria_Alphaproteobacteria_CLASS_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")



Micro <- otutab5_melt[otutab5_melt$Class=="Deinococci",]
Rhodo <- otutab5_melt[otutab5_melt$Class=="Gammaproteobacteria",]

nrow(Micro)
nrow(Rhodo)

Micro <- subset(Micro, select = -c(Class))
Rhodo <- subset(Rhodo, select = -c(Class))

MicroRhodo <- merge(Micro,Rhodo,by="variable")

names(MicroRhodo)[names(MicroRhodo) == 'value.x'] <- 'Deinococci'
names(MicroRhodo)[names(MicroRhodo) == 'value.y'] <- 'Gammaproteobacteria'


ggplot(MicroRhodo, aes(x=Deinococci,y=Gammaproteobacteria)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Deinococci transformed abundance") + ylab("Gammaproteobacteria transformed abundance") +scale_x_continuous(limits=c(-4,12),breaks=c(seq(-4,12,1))) +scale_y_continuous(limits=c(-4,12),breaks=c(seq(-4,12,1)))


ggsave("Deinococci_Gammaproteobacteria_Class_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")


#phylum



ps3phylum <- tax_glom(ps3, "Phylum", NArm = TRUE)


taxtab <- as.data.frame(tax_table(ps3phylum))

otutab <- as.data.frame(ps3phylum@otu_table)

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

# use clone to get a deep copy of soil_amp (R6 object)
tmp <- clone(mt)

# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
    #the pairwise correlation coefficients are in tmp$res_cor_p as a list! (see manual entry for trans_network())


# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_all_phylum.gexf")
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
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_surface_phylum.gexf")

# put the network into the list
mt_network$surface <- tmp
# select samples of "TW" group
tmp <- clone(mt)
tmp$sample_table %<>% subset(surface_or_interior == "interior")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.001)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
tmp$cal_module(method = "cluster_fast_greedy")
tmp$save_network(filepath = "network_interior_phylum.gexf")

mt_network$interior <- tmp



str(mt_network$all$res_cor_p)

all_cor_df <- as.data.frame(mt_network$all$res_cor_p$cor)

write.table(all_cor_df,"all_cor_df_phylum.tsv",sep="\t")

all_p_df <- as.data.frame(mt_network$all$res_cor_p$p)

all_cor_df$OTU1 <- rownames(all_cor_df)
all_p_df$OTU1 <- rownames(all_p_df)

library(reshape2)

all_cor_df_melt <- reshape2::melt(all_cor_df)
all_p_df_melt <- reshape2::melt(all_p_df)

names(all_cor_df_melt)[names(all_cor_df_melt) == 'variable'] <- 'OTU2'
names(all_p_df_melt)[names(all_p_df_melt) == 'variable'] <- 'OTU2'

names(all_cor_df_melt)[names(all_cor_df_melt) == 'value'] <- 'cor'
names(all_p_df_melt)[names(all_p_df_melt) == 'value'] <- 'p'

all_cor_df_melt$OTU1_OTU2 <- paste(all_cor_df_melt$OTU1,all_cor_df_melt$OTU2)

all_p_df_melt$OTU1_OTU2 <- paste(all_p_df_melt$OTU1,all_p_df_melt$OTU2)

all_cor_df_melt_ii <- data.frame(OTU1_OTU2=all_cor_df_melt$OTU1_OTU2,cor=all_cor_df_melt$cor)
all_p_df_melt_ii <- data.frame(OTU1_OTU2=all_p_df_melt$OTU1_OTU2,p=all_p_df_melt$p)

p_cor_df_melt_merge <- merge(all_cor_df_melt_ii,all_p_df_melt_ii)


library(tidyr)

#split into two
p_cor_df_melt_merge_sep <- separate(p_cor_df_melt_merge, OTU1_OTU2, into = c("OTU1","OTU2"), sep = " (?=[^ ]+$)")

head(p_cor_df_melt_merge_sep)

p_cor_df_melt_merge_sep_sort <- p_cor_df_melt_merge_sep[order(p_cor_df_melt_merge_sep$cor, decreasing = TRUE),]  



library(ggplot2)
library(cowplot)
ggplot(p_cor_df_melt_merge_sep_sort, aes(x=OTU1, y=OTU2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal()

taxonomy_df = read.csv("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/05_qiime_feature-classifier_classify-consensus-vsearch/export/taxonomy_revised_fig_samples_only_5-2025.tsv", header=TRUE)


#okay, replace ASV id with Phylum id
taxonomy_df_OTU_cor_OTU1 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU1),]

taxonomy_df_OTU_cor_OTU2 <- taxonomy_df[taxonomy_df$ASV %in% unique(p_cor_df_melt_merge_sep_sort$OTU2),]


p_cor_df_melt_merge_sep_sort$Phylum_1 <- "a"

p_cor_df_melt_merge_sep_sort$Phylum_2 <- "b"
p_cor_df_melt_merge_sep_sort$OTU1 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU1)
p_cor_df_melt_merge_sep_sort$OTU2 <- as.factor(p_cor_df_melt_merge_sep_sort$OTU2)

p_cor_df_melt_merge_sep_sort2 <- p_cor_df_melt_merge_sep_sort


for(i in unique(p_cor_df_melt_merge_sep_sort2$OTU1)){
  the_otu <- i
  the_otu_tax <- taxonomy_df_OTU_cor_OTU1[taxonomy_df_OTU_cor_OTU1$ASV == the_otu,]
  p_cor_df_melt_merge_sep_sort2$Phylum_1[p_cor_df_melt_merge_sep_sort2$OTU1 == the_otu] <- the_otu_tax$Phylum
  p_cor_df_melt_merge_sep_sort2$Phylum_2[p_cor_df_melt_merge_sep_sort2$OTU2 == the_otu] <- the_otu_tax$Phylum

}


ggplot(p_cor_df_melt_merge_sep_sort2, aes(x=Phylum_1, y=Phylum_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))

#can i even try to do a upper triangle heat map...  https://albert-rapp.de/posts/ggplot2-tips/24_correlation_heat_map/24_correlation_heat_map.html

p_cor_df_melt_merge_sep_sort2$Phylum_1 <- as.factor(p_cor_df_melt_merge_sep_sort2$Phylum_1)
p_cor_df_melt_merge_sep_sort2$Phylum_2 <- as.factor(p_cor_df_melt_merge_sep_sort2$Phylum_2)


p_cor_df_melt_merge_sep_sort3 <- p_cor_df_melt_merge_sep_sort2 |> 
  mutate(
    lvl_a = as.numeric(Phylum_1),
    lvl_b = as.numeric(Phylum_2),
    cor = if_else(lvl_a < lvl_b, cor, NA)
  ) 



ggplot(p_cor_df_melt_merge_sep_sort3, aes(x=Phylum_1, y=Phylum_2, fill=cor)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "darkgreen", mid = "white", midpoint = 0, limit = c(-1,1), name="Spearman\nCorrelation",na.value = 'white') + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1))


ggsave("all_Phylum_heatmap.pdf",useDingbats=FALSE,bg="white")
#omit redundant rows (ie, Genus 2:Genus 1 pair is the same as Genus 1: Genus 2 pair)
p_cor_df_melt_merge_sep_sort6 <- na.omit(p_cor_df_melt_merge_sep_sort3)
#get only significant correlations
p_cor_df_melt_merge_sep_sort7 <- p_cor_df_melt_merge_sep_sort6[p_cor_df_melt_merge_sep_sort6$p < 0.05,]

#get top ten
top10corgenus <- head(p_cor_df_melt_merge_sep_sort7,10)
#get bottom ten
bottom10corgenus <- tail(p_cor_df_melt_merge_sep_sort7,10)
#bind
corpairs <- rbind(top10corgenus,bottom10corgenus)

#get pair for plotting
corpairs$genpair <- paste(corpairs$Phylum_1,corpairs$Phylum_2,sep="/")



ggplot(corpairs, aes(x=cor,y=reorder(genpair, cor))) + geom_col(fill="lightblue") + theme_cowplot() + background_grid() + scale_x_continuous(limits=c(-1,1),breaks=c(seq(-1,1,by=0.1))) +xlab("Spearman correlation coefficient") + ylab("Phylum pair")

ggsave("Phylum_correlation_top_bottom_ten.pdf",useDingbats=FALSE,bg="white")



#okay, just trying to gut check this by looking at some scatterplots


#wait, let's transform first?

ps_clr <- microbiome::transform(ps3phylum, "clr")

otutab3 <- as.data.frame(ps_clr@otu_table)

otutab3$OTU <- as.factor(rownames(otutab))
otutab3$Phylum <- "a"

otutab4 <- otutab3

#add Phyla
for(i in unique(otutab4$OTU)){
  the_otu <- i
  the_otu_tax <- taxonomy_df[taxonomy_df$ASV == the_otu,]
  otutab4$Phylum[otutab4$OTU == the_otu] <- the_otu_tax$Phylum
}

otutab5 <- subset(otutab4, select = -c(OTU))
otutab5_melt <- reshape2::melt(otutab5,id.vars="Phylum")

#top Abditibacteriota/Bdellovibrionota


Frank <- otutab5_melt[otutab5_melt$Phylum=="Abditibacteriota",]
Kineo <- otutab5_melt[otutab5_melt$Phylum=="Bdellovibrionota",]

nrow(Frank)
nrow(Kineo)

Frank <- subset(Frank, select = -c(Phylum))
Kineo <- subset(Kineo, select = -c(Phylum))

FrankKineo <- merge(Frank,Kineo,by="variable")

names(FrankKineo)[names(FrankKineo) == 'value.x'] <- 'Abditibacteriota'
names(FrankKineo)[names(FrankKineo) == 'value.y'] <- 'Bdellovibrionota'


ggplot(FrankKineo, aes(x=Abditibacteriota,y=Bdellovibrionota)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Abditibacteriota transformed abundance") + ylab("Bdellovibrionota transformed abundance") +scale_x_continuous(limits=c(-4,11),breaks=c(seq(-4,11,1))) +scale_y_continuous(limits=c(-4,11),breaks=c(seq(-4,11,1)))


ggsave("Abditibacteriota_Bdellovibrionota_PHYLUM_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")



Micro <- otutab5_melt[otutab5_melt$Class=="Deinococci",]
Rhodo <- otutab5_melt[otutab5_melt$Class=="Gammaproteobacteria",]

nrow(Micro)
nrow(Rhodo)

Micro <- subset(Micro, select = -c(Class))
Rhodo <- subset(Rhodo, select = -c(Class))

MicroRhodo <- merge(Micro,Rhodo,by="variable")

names(MicroRhodo)[names(MicroRhodo) == 'value.x'] <- 'Deinococci'
names(MicroRhodo)[names(MicroRhodo) == 'value.y'] <- 'Gammaproteobacteria'


ggplot(MicroRhodo, aes(x=Deinococci,y=Gammaproteobacteria)) + geom_point() + geom_smooth(method="lm",se=FALSE,linetype="dotted") + theme_cowplot() +xlab("Deinococci transformed abundance") + ylab("Gammaproteobacteria transformed abundance") +scale_x_continuous(limits=c(-4,12),breaks=c(seq(-4,12,1))) +scale_y_continuous(limits=c(-4,12),breaks=c(seq(-4,12,1)))


ggsave("Deinococci_Gammaproteobacteria_Class_transformed_abundance_scatterplot_genus.pdf",useDingbats=FALSE,bg="white")

