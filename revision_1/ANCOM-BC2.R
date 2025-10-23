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


#just fig samples

#ps3_just_figs <- subset_samples(ps3, nematode_species == "C. inopinata")

#just fig interiors (suspensions)
#ps4 <- subset_samples(ps3_just_figs, category == "woodruff_fig_interior")


#https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html

library(ANCOMBC)

#ASV | All samples
set.seed(666)
?ancombc2
#what is *needed* here?
    #fix_formula="category"
    #do p_adj_method="BH"
    #group="category" 
    #global=TRUE
    #pairwise=TRUE
    #dunnet = TRUE

head(tax_table(ps3))
#ASV
ABC2_ps3_otu_res <- ancombc2(data = ps3,tax_level = "OTU",fix_formula="category",group="category",p_adj_method="BH",global=TRUE,pairwise=TRUE,dunnet = TRUE)

write.table(ABC2_ps3_otu_res$ss_tab,"ABC2_CATEGORY_ASV_res_ss_tab.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_otu_res$res,"ABC2_CATEGORY_ASV_res_res.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_otu_res$res_global,"ABC2_CATEGORY_ASV_res_res_global.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_otu_res$res_pair,"ABC2_CATEGORY_ASV_res_res_pair.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_otu_res$res_dunn,"ABC2_CATEGORY_ASV_res_res_dunn.tsv",quote=FALSE,sep="\t",row.names = FALSE)

ABC2_ps3_genus_res <- ancombc2(data = ps3,tax_level = "Genus",fix_formula="category",group="category",p_adj_method="BH",global=TRUE,pairwise=TRUE,dunnet = TRUE)

write.table(ABC2_ps3_genus_res$ss_tab,"ABC2_CATEGORY_genus_res_ss_tab.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_genus_res$res,"ABC2_CATEGORY_genus_res_res.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_genus_res$res_global,"ABC2_CATEGORY_genus_res_res_global.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_genus_res$res_pair,"ABC2_CATEGORY_genus_res_res_pair.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_genus_res$res_dunn,"ABC2_CATEGORY_genus_res_res_dunn.tsv",quote=FALSE,sep="\t",row.names = FALSE)


ABC2_ps3_Family_res <- ancombc2(data = ps3,tax_level = "Family",fix_formula="category",group="category",p_adj_method="BH",global=TRUE,pairwise=TRUE,dunnet = TRUE)

write.table(ABC2_ps3_Family_res$ss_tab,"ABC2_CATEGORY_Family_res_ss_tab.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_Family_res$res,"ABC2_CATEGORY_Family_res_res.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_Family_res$res_global,"ABC2_CATEGORY_Family_res_res_global.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_Family_res$res_pair,"ABC2_CATEGORY_Family_res_res_pair.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_Family_res$res_dunn,"ABC2_CATEGORY_Family_res_res_dunn.tsv",quote=FALSE,sep="\t",row.names = FALSE)

#you know, let's just try a for loop

the_ranks <- c("Domain","Phylum","Class","Order","Family","Genus","Species","OTU")

big_ss_tab <- NULL
big_res <- NULL
big_res_global <- NULL
big_res_pair <- NULL
big_res_dunn <- NULL

ABC2_results <- NULL

for (i in the_ranks){
    ABC2_results <- ancombc2(data = ps3,tax_level = i,fix_formula="category",group="category",p_adj_method="BH",global=TRUE,pairwise=TRUE,dunnet = TRUE)
    the_ss_tab <- as.data.frame(ABC2_results$ss_tab)
    the_res <- as.data.frame(ABC2_results$res)
    the_res_global <- as.data.frame(ABC2_results$res_global)
    the_res_pair <- as.data.frame(ABC2_results$res_pair)
    the_res_dunn <- as.data.frame(ABC2_results$res_dunn)
    the_ss_tab$rank <- i
    the_res$rank <- i
    the_res_global$rank <- i
    the_res_pair$rank <- i
    the_res_dunn$rank <- i
    big_ss_tab <- rbind(big_ss_tab,the_ss_tab)
    big_res <- rbind(big_res,the_res)
    big_res_global <- rbind(big_res_global,the_res_global)
    big_res_pair <- rbind(big_res_pair,the_res_pair)
    big_res_dunn <- rbind(big_res_dunn,the_res_dunn)
}

write.table(big_ss_tab,"ABC2_CATEGORY_all_ranks_results_ss_tab.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(big_res,"ABC2_CATEGORY_all_ranks_results_res.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(big_res_global,"ABC2_CATEGORY_all_ranks_results_res_global.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(big_res_pair,"ABC2_CATEGORY_all_ranks_results_res_pair.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(big_res_dunn,"ABC2_CATEGORY_all_ranks_results_res_dunn.tsv",quote=FALSE,sep="\t",row.names = FALSE)

#now, the elegans-inopinata species comparison


ABC2_ps3_otu_INOELE_res <- ancombc2(data = ps3,tax_level = "OTU",fix_formula="nematode_species",group="nematode_species",p_adj_method="BH")

write.table(ABC2_ps3_otu_INOELE_res$ss_tab,"ABC2_INOELE_ASV_res_ss_tab.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_otu_INOELE_res$res,"ABC2_INOELE_ASV_res_res.tsv",quote=FALSE,sep="\t",row.names = FALSE)

ABC2_ps3_otu_INOELE_res_sig <- ABC2_ps3_otu_INOELE_res$res[ABC2_ps3_otu_INOELE_res$res$`diff_nematode_speciesC. inopinata` == TRUE,]
ABC2_ps3_otu_INOELE_res_sig_pass_ss <- ABC2_ps3_otu_INOELE_res_sig[ABC2_ps3_otu_INOELE_res_sig$`passed_ss_nematode_speciesC. inopinata` == TRUE,]

write.table(ABC2_ps3_otu_INOELE_res_sig_pass_ss,"ABC2_INOELE_ASV_res_res_sig_pass_sensitivity.tsv",quote=FALSE,sep="\t",row.names = FALSE)


#genus


ABC2_ps3_genus_INOELE_res <- ancombc2(data = ps3,tax_level = "Genus",fix_formula="nematode_species",p_adj_method="BH")

write.table(ABC2_ps3_genus_INOELE_res$ss_tab,"ABC2_INOELE_genus_res_ss_tab.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(ABC2_ps3_genus_INOELE_res$res,"ABC2_INOELE_genus_res_res.tsv",quote=FALSE,sep="\t",row.names = FALSE)

ABC2_ps3_genus_INOELE_res_sig <- ABC2_ps3_genus_INOELE_res$res[ABC2_ps3_genus_INOELE_res$res$`diff_nematode_speciesC. inopinata` == TRUE,]
ABC2_ps3_genus_INOELE_res_sig_pass_ss <- ABC2_ps3_genus_INOELE_res_sig[ABC2_ps3_genus_INOELE_res_sig$`passed_ss_nematode_speciesC. inopinata` == TRUE,]

write.table(ABC2_ps3_genus_INOELE_res_sig_pass_ss,"ABC2_INOELE_genus_res_res_sig_pass_sensitivity.tsv",quote=FALSE,sep="\t",row.names = FALSE)

#just do them all...
big_ss_tab <- NULL
big_res <- NULL
big_res_global <- NULL
big_res_pair <- NULL
big_res_dunn <- NULL

ABC2_results <- NULL

for (i in the_ranks){
    ABC2_results <- ancombc2(data = ps3,tax_level = i,fix_formula="nematode_species",group="nematode_species",p_adj_method="BH")
    the_ss_tab <- as.data.frame(ABC2_results$ss_tab)
    the_res <- as.data.frame(ABC2_results$res)
    the_ss_tab$rank <- i
    the_res$rank <- i
    big_ss_tab <- rbind(big_ss_tab,the_ss_tab)
    big_res <- rbind(big_res,the_res)
}

write.table(big_ss_tab,"ABC2_nematode_species_all_ranks_results_ss_tab.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(big_res,"ABC2_nematode_species_all_ranks_results_res.tsv",quote=FALSE,sep="\t",row.names = FALSE)

big_res_sig_pass_ss <- big_res[big_res$`diff_robust_nematode_speciesC. inopinata` == TRUE,]

write.table(big_res_sig_pass_ss,"ABC2_nematode_species_all_ranks_results_sig_pass_sensitivity.tsv",quote=FALSE,sep="\t",row.names = FALSE)


big_res_sig_pass_ss <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/ANCOMBC/ABC2_nematode_species_all_ranks_results_sig_pass_sensitivity.tsv",header=TRUE,sep="\t")

unique(big_res_sig_pass_ss$rank)

big_res_sig_pass_ss_OTU <- big_res_sig_pass_ss[big_res_sig_pass_ss$rank == "OTU",]
big_res_sig_pass_ss_genus <- big_res_sig_pass_ss[big_res_sig_pass_ss$rank == "Genus",]
big_res_sig_pass_ss_family <- big_res_sig_pass_ss[big_res_sig_pass_ss$rank == "Family",]


nrow(big_res_sig_pass_ss_OTU)
#[1] 97

nrow(big_res_sig_pass_ss_genus)
#[1] 210

nrow(big_res_sig_pass_ss_family)
#[1] 133

sigdf <- big_res_sig_pass_ss


#figs only

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
ps3_figs_only <- subset_taxa(no_controls, !OTU %in% control_OTU$OTU)

#interior samples only, worm occupancy

ps4_fig_interiors <- subset_samples(ps3_figs_only, surface_or_interior == "interior")



#inside v outside, figs only


#just do them all...
big_ss_tab <- NULL
big_res <- NULL
big_res_global <- NULL
big_res_pair <- NULL
big_res_dunn <- NULL

ABC2_results <- NULL

for (i in the_ranks){
    ABC2_results <- ancombc2(data = ps3_figs_only,tax_level = i,fix_formula="surface_or_interior",group="surface_or_interior",p_adj_method="BH")
    the_ss_tab <- as.data.frame(ABC2_results$ss_tab)
    the_res <- as.data.frame(ABC2_results$res)
    the_ss_tab$rank <- i
    the_res$rank <- i
    big_ss_tab <- rbind(big_ss_tab,the_ss_tab)
    big_res <- rbind(big_res,the_res)
}


write.table(big_ss_tab,"ABC2_figs_only_surface_or_interior_all_ranks_results_ss_tab.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(big_res,"ABC2_figs_only_surface_or_interior_all_ranks_results_res.tsv",quote=FALSE,sep="\t",row.names = FALSE)

big_res_sig_pass_ss <- big_res[big_res$diff_robust_surface_or_interiorsurface == TRUE,]

write.table(big_res_sig_pass_ss,"ABC2_figs_only_surface_or_interior_all_ranks_results_sig_pass_sensitivity.tsv",quote=FALSE,sep="\t",row.names = FALSE)





#just do them all...
big_ss_tab <- NULL
big_res <- NULL
big_res_global <- NULL
big_res_pair <- NULL
big_res_dunn <- NULL

ABC2_results <- NULL

for (i in the_ranks){
    ABC2_results <- ancombc2(data = ps4_fig_interiors,tax_level = i,fix_formula="worms_present",group="worms_present",p_adj_method="BH")
    the_ss_tab <- as.data.frame(ABC2_results$ss_tab)
    the_res <- as.data.frame(ABC2_results$res)
    the_ss_tab$rank <- i
    the_res$rank <- i
    big_ss_tab <- rbind(big_ss_tab,the_ss_tab)
    big_res <- rbind(big_res,the_res)
}


write.table(big_ss_tab,"ABC2_fig_interior_only_worms_present_all_ranks_results_ss_tab.tsv",quote=FALSE,sep="\t",row.names = FALSE)
write.table(big_res,"ABC2_fig_interior_only_worms_present_all_ranks_results_res.tsv",quote=FALSE,sep="\t",row.names = FALSE)

big_res_sig_pass_ss <- big_res[big_res$diff_robust_worms_presentyes == TRUE,]

write.table(big_res_sig_pass_ss,"ABC2_fig_interior_only_worms_present_all_ranks_results_sig_pass_sensitivity.tsv",quote=FALSE,sep="\t",row.names = FALSE)




#the supplemental figures

#inop v eleg substrates


#ASV-- let's top and bottom five in terms of effect size, among those robust, differentially abundant ASV's

sigdf <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/ANCOMBC/ABC2_nematode_species_all_ranks_results_sig_pass_sensitivity.tsv",header=TRUE,sep="\t")

sigdf_ASV <- sigdf[sigdf$rank == "OTU",]

#how many?
nrow(sigdf_ASV)
#[1] 97

sigdf_sort <- sigdf_ASV[order(sigdf_ASV$`lfc_nematode_speciesC..inopinata`),]

#how many diff abundant in elegans

nrow(sigdf_ASV[sigdf_ASV$`lfc_nematode_speciesC..inopinata` < 0,])
#[1] 34

#how many diff abundant in inopinata
nrow(sigdf_ASV[sigdf_ASV$`lfc_nematode_speciesC..inopinata` > 0,])
#[1] 63

top_five <- head(sigdf_sort,5)

bottom_five <- tail(sigdf_sort,5)

asv_to_get <- c(top_five$taxon,bottom_five$taxon)

taxdf <- as.data.frame(tax_table(ps3))

asv_to_get_taxonomy <- taxdf[taxdf$OTU %in% asv_to_get, ]

#okay prep df


#ps_family <- tax_glom(ps3, "Family", NArm = TRUE)


ps3_clr <- microbiome::transform(ps3, "clr")



mphyseq = psmelt(ps3_clr)

mphyseq$taxa_OTU <- as.factor(mphyseq$taxa_OTU)

mphyseqbigeff <- mphyseq[mphyseq$taxa_OTU %in% asv_to_get, ]

mphyseqbigeff$taxa_OTU <- factor(mphyseqbigeff$taxa_OTU, levels=asv_to_get)

#                                                          Genus
#1f23f702c272546c7a741f770101aba9                    Luteibacter
#90c8b1fafd90f04d6364c7620d422945       Candidatus_Nitrocosmicus
#8db7958e664794610207ab124860e386
#7769cac4c6bfbe7e3cddd33636c22519 Methylobacterium-Methylorubrum
#64985768bd5c61edeb534b014769186f                    Pseudomonas
#313ff9f5a805c6777206081ea4910b17                    Pseudomonas
#ad0afd7293b32bf565de97915ef4c20a                    Pseudomonas
#293f4657f2ecfd93549cdf706d420f00
#af0fca738bb4eea2a617983a19b1c74e                      Kosakonia
#5903d6eeb80fd89396b445f2c4874865                 Chthoniobacter

levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='1f23f702c272546c7a741f770101aba9'] <- 'Luteibacter ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='90c8b1fafd90f04d6364c7620d422945'] <- 'Ca. Nitrocosmicus ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='8db7958e664794610207ab124860e386'] <- 'Rhodobacteraceae ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='7769cac4c6bfbe7e3cddd33636c22519'] <- 'Methylobacterium-Methylorubrum ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='64985768bd5c61edeb534b014769186f'] <- 'Pseudomonas ASV 1'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='313ff9f5a805c6777206081ea4910b17'] <- 'Pseudomonas ASV 2'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='ad0afd7293b32bf565de97915ef4c20a'] <- 'Pseudomonas ASV 3'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='293f4657f2ecfd93549cdf706d420f00'] <- 'Stenotrophomonas ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='af0fca738bb4eea2a617983a19b1c74e'] <- 'Enterobacterales ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='5903d6eeb80fd89396b445f2c4874865'] <- 'Chthoniobacter ASV'


ggplot(mphyseqbigeff, aes(x=taxa_OTU,y=Abundance)) + geom_sina(aes(shape=category,colour=nematode_species),scale="width",alpha=0.5) + stat_summary(aes(group=category,colour=nematode_species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1,position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic"),legend.title=element_blank()) + scale_colour_brewer(palette="Set1", guide =guide_legend(label.theme = element_text(face = "italic"))) + xlab("ASV") + ylab("Transformed abundance")

ggsave("supplemental_figure_ASV_elegans_inopinata_8-4-2025.pdf",bg="white",height=8,width=15,units="in",useDingbats=FALSE)


#genus



sigdf_genus <- sigdf[sigdf$rank == "Genus",]

sigdf_sort <- sigdf_genus[order(sigdf_genus$`lfc_nematode_speciesC..inopinata`),]

top_five <- head(sigdf_sort,5)

bottom_five <- tail(sigdf_sort,5)

#a bunch of weird (read: not valid) genera here
genera_to_get <- c("Luteitalea","Gluconobacter","Planctomicrobium","Nitrospira","Pirellula",bottom_five$taxon)



#how many?
nrow(sigdf_genus)
#[1] 210
#how many diff abundant in elegans

nrow(sigdf_genus[sigdf_genus$`lfc_nematode_speciesC..inopinata` < 0,])
#[1] 142

#how many diff abundant in inopinata
nrow(sigdf_genus[sigdf_genus$`lfc_nematode_speciesC..inopinata` > 0,])
#[1] 68


#taxdf <- as.data.frame(tax_table(ps3))

genus_to_get_taxonomy <- taxdf[taxdf$Genus %in% genera_to_get, ]

#okay prep df


ps_genus <- tax_glom(ps3, "Genus", NArm = TRUE)


ps3_clr <- microbiome::transform(ps_genus, "clr")



mphyseq = psmelt(ps3_clr)

mphyseq$Genus <- as.factor(mphyseq$Genus)

mphyseqbigeff <- mphyseq[mphyseq$Genus %in% genera_to_get, ]

mphyseqbigeff$Genus <- factor(mphyseqbigeff$Genus, levels=genera_to_get)



ggplot(mphyseqbigeff, aes(x=Genus,y=Abundance)) + geom_sina(aes(shape=category,colour=nematode_species),scale="width",alpha=0.5) + stat_summary(aes(group=category,colour=nematode_species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1,position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic"),legend.title=element_blank()) + scale_colour_brewer(palette="Set1", guide =guide_legend(label.theme = element_text(face = "italic"))) + xlab("Genus") + ylab("Transformed abundance")

ggsave("supplemental_figure_genus_elegans_inopinata_8-4-2025.pdf",bg="white",height=8,width=15,units="in",useDingbats=FALSE)




#family




sigdf_family <- sigdf[sigdf$rank == "Family",]

sigdf_sort <- sigdf_family[order(sigdf_family$`lfc_nematode_speciesC..inopinata`),]

nrow(sigdf_family)
#[1] 133
nrow(sigdf_family[sigdf_family$`lfc_nematode_speciesC..inopinata` > 0,])
#[1] 30

nrow(sigdf_family[sigdf_family$`lfc_nematode_speciesC..inopinata` < 0,])
#[1] 103



sigdf_asv <- sigdf[sigdf$rank == "OTU",]

nrow(sigdf_asv)
#[1] 97
nrow(sigdf_asv[sigdf_asv$`lfc_nematode_speciesC..inopinata` > 0,])
#[1] 63

nrow(sigdf_asv[sigdf_asv$`lfc_nematode_speciesC..inopinata` < 0,])
#[1] 34



sigdf_genus <- sigdf[sigdf$rank == "Genus",]

nrow(sigdf_genus)
#[1] 210
nrow(sigdf_genus[sigdf_genus$`lfc_nematode_speciesC..inopinata` > 0,])
#[1] 68

nrow(sigdf_genus[sigdf_genus$`lfc_nematode_speciesC..inopinata` < 0,])
#[1] 142


head(sigdf_sort,10)



top_five <- head(sigdf_sort,5)

bottom_five <- tail(sigdf_sort,5)

#a bunch of weird (read: not valid) families here
families_to_get <- c("Saprospiraceae","Pirellulaceae","Rhodanobacteraceae","Vicinamibacteraceae","Phycisphaeraceae",bottom_five$taxon)

#taxdf <- as.data.frame(tax_table(ps3))

family_to_get_taxonomy <- taxdf[taxdf$Genus %in% families_to_get, ]

#okay prep df


ps_family <- tax_glom(ps3, "Family", NArm = TRUE)


ps3_clr <- microbiome::transform(ps_family, "clr")



mphyseq = psmelt(ps3_clr)

mphyseq$Family <- as.factor(mphyseq$Family)

mphyseqbigeff <- mphyseq[mphyseq$Family %in% families_to_get, ]

mphyseqbigeff$Family <- factor(mphyseqbigeff$Family, levels=families_to_get)


ggplot(mphyseqbigeff, aes(x=Family,y=Abundance)) + geom_sina(aes(shape=category,colour=nematode_species),scale="width",alpha=0.5) + stat_summary(aes(group=category,colour=nematode_species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1,position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic"),legend.title=element_blank()) + scale_colour_brewer(palette="Set1", guide =guide_legend(label.theme = element_text(face = "italic"))) + xlab("Family") + ylab("Transformed abundance")

ggsave("supplemental_figure_family_elegans_inopinata_8-4-2025.pdf",bg="white",height=8,width=15,units="in",useDingbats=FALSE)



#in figs only, interior/exterior















sigdf <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/ANCOMBC/ABC2_figs_only_surface_or_interior_all_ranks_results_sig_pass_sensitivity.tsv",header=TRUE,sep="\t")

sigdf_ASV <- sigdf[sigdf$rank == "OTU",]

#how many?
nrow(sigdf_ASV)
#[1] 11
#just plot them all!

sigdf_sort <- sigdf_ASV[order(sigdf_ASV$lfc_surface_or_interiorsurface),]

#how many diff abundant in interior

nrow(sigdf_ASV[sigdf_ASV$lfc_surface_or_interiorsurface < 0,])
#[1] 10

#how many diff abundant in surface
nrow(sigdf_ASV[sigdf_ASV$lfc_surface_or_interiorsurface > 0,])
#[1] 1

#top_five <- head(sigdf_sort,5)
#
#bottom_five <- tail(sigdf_sort,5)
#
#asv_to_get <- c(top_five$taxon,bottom_five$taxon)

taxdf <- as.data.frame(tax_table(ps3_figs_only))

asv_to_get_taxonomy <- taxdf[taxdf$OTU %in% sigdf_sort$taxon, ]

#okay prep df


#ps_family <- tax_glom(ps3, "Family", NArm = TRUE)


ps3_clr <- microbiome::transform(ps3_figs_only, "clr")



mphyseq = psmelt(ps3_clr)

mphyseq$taxa_OTU <- as.factor(mphyseq$taxa_OTU)

mphyseqbigeff <- mphyseq[mphyseq$taxa_OTU %in% sigdf_sort$taxon, ]

mphyseqbigeff$taxa_OTU <- factor(mphyseqbigeff$taxa_OTU, levels=sigdf_sort$taxon)

                                                                              #Genus
#7521edd61e578d0d49cfd01fbf8cc2f5                                        Xanthomonas
#40ffac9cd8b8ef62bba893e2223310f0
#7ce647f7d457226fc3950cd5ec4eefdf
#34ed5dee580cb7e4333f93fa7bbb2354                                          Kosakonia
#4552af3e905206840cf8b67bcaf2693f                                         Raoultella
#3d3f89fadf28973358458a50d436e91e                     Methylobacterium-Methylorubrum
#b37bf03859da8d004dbb9f4f4a220214
#718981e648d6a3a51625c565106f6881                                       Ochrobactrum
#fc6622c636a5210293fb2873fc4761d9 Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium
#607aee8b5f8abc3baf12a67738ff12d1
#42c4ad8258e5ae1d878307dd458b5990                                         Paracoccus

levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='7521edd61e578d0d49cfd01fbf8cc2f5'] <- 'Xanthomonas ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='40ffac9cd8b8ef62bba893e2223310f0'] <- 'Comamonadaceae ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='7ce647f7d457226fc3950cd5ec4eefdf'] <- 'Enterobacterales ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='34ed5dee580cb7e4333f93fa7bbb2354'] <- 'Kosakonia ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='4552af3e905206840cf8b67bcaf2693f'] <- 'Raoultella ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='3d3f89fadf28973358458a50d436e91e'] <- 'Methylobacterium-Methylorubrum ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='b37bf03859da8d004dbb9f4f4a220214'] <- 'Xanthobacteraceae ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='718981e648d6a3a51625c565106f6881'] <- 'Ochrobactrum ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='fc6622c636a5210293fb2873fc4761d9'] <- 'Allo.-Neo.-Para.-Rhizobium ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='607aee8b5f8abc3baf12a67738ff12d1'] <- 'Rhodobacteraceae ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='42c4ad8258e5ae1d878307dd458b5990'] <- 'Paracoccus ASV'



ggplot(mphyseqbigeff, aes(x=taxa_OTU,y=Abundance)) + geom_sina(aes(colour=surface_or_interior),scale="width",alpha=0.5) + stat_summary(aes(group=surface_or_interior,colour=surface_or_interior),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1,position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none"))+ scale_color_manual(labels = c("Fig supsensions", "Fig surface\nwashes"),values = c("#D5770C","#760D7D")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic",size=11),legend.title=element_blank()) + scale_y_continuous(breaks=c(-2:13)) + xlab("ASV") + ylab("Transformed abundance")


ggsave("supplemental_figure_ASV_surface_interior_8-4-2025.pdf",bg="white",height=8,width=15,units="in",useDingbats=FALSE)




#genus

sigdf_genus <- sigdf[sigdf$rank == "Genus",]

#how many?
nrow(sigdf_genus)
#[1] 12
#just plot them all!

sigdf_sort <- sigdf_genus[order(sigdf_genus$lfc_surface_or_interiorsurface),]

#how many diff abundant in interior

nrow(sigdf_genus[sigdf_genus$lfc_surface_or_interiorsurface < 0,])
#[1] 9

#how many diff abundant in surface
nrow(sigdf_genus[sigdf_genus$lfc_surface_or_interiorsurface > 0,])
#[1] 3

#top_five <- head(sigdf_sort,5)
#
#bottom_five <- tail(sigdf_sort,5)
#
genera_to_get <- c(top_five$taxon,bottom_five$taxon)

taxdf <- as.data.frame(tax_table(ps3_figs_only))

asv_to_get_taxonomy <- taxdf[taxdf$OTU %in% sigdf_genus$taxon, ]

#okay prep df


ps_genus <- tax_glom(ps3_figs_only, "Genus", NArm = TRUE)


ps3_clr <- microbiome::transform(ps_genus, "clr")



mphyseq = psmelt(ps3_clr)

mphyseq$Genus <- as.factor(mphyseq$Genus)

mphyseqbigeff <- mphyseq[mphyseq$Genus %in% sigdf_sort$taxon, ]

mphyseqbigeff$Genus <- factor(mphyseqbigeff$Genus, levels=sigdf_sort$taxon)
#                                                                           taxon
#13                                                                    Taibaiella
#17 Bacteria_Proteobacteria_Alphaproteobacteria_Rhodobacterales_Rhodobacteraceae_
#11                                                                      Gordonia
#21                                                                    Raoultella
#15                                                                  Ochrobactrum
#16    Bacteria_Proteobacteria_Alphaproteobacteria_Rhizobiales_Xanthobacteraceae_
#19                Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacterales__
#20                                                                     Kosakonia
#22                                                                   Xanthomonas
#12     Bacteria_Actinobacteriota_Actinobacteria_Micrococcales_Microbacteriaceae_
#18                                                                  Sphingomonas
#14                                                Methylobacterium-Methylorubrum




ggplot(mphyseqbigeff, aes(x=Genus,y=Abundance)) + geom_sina(aes(colour=surface_or_interior),scale="width",alpha=0.5) + stat_summary(aes(group=surface_or_interior,colour=surface_or_interior),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1,position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none"))+ scale_color_manual(labels = c("Fig supsensions", "Fig surface\nwashes"),values = c("#D5770C","#760D7D")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic",size=11),legend.title=element_blank()) + scale_y_continuous(breaks=c(-2:13)) + xlab("Genus") + ylab("Transformed abundance")


ggsave("supplemental_figure_Genus_surface_interior_8-4-2025.pdf",bg="white",height=8,width=15,units="in",useDingbats=FALSE)











#family

sigdf_family <- sigdf[sigdf$rank == "Family",]

#how many?
nrow(sigdf_family)
#[1] 12
#just plot them all!

sigdf_sort <- sigdf_family[order(sigdf_family$lfc_surface_or_interiorsurface),]

#how many diff abundant in interior

nrow(sigdf_family[sigdf_family$lfc_surface_or_interiorsurface < 0,])
#[1] 4

#how many diff abundant in surface
nrow(sigdf_family[sigdf_family$lfc_surface_or_interiorsurface > 0,])
#[1] 2

#top_five <- head(sigdf_sort,5)
#
#bottom_five <- tail(sigdf_sort,5)
#
#genera_to_get <- c(top_five$taxon,bottom_five$taxon)

#taxdf <- as.data.frame(tax_table(ps3_figs_only))

#asv_to_get_taxonomy <- taxdf[taxdf$OTU %in% sigdf_sort$taxon, ]

#okay prep df


ps_family <- tax_glom(ps3_figs_only, "Family", NArm = TRUE)


ps3_clr <- microbiome::transform(ps_family, "clr")



mphyseq = psmelt(ps3_clr)

mphyseq$Family <- as.factor(mphyseq$Family)

mphyseqbigeff <- mphyseq[mphyseq$Family %in% sigdf_sort$taxon, ]

mphyseqbigeff$Family <- factor(mphyseqbigeff$Family, levels=sigdf_sort$taxon)




ggplot(mphyseqbigeff, aes(x=Family,y=Abundance)) + geom_sina(aes(colour=surface_or_interior),scale="width",alpha=0.5) + stat_summary(aes(group=surface_or_interior,colour=surface_or_interior),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1,position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none"))+ scale_color_manual(labels = c("Fig supsensions", "Fig surface\nwashes"),values = c("#D5770C","#760D7D")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic",size=11),legend.title=element_blank()) + scale_y_continuous(breaks=c(-2:13)) + xlab("Family") + ylab("Transformed abundance")


ggsave("supplemental_figure_Family_surface_interior_8-4-2025.pdf",bg="white",height=8,width=15,units="in",useDingbats=FALSE)





#okay, worms....








resdf <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/ANCOMBC/ABC2_fig_interior_only_worms_present_all_ranks_results_res.tsv",header=TRUE,sep="\t")
resdf_sig <- resdf[resdf$diff_worms_presentyes == TRUE,]

resdf_ASV <- resdf_sig[resdf_sig$rank == "OTU",]

    #only ONE taxon robus (sig and passed sensitivity test)
#how many?
nrow(resdf_ASV)
#[1] 70


#how many diff abundant in figs without worms

nrow(resdf_ASV[resdf_ASV$lfc_worms_presentyes < 0,])
#[1] 35

#how many diff abundant in figs with worms
nrow(resdf_ASV[resdf_ASV$lfc_worms_presentyes > 0,])
#[1] 35
    #huh

sigdf_sort <- resdf_ASV[order(resdf_ASV$lfc_worms_presentyes),]

top_five <- head(sigdf_sort,5)
#
bottom_five <- tail(sigdf_sort,5)
#
asv_to_get <- c(top_five$taxon,bottom_five$taxon,"51d810a3fed597ec8151df6d57268336")




taxdf <- as.data.frame(tax_table(ps4_fig_interiors))

asv_to_get_taxonomy <- taxdf[taxdf$OTU %in% asv_to_get, ]

#okay prep df


#ps_family <- tax_glom(ps3, "Family", NArm = TRUE)


ps4_clr <- microbiome::transform(ps4_fig_interiors, "clr")



mphyseq = psmelt(ps4_clr)

mphyseq$taxa_OTU <- as.factor(mphyseq$taxa_OTU)

mphyseqbigeff <- mphyseq[mphyseq$taxa_OTU %in% asv_to_get, ]

mphyseqbigeff$taxa_OTU <- factor(mphyseqbigeff$taxa_OTU, levels=asv_to_get)

#                                                                              Genus
#a36a2714809a2e0c8f45edc9f26c67e6
#59a351402e79eb62104b5dbaf7c7af68                                          Spirosoma
#12dd6509c32da6a9937d21f5b52f07dc                                         env.OPS_17
#2d77cf80c9b03b82ce7ed6eacfa04e12                                   Chryseobacterium
#912f795ef3617aeb65ff9ad69b9567a8                                   Stenotrophomonas
#dee42ce3aa52ffbe2c350ee23b491242                                   Stenotrophomonas
#f50f7b988cd2a4793882818557fb0e69                                        Cystobacter
#749906e6079c81c5b2979a147e503684                                         Williamsia
#c50d36a6e593b5147afb76408cf3dc38
#4ba6bad752253e930876e8c688e4d9fb Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium

levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='a36a2714809a2e0c8f45edc9f26c67e6'] <- 'Unassigned ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='59a351402e79eb62104b5dbaf7c7af68'] <- 'Spirosoma ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='12dd6509c32da6a9937d21f5b52f07dc'] <- 'Sphingobacteriales ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='2d77cf80c9b03b82ce7ed6eacfa04e12'] <- 'Chryseobacterium'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='912f795ef3617aeb65ff9ad69b9567a8'] <- 'Stenotrophomonas ASV 1'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='dee42ce3aa52ffbe2c350ee23b491242'] <- 'Stenotrophomonas ASV 2'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='f50f7b988cd2a4793882818557fb0e69'] <- 'Cystobacter ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='749906e6079c81c5b2979a147e503684'] <- 'Williamsia ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='c50d36a6e593b5147afb76408cf3dc38'] <- 'Rickettsiaceae ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='4ba6bad752253e930876e8c688e4d9fb'] <- 'Allo.-Neo.-Para.-Rhizobium ASV'
levels(mphyseqbigeff$taxa_OTU)[levels(mphyseqbigeff$taxa_OTU)=='51d810a3fed597ec8151df6d57268336'] <- 'Ochrobactrum ASV'

ggplot(mphyseqbigeff, aes(x=taxa_OTU,y=Abundance)) + geom_sina(aes(colour=worms_present),scale="width") + stat_summary(aes(group=worms_present),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open(font_size =16) + background_grid(major = c("x"),minor = c("none")) + scale_color_manual(values = c("#fabb01", "#0081c7")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic",size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15)) + scale_y_continuous(breaks=c(-2:13)) + xlab("ASV") + ylab("Transformed abundance") +labs(colour="Worms\npresent?")


ggsave("supplemental_figure_ASV_fig_suspension_nematode_occupancy_9-10-2025.pdf",bg="white",height=8,width=15,units="in",useDingbats=FALSE)

