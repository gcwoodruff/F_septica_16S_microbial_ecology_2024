#Code for figures for associated with Woodruff et al. 2024, "16S metabarcoding reveals bacteria associated with a fig microcommunity"

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
library(ggmap)

#figure 1, map

#put google maps key and api secret in quotes
ggmap::register_google(key = "")
api_secret <- ''

	#"location_data.csv" is a truncated version of the sample metadata file ('fig_microbe_sample_taiwan_2019_metadata_d_4-23-22.csv')
loc_dat <- read.csv("location_data.csv", header=TRUE)
loc_dat$Caenorhabditis <- factor(loc_dat$Caenorhabditis, levels=c("yes", "no"))

#set styles for map
style1 <- c(feature = "poi", element = "labels", visibility = "off")
style2 <- c("&style=", feature = "administrative.locality", element = "labels", visibility = "off")

style <- c(style1, style2)

#get map with style and location

amap <- get_googlemap(center = c(121.543720, 25.011385), zoom = 13,style = style) 

#put data on map
ggmap(amap) + geom_point(data = loc_dat, mapping = aes(x = longitude, y = latitude, colour=Caenorhabditis), size=2) + scale_colour_manual(values=c("#0571b0","#ca0020")) + theme_void() + theme(legend.position="none")
#save map, this is figure 1
ggsave("taipei_2.png", height=4.2,width=4.2, units="in")


#figure 2


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


#remove mitochondrial and chloroplast reads

nochmi <- subset_taxa(ps.ng.tax, !Family %in% "Mitochondria" & !Order %in% "Chloroplast")


#subset the controls
controls <- prune_samples((sample_names(nochmi) %in% c("EXTRNEG","ExtrNeg1","EXTRPOS","GW10","GW11","GW31","GW9","PCRNeg1","PCRNeg2","PCRPos1","PCRPos2")), nochmi)


#remove control samples

no_controls <- prune_samples(!(sample_names(nochmi) %in% c("EXTRNEG","ExtrNeg1","EXTRPOS","GW10","GW11","GW31","GW9","PCRNeg1","PCRNeg2","PCRPos1","PCRPos2")), nochmi)



prevdf = apply(X = otu_table(controls),
               MARGIN = ifelse(taxa_are_rows(controls), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(controls),
                    tax_table(controls))


control_OTU <- prevdf[prevdf$Prevalence > 0,]

#remove OTU's found in controls and controls; remove control samples
ps3 <- subset_taxa(no_controls, !OTU %in% control_OTU$OTU)




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

#get factor levels right


plot2df$Family <- as.factor(plot2df$Family)

plot2df$Family <- factor(plot2df$Family, levels=c("Geodermatophilaceae", "Kineosporiaceae", "Microbacteriaceae", "Nocardiaceae", "Chitinophagaceae", "Hymenobacteraceae", "Sphingobacteriaceae", "Spirosomaceae", "Weeksellaceae", "Acetobacteraceae", "Beijerinckiaceae", "Rhizobiaceae", "Rhodobacteraceae", "Sphingomonadaceae", "Xanthobacteraceae", "Comamonadaceae", "Enterobacteriaceae", "Erwiniaceae", "Moraxellaceae", "Xanthomonadaceae"))


plot2df$Class <- as.factor(plot2df$Class)

plot2df$Class <- factor(plot2df$Class, levels=c("Actinobacteria","Bacteroidia","Alphaproteobacteria","Gammaproteobacteria"))


ggplot(plot2df,aes(y=Percent_prevalence,x=Family)) + geom_hline(yintercept=80,linetype="dotted") + geom_point(aes(shape=Group, colour=Class),size=2)  +theme_half_open(font_size =14) + background_grid(major = c("x"),minor = c("none")) + theme(axis.text.x = element_text(angle = 45, hjust=1,face="italic")) + scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) + scale_colour_manual("Class",values=c("yellow3","#7b3294","#bd0026","#3182bd")) + ylab("Prevalence (%)")


#ggplot(plot2df,aes(y=Percent_prevalence,x=Family)) + geom_hline(yintercept=80,linetype="dotted") + geom_point(aes(shape=Group, colour=Class),size=2)  + theme_cowplot(font_size=14) + theme(axis.text.x = element_text(angle = 45, hjust=1,face="italic")) + scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) + scale_colour_manual("Class",values=c("yellow3","#7b3294","#bd0026","#3182bd")) + ylab("Prevalence (%)")




#this is figure 2
#ggsave("prevalence_figure_v2.pdf",height=5.6,width=10,units="in")
ggsave("prevalence_figure_v3.pdf",height=5.6,width=10,units="in")








#Ordination with all samples , no chloro/mito

pslog <- transform_sample_counts(nochmi, function(x) log(1 + x))
#PCoA using Bray-Curtis dissimilarity
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
pslog <- transform_sample_counts(nochmi, function(x) log(1 + x))
#PCoA using Bray-Curtis dissimilarity
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]

ord_df <- plot_ordination(pslog, out.pcoa.log, color = "surface_or_interior", justDF = TRUE)

orddf1 <- ggplot(ord_df, aes(x=Axis.1,y=Axis.2)) + geom_point(aes(color = surface_or_interior),size=1) + coord_fixed(sqrt(evals[2] / evals[1])) + scale_color_manual(labels = c("Control", "Fig supsensions", "Fig surface washes"), values = c("#fde725", "#21918c","#440154")) +labs(colour="Sample Type") + theme_cowplot()  + xlab("Axis 1 (12%)") + ylab("Axis 2 (10.1%)")
    #this is figure 3B  

#making the composition plot for figure 3
ps1.com <- ps3

taxa_names(ps1.com) <- paste0("OTU_", rownames(tax_table(ps1.com)))

taxic <- as.data.frame(ps1.com@tax_table) 

taxic$OTU <- rownames(taxic) 
colnames(taxic)

taxmat <- as.matrix(taxic) 
new.tax <- tax_table(taxmat) 
tax_table(ps1.com) <- new.tax 

tax_table(ps1.com)[tax_table(ps1.com)[, "Family"] == "", "Family"] <- "Unclassified family"

guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 11,
  face = "italic", colour = "Black", angle = 0
)))
# first remove the phy_tree
ps1.com@phy_tree <- NULL
# Second merge at family level
ps1.com.fam <- aggregate_top_taxa2(ps1.com, "Family", top = 12)
#transform count
ps1.com.fam.ra = transform_sample_counts(ps1.com.fam, function(x){x / sum(x)})
#this is figure 3A
compplot1 <- plot_composition(ps1.com.fam.ra,sample.sort="surface_interior_worms",otu.sort="abundance") + theme(legend.position = "bottom") + scale_fill_manual("Family",values=c("gray36","#3182bd","gray69","#bdd7e7","#b30000","#fc8d59","#6baed6","#fef0d9","#7b3294","#eff3ff","#08519c","#ffff99","grey")) +  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + guide_italics + theme(legend.title = element_text(size = 13), axis.text.x = element_text(size=8)) + ylab("Relative Abundance")
#this is figure 3
compplot1/orddf1 


#figure 4a


#just get interior samples
ps4 <- subset_samples(ps3, surface_or_interior == "interior")


ps1.com <- ps4

taxa_names(ps1.com) <- paste0("OTU_", rownames(tax_table(ps1.com)))

taxic <- as.data.frame(ps1.com@tax_table) # this will help in setting large color options


taxic$OTU <- rownames(taxic) 
colnames(taxic) 

taxmat <- as.matrix(taxic)
new.tax <- tax_table(taxmat) 
tax_table(ps1.com) <- new.tax

# now edit the unclassified taxa
tax_table(ps1.com)[tax_table(ps1.com)[, "Family"] == "", "Family"] <- "Unclassified family"


guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 11,
  face = "italic", colour = "Black", angle = 0
)))

ps1.com@phy_tree <- NULL


ps1.com.fam <- aggregate_top_taxa2(ps1.com, "Family", top = 12)


ps1.com.fam.ra = transform_sample_counts(ps1.com.fam, function(x){x / sum(x)})

plot_composition(ps1.com.fam.ra,sample.sort="surface_interior_worms",otu.sort="abundance") + theme(legend.position = "bottom") + scale_fill_manual("Family",values=c("#3182bd","gray36","gray69","#bdd7e7","#b30000","#fc8d59","#eff3ff","#7b3294","#08519c","#6baed6","#fdcc8a","#fef0d9","#d7b5d8","grey")) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + guide_italics + theme(legend.title = element_text(size = 15)) + ylab("Relative Abundance")

#this is figure 4a
ggsave("worms_present_family_composition_2024.pdf", height=6,width=6, units="in")


#figure 4b
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

  #MDS/PCoA Performs principal coordinate analysis (also called
              #principle coordinate decomposition, multidimensional
              #scaling (MDS), or classical scaling) of a distance matrix
              #(Gower 1966), including two correction methods for
              #negative eigenvalues.  See ‘pcoa’ for further details.

evals <- out.pcoa.log$values[,1]

ord_df <-plot_ordination(pslog, out.pcoa.log, color = "worms_present", justDF = TRUE)

ggplot(ord_df, aes(x=Axis.1,y=Axis.2)) + geom_point(aes(color = worms_present),size=1) + coord_fixed(sqrt(evals[2] / evals[1])) + scale_color_manual(values = c("#fabb01", "#0081c7")) +labs(colour="Worms present?") + theme_cowplot()  + xlab("Axis 1 (19.4%)") + ylab("Axis 2 (8.1%)")
#this is figure 4b


#figure 5


#diversity metrics
#worms present/absent

    #remove samples with very low read counts
ps3_no_small <- prune_samples((!sample_names(ps4) %in% c("GW20","GW34")), ps4)

    #rarefy
ps_rare <- phyloseq::rarefy_even_depth(ps3_no_small, rngseed = 123, replace = FALSE)

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

#this is figure 5
ggsave("worms_present_diversity_sina_2024.pdf", height=4,width=7, units="in",bg="white",useDingbats=FALSE)




#figure 6


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

#figure 6b
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
