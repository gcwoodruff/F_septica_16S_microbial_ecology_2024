
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
setwd("/Users/gavin/genome/genome_old/16S_enrivonmental_fig_microbe_Illumina_metabarcoding_2022/phyloseq_attempt_4-23-22/pub_prep/")

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


#remove mitochondrial and chloroplast reads

nochmi <- subset_taxa(ps.ng.tax, !Family %in% "Mitochondria" & !Order %in% "Chloroplast")


#remove control OTUs


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


#remove control samples

no_controls <- prune_samples(!(sample_names(nochmi) %in% c("EXTRNEG","ExtrNeg1","EXTRPOS","GW10","GW11","GW31","GW9","PCRNeg1","PCRNeg2","PCRPos1","PCRPos2")), nochmi)

#remove OTU's found in controls and controls; remove control samples
ps3 <- subset_taxa(no_controls, !OTU %in% control_OTU$OTU)


#just get fig interiors (suspensions)


ps4 <- subset_samples(ps3, surface_or_interior == "interior")

adiv_no_rare <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps4, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps4, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps4)))), tree = phyloseq::phy_tree(ps4))[, 1],
  "foundress_number" = phyloseq::sample_data(ps4)$foundress_number)


ps_no_small <- prune_samples((!sample_names(ps4) %in% c("GW20","GW34")), ps4)


ps_rare <- phyloseq::rarefy_even_depth(ps_no_small, rngseed = 123, replace = FALSE)


adiv_with_rare <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps_rare)))), tree = phyloseq::phy_tree(ps_rare))[, 1],
  "foundress_number" = phyloseq::sample_data(ps_rare)$foundress_number)

summary(adiv_with_rare$Observed)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   37.0    72.0   102.0   158.4   201.0   651.0

summary(adiv_with_rare$Shannon)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.002   2.227   2.914   2.938   3.563   4.828
summary(adiv_with_rare$Phylogenetic)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  2.438   5.033   6.941   8.371  10.504  24.756




summary(adiv_no_rare$Observed)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  37.00   73.25  104.00  167.37  209.00  704.00
summary(adiv_no_rare$Shannon)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.008   2.274   2.895   2.933   3.567   4.858

summary(adiv_no_rare$Phylogenetic)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  2.438   4.897   6.770   8.774  11.178  27.317






#figure 5, no rarefaction


adiv <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps4, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps4, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps4)))), tree = phyloseq::phy_tree(ps4))[, 1],
  "worms_present" = phyloseq::sample_data(ps4)$worms_present)

adiv$sample_id <- rownames(adiv)
names(adiv)[names(adiv) == "Observed"] <- "Number_of_OTU"

adiv_melt <- reshape2::melt(adiv, id.vars= c("sample_id","worms_present"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic"))

levels(adiv_melt$variable)[levels(adiv_melt$variable)=="Number_of_OTU"] <- "Number of ASV"

worm_occupancy_plot_df <- adiv_melt

ggplot(adiv_melt, aes(x=worms_present,y=value)) + geom_sina(scale="width") + stat_summary(aes(group=worms_present),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + facet_wrap(~variable,nrow=1,scales="free") +theme_cowplot() + scale_y_continuous(limits = c(0,NA)) + xlab("Caenorhabditis nematodes present?") +ylab("Diversity measure") + theme(strip.background = element_rect(colour="white", fill="white"))


a <- ggplot(worm_occupancy_plot_df, aes(x=worms_present,y=value)) + geom_sina(scale="width") + stat_summary(aes(group=worms_present),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + facet_wrap(~variable,nrow=1,scales="free") +theme_cowplot() + scale_y_continuous(limits = c(0,NA)) + xlab("Caenorhabditis nematodes present?") +ylab("Diversity measure") + theme(strip.background = element_rect(colour="white", fill="white"))


#this is figure 5
ggsave("/Users/gavin/genome/septica_microbiome_revisions/5-2025/new_figures_7-29-25/worms_present_diversity_sina_all_reads_2025.pdf", height=4,width=7, units="in",bg="white",useDingbats=FALSE)






adiv <- data.frame(
  "Number_of_OTU" = phyloseq::estimate_richness(ps4, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps4, measures = "Shannon"),
  "Phylogenetic" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps4)))), tree = phyloseq::phy_tree(ps4))[, 1],
  "foundress_number" = phyloseq::sample_data(ps4)$foundress_number)

adiv$sample_id <- rownames(adiv)
names(adiv)[names(adiv) == "Observed"] <- "Number_of_OTU"

adiv_melt <- reshape2::melt(adiv, id.vars= c("sample_id","foundress_number"),measure.vars = c("Number_of_OTU", "Shannon","Phylogenetic"))

levels(adiv_melt$variable)[levels(adiv_melt$variable)=="Number_of_OTU"] <- "Number of ASV"

str(adiv_melt)

adiv_melt$foundress_number <- as.numeric(adiv_melt$foundress_number)

foundress_number_plot_df <- adiv_melt

ggplot(adiv_melt, aes(x=foundress_number,y=value)) + geom_point() + facet_wrap(~variable,nrow=1,scales="free") +theme_cowplot() + stat_smooth() + scale_y_continuous(limits = c(0,NA)) + xlab("Foundress number") +ylab("Diversity measure") + theme(strip.background = element_rect(colour="white", fill="white")) + scale_y_continuous(oob=scales::rescale_none)
    #this is figure 6

b <- ggplot(adiv_melt, aes(x=foundress_number,y=value)) + geom_point() + facet_wrap(~variable,nrow=1,scales="free") +theme_cowplot() + stat_smooth() + scale_y_continuous(limits = c(0,NA)) + xlab("Foundress number") +ylab("Diversity measure") + theme(strip.background = element_rect(colour="white", fill="white")) + scale_y_continuous(oob=scales::rescale_none)

ggsave("/Users/gavin/genome/septica_microbiome_revisions/5-2025/new_figures_7-29-25/foundress_number_diversity_sina_all_reads_2025.pdf", height=4,width=7, units="in",bg="white",useDingbats=FALSE)



#nonlinear models, no rarefaction

adiv$foundress_number <- as.numeric(adiv$foundress_number)


asymmodNumber_of_OTU <- nls(Number_of_OTU ~ SSasymp(foundress_number, Asym, R0, lrc), data = adiv)

summary(asymmodNumber_of_OTU)

#Parameters:
#      Estimate Std. Error t value Pr(>|t|)
#Asym  95.36164   41.91796   2.275   0.0291 *
#R0   285.24411   39.47957   7.225 1.96e-08 ***
#lrc   -0.04202    0.86321  -0.049   0.9615
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 131.5 on 35 degrees of freedom
#
#Number of iterations to convergence: 14
#Achieved convergence tolerance: 9.115e-06





adiv_melt_just_OTU <- adiv_melt[adiv_melt$variable == "Number of ASV",]

ggplot(adiv_melt_just_OTU, aes(x=foundress_number,y=value)) + geom_point() +theme_cowplot() + stat_smooth() + scale_y_continuous(limits = c(0,NA)) + xlab("Foundress number") +ylab("Diversity measure") + theme(strip.background = element_rect(colour="white", fill="white")) + scale_y_continuous(oob=scales::rescale_none)

pred_x <- seq(0, 20, length.out = 37)

#pred_y <- predict(asymmodNumber_of_OTU, pred_x)
    #not sure what happened here
    #Asym+(R0-Asym)*exp(-exp(lrc)*input)
pred_y2 <- (95.36164+(285.24411-95.36164)*exp(-exp(-0.04202)*pred_x))
    #this is the function R says SSasymp is trying to fit, I just plugged in the fitted parameters

predicted_data <- data.frame(x = pred_x,y=pred_y2)
observed_data <- data.frame(x=adiv_melt_just_OTU$foundress_number, y=adiv_melt_just_OTU$value)

se_upper <- ((95.36164+41.91796)+((285.24411+39.47957)-(95.36164+41.91796))*exp(-exp((-0.04202))*pred_x))

se_lower <- ((95.36164-41.91796)+((285.24411-39.47957)-(95.36164-41.91796))*exp(-exp((-0.04202))*pred_x))

se_high_df <- data.frame(x = pred_x,y=se_upper)
se_low_df <- data.frame(x = pred_x,y=se_lower)

ggplot(predicted_data, aes(x=x,y=y)) + geom_line(colour="blue") + geom_line(data=se_high_df,colour="blue",linetype="dotted") + geom_line(data=se_low_df,colour="blue",linetype="dotted") + geom_point(data = observed_data) +theme_cowplot() + xlab("Foundress number") + ylab("Number of ASV")

c <- ggplot(predicted_data, aes(x=x,y=y)) + geom_line(colour="blue") + geom_line(data=se_high_df,colour="blue",linetype="dotted") + geom_line(data=se_low_df,colour="blue",linetype="dotted") + geom_point(data = observed_data) +theme_cowplot() + xlab("Foundress number") + ylab("Number of ASV")


    #okay that *looks* reasonable/not insane



asymmodShannon <- nls(Shannon ~ SSasymp(foundress_number, Asym, R0, lrc), data = adiv)

summary(asymmodShannon)
#Error in nls(y ~ cbind(1 - exp(-exp(lrc) * x), exp(-exp(lrc) * x)), data = xy,  :
#  step factor 0.000488281 reduced below 'minFactor' of 0.000976562



adiv_melt_just_Shannon <- adiv_melt[adiv_melt$variable == "Shannon",]


ggplot(adiv_melt_just_Shannon, aes(x=foundress_number,y=value)) + geom_point() +theme_cowplot() + xlab("Foundress number") + ylab("Shannon diversity")


d <- ggplot(adiv_melt_just_Shannon, aes(x=foundress_number,y=value)) + geom_point() +theme_cowplot() + xlab("Foundress number") + ylab("Shannon diversity") + scale_x_continuous(limits=c(0,20),breaks=c(0,5,10,15,20))


asymmodPhylogenetic <- nls(Phylogenetic ~ SSasymp(foundress_number, Asym, R0, lrc), data = adiv)

summary(asymmodPhylogenetic)

#Formula: Phylogenetic ~ SSasymp(foundress_number, Asym, R0, lrc)
#
#Parameters:
#     Estimate Std. Error t value Pr(>|t|)
#Asym   5.3571     1.7643   3.036   0.0045 **
#R0    13.0918     1.3641   9.597 2.46e-11 ***
#lrc   -0.4439     0.6928  -0.641   0.5259
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 4.583 on 35 degrees of freedom
#
#Number of iterations to convergence: 4
#Achieved convergence tolerance: 3.559e-06
adiv_melt_just_Phylogenetic <- adiv_melt[adiv_melt$variable == "Phylogenetic",]

ggplot(adiv_melt_just_Phylogenetic, aes(x=foundress_number,y=value)) + geom_point() +theme_cowplot() + stat_smooth() + scale_y_continuous(limits = c(0,NA)) + xlab("Foundress number") +ylab("Diversity measure") + theme(strip.background = element_rect(colour="white", fill="white")) + scale_y_continuous(oob=scales::rescale_none)

pred_x <- seq(0, 20, length.out = 37)

#pred_y <- predict(asymmodNumber_of_OTU, pred_x)
    #not sure what happened here

pred_y2 <- (5.3571+(13.0918-5.3571)*exp(-exp(-0.4439)*pred_x))
    #this is the function R says SSasymp is trying to fit, I just plugged in the fitted parameters

predicted_data <- data.frame(x = pred_x,y=pred_y2)
observed_data <- data.frame(x=adiv_melt_just_Phylogenetic$foundress_number, y=adiv_melt_just_Phylogenetic$value)


#Parameters:
#     Estimate Std. Error t value Pr(>|t|)
#Asym   5.3571     1.7643   3.036   0.0045 **
#R0    13.0918     1.3641   9.597 2.46e-11 ***
#lrc   -0.4439     0.6928  -0.641   0.5259
#---

se_upper <- ((5.3571+1.7643)+((13.0918+1.3641)-(5.3571+1.7643))*exp(-exp(-0.4439)*pred_x))

se_lower <- ((5.3571-1.7643)+((13.0918-1.3641)-(5.3571-1.7643))*exp(-exp(-0.4439)*pred_x))


se_high_df <- data.frame(x = pred_x,y=se_upper)
se_low_df <- data.frame(x = pred_x,y=se_lower)



ggplot(predicted_data, aes(x=x,y=y)) + geom_line(colour="blue") + geom_line(data=se_high_df,colour="blue",linetype="dotted") + geom_line(data=se_low_df,colour="blue",linetype="dotted") + geom_point(data = observed_data) +theme_cowplot() + xlab("Foundress number") + ylab("Phylogenetic diversity")

e <- ggplot(predicted_data, aes(x=x,y=y)) + geom_line(colour="blue") + geom_line(data=se_high_df,colour="blue",linetype="dotted") + geom_line(data=se_low_df,colour="blue",linetype="dotted") + geom_point(data = observed_data) +theme_cowplot() + xlab("Foundress number") + ylab("Phylogenetic diversity") +scale_y_continuous()

    #okay that *looks* reasonable/not insane

a/(c+d+e)

ggsave("/Users/gavin/genome/septica_microbiome_revisions/5-2025/new_figures_7-29-25/nematodes_foundress_number_diversity_composite_figure_5_7-31-2025.pdf", height=6,width=8, units="in",bg="white",useDingbats=FALSE)

summary(lm(Number_of_OTU ~ foundress_number, data=adiv))

#Call:
#lm(formula = Number_of_OTU ~ foundress_number, data = adiv)
#
#Residuals:
#    Min      1Q  Median      3Q     Max
#-154.68  -94.43  -39.69   53.45  497.82
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)       206.179     28.544   7.223 1.68e-08 ***
#foundress_number  -12.498      5.405  -2.313   0.0266 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 142.3 on 36 degrees of freedom
#Multiple R-squared:  0.1293,    Adjusted R-squared:  0.1052
#F-statistic: 5.348 on 1 and 36 DF,  p-value: 0.02658

summary(lm(Shannon ~ foundress_number, data=adiv))

#summary(lm(Shannon ~ foundress_number, data=adiv))
#
#Call:
#lm(formula = Shannon ~ foundress_number, data = adiv)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-1.99060 -0.61015 -0.00647  0.45610  1.74162
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)       3.11679    0.19306  16.144   <2e-16 ***
#foundress_number -0.05914    0.03655  -1.618    0.114
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.9627 on 36 degrees of freedom
#Multiple R-squared:  0.06777,   Adjusted R-squared:  0.04188
#F-statistic: 2.617 on 1 and 36 DF,  p-value: 0.1144

summary(lm(Phylogenetic ~ foundress_number, data=adiv))


#summary(lm(Phylogenetic ~ foundress_number, data=adiv))
#
#Call:
#lm(formula = Phylogenetic ~ foundress_number, data = adiv)
#
#Residuals:
#    Min      1Q  Median      3Q     Max
#-6.4497 -3.6458 -0.8723  2.5423 17.0741
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)       10.2428     1.0141  10.100 4.75e-12 ***
#foundress_number  -0.4732     0.1920  -2.464   0.0186 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 5.057 on 36 degrees of freedom
#Multiple R-squared:  0.1443,    Adjusted R-squared:  0.1206
#F-statistic: 6.073 on 1 and 36 DF,  p-value: 0.01864

