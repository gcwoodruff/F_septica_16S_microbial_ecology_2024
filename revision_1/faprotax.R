
library(ggplot2)
library(ggforce)
library(cowplot)

woodruff_fig_surface_samples <- c("GW12", "GW16", "GW18", "GW20", "GW24", "GW26", "GW28", "GW30", "GW33", "GW35", "GW37", "GW39", "GW40", "GW41", "GW44", "GW46", "GW49", "GW5", "GW53", "GW54", "GW56", "GW57", "GW59", "GW6", "GW60", "GW62", "GW67", "GW69", "GW7", "GW71", "GW72", "GW8")

woodruff_fig_interior_samples <- c("GW1", "GW13", "GW14", "GW15", "GW17", "GW19", "GW2", "GW21", "GW22", "GW23", "GW25", "GW27", "GW29", "GW3", "GW32", "GW34", "GW36", "GW38", "GW4", "GW42", "GW43", "GW45", "GW47", "GW48", "GW50", "GW51", "GW52", "GW55", "GW58", "GW61", "GW63", "GW64", "GW65", "GW66", "GW68", "GW70", "GW73", "GW74")

samuel_substrate_samples <- c("SRR5094338", "SRR5094346", "SRR5094349", "SRR5094351", "SRR5094352", "SRR5094356", "SRR5094361", "SRR5094362", "SRR5094366", "SRR5094370", "SRR5094374", "SRR5094376", "SRR5094384", "SRR5094388", "SRR5094389", "SRR5094400", "SRR5094401", "SRR5094404", "SRR5094405", "SRR5094408", "SRR5094415", "SRR5094419", "SRR5094420", "SRR5094421", "SRR5094424", "SRR5094428", "SRR5094432", "SRR5094436", "SRR5094437", "SRR5094438", "SRR5094440", "SRR5094443", "SRR5094445", "SRR5094449", "SRR5094460", "SRR5094469", "SRR5094471", "SRR5094472", "SRR5094475", "SRR5094478", "SRR5094479", "SRR5094480", "SRR5094485", "SRR5094489", "SRR5094515", "SRR5094517", "SRR5094525", "SRR5094528", "SRR5094531", "SRR5094532", "SRR5094534", "SRR5094536", "SRR5094551", "SRR5094552", "SRR5094554", "SRR5094565", "SRR5094566", "SRR5094569", "SRR5094574", "SRR5094576")

dirksen_substrate_samples <- c("ERR1307271", "ERR1307272", "ERR1307273", "ERR1307274", "ERR1307275", "ERR1307276", "ERR1307277", "ERR1307278", "ERR1307279", "ERR1307280", "ERR1307281", "ERR1307282", "ERR1307283", "ERR1307284", "ERR1307285", "ERR1307286", "ERR1307287", "ERR1307288", "ERR1307290", "ERR1307291", "ERR1307292", "ERR1307293", "ERR1307294", "ERR1307295", "ERR1307296", "ERR1307297", "ERR1307298", "ERR1307299", "ERR1307300", "ERR1307301", "ERR1307302", "ERR1307303", "ERR1307304", "ERR1307305", "ERR1307306", "ERR1307307", "ERR1307308", "ERR1307309", "ERR1307310", "ERR1307311", "ERR1307312", "ERR1307313", "ERR1307314", "ERR1307315", "ERR1307316", "ERR1307317", "ERR1307318", "ERR1307319", "ERR1307320", "ERR1307321", "ERR1307322", "ERR1307323", "ERR1307324", "ERR1307325", "ERR1307326", "ERR1307327", "ERR1307328", "ERR1307329", "ERR1307330", "ERR1307331", "ERR1307332", "ERR1307333", "ERR1307334", "ERR1307335", "ERR1307336", "ERR1307337", "ERR1307338", "ERR1307339", "ERR1307340", "ERR1307341", "ERR1307342", "ERR1307343", "ERR1307344", "ERR1307345", "ERR1307346", "ERR1307347", "ERR1307348", "ERR1307349", "ERR1307350", "ERR1307351", "ERR1307352", "ERR1307353", "ERR1307354", "ERR1307355", "ERR1307356", "ERR1307357", "ERR1307358", "ERR1307359", "ERR1307360", "ERR1307361", "ERR1307362", "ERR1307363", "ERR1307364", "ERR1307365", "ERR1307366", "ERR1307367", "ERR1307368", "ERR1307369", "ERR1307370", "ERR1307371", "ERR1307372", "ERR1307373", "ERR1307374", "ERR1307375", "ERR1307376", "ERR1307377", "ERR1307378", "ERR1307379", "ERR1307380", "ERR1307381", "ERR1307382", "ERR1307383", "ERR1307384", "ERR1307385", "ERR1307386", "ERR1307387", "ERR1307388", "ERR1307389", "ERR1307390", "ERR1307391", "ERR1307392", "ERR1307393", "ERR1307394", "ERR1307395", "ERR1307396", "ERR1307397", "ERR1307398", "ERR1307399", "ERR1307400", "ERR1307401", "ERR1307402", "ERR1307403", "ERR1307404", "ERR1307405", "ERR1307406", "ERR1307407", "ERR1307408", "ERR1307409", "ERR1307410", "ERR1307411", "ERR1307412", "ERR1307413", "ERR1307414", "ERR1307415", "ERR1307416", "ERR1307417", "ERR1307418", "ERR1307419", "ERR1307420", "ERR1307421", "ERR1307422", "ERR1307423", "ERR1307424", "ERR1307425", "ERR1307426", "ERR1307427", "ERR1307428", "ERR1307429", "ERR1307430", "ERR1307431", "ERR1307432", "ERR1307433", "ERR1307434", "ERR1307435", "ERR1307436", "ERR1307437", "ERR1307438", "ERR1307439", "ERR1307440", "ERR1307441", "ERR1307442", "ERR1307443", "ERR1307444", "ERR1307445", "ERR1307446", "ERR1307447", "ERR1307448", "ERR1307449", "ERR1307450", "ERR1307451", "ERR1354120")

interior_no_caeno_samples <-c("GW1", "GW14", "GW15", "GW22", "GW32", "GW34", "GW36", "GW4", "GW43", "GW48", "GW50", "GW51", "GW52", "GW55", "GW64", "GW65", "GW66", "GW73", "GW74")

interior_yes_caeno_samples <- c("GW38", "GW45", "GW13", "GW17", "GW19", "GW2", "GW21", "GW23", "GW25", "GW27", "GW29", "GW3", "GW42", "GW47", "GW58", "GW61", "GW63", "GW68", "GW70")

dat <- read.table("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/faprotax_functional_table_w_elegans.tsv",sep="\t",header=TRUE)

dat_melt <- reshape2::melt(dat)

fig_surface_df <- dat_melt[dat_melt$variable %in% woodruff_fig_surface_samples, ]
fig_interior_df <- dat_melt[dat_melt$variable %in% woodruff_fig_interior_samples, ]
samuel_df <- dat_melt[dat_melt$variable %in% samuel_substrate_samples, ]
dirksen_df <- dat_melt[dat_melt$variable %in% dirksen_substrate_samples, ]

fig_surface_df$species <- "C. inopinata"
fig_interior_df$species <- "C. inopinata"
samuel_df$species <- "C. elegans"
dirksen_df$species <- "C. elegans"


fig_surface_df$category <- "Fig surface wash"
fig_interior_df$category <- "Fig suspension"
samuel_df$category <- "Samuel et al. substrate"
dirksen_df$category <- "Dirksen et al. substrate"

plotdf <- rbind(fig_surface_df,fig_interior_df,samuel_df,dirksen_df)


ggplot(plotdf, aes(x=value,y=group)) + geom_sina(aes(colour=species),scale="width",alpha=0.5) + stat_summary(aes(group=species),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open() + background_grid(major = c("y"),minor = c("none"))

gr1 <- c("acetoclastic_methanogenesis", "aerobic_ammonia_oxidation", "aerobic_anoxygenic_phototrophy", "aerobic_chemoheterotrophy", "aerobic_nitrite_oxidation", "aliphatic_non_methane_hydrocarbon_degradation", "anammox", "animal_parasites_or_symbionts", "anoxygenic_photoautotrophy", "anoxygenic_photoautotrophy_Fe_oxidizing", "anoxygenic_photoautotrophy_H2_oxidizing", "anoxygenic_photoautotrophy_S_oxidizing", "aromatic_compound_degradation", "aromatic_hydrocarbon_degradation", "arsenate_detoxification", "arsenate_respiration", "arsenite_oxidation_detoxification", "arsenite_oxidation_energy_yielding", "cellulolysis", "chemoheterotrophy")

gr2 <- c("chitinolysis", "chlorate_reducers", "chloroplasts", "dark_hydrogen_oxidation", "dark_iron_oxidation", "dark_oxidation_of_sulfur_compounds", "dark_sulfide_oxidation", "dark_sulfite_oxidation", "dark_sulfur_oxidation", "dark_thiosulfate_oxidation", "denitrification", "dissimilatory_arsenate_reduction", "dissimilatory_arsenite_oxidation", "fermentation", "fish_parasites", "fumarate_respiration", "human_associated", "human_gut", "human_pathogens_all", "human_pathogens_diarrhea")

gr3 <- c("human_pathogens_gastroenteritis", "human_pathogens_meningitis", "human_pathogens_nosocomia", "human_pathogens_pneumonia", "human_pathogens_septicemia", "hydrocarbon_degradation", "hydrogenotrophic_methanogenesis", "intracellular_parasites", "invertebrate_parasites", "iron_respiration", "knallgas_bacteria", "ligninolysis", "mammal_gut", "manganese_oxidation", "manganese_respiration", "methanogenesis", "methanogenesis_by_CO2_reduction_with_H2", "methanogenesis_by_disproportionation_of_methyl_groups", "methanogenesis_by_reduction_of_methyl_compounds_with_H2", "methanogenesis_using_formate")

gr4 <- c("methanol_oxidation", "methanotrophy", "methylotrophy", "nitrate_ammonification", "nitrate_denitrification", "nitrate_reduction", "nitrate_respiration", "nitrification", "nitrite_ammonification", "nitrite_denitrification", "nitrite_respiration", "nitrogen_fixation", "nitrogen_respiration", "nitrous_oxide_denitrification", "nonphotosynthetic_cyanobacteria", "oil_bioremediation", "oxygenic_photoautotrophy", "photoautotrophy", "photoheterotrophy", "photosynthetic_cyanobacteria")

gr5 <- c("phototrophy", "plant_pathogen", "plastic_degradation", "predatory_or_exoparasitic", "reductive_acetogenesis", "respiration_of_sulfur_compounds", "sulfate_respiration", "sulfite_respiration", "sulfur_respiration", "thiosulfate_respiration", "ureolysis", "xylanolysis")

plotdf1 <- plotdf[plotdf$group %in% gr1, ]
plotdf2 <- plotdf[plotdf$group %in% gr2, ]
plotdf3 <- plotdf[plotdf$group %in% gr3, ]
plotdf4 <- plotdf[plotdf$group %in% gr4, ]
plotdf5 <- plotdf[plotdf$group %in% gr5, ]


ggplot(plotdf1, aes(x=value,y=group)) + geom_sina(aes(colour=species),scale="width",alpha=0.5) + stat_summary(aes(group=species),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open() + background_grid(major = c("y"),minor = c("none"))


ggplot(plotdf2, aes(x=value,y=group)) + geom_sina(aes(colour=species),scale="width",alpha=0.5) + stat_summary(aes(group=species),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open() + background_grid(major = c("y"),minor = c("none"))


ggplot(plotdf3, aes(x=value,y=group)) + geom_sina(aes(colour=species),scale="width",alpha=0.5) + stat_summary(aes(group=species),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open() + background_grid(major = c("y"),minor = c("none"))


ggplot(plotdf4, aes(x=value,y=group)) + geom_sina(aes(colour=species),scale="width",alpha=0.5) + stat_summary(aes(group=species),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open() + background_grid(major = c("y"),minor = c("none"))


ggplot(plotdf5, aes(x=value,y=group)) + geom_sina(aes(colour=species),scale="width",alpha=0.5) + stat_summary(aes(group=species),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open() + background_grid(major = c("y"),minor = c("none"))


plotdf$group <-as.factor(plotdf$group)
plotdf$category <-as.factor(plotdf$category)
plotdf$species <-as.factor(plotdf$species)

effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL
group_id <- NULL

library(effsize)

for (i in levels(plotdf$group)){
 dat <- plotdf[plotdf$group == i,]
 inop <- dat[dat$species == "C. inopinata",]
 eleg <- dat[dat$species == "C. elegans",]
 group_id <- rbind(group_id, i)
 effect_size_lower <- rbind(effect_size_lower, cohen.d(inop$value,eleg$value)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(inop$value,eleg$value)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(inop$value,eleg$value)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(inop$value,eleg$value)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(inop$value,eleg$value)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(inop$value,eleg$value)$p.value)
}

stat_df <- as.data.frame(cbind(group_id,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p))
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("group","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)
#correct for multiple tests
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

write.table(stat_df, "/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/wilcox_tests_faprotax_groups_species.tsv",sep="\t",row.names=FALSE, quote=FALSE)


sigcat <- subset(stat_df, p.adj < 0.05)

plotsig <- plotdf[plotdf$group %in% sigcat$group, ]

plotsig2 <- plotsig[plotsig$group != "chemoheterotrophy",]

ggplot(plotsig2, aes(x=value,y=group)) + geom_sina(aes(colour=species),scale="width",alpha=0.25) + stat_summary(aes(colour=species),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1,position = position_dodge(width = 0.9)) +theme_half_open() + background_grid(major = c("y"),minor = c("none")) + xlim(0,0.1) + scale_colour_brewer(palette="Set1") + theme(legend.text = element_text(face ="italic")) + xlab("ASV Proportion")

ggsave("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/faprotax_sig_eleg_inop.pdf",useDingbats=FALSE,unit="in",height=7,width=10)


write.table(sigcat, "/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/wilcox_tests_faprotax_groups_species_significant.tsv",sep="\t",row.names=FALSE, quote=FALSE)


intextplotdf <- rbind(fig_surface_df,fig_interior_df)














intextplotdf$group <-as.factor(intextplotdf$group)
intextplotdf$category <-as.factor(intextplotdf$category)
intextplotdf$species <-as.factor(intextplotdf$species)

effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL
group_id <- NULL

library(effsize)

for (i in levels(intextplotdf$group)){
 dat <- intextplotdf[intextplotdf$group == i,]
 int <- dat[dat$category == "Fig suspension",]
 ext <- dat[dat$category == "Fig surface wash",]
 group_id <- rbind(group_id, i)
 effect_size_lower <- rbind(effect_size_lower, cohen.d(int$value,ext$value)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(int$value,ext$value)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(int$value,ext$value)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(int$value,ext$value)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(int$value,ext$value)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(int$value,ext$value)$p.value)
}

stat_df <- as.data.frame(cbind(group_id,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p))
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("group","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)
#correct for multiple tests
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

write.table(stat_df, "/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/wilcox_tests_faprotax_groups_figs_only_interior_v_exterior.tsv",sep="\t",row.names=FALSE, quote=FALSE)




sigcat <- subset(stat_df, p.adj < 0.05)

plotsig <- intextplotdf[intextplotdf$group %in% sigcat$group, ]

ggplot(plotsig, aes(x=value,y=group)) + geom_sina(aes(colour=category),scale="width",alpha=0.5) + stat_summary(aes(group=category),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open() + background_grid(major = c("y"),minor = c("none")) + scale_color_manual(values = c("#D5770C","#760D7D")) + xlab("ASV Proportion")


ggsave("faprotax_sig_figs_only_interior_exterior.pdf",useDingbats=FALSE,unit="in",height=7,width=10)




with_worms <- dat_melt[dat_melt$variable %in% interior_yes_caeno_samples, ]
without_worms <- dat_melt[dat_melt$variable %in% interior_no_caeno_samples, ]

with_worms$worms <- "Yes"
without_worms$worms <- "No"

wormsdf <- rbind(with_worms,without_worms)








wormsdf$group <-as.factor(wormsdf$group)
wormsdf$worms <- as.factor(wormsdf$worms)

effect_size_lower <- NULL
effect_size <- NULL
effect_size_upper <- NULL
magn <- NULL
wilcox_stat <- NULL
wilcox_p <- NULL
group_id <- NULL


for (i in levels(wormsdf$group)){
 dat <- wormsdf[wormsdf$group == i,]
 withworms <- dat[dat$worms == "Yes",]
 withoutworms <- dat[dat$worms == "No",]
 group_id <- rbind(group_id, i)
 effect_size_lower <- rbind(effect_size_lower, cohen.d(withworms$value,withoutworms$value)$conf.int[1]) 
 effect_size <- rbind(effect_size, cohen.d(withworms$value,withoutworms$value)$estimate)
 effect_size_upper <- rbind(effect_size_upper, cohen.d(withworms$value,withoutworms$value)$conf.int[2]) 
 magn <- rbind(magn, cohen.d(withworms$value,withoutworms$value)$magnitude) 
 wilcox_stat <- rbind(wilcox_stat, wilcox.test(withworms$value,withoutworms$value)$statistic) 
 wilcox_p <- rbind(wilcox_p, wilcox.test(withworms$value,withoutworms$value)$p.value)
}

stat_df <- as.data.frame(cbind(group_id,effect_size_lower,effect_size,effect_size_upper,magn,wilcox_stat,wilcox_p))
  #getting data right, column names
rownames(stat_df) <- NULL
colnames(stat_df) <- c("group","effect_size_lower","effect_size","effect_size_upper","effect_size_magnitude","wilcox_stat","wilcox_p")
  #getting data structure right
stat_df$effect_size_lower <- as.numeric(stat_df$effect_size_lower)
stat_df$effect_size <- as.numeric(stat_df$effect_size)
stat_df$effect_size_upper <- as.numeric(stat_df$effect_size_upper)
stat_df$effect_size_magnitude <- as.numeric(stat_df$effect_size_magnitude)
stat_df$wilcox_stat <- as.numeric(stat_df$wilcox_stat)
stat_df$wilcox_p <- as.numeric(stat_df$wilcox_p)
#correct for multiple tests
stat_df$p.adj <- p.adjust(stat_df$wilcox_p, method="BH")

write.table(stat_df, "/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/okay/wilcox_tests_faprotax_groups_fig_interios_only_worms_v_no_worms.tsv",sep="\t",row.names=FALSE, quote=FALSE)




sigcat <- subset(stat_df, p.adj < 0.05)

plotsig <- wormsdf[wormsdf$group %in% sigcat$group, ]

ggplot(plotsig, aes(x=value,y=group)) + geom_sina(aes(colour=worms),scale="width",alpha=0.5) + stat_summary(aes(group=worms),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1, colour="black",position = position_dodge(width = 0.9)) +theme_half_open() + background_grid(major = c("y"),minor = c("none")) + scale_color_manual(values = c("#D5770C","#760D7D")) + xlab("ASV Proportion")

summary(stat_df$p.adj)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.5413  0.5413  0.7427  0.7557  0.9465  1.0000      48

#alright