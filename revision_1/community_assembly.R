
##okay try community assembly, neutral/null models with iCAMP

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

#install.packages("iCAMP")
library(iCAMP)


ps3 <- readRDS("ps3_figs_only.rds")



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

env_num_only <- data.frame(latitude=env$latitude,longitude=env$longitude,day_of_aug_2019_picked=env$day_of_aug_2019_picked,day_of_aug_2019_dissected=env$day_of_aug_2019_dissected,foundress_number=env$foundress_number,row.names = rownames(env),control_interior_surface=as.numeric(as.factor(env$control_interior_surface)),caeno_present=as.numeric(as.factor(env$caeno_present)))

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
bin.size.limit = 12 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
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

write.table(icres$CbMPDiCBraya,"/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/icamp_out_figs_only/icamp_big_relative_importance_pairwise.tsv",sep="\t",row.names=FALSE)


detail.null=TRUE
bin.size.limit = 12 
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

write.table(icres2$CbMPDiCBraya,"/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/icamp_out_figs_only/icamp_big_SES.RC_relative_importance_pairwise.tsv",sep="\t",row.names=FALSE)




# 9.2.2 # normality test
nntest=iCAMP::null.norm(icamp.output=icres2, p.norm.cut=0.05, detail.out=FALSE)
# output shows non-normal distribution ratio in each bin, i.e. the proportion of turnovers which have null values significantly deviated from normal distribution.
# if some ratio values are very high, may need to change to use "Confidence" as sig.index.
  #1: In nortest::cvm.test(x) :
  #p-value is smaller than 7.37e-10, cannot be computed more accurately

write.table(nntest,"/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/icamp_out_figs_only/icamp_big_SES.RC_nntest.tsv",sep="\t",row.names=FALSE)

# 9.2.3 # change sig.index to "Confidence".
icres3=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "Confidence", detail.save = TRUE, detail.null = FALSE, conf.cut = 0.975)
head(icres3$CbMPDiCBraya)


write.table(icres3$CbMPDiCBraya,"/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/icamp_out_figs_only/icamp_big_icres3.tsv",sep="\t",row.names=FALSE)


# 9.2.4 # change sig.index to "RC" for both phylogenetic and taxonomic metrics.
icres4=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "RC", detail.save = TRUE, detail.null = FALSE, rc.cut = 0.95)
head(icres4$RCbMPDiRCbraya)


write.table(icres4$CbMPDiCBraya,"/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/icamp_out_figs_only/icamp_big_icres4.tsv",sep="\t",row.names=FALSE)


# 9.2.5 # the function can also change the significance threshold.
icres5=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "SES.RC", detail.save = TRUE, detail.null = FALSE, ses.cut = 1.64, rc.cut = 0.9)
head(icres5$bNRIiRCbraya)

write.table(icres5$bNRIiRCbraya,"/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/icamp_out_figs_only/icamp_big_icres5.tsv",sep="\t",row.names=FALSE)



# 9.5 # input community matrix as relative abundances (values < 1) rather than counts
comra=comm/rowSums(comm)
prefixra=paste0(prefix,"RA")
bin.size.limit = 80 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
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


save.image(file = "community_assembly_workspace_8-11-2025.RData")





#following Johnke et al. 2025...


#icres = iCAMP::icamp.big (ds = 0.2, pd.cut = NA, sp.check = TRUE, phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all", phylo.metric = "bMPD", bin.size.limit = 24, detail.null = FALSE, ignore.zero = TRUE, correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend", ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975, omit.option = "no",meta.ab =NULL)

#icbin = iCAMP::icamp.bins(icamp.detail = icres$detail,treat = treatment, clas = clas,silent = FALSE, boot = TRUE, rand.time = rand.time,between.group =TRUE)

load("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/icamp_out_figs_only/community_assembly_workspace_8-11-2025.RData")


icresJohnke = iCAMP::icamp.big(comm=comm,tree=tree,ds = 0.2, pd.cut = NA, sp.check = TRUE, phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all", phylo.metric = "bMPD", bin.size.limit = 24, detail.null = FALSE, ignore.zero = TRUE, correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend", ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975, omit.option = "no",meta.ab =NULL)

newtreat <- data.frame(surface_or_interior=treat$surface_or_interior,row.names=rownames(treat))

icbin = iCAMP::icamp.bins(icamp.detail = icresJohnke$detail,treat = newtreat, clas = clas,silent = FALSE, boot = TRUE, rand.time = rand.time,between.group =TRUE)


save.image(file = "community_assembly_workspace_9-29-2025.RData")

plotdat <- icbin$Pt

plotdat1 <- plotdat[c(1:2),c(3:8)]

plotdat1melt <- reshape2::melt(plotdat1,id.vars="Group")

levels(plotdat1melt$variable)[match("HeS",levels(plotdat1melt$variable))] <- "Heterogeneous Selection"
levels(plotdat1melt$variable)[match("HoS",levels(plotdat1melt$variable))] <- "Homogeneous Selection"
levels(plotdat1melt$variable)[match("DL",levels(plotdat1melt$variable))] <- "Dispersal Limitation"
levels(plotdat1melt$variable)[match("HD",levels(plotdat1melt$variable))] <- "Homogenizing Dispersal"
levels(plotdat1melt$variable)[match("DR",levels(plotdat1melt$variable))] <- "Drift and Others"

plotdat1melt$Group <- as.factor(plotdat1melt$Group)

levels(plotdat1melt$Group)[match("interior",levels(plotdat1melt$Group))] <- "Fig suspension"
levels(plotdat1melt$Group)[match("surface",levels(plotdat1melt$Group))] <- "Fig surface wash"

plotdat1melt$value <- as.numeric(plotdat1melt$value)


ggplot(plotdat1melt, aes(x = Group, y = value, fill = variable)) + geom_col() + scale_fill_brewer(palette='Set1') + theme_cowplot() +ylab("Proportion")




write.table(icbin$Pt,"/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/icamp_out_figs_only/icbin_relative_importance_sufrace_interior_9-29-25.tsv",sep="\t",row.names=FALSE)

