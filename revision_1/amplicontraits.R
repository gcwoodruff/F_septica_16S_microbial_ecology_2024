setwd("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/amplicontraits/01_elegans/00_usearch_out/") 
d1_up10 <- read.table("d1_up10.txt",header=FALSE,sep='\t')
cellshape <- read.table("cellshape.txt",header=FALSE,sep='\t')
d1_lo5 <- read.table("d1_lo5.txt",header=FALSE,sep='\t')
d1_lo10 <- read.table("d1_lo10.txt",header=FALSE,sep='\t')
d1_lo20 <- read.table("d1_lo20.txt",header=FALSE,sep='\t')
d1_lo30 <- read.table("d1_lo30.txt",header=FALSE,sep='\t')
d1_up5 <- read.table("d1_up5.txt",header=FALSE,sep='\t')
d1_up20 <- read.table("d1_up20.txt",header=FALSE,sep='\t')
d1_up30 <- read.table("d1_up30.txt",header=FALSE,sep='\t')
d2_lo5 <- read.table("d2_lo5.txt",header=FALSE,sep='\t')
d2_lo10 <- read.table("d2_lo10.txt",header=FALSE,sep='\t')
d2_lo20 <- read.table("d2_lo20.txt",header=FALSE,sep='\t')
d2_lo30 <- read.table("d2_lo30.txt",header=FALSE,sep='\t')
d2_up5 <- read.table("d2_up5.txt",header=FALSE,sep='\t')
d2_up10 <- read.table("d2_up10.txt",header=FALSE,sep='\t')
d2_up20 <- read.table("d2_up20.txt",header=FALSE,sep='\t')
d2_up30 <- read.table("d2_up30.txt",header=FALSE,sep='\t')
doubling_h5 <- read.table("doubling_h5.txt",header=FALSE,sep='\t')
doubling_h10 <- read.table("doubling_h10.txt",header=FALSE,sep='\t')
doubling_h20 <- read.table("doubling_h20.txt",header=FALSE,sep='\t')
doubling_h30 <- read.table("doubling_h30.txt",header=FALSE,sep='\t')
doubling_h40 <- read.table("doubling_h40.txt",header=FALSE,sep='\t')
doubling_h50 <- read.table("doubling_h50.txt",header=FALSE,sep='\t')
genome_size5 <- read.table("genome_size5.txt",header=FALSE,sep='\t')
genome_size10 <- read.table("genome_size10.txt",header=FALSE,sep='\t')
genome_size20 <- read.table("genome_size20.txt",header=FALSE,sep='\t')
gramstain <- read.table("gramstain.txt",header=FALSE,sep='\t')
metabolism <- read.table("metabolism.txt",header=FALSE,sep='\t')
motility <- read.table("motility.txt",header=FALSE,sep='\t')
optimum_ph5 <- read.table("optimum_ph5.txt",header=FALSE,sep='\t')
optimum_ph10 <- read.table("optimum_ph10.txt",header=FALSE,sep='\t')
optimum_ph20 <- read.table("optimum_ph20.txt",header=FALSE,sep='\t')
optimum_tmp10 <- read.table("optimum_tmp10.txt",header=FALSE,sep='\t')
optimum_tmp20 <- read.table("optimum_tmp20.txt",header=FALSE,sep='\t')
rangesalinity <- read.table("rangesalinity.txt",header=FALSE,sep='\t')
rangetmp <- read.table("rangetmp.txt",header=FALSE,sep='\t')
rRNA16S_genes_exact_int <- read.table("rRNA16S_genes_exact_int.txt",header=FALSE,sep='\t')
rRNA16S_genes_exact <- read.table("rRNA16S_genes_exact.txt",header=FALSE,sep='\t')
rRNA16S_genes5 <- read.table("rRNA16S_genes5.txt",header=FALSE,sep='\t')
rRNA16S_genes10 <- read.table("rRNA16S_genes10.txt",header=FALSE,sep='\t')
sporulation <- read.table("sporulation.txt",header=FALSE,sep='\t')


d1_up10_bs75 <- d1_up10[d1_up10$V3 >74,]
cellshape_bs75 <- cellshape[cellshape$V3 >74,]
d1_lo5_bs75 <- d1_lo5[d1_lo5$V3 >74,]
d1_lo10_bs75 <- d1_lo10[d1_lo10$V3 >74,]
d1_lo20_bs75 <- d1_lo20[d1_lo20$V3 >74,]
d1_lo30_bs75 <- d1_lo30[d1_lo30$V3 >74,]
d1_up5_bs75 <- d1_up5[d1_up5$V3 >74,]
d1_up20_bs75 <- d1_up20[d1_up20$V3 >74,]
d1_up30_bs75 <- d1_up30[d1_up30$V3 >74,]
d2_lo5_bs75 <- d2_lo5[d2_lo5$V3 >74,]
d2_lo10_bs75 <- d2_lo10[d2_lo10$V3 >74,]
d2_lo20_bs75 <- d2_lo20[d2_lo20$V3 >74,]
d2_lo30_bs75 <- d2_lo30[d2_lo30$V3 >74,]
d2_up5_bs75 <- d2_up5[d2_up5$V3 >74,]
d2_up10_bs75 <- d2_up10[d2_up10$V3 >74,]
d2_up20_bs75 <- d2_up20[d2_up20$V3 >74,]
d2_up30_bs75 <- d2_up30[d2_up30$V3 >74,]
doubling_h5_bs75 <- doubling_h5[doubling_h5$V3 >74,]
doubling_h10_bs75 <- doubling_h10[doubling_h10$V3 >74,]
doubling_h20_bs75 <- doubling_h20[doubling_h20$V3 >74,]
doubling_h30_bs75 <- doubling_h30[doubling_h30$V3 >74,]
doubling_h40_bs75 <- doubling_h40[doubling_h40$V3 >74,]
doubling_h50_bs75 <- doubling_h50[doubling_h50$V3 >74,]
genome_size5_bs75 <- genome_size5[genome_size5$V3 >74,]
genome_size10_bs75 <- genome_size10[genome_size10$V3 >74,]
genome_size20_bs75 <- genome_size20[genome_size20$V3 >74,]
gramstain_bs75 <- gramstain[gramstain$V3 >74,]
metabolism_bs75 <- metabolism[metabolism$V3 >74,]
motility_bs75 <- motility[motility$V3 >74,]
optimum_ph5_bs75 <- optimum_ph5[optimum_ph5$V3 >74,]
optimum_ph10_bs75 <- optimum_ph10[optimum_ph10$V3 >74,]
optimum_ph20_bs75 <- optimum_ph20[optimum_ph20$V3 >74,]
optimum_tmp10_bs75 <- optimum_tmp10[optimum_tmp10$V3 >74,]
optimum_tmp20_bs75 <- optimum_tmp20[optimum_tmp20$V3 >74,]
rangesalinity_bs75 <- rangesalinity[rangesalinity$V3 >74,]
rangetmp_bs75 <- rangetmp[rangetmp$V3 >74,]
rRNA16S_genes_exact_int_bs75 <- rRNA16S_genes_exact_int[rRNA16S_genes_exact_int$V3 >74,]
rRNA16S_genes_exact_bs75 <- rRNA16S_genes_exact[rRNA16S_genes_exact$V3 >74,]
rRNA16S_genes5_bs75 <- rRNA16S_genes5[rRNA16S_genes5$V3 >74,]
rRNA16S_genes10_bs75 <- rRNA16S_genes10[rRNA16S_genes10$V3 >74,]
sporulation_bs75 <- sporulation[sporulation$V3 >74,]

merge1 <- merge(d1_up10_bs75,cellshape_bs75,by="V1",all=TRUE)
merge1b <- data.frame(V1=merge1$V1,d1_up10=merge1$V2.x,cellshape=merge1$V2.y)

merge2<- merge(merge1b,d1_lo5_bs75,by="V1",all=TRUE)
merge2b <- merge2[ , !(names(merge2) %in% c("V3","V4"))]
names(merge2b)[names(merge2b) == 'V2'] <- 'd1_lo5'

merge3<- merge(merge2b,d1_lo10_bs75,by="V1",all=TRUE)
merge3b <- merge3[ , !(names(merge3) %in% c("V3","V4"))]
names(merge3b)[names(merge3b) == 'V2'] <- 'd1_lo10'

merge4<- merge(merge3b,d1_lo20_bs75,by="V1",all=TRUE)
merge4b <- merge4[ , !(names(merge4) %in% c("V3","V4"))]
names(merge4b)[names(merge4b) == 'V2'] <- 'd1_lo20'

#d1_up10==
#cellshape==
#d1_lo5==
#d1_lo10==
#d1_lo20==
#d1_lo30==
#d1_up5==
#d1_up20==
#d1_up30==
#d2_lo5==
#d2_lo10==
#d2_lo20==
#d2_lo30==
#d2_up5==
#d2_up10==
#d2_up20==
#d2_up30==
#doubling_h5==
#doubling_h10==
#doubling_h20==
#doubling_h30==
#doubling_h40==
#doubling_h50==
#genome_size5==
#genome_size10==
#genome_size20==
#gramstain==
#metabolism==
#motility==
#optimum_ph5==
#optimum_ph10==
#optimum_ph20==
#optimum_tmp10==
#optimum_tmp20==
#rangesalinity==
#rangetmp==
#rRNA16S_genes_exact_int==
#rRNA16S_genes_exact==
#rRNA16S_genes5==
#rRNA16S_genes10==
#sporulation

merge5<- merge(merge4b,d1_lo30_bs75,by="V1",all=TRUE)
merge5b <- merge5[ , !(names(merge5) %in% c("V3","V4"))]
names(merge5b)[names(merge5b) == 'V2'] <- 'd1_lo30_bs75'
names(merge5b)[names(merge5b) == 'd1_lo30_bs75'] <- 'd1_lo30'


merge6<- merge(merge5b,d1_up5_bs75,by="V1",all=TRUE)
merge6b <- merge6[ , !(names(merge6) %in% c("V3","V4"))]
names(merge6b)[names(merge6b) == 'V2'] <- 'd1_up5'

merge7<- merge(merge6b,d1_up20_bs75,by="V1",all=TRUE)
merge7b <- merge7[ , !(names(merge7) %in% c("V3","V4"))]
names(merge7b)[names(merge7b) == 'V2'] <- 'd1_up20'


merge8<- merge(merge7b,d1_up30_bs75,by="V1",all=TRUE)
merge8b <- merge8[ , !(names(merge8) %in% c("V3","V4"))]
names(merge8b)[names(merge8b) == 'V2'] <- 'd1_up30'

merge9<- merge(merge8b,d2_lo5_bs75,by="V1",all=TRUE)
merge9b <- merge9[ , !(names(merge9) %in% c("V3","V4"))]
names(merge9b)[names(merge9b) == 'V2'] <- 'd2_lo5'

merge10<- merge(merge9b,d2_lo10_bs75,by="V1",all=TRUE)
merge10b <- merge10[ , !(names(merge10) %in% c("V3","V4"))]
names(merge10b)[names(merge10b) == 'V2'] <- 'd2_lo10'

merge11<- merge(merge10b,d2_lo20_bs75,by="V1",all=TRUE)
merge11b <- merge11[ , !(names(merge11) %in% c("V3","V4"))]
names(merge11b)[names(merge11b) == 'V2'] <- 'd2_lo20'


merge12<- merge(merge11b,d2_lo30_bs75,by="V1",all=TRUE)
merge12b <- merge12[ , !(names(merge12) %in% c("V3","V4"))]
names(merge12b)[names(merge12b) == 'V2'] <- 'd2_lo30'


merge13<- merge(merge12b,d2_up5_bs75,by="V1",all=TRUE)
merge13b <- merge13[ , !(names(merge13) %in% c("V3","V4"))]
names(merge13b)[names(merge13b) == 'V2'] <- 'd2_up5'

merge14<- merge(merge13b,d2_up10_bs75,by="V1",all=TRUE)
merge14b <- merge14[ , !(names(merge14) %in% c("V3","V4"))]
names(merge14b)[names(merge14b) == 'V2'] <- 'd2_up10'


merge15<- merge(merge14b,d2_up20_bs75,by="V1",all=TRUE)
merge15b <- merge15[ , !(names(merge15) %in% c("V3","V4"))]
names(merge15b)[names(merge15b) == 'V2'] <- 'd2_up20'


merge16<- merge(merge15b,d2_up30_bs75,by="V1",all=TRUE)
merge16b <- merge16[ , !(names(merge16) %in% c("V3","V4"))]
names(merge16b)[names(merge16b) == 'V2'] <- 'd2_up30'


merge17<- merge(merge16b,doubling_h5_bs75,by="V1",all=TRUE)
merge17b <- merge17[ , !(names(merge17) %in% c("V3","V4"))]
names(merge17b)[names(merge17b) == 'V2'] <- 'doubling_h5'



merge18<- merge(merge17b,doubling_h10_bs75,by="V1",all=TRUE)
merge18b <- merge18[ , !(names(merge18) %in% c("V3","V4"))]
names(merge18b)[names(merge18b) == 'V2'] <- 'doubling_h10'

merge19<- merge(merge18b,doubling_h20_bs75,by="V1",all=TRUE)
merge19b <- merge19[ , !(names(merge19) %in% c("V3","V4"))]
names(merge19b)[names(merge19b) == 'V2'] <- 'doubling_h20'



merge20<- merge(merge19b,doubling_h30_bs75,by="V1",all=TRUE)
merge20b <- merge20[ , !(names(merge20) %in% c("V3","V4"))]
names(merge20b)[names(merge20b) == 'V2'] <- 'doubling_h30'




merge21<- merge(merge20b,doubling_h40_bs75,by="V1",all=TRUE)
merge21b <- merge21[ , !(names(merge21) %in% c("V3","V4"))]
names(merge21b)[names(merge21b) == 'V2'] <- 'doubling_h40'



merge22<- merge(merge21b,doubling_h50_bs75,by="V1",all=TRUE)
merge22b <- merge22[ , !(names(merge22) %in% c("V3","V4"))]
names(merge22b)[names(merge22b) == 'V2'] <- 'doubling_h50'


merge23<- merge(merge22b,genome_size5_bs75,by="V1",all=TRUE)
merge23b <- merge23[ , !(names(merge23) %in% c("V3","V4"))]
names(merge23b)[names(merge23b) == 'V2'] <- 'genome_size5'



merge24<- merge(merge23b,genome_size10_bs75,by="V1",all=TRUE)
merge24b <- merge24[ , !(names(merge24) %in% c("V3","V4"))]
names(merge24b)[names(merge24b) == 'V2'] <- 'genome_size10'


merge25<- merge(merge24b,genome_size20_bs75,by="V1",all=TRUE)
merge25b <- merge25[ , !(names(merge25) %in% c("V3","V4"))]
names(merge25b)[names(merge25b) == 'V2'] <- 'genome_size20'


merge26<- merge(merge25b,gramstain_bs75,by="V1",all=TRUE)
merge26b <- merge26[ , !(names(merge26) %in% c("V3","V4"))]
names(merge26b)[names(merge26b) == 'V2'] <- 'gramstain'



merge27<- merge(merge26b,metabolism_bs75,by="V1",all=TRUE)
merge27b <- merge27[ , !(names(merge27) %in% c("V3","V4"))]
names(merge27b)[names(merge27b) == 'V2'] <- 'metabolism'

merge28<- merge(merge27b,motility_bs75,by="V1",all=TRUE)
merge28b <- merge28[ , !(names(merge28) %in% c("V3","V4"))]
names(merge28b)[names(merge28b) == 'V2'] <- 'motility'



merge29<- merge(merge28b,optimum_ph5_bs75,by="V1",all=TRUE)
merge29b <- merge29[ , !(names(merge29) %in% c("V3","V4"))]
names(merge29b)[names(merge29b) == 'V2'] <- 'optimum_ph5'




merge30<- merge(merge29b,optimum_ph10_bs75,by="V1",all=TRUE)
merge30b <- merge30[ , !(names(merge30) %in% c("V3","V4"))]
names(merge30b)[names(merge30b) == 'V2'] <- 'optimum_ph10'




merge31<- merge(merge30b,optimum_ph20_bs75,by="V1",all=TRUE)
merge31b <- merge31[ , !(names(merge31) %in% c("V3","V4"))]
names(merge31b)[names(merge31b) == 'V2'] <- 'optimum_ph20'


merge32<- merge(merge31b,optimum_tmp10_bs75,by="V1",all=TRUE)
merge32b <- merge32[ , !(names(merge32) %in% c("V3","V4"))]
names(merge32b)[names(merge32b) == 'V2'] <- 'optimum_tmp10'


merge33<- merge(merge32b,optimum_tmp20_bs75,by="V1",all=TRUE)
merge33b <- merge33[ , !(names(merge33) %in% c("V3","V4"))]
names(merge33b)[names(merge33b) == 'V2'] <- 'optimum_tmp20'

merge34<- merge(merge33b,rangesalinity_bs75,by="V1",all=TRUE)
merge34b <- merge34[ , !(names(merge34) %in% c("V3","V4"))]
names(merge34b)[names(merge34b) == 'V2'] <- 'rangesalinity'

merge35<- merge(merge34b,rangetmp_bs75,by="V1",all=TRUE)
merge35b <- merge35[ , !(names(merge35) %in% c("V3","V4"))]
names(merge35b)[names(merge35b) == 'V2'] <- 'rangetmp'

merge36<- merge(merge35b,rRNA16S_genes_exact_int_bs75,by="V1",all=TRUE)
merge36b <- merge36[ , !(names(merge36) %in% c("V3","V4"))]
names(merge36b)[names(merge36b) == 'V2'] <- 'rRNA16S_genes_exact_int'


merge37<- merge(merge36b,rRNA16S_genes_exact_bs75,by="V1",all=TRUE)
merge37b <- merge37[ , !(names(merge37) %in% c("V3","V4"))]
names(merge37b)[names(merge37b) == 'V2'] <- 'rRNA16S_genes_exact'


merge38<- merge(merge37b,rRNA16S_genes5_bs75,by="V1",all=TRUE)
merge38b <- merge38[ , !(names(merge38) %in% c("V3","V4"))]
names(merge38b)[names(merge38b) == 'V2'] <- 'rRNA16S_genes5'


merge39<- merge(merge38b,rRNA16S_genes10_bs75,by="V1",all=TRUE)
merge39b <- merge39[ , !(names(merge39) %in% c("V3","V4"))]
names(merge39b)[names(merge39b) == 'V2'] <- 'rRNA16S_genes10'


merge40<- merge(merge39b,sporulation_bs75,by="V1",all=TRUE)
merge40b <- merge40[ , !(names(merge40) %in% c("V3","V4"))]
names(merge40b)[names(merge40b) == 'V2'] <- 'sporulation'

#woo!

write.table(merge40b,"function_table.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#just get some sense of trait distribution...

as.data.frame(table(merge40b$metabolism))
#                Var1  Freq
#1            aerobic 13287
#2          anaerobic  1979
#3        facultative   934
#4    microaerophilic  3177
#5   obligate_aerobic  5561
#6 obligate_anaerobic   468

as.data.frame(table(merge40b$cellshape))
#            Var1  Freq
#1       bacillus 43915
#2        branced     4
#3  coccobacillus   549
#4         coccus  2470
#5       filament   432
#6       fusiform     4
#7    pleomorphic    10
#8        spindle     2
#9         spiral   118
#10          star    18
#11        vibrio   255

as.data.frame(table(merge40b$motility))
#            Var1  Freq
#1 axial_filament     5
#2       flagella  1259
#3        gliding   436
#4             no 17730
#5            yes  8779





###figs only

setwd("/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/amplicontraits/00_usearch_out/") 
d1_up10 <- read.table("d1_up10.txt",header=FALSE,sep='\t')
cellshape <- read.table("cellshape.txt",header=FALSE,sep='\t')
d1_lo5 <- read.table("d1_lo5.txt",header=FALSE,sep='\t')
d1_lo10 <- read.table("d1_lo10.txt",header=FALSE,sep='\t')
d1_lo20 <- read.table("d1_lo20.txt",header=FALSE,sep='\t')
d1_lo30 <- read.table("d1_lo30.txt",header=FALSE,sep='\t')
d1_up5 <- read.table("d1_up5.txt",header=FALSE,sep='\t')
d1_up20 <- read.table("d1_up20.txt",header=FALSE,sep='\t')
d1_up30 <- read.table("d1_up30.txt",header=FALSE,sep='\t')
d2_lo5 <- read.table("d2_lo5.txt",header=FALSE,sep='\t')
d2_lo10 <- read.table("d2_lo10.txt",header=FALSE,sep='\t')
d2_lo20 <- read.table("d2_lo20.txt",header=FALSE,sep='\t')
d2_lo30 <- read.table("d2_lo30.txt",header=FALSE,sep='\t')
d2_up5 <- read.table("d2_up5.txt",header=FALSE,sep='\t')
d2_up10 <- read.table("d2_up10.txt",header=FALSE,sep='\t')
d2_up20 <- read.table("d2_up20.txt",header=FALSE,sep='\t')
d2_up30 <- read.table("d2_up30.txt",header=FALSE,sep='\t')
doubling_h5 <- read.table("doubling_h5.txt",header=FALSE,sep='\t')
doubling_h10 <- read.table("doubling_h10.txt",header=FALSE,sep='\t')
doubling_h20 <- read.table("doubling_h20.txt",header=FALSE,sep='\t')
doubling_h30 <- read.table("doubling_h30.txt",header=FALSE,sep='\t')
doubling_h40 <- read.table("doubling_h40.txt",header=FALSE,sep='\t')
doubling_h50 <- read.table("doubling_h50.txt",header=FALSE,sep='\t')
genome_size5 <- read.table("genome_size5.txt",header=FALSE,sep='\t')
genome_size10 <- read.table("genome_size10.txt",header=FALSE,sep='\t')
genome_size20 <- read.table("genome_size20.txt",header=FALSE,sep='\t')
gramstain <- read.table("gramstain.txt",header=FALSE,sep='\t')
metabolism <- read.table("metabolism.txt",header=FALSE,sep='\t')
motility <- read.table("motility.txt",header=FALSE,sep='\t')
optimum_ph5 <- read.table("optimum_ph5.txt",header=FALSE,sep='\t')
optimum_ph10 <- read.table("optimum_ph10.txt",header=FALSE,sep='\t')
optimum_ph20 <- read.table("optimum_ph20.txt",header=FALSE,sep='\t')
optimum_tmp10 <- read.table("optimum_tmp10.txt",header=FALSE,sep='\t')
optimum_tmp20 <- read.table("optimum_tmp20.txt",header=FALSE,sep='\t')
rangesalinity <- read.table("rangesalinity.txt",header=FALSE,sep='\t')
rangetmp <- read.table("rangetmp.txt",header=FALSE,sep='\t')
rRNA16S_genes_exact_int <- read.table("rRNA16S_genes_exact_int.txt",header=FALSE,sep='\t')
rRNA16S_genes_exact <- read.table("rRNA16S_genes_exact.txt",header=FALSE,sep='\t')
rRNA16S_genes5 <- read.table("rRNA16S_genes5.txt",header=FALSE,sep='\t')
rRNA16S_genes10 <- read.table("rRNA16S_genes10.txt",header=FALSE,sep='\t')
sporulation <- read.table("sporulation.txt",header=FALSE,sep='\t')


d1_up10_bs75 <- d1_up10[d1_up10$V3 >74,]
cellshape_bs75 <- cellshape[cellshape$V3 >74,]
d1_lo5_bs75 <- d1_lo5[d1_lo5$V3 >74,]
d1_lo10_bs75 <- d1_lo10[d1_lo10$V3 >74,]
d1_lo20_bs75 <- d1_lo20[d1_lo20$V3 >74,]
d1_lo30_bs75 <- d1_lo30[d1_lo30$V3 >74,]
d1_up5_bs75 <- d1_up5[d1_up5$V3 >74,]
d1_up20_bs75 <- d1_up20[d1_up20$V3 >74,]
d1_up30_bs75 <- d1_up30[d1_up30$V3 >74,]
d2_lo5_bs75 <- d2_lo5[d2_lo5$V3 >74,]
d2_lo10_bs75 <- d2_lo10[d2_lo10$V3 >74,]
d2_lo20_bs75 <- d2_lo20[d2_lo20$V3 >74,]
d2_lo30_bs75 <- d2_lo30[d2_lo30$V3 >74,]
d2_up5_bs75 <- d2_up5[d2_up5$V3 >74,]
d2_up10_bs75 <- d2_up10[d2_up10$V3 >74,]
d2_up20_bs75 <- d2_up20[d2_up20$V3 >74,]
d2_up30_bs75 <- d2_up30[d2_up30$V3 >74,]
doubling_h5_bs75 <- doubling_h5[doubling_h5$V3 >74,]
doubling_h10_bs75 <- doubling_h10[doubling_h10$V3 >74,]
doubling_h20_bs75 <- doubling_h20[doubling_h20$V3 >74,]
doubling_h30_bs75 <- doubling_h30[doubling_h30$V3 >74,]
doubling_h40_bs75 <- doubling_h40[doubling_h40$V3 >74,]
doubling_h50_bs75 <- doubling_h50[doubling_h50$V3 >74,]
genome_size5_bs75 <- genome_size5[genome_size5$V3 >74,]
genome_size10_bs75 <- genome_size10[genome_size10$V3 >74,]
genome_size20_bs75 <- genome_size20[genome_size20$V3 >74,]
gramstain_bs75 <- gramstain[gramstain$V3 >74,]
metabolism_bs75 <- metabolism[metabolism$V3 >74,]
motility_bs75 <- motility[motility$V3 >74,]
optimum_ph5_bs75 <- optimum_ph5[optimum_ph5$V3 >74,]
optimum_ph10_bs75 <- optimum_ph10[optimum_ph10$V3 >74,]
optimum_ph20_bs75 <- optimum_ph20[optimum_ph20$V3 >74,]
optimum_tmp10_bs75 <- optimum_tmp10[optimum_tmp10$V3 >74,]
optimum_tmp20_bs75 <- optimum_tmp20[optimum_tmp20$V3 >74,]
rangesalinity_bs75 <- rangesalinity[rangesalinity$V3 >74,]
rangetmp_bs75 <- rangetmp[rangetmp$V3 >74,]
rRNA16S_genes_exact_int_bs75 <- rRNA16S_genes_exact_int[rRNA16S_genes_exact_int$V3 >74,]
rRNA16S_genes_exact_bs75 <- rRNA16S_genes_exact[rRNA16S_genes_exact$V3 >74,]
rRNA16S_genes5_bs75 <- rRNA16S_genes5[rRNA16S_genes5$V3 >74,]
rRNA16S_genes10_bs75 <- rRNA16S_genes10[rRNA16S_genes10$V3 >74,]
sporulation_bs75 <- sporulation[sporulation$V3 >74,]

merge1 <- merge(d1_up10_bs75,cellshape_bs75,by="V1",all=TRUE)
merge1b <- data.frame(V1=merge1$V1,d1_up10=merge1$V2.x,cellshape=merge1$V2.y)

merge2<- merge(merge1b,d1_lo5_bs75,by="V1",all=TRUE)
merge2b <- merge2[ , !(names(merge2) %in% c("V3","V4"))]
names(merge2b)[names(merge2b) == 'V2'] <- 'd1_lo5'

merge3<- merge(merge2b,d1_lo10_bs75,by="V1",all=TRUE)
merge3b <- merge3[ , !(names(merge3) %in% c("V3","V4"))]
names(merge3b)[names(merge3b) == 'V2'] <- 'd1_lo10'

merge4<- merge(merge3b,d1_lo20_bs75,by="V1",all=TRUE)
merge4b <- merge4[ , !(names(merge4) %in% c("V3","V4"))]
names(merge4b)[names(merge4b) == 'V2'] <- 'd1_lo20'

#d1_up10==
#cellshape==
#d1_lo5==
#d1_lo10==
#d1_lo20==
#d1_lo30==
#d1_up5==
#d1_up20==
#d1_up30==
#d2_lo5==
#d2_lo10==
#d2_lo20==
#d2_lo30==
#d2_up5==
#d2_up10==
#d2_up20==
#d2_up30==
#doubling_h5==
#doubling_h10==
#doubling_h20==
#doubling_h30==
#doubling_h40==
#doubling_h50==
#genome_size5==
#genome_size10==
#genome_size20==
#gramstain==
#metabolism==
#motility==
#optimum_ph5==
#optimum_ph10==
#optimum_ph20==
#optimum_tmp10==
#optimum_tmp20==
#rangesalinity==
#rangetmp==
#rRNA16S_genes_exact_int==
#rRNA16S_genes_exact==
#rRNA16S_genes5==
#rRNA16S_genes10==
#sporulation

merge5<- merge(merge4b,d1_lo30_bs75,by="V1",all=TRUE)
merge5b <- merge5[ , !(names(merge5) %in% c("V3","V4"))]
names(merge5b)[names(merge5b) == 'V2'] <- 'd1_lo30_bs75'
names(merge5b)[names(merge5b) == 'd1_lo30_bs75'] <- 'd1_lo30'


merge6<- merge(merge5b,d1_up5_bs75,by="V1",all=TRUE)
merge6b <- merge6[ , !(names(merge6) %in% c("V3","V4"))]
names(merge6b)[names(merge6b) == 'V2'] <- 'd1_up5'

merge7<- merge(merge6b,d1_up20_bs75,by="V1",all=TRUE)
merge7b <- merge7[ , !(names(merge7) %in% c("V3","V4"))]
names(merge7b)[names(merge7b) == 'V2'] <- 'd1_up20'


merge8<- merge(merge7b,d1_up30_bs75,by="V1",all=TRUE)
merge8b <- merge8[ , !(names(merge8) %in% c("V3","V4"))]
names(merge8b)[names(merge8b) == 'V2'] <- 'd1_up30'

merge9<- merge(merge8b,d2_lo5_bs75,by="V1",all=TRUE)
merge9b <- merge9[ , !(names(merge9) %in% c("V3","V4"))]
names(merge9b)[names(merge9b) == 'V2'] <- 'd2_lo5'

merge10<- merge(merge9b,d2_lo10_bs75,by="V1",all=TRUE)
merge10b <- merge10[ , !(names(merge10) %in% c("V3","V4"))]
names(merge10b)[names(merge10b) == 'V2'] <- 'd2_lo10'

merge11<- merge(merge10b,d2_lo20_bs75,by="V1",all=TRUE)
merge11b <- merge11[ , !(names(merge11) %in% c("V3","V4"))]
names(merge11b)[names(merge11b) == 'V2'] <- 'd2_lo20'


merge12<- merge(merge11b,d2_lo30_bs75,by="V1",all=TRUE)
merge12b <- merge12[ , !(names(merge12) %in% c("V3","V4"))]
names(merge12b)[names(merge12b) == 'V2'] <- 'd2_lo30'


merge13<- merge(merge12b,d2_up5_bs75,by="V1",all=TRUE)
merge13b <- merge13[ , !(names(merge13) %in% c("V3","V4"))]
names(merge13b)[names(merge13b) == 'V2'] <- 'd2_up5'

merge14<- merge(merge13b,d2_up10_bs75,by="V1",all=TRUE)
merge14b <- merge14[ , !(names(merge14) %in% c("V3","V4"))]
names(merge14b)[names(merge14b) == 'V2'] <- 'd2_up10'


merge15<- merge(merge14b,d2_up20_bs75,by="V1",all=TRUE)
merge15b <- merge15[ , !(names(merge15) %in% c("V3","V4"))]
names(merge15b)[names(merge15b) == 'V2'] <- 'd2_up20'


merge16<- merge(merge15b,d2_up30_bs75,by="V1",all=TRUE)
merge16b <- merge16[ , !(names(merge16) %in% c("V3","V4"))]
names(merge16b)[names(merge16b) == 'V2'] <- 'd2_up30'


merge17<- merge(merge16b,doubling_h5_bs75,by="V1",all=TRUE)
merge17b <- merge17[ , !(names(merge17) %in% c("V3","V4"))]
names(merge17b)[names(merge17b) == 'V2'] <- 'doubling_h5'



merge18<- merge(merge17b,doubling_h10_bs75,by="V1",all=TRUE)
merge18b <- merge18[ , !(names(merge18) %in% c("V3","V4"))]
names(merge18b)[names(merge18b) == 'V2'] <- 'doubling_h10'

merge19<- merge(merge18b,doubling_h20_bs75,by="V1",all=TRUE)
merge19b <- merge19[ , !(names(merge19) %in% c("V3","V4"))]
names(merge19b)[names(merge19b) == 'V2'] <- 'doubling_h20'



merge20<- merge(merge19b,doubling_h30_bs75,by="V1",all=TRUE)
merge20b <- merge20[ , !(names(merge20) %in% c("V3","V4"))]
names(merge20b)[names(merge20b) == 'V2'] <- 'doubling_h30'




merge21<- merge(merge20b,doubling_h40_bs75,by="V1",all=TRUE)
merge21b <- merge21[ , !(names(merge21) %in% c("V3","V4"))]
names(merge21b)[names(merge21b) == 'V2'] <- 'doubling_h40'



merge22<- merge(merge21b,doubling_h50_bs75,by="V1",all=TRUE)
merge22b <- merge22[ , !(names(merge22) %in% c("V3","V4"))]
names(merge22b)[names(merge22b) == 'V2'] <- 'doubling_h50'


merge23<- merge(merge22b,genome_size5_bs75,by="V1",all=TRUE)
merge23b <- merge23[ , !(names(merge23) %in% c("V3","V4"))]
names(merge23b)[names(merge23b) == 'V2'] <- 'genome_size5'



merge24<- merge(merge23b,genome_size10_bs75,by="V1",all=TRUE)
merge24b <- merge24[ , !(names(merge24) %in% c("V3","V4"))]
names(merge24b)[names(merge24b) == 'V2'] <- 'genome_size10'


merge25<- merge(merge24b,genome_size20_bs75,by="V1",all=TRUE)
merge25b <- merge25[ , !(names(merge25) %in% c("V3","V4"))]
names(merge25b)[names(merge25b) == 'V2'] <- 'genome_size20'


merge26<- merge(merge25b,gramstain_bs75,by="V1",all=TRUE)
merge26b <- merge26[ , !(names(merge26) %in% c("V3","V4"))]
names(merge26b)[names(merge26b) == 'V2'] <- 'gramstain'



merge27<- merge(merge26b,metabolism_bs75,by="V1",all=TRUE)
merge27b <- merge27[ , !(names(merge27) %in% c("V3","V4"))]
names(merge27b)[names(merge27b) == 'V2'] <- 'metabolism'

merge28<- merge(merge27b,motility_bs75,by="V1",all=TRUE)
merge28b <- merge28[ , !(names(merge28) %in% c("V3","V4"))]
names(merge28b)[names(merge28b) == 'V2'] <- 'motility'



merge29<- merge(merge28b,optimum_ph5_bs75,by="V1",all=TRUE)
merge29b <- merge29[ , !(names(merge29) %in% c("V3","V4"))]
names(merge29b)[names(merge29b) == 'V2'] <- 'optimum_ph5'




merge30<- merge(merge29b,optimum_ph10_bs75,by="V1",all=TRUE)
merge30b <- merge30[ , !(names(merge30) %in% c("V3","V4"))]
names(merge30b)[names(merge30b) == 'V2'] <- 'optimum_ph10'




merge31<- merge(merge30b,optimum_ph20_bs75,by="V1",all=TRUE)
merge31b <- merge31[ , !(names(merge31) %in% c("V3","V4"))]
names(merge31b)[names(merge31b) == 'V2'] <- 'optimum_ph20'


merge32<- merge(merge31b,optimum_tmp10_bs75,by="V1",all=TRUE)
merge32b <- merge32[ , !(names(merge32) %in% c("V3","V4"))]
names(merge32b)[names(merge32b) == 'V2'] <- 'optimum_tmp10'


merge33<- merge(merge32b,optimum_tmp20_bs75,by="V1",all=TRUE)
merge33b <- merge33[ , !(names(merge33) %in% c("V3","V4"))]
names(merge33b)[names(merge33b) == 'V2'] <- 'optimum_tmp20'

merge34<- merge(merge33b,rangesalinity_bs75,by="V1",all=TRUE)
merge34b <- merge34[ , !(names(merge34) %in% c("V3","V4"))]
names(merge34b)[names(merge34b) == 'V2'] <- 'rangesalinity'

merge35<- merge(merge34b,rangetmp_bs75,by="V1",all=TRUE)
merge35b <- merge35[ , !(names(merge35) %in% c("V3","V4"))]
names(merge35b)[names(merge35b) == 'V2'] <- 'rangetmp'

merge36<- merge(merge35b,rRNA16S_genes_exact_int_bs75,by="V1",all=TRUE)
merge36b <- merge36[ , !(names(merge36) %in% c("V3","V4"))]
names(merge36b)[names(merge36b) == 'V2'] <- 'rRNA16S_genes_exact_int'


merge37<- merge(merge36b,rRNA16S_genes_exact_bs75,by="V1",all=TRUE)
merge37b <- merge37[ , !(names(merge37) %in% c("V3","V4"))]
names(merge37b)[names(merge37b) == 'V2'] <- 'rRNA16S_genes_exact'


merge38<- merge(merge37b,rRNA16S_genes5_bs75,by="V1",all=TRUE)
merge38b <- merge38[ , !(names(merge38) %in% c("V3","V4"))]
names(merge38b)[names(merge38b) == 'V2'] <- 'rRNA16S_genes5'


merge39<- merge(merge38b,rRNA16S_genes10_bs75,by="V1",all=TRUE)
merge39b <- merge39[ , !(names(merge39) %in% c("V3","V4"))]
names(merge39b)[names(merge39b) == 'V2'] <- 'rRNA16S_genes10'


merge40<- merge(merge39b,sporulation_bs75,by="V1",all=TRUE)
merge40b <- merge40[ , !(names(merge40) %in% c("V3","V4"))]
names(merge40b)[names(merge40b) == 'V2'] <- 'sporulation'

#woo!

write.table(merge40b,"function_table_figs_only.tsv",sep="\t",quote=FALSE,row.names=FALSE)
