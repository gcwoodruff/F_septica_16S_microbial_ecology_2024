#Code for read processing, OTU assignment, and OTU counts in Woodruff et al. 2024, "The bacteria of a fig microcommunity"


###download reads
  #now on SRA, BioProject ID PRJNA1170329


#make working directory
mkdir 16S_fig_metabarcoding_3-2022
#navigate to working directory
cd 16S_fig_metabarcoding_3-2022
#make directory for raw reads
mkdir 00_raw_fastq
#make directory for slurm scripts
mkdir slurm_scripts

#navigate to raw reads directory
cd 00_raw_fastq

#download reads
	#this command was put into a slurm script
wget -c -r -l 10 --no-parent --reject '*.html*' -nH https://genomics-core-public.omrf.org/cclwood-gwcclwoodraff-kio9PyMHRLgMD8Ps/


###FastQC


#set working directory

wkdir='/scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/'

#load FastQC module for read quality
module load FastQC/0.11.5-Java-1.8.0_91

#make directory for FastQC output 
mkdir $wkdir/01_FastQC
#navigate to reads
cd $wkdir/00_raw_fastq/cclwood-gwcclwoodraff-kio9PyMHRLgMD8Ps
#run FastQC for all fastq files
for i in *fastq.gz; do fastqc -o $wkdir/01_FastQC $i; done


###read counts per sample


#get number of reads per fastqc
#make output directory
mkdir $wkdir/02_read_counts
#navigate to reads
cd $wkdir/00_raw_fastq/cclwood-gwcclwoodraff-kio9PyMHRLgMD8Ps
#count reads
for i in *fastq.gz; do zcat $i | echo $((`wc -l`/4)) > $wkdir/02_read_counts/${i%.fastq.gz}.txt; done
#make sample id file
cd $wkdir/02_read_counts

ls > $wkdir/fastq_ids.txt
#put all read counts in single file
cat * > read_counts.tmp
#combine read ids with read counts
paste  $wkdir/fastq_ids.txt read_counts.tmp > read_counts.tsv
#remove temp files
rm *.txt
rm read_counts.tmp
#make read count table; these data are in supplemental table 3
cat read_counts.tsv
  #aka read_counts_one_per_sample.tsv

###renaming fastq files so they can be imported into Qiime2


#copy fastq in case the renaming goes wrong
cd  $wkdir/00_raw_fastq/cclwood-gwcclwoodraff-kio9PyMHRLgMD8Ps

mkdir $wkdir/00_raw_fastq/00_cp/

for i in *.fastq.gz; do cp $i /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/03_qiime_import_fastq/00_raw_fastq/00_cp/$i; done

cd $wkdir/00_raw_fastq/00_cp/

#rename fastq files for Qiime2-- this software wants a very specific file name structure

mv Extr_Neg_1_S61_R1_001.fastq.gz ExtrNeg1_S61_L001_R1_001.fastq.gz
mv Extr_Neg_1_S61_R2_001.fastq.gz ExtrNeg1_S61_L001_R2_001.fastq.gz
mv EXTR_NEG_S73_R1_001.fastq.gz EXTRNEG_S73_L001_R1_001.fastq.gz
mv EXTR_NEG_S73_R2_001.fastq.gz EXTRNEG_S73_L001_R2_001.fastq.gz
mv EXTR_POS_S72_R1_001.fastq.gz EXTRPOS_S72_L001_R1_001.fastq.gz
mv EXTR_POS_S72_R2_001.fastq.gz EXTRPOS_S72_L001_R2_001.fastq.gz
mv GW_10_S2_R1_001.fastq.gz GW10_S2_L001_R1_001.fastq.gz
mv GW_10_S2_R2_001.fastq.gz GW10_S2_L001_R2_001.fastq.gz
mv GW_11_S3_R1_001.fastq.gz GW11_S3_L001_R1_001.fastq.gz
mv GW_11_S3_R2_001.fastq.gz GW11_S3_L001_R2_001.fastq.gz
mv GW_12_S4_R1_001.fastq.gz GW12_S4_L001_R1_001.fastq.gz
mv GW_12_S4_R2_001.fastq.gz GW12_S4_L001_R2_001.fastq.gz
mv GW_13_S5_R1_001.fastq.gz GW13_S5_L001_R1_001.fastq.gz
mv GW_13_S5_R2_001.fastq.gz GW13_S5_L001_R2_001.fastq.gz
mv GW_14_S6_R1_001.fastq.gz GW14_S6_L001_R1_001.fastq.gz
mv GW_14_S6_R2_001.fastq.gz GW14_S6_L001_R2_001.fastq.gz
mv GW_15_S7_R1_001.fastq.gz GW15_S7_L001_R1_001.fastq.gz
mv GW_15_S7_R2_001.fastq.gz GW15_S7_L001_R2_001.fastq.gz
mv GW_16_S74_R1_001.fastq.gz GW16_S74_L001_R1_001.fastq.gz
mv GW_16_S74_R2_001.fastq.gz GW16_S74_L001_R2_001.fastq.gz
mv GW_17_S8_R1_001.fastq.gz GW17_S8_L001_R1_001.fastq.gz
mv GW_17_S8_R2_001.fastq.gz GW17_S8_L001_R2_001.fastq.gz
mv GW_18_S9_R1_001.fastq.gz GW18_S9_L001_R1_001.fastq.gz
mv GW_18_S9_R2_001.fastq.gz GW18_S9_L001_R2_001.fastq.gz
mv GW_19_S10_R1_001.fastq.gz GW19_S10_L001_R1_001.fastq.gz
mv GW_19_S10_R2_001.fastq.gz GW19_S10_L001_R2_001.fastq.gz
mv GW_1_S64_R1_001.fastq.gz GW1_S64_L001_R1_001.fastq.gz
mv GW_1_S64_R2_001.fastq.gz GW1_S64_L001_R2_001.fastq.gz
mv GW_20_S11_R1_001.fastq.gz GW20_S11_L001_R1_001.fastq.gz
mv GW_20_S11_R2_001.fastq.gz GW20_S11_L001_R2_001.fastq.gz
mv GW_21_S12_R1_001.fastq.gz GW21_S12_L001_R1_001.fastq.gz
mv GW_21_S12_R2_001.fastq.gz GW21_S12_L001_R2_001.fastq.gz
mv GW_22_S13_R1_001.fastq.gz GW22_S13_L001_R1_001.fastq.gz
mv GW_22_S13_R2_001.fastq.gz GW22_S13_L001_R2_001.fastq.gz
mv GW_23_S14_R1_001.fastq.gz GW23_S14_L001_R1_001.fastq.gz
mv GW_23_S14_R2_001.fastq.gz GW23_S14_L001_R2_001.fastq.gz
mv GW_24_S15_R1_001.fastq.gz GW24_S15_L001_R1_001.fastq.gz
mv GW_24_S15_R2_001.fastq.gz GW24_S15_L001_R2_001.fastq.gz
mv GW_25_S16_R1_001.fastq.gz GW25_S16_L001_R1_001.fastq.gz
mv GW_25_S16_R2_001.fastq.gz GW25_S16_L001_R2_001.fastq.gz
mv GW_26_S17_R1_001.fastq.gz GW26_S17_L001_R1_001.fastq.gz
mv GW_26_S17_R2_001.fastq.gz GW26_S17_L001_R2_001.fastq.gz
mv GW_27_S18_R1_001.fastq.gz GW27_S18_L001_R1_001.fastq.gz
mv GW_27_S18_R2_001.fastq.gz GW27_S18_L001_R2_001.fastq.gz
mv GW_28_S75_R1_001.fastq.gz GW28_S75_L001_R1_001.fastq.gz
mv GW_28_S75_R2_001.fastq.gz GW28_S75_L001_R2_001.fastq.gz
mv GW_29_S19_R1_001.fastq.gz GW29_S19_L001_R1_001.fastq.gz
mv GW_29_S19_R2_001.fastq.gz GW29_S19_L001_R2_001.fastq.gz
mv GW_2_S65_R1_001.fastq.gz GW2_S65_L001_R1_001.fastq.gz
mv GW_2_S65_R2_001.fastq.gz GW2_S65_L001_R2_001.fastq.gz
mv GW_30_S20_R1_001.fastq.gz GW30_S20_L001_R1_001.fastq.gz
mv GW_30_S20_R2_001.fastq.gz GW30_S20_L001_R2_001.fastq.gz
mv GW_31_S21_R1_001.fastq.gz GW31_S21_L001_R1_001.fastq.gz
mv GW_31_S21_R2_001.fastq.gz GW31_S21_L001_R2_001.fastq.gz
mv GW_32_S22_R1_001.fastq.gz GW32_S22_L001_R1_001.fastq.gz
mv GW_32_S22_R2_001.fastq.gz GW32_S22_L001_R2_001.fastq.gz
mv GW_33_S76_R1_001.fastq.gz GW33_S76_L001_R1_001.fastq.gz
mv GW_33_S76_R2_001.fastq.gz GW33_S76_L001_R2_001.fastq.gz
mv GW_34_S23_R1_001.fastq.gz GW34_S23_L001_R1_001.fastq.gz
mv GW_34_S23_R2_001.fastq.gz GW34_S23_L001_R2_001.fastq.gz
mv GW_35_S77_R1_001.fastq.gz GW35_S77_L001_R1_001.fastq.gz
mv GW_35_S77_R2_001.fastq.gz GW35_S77_L001_R2_001.fastq.gz
mv GW_36_S24_R1_001.fastq.gz GW36_S24_L001_R1_001.fastq.gz
mv GW_36_S24_R2_001.fastq.gz GW36_S24_L001_R2_001.fastq.gz
mv GW_37_S25_R1_001.fastq.gz GW37_S25_L001_R1_001.fastq.gz
mv GW_37_S25_R2_001.fastq.gz GW37_S25_L001_R2_001.fastq.gz
mv GW_38_S26_R1_001.fastq.gz GW38_S26_L001_R1_001.fastq.gz
mv GW_38_S26_R2_001.fastq.gz GW38_S26_L001_R2_001.fastq.gz
mv GW_39_S27_R1_001.fastq.gz GW39_S27_L001_R1_001.fastq.gz
mv GW_39_S27_R2_001.fastq.gz GW39_S27_L001_R2_001.fastq.gz
mv GW_3_S66_R1_001.fastq.gz GW3_S66_L001_R1_001.fastq.gz
mv GW_3_S66_R2_001.fastq.gz GW3_S66_L001_R2_001.fastq.gz
mv GW_40_S28_R1_001.fastq.gz GW40_S28_L001_R1_001.fastq.gz
mv GW_40_S28_R2_001.fastq.gz GW40_S28_L001_R2_001.fastq.gz
mv GW_41_S29_R1_001.fastq.gz GW41_S29_L001_R1_001.fastq.gz
mv GW_41_S29_R2_001.fastq.gz GW41_S29_L001_R2_001.fastq.gz
mv GW_42_S30_R1_001.fastq.gz GW42_S30_L001_R1_001.fastq.gz
mv GW_42_S30_R2_001.fastq.gz GW42_S30_L001_R2_001.fastq.gz
mv GW_43_S31_R1_001.fastq.gz GW43_S31_L001_R1_001.fastq.gz
mv GW_43_S31_R2_001.fastq.gz GW43_S31_L001_R2_001.fastq.gz
mv GW_44_S32_R1_001.fastq.gz GW44_S32_L001_R1_001.fastq.gz
mv GW_44_S32_R2_001.fastq.gz GW44_S32_L001_R2_001.fastq.gz
mv GW_45_S33_R1_001.fastq.gz GW45_S33_L001_R1_001.fastq.gz
mv GW_45_S33_R2_001.fastq.gz GW45_S33_L001_R2_001.fastq.gz
mv GW_46_S34_R1_001.fastq.gz GW46_S34_L001_R1_001.fastq.gz
mv GW_46_S34_R2_001.fastq.gz GW46_S34_L001_R2_001.fastq.gz
mv GW_47_S35_R1_001.fastq.gz GW47_S35_L001_R1_001.fastq.gz
mv GW_47_S35_R2_001.fastq.gz GW47_S35_L001_R2_001.fastq.gz
mv GW_48_S36_R1_001.fastq.gz GW48_S36_L001_R1_001.fastq.gz
mv GW_48_S36_R2_001.fastq.gz GW48_S36_L001_R2_001.fastq.gz
mv GW_49_S37_R1_001.fastq.gz GW49_S37_L001_R1_001.fastq.gz
mv GW_49_S37_R2_001.fastq.gz GW49_S37_L001_R2_001.fastq.gz
mv GW_4_S67_R1_001.fastq.gz GW4_S67_L001_R1_001.fastq.gz
mv GW_4_S67_R2_001.fastq.gz GW4_S67_L001_R2_001.fastq.gz
mv GW_50_S38_R1_001.fastq.gz GW50_S38_L001_R1_001.fastq.gz
mv GW_50_S38_R2_001.fastq.gz GW50_S38_L001_R2_001.fastq.gz
mv GW_51_S39_R1_001.fastq.gz GW51_S39_L001_R1_001.fastq.gz
mv GW_51_S39_R2_001.fastq.gz GW51_S39_L001_R2_001.fastq.gz
mv GW_52_S40_R1_001.fastq.gz GW52_S40_L001_R1_001.fastq.gz
mv GW_52_S40_R2_001.fastq.gz GW52_S40_L001_R2_001.fastq.gz
mv GW_53_S41_R1_001.fastq.gz GW53_S41_L001_R1_001.fastq.gz
mv GW_53_S41_R2_001.fastq.gz GW53_S41_L001_R2_001.fastq.gz
mv GW_54_S42_R1_001.fastq.gz GW54_S42_L001_R1_001.fastq.gz
mv GW_54_S42_R2_001.fastq.gz GW54_S42_L001_R2_001.fastq.gz
mv GW_55_S43_R1_001.fastq.gz GW55_S43_L001_R1_001.fastq.gz
mv GW_55_S43_R2_001.fastq.gz GW55_S43_L001_R2_001.fastq.gz
mv GW_56_S44_R1_001.fastq.gz GW56_S44_L001_R1_001.fastq.gz
mv GW_56_S44_R2_001.fastq.gz GW56_S44_L001_R2_001.fastq.gz
mv GW_57_S45_R1_001.fastq.gz GW57_S45_L001_R1_001.fastq.gz
mv GW_57_S45_R2_001.fastq.gz GW57_S45_L001_R2_001.fastq.gz
mv GW_58_S46_R1_001.fastq.gz GW58_S46_L001_R1_001.fastq.gz
mv GW_58_S46_R2_001.fastq.gz GW58_S46_L001_R2_001.fastq.gz
mv GW_59_S78_R1_001.fastq.gz GW59_S78_L001_R1_001.fastq.gz
mv GW_59_S78_R2_001.fastq.gz GW59_S78_L001_R2_001.fastq.gz
mv GW_5_S68_R1_001.fastq.gz GW5_S68_L001_R1_001.fastq.gz
mv GW_5_S68_R2_001.fastq.gz GW5_S68_L001_R2_001.fastq.gz
mv GW_60_S47_R1_001.fastq.gz GW60_S47_L001_R1_001.fastq.gz
mv GW_60_S47_R2_001.fastq.gz GW60_S47_L001_R2_001.fastq.gz
mv GW_61_S48_R1_001.fastq.gz GW61_S48_L001_R1_001.fastq.gz
mv GW_61_S48_R2_001.fastq.gz GW61_S48_L001_R2_001.fastq.gz
mv GW_62_S79_R1_001.fastq.gz GW62_S79_L001_R1_001.fastq.gz
mv GW_62_S79_R2_001.fastq.gz GW62_S79_L001_R2_001.fastq.gz
mv GW_63_S49_R1_001.fastq.gz GW63_S49_L001_R1_001.fastq.gz
mv GW_63_S49_R2_001.fastq.gz GW63_S49_L001_R2_001.fastq.gz
mv GW_64_S50_R1_001.fastq.gz GW64_S50_L001_R1_001.fastq.gz
mv GW_64_S50_R2_001.fastq.gz GW64_S50_L001_R2_001.fastq.gz
mv GW_65_S51_R1_001.fastq.gz GW65_S51_L001_R1_001.fastq.gz
mv GW_65_S51_R2_001.fastq.gz GW65_S51_L001_R2_001.fastq.gz
mv GW_66_S52_R1_001.fastq.gz GW66_S52_L001_R1_001.fastq.gz
mv GW_66_S52_R2_001.fastq.gz GW66_S52_L001_R2_001.fastq.gz
mv GW_67_S53_R1_001.fastq.gz GW67_S53_L001_R1_001.fastq.gz
mv GW_67_S53_R2_001.fastq.gz GW67_S53_L001_R2_001.fastq.gz
mv GW_68_S54_R1_001.fastq.gz GW68_S54_L001_R1_001.fastq.gz
mv GW_68_S54_R2_001.fastq.gz GW68_S54_L001_R2_001.fastq.gz
mv GW_69_S55_R1_001.fastq.gz GW69_S55_L001_R1_001.fastq.gz
mv GW_69_S55_R2_001.fastq.gz GW69_S55_L001_R2_001.fastq.gz
mv GW_6_S69_R1_001.fastq.gz GW6_S69_L001_R1_001.fastq.gz
mv GW_6_S69_R2_001.fastq.gz GW6_S69_L001_R2_001.fastq.gz
mv GW_70_S56_R1_001.fastq.gz GW70_S56_L001_R1_001.fastq.gz
mv GW_70_S56_R2_001.fastq.gz GW70_S56_L001_R2_001.fastq.gz
mv GW_71_S57_R1_001.fastq.gz GW71_S57_L001_R1_001.fastq.gz
mv GW_71_S57_R2_001.fastq.gz GW71_S57_L001_R2_001.fastq.gz
mv GW_72_S58_R1_001.fastq.gz GW72_S58_L001_R1_001.fastq.gz
mv GW_72_S58_R2_001.fastq.gz GW72_S58_L001_R2_001.fastq.gz
mv GW_73_S59_R1_001.fastq.gz GW73_S59_L001_R1_001.fastq.gz
mv GW_73_S59_R2_001.fastq.gz GW73_S59_L001_R2_001.fastq.gz
mv GW_74_S60_R1_001.fastq.gz GW74_S60_L001_R1_001.fastq.gz
mv GW_74_S60_R2_001.fastq.gz GW74_S60_L001_R2_001.fastq.gz
mv GW_7_S70_R1_001.fastq.gz GW7_S70_L001_R1_001.fastq.gz
mv GW_7_S70_R2_001.fastq.gz GW7_S70_L001_R2_001.fastq.gz
mv GW_8_S71_R1_001.fastq.gz GW8_S71_L001_R1_001.fastq.gz
mv GW_8_S71_R2_001.fastq.gz GW8_S71_L001_R2_001.fastq.gz
mv GW_9_S1_R1_001.fastq.gz GW9_S1_L001_R1_001.fastq.gz
mv GW_9_S1_R2_001.fastq.gz GW9_S1_L001_R2_001.fastq.gz
mv PCR_Neg_1_S63_R1_001.fastq.gz PCRNeg1_S63_L001_R1_001.fastq.gz
mv PCR_Neg_1_S63_R2_001.fastq.gz PCRNeg1_S63_L001_R2_001.fastq.gz
mv PCR_Neg_2_S81_R1_001.fastq.gz PCRNeg2_S81_L001_R1_001.fastq.gz
mv PCR_Neg_2_S81_R2_001.fastq.gz PCRNeg2_S81_L001_R2_001.fastq.gz
mv PCR_Pos_1_S62_R1_001.fastq.gz PCRPos1_S62_L001_R1_001.fastq.gz
mv PCR_Pos_1_S62_R2_001.fastq.gz PCRPos1_S62_L001_R2_001.fastq.gz
mv PCR_Pos_2_S80_R1_001.fastq.gz PCRPos2_S80_L001_R1_001.fastq.gz
mv PCR_Pos_2_S80_R2_001.fastq.gz PCRPos2_S80_L001_R2_001.fastq.gz


##import reads into Qiime2


#load module on HPC
#module load QIIME2/2019.4-Miniconda3

#make output directory
mkdir $wkdir/03_qiime_import_fastq/

#import reads into Qiime2
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/00_raw_fastq/00_cp/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/03_qiime_import_fastq/demux-paired-end.qza


##run dada2 for OTU counts


mkdir $wkdir/04_qiime_dada2/

#run dada2
qiime dada2 denoise-paired --i-demultiplexed-seqs /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/03_qiime_import_fastq/demux-paired-end.qza --o-table /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/04_qiime_dada2/table.qza --o-representative-sequences /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/04_qiime_dada2/rep-seqs.qza --o-denoising-stats /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/04_qiime_dada2/denoising-stats.qza --p-trunc-len-f 180 --p-trunc-len-r 115


#download reference taxonomy for taxon assignment


mkdir $wkdir/08_qiime_reference_taxonomy/
cd $wkdir/08_qiime_reference_taxonomy/
wget https://data.qiime2.org/2022.2/common/silva-138-99-seqs-515-806.qza
wget https://data.qiime2.org/2022.2/common/silva-138-99-tax-515-806.qza


##run feature classifier to assign taxonomy to OTU's


mkdir $wkdir/09_qiime_feature-classifier_classify-consensus-vsearch/

qiime feature-classifier classify-consensus-vsearch \
  --i-query /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/04_qiime_dada2/rep-seqs.qza \
  --i-reference-reads /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/08_qiime_reference_taxonomy/silva-138-99-seqs-515-806.qza \
  --i-reference-taxonomy /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/08_qiime_reference_taxonomy/silva-138-99-tax-515-806.qza \
  --o-classification /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/09_qiime_feature-classifier_classify-consensus-vsearch/taxonomy.qza


##Generate a tree for phylogenetic diversity analyses


mkdir $wkdir/11_B_qiime_phylogeny_align-to-tree-mafft-fasttree/


qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/04_qiime_dada2/rep-seqs.qza \
  --o-alignment /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/11_B_qiime_phylogeny_align-to-tree-mafft-fasttree/aligned-rep-seqs.qza \
  --o-masked-alignment /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/11_B_qiime_phylogeny_align-to-tree-mafft-fasttree/masked-aligned-rep-seqs.qza \
  --o-tree /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/11_B_qiime_phylogeny_align-to-tree-mafft-fasttree/unrooted-tree.qza \
  --o-rooted-tree /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/11_B_qiime_phylogeny_align-to-tree-mafft-fasttree/rooted-tree.qza



##Convert .qza files to file formats that I can actually use (csv, tsv, nwk)

mkdir $wkdir/05_qiime_export_dada2/

cd $wkdir/04_qiime_dada2/

qiime tools export \
  --input-path table.qza \
  --output-path /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/05_qiime_export_dada2/exported-table

#make the feature table readable
cd $wkdir/05_qiime_export_dada2/exported-table
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv

    #feature-table.tsv was converted to feature-table.csv ; this is the feature table used for downstream analyses.

#get usable nwk files
mkdir $wkdir/11_B_qiime_phylogeny_align-to-tree-mafft-fasttree/export/
cd $wkdir/11_B_qiime_phylogeny_align-to-tree-mafft-fasttree/
for i in *; do qiime tools export --input-path $i --output-path /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/11_B_qiime_phylogeny_align-to-tree-mafft-fasttree/export/export/$i; done
	#the rooted tree tree.nwk was renamed to tree.tree for phyloseq and used for downstream analyses


#export taxonomy
mkdir $wkdir/09_qiime_feature-classifier_classify-consensus-vsearch/export/
cd $wkdir/09_qiime_feature-classifier_classify-consensus-vsearch/
qiime tools export --input-path taxonomy.qza --output-path /scratch/gcwoodruff/16S_fig_metabarcoding_3-2022/09_qiime_feature-classifier_classify-consensus-vsearch/export/
  #taxonomy.tsv was used for downstream analyses


#further analyses can be found in workflow.R





