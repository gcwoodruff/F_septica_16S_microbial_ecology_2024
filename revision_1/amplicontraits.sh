cd /Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/amplicontraits

mkdir traitdatabases
cd traitdatabases

wget https://erda.ku.dk/archives/f5d4b1d41f74ba3d6f73b212dbb11591/traitDatabases.tar

tar -xvf traitDatabases.tar

for i in *; do gunzip $i; done


#downloaded usearch5.2.236_i86osx32 from https://drive5.com/usearch/download.html
#renamed usearch
 
cd /Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/amplicontraits/

chmod +x usearch

export PATH=$PATH:/Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/amplicontraits/


arch -x86_64 usearch
#this dont work on my mac
#okay... again, need linux and oscer......

scp -r /Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/amplicontraits/ gcwoodruff@schooner.oscer.ou.edu:/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/amplicontraits_8-8-2025


cd /home/gcwoodruff/download/

wget https://drive5.com/downloads/usearch11.0.667_i86linux32.gz

gunzip usearch11.0.667_i86linux32.gz
mv usearch11.0.667_i86linux32 usearch
chmod +x usearch
export PATH=$PATH:/home/gcwoodruff/download/
#okay it works, great


cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/amplicontraits_8-8-2025/amplicontraits/


scp /Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/04_qiime_dada2/rep-seqs_fa/dna-sequences.fasta gcwoodruff@schooner.oscer.ou.edu:/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/amplicontraits_8-8-2025/amplicontraits/dna-sequences.fasta


#let's try one...
mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/amplicontraits_8-8-2025/amplicontraits/00_usearch_out/

wkdir='/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/amplicontraits_8-8-2025/amplicontraits/'
cd $wkdir/


usearch -sinaps dna-sequences.fasta -db $wkdir/traitdatabases/genome_size10.fasta -attr genome_size10 -tabbedout $wkdir/00_usearch_out/genome_size10.txt -strand plus


#okay, cool, just run it!
.txt
cd $wkdir/traitdatabases/

for i in *; do usearch -sinaps $wkdir/dna-sequences.fasta -db $i -attr ${i%.fasta} -tabbedout $wkdir/00_usearch_out/${i%.fasta}.txt -strand plus; done

usearch -sinaps $wkdir/dna-sequences.fasta -db $wkdir/traitdatabases/cellShape.fasta -attr cellshape -tabbedout $wkdir/00_usearch_out/cellshape.txt -strand plus


#woo, we did it!


scp -r gcwoodruff@schooner.oscer.ou.edu:/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/amplicontraits_8-8-2025/amplicontraits/00_usearch_out/ /Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/amplicontraits/



#the elegans asv's..........

wkdir='/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/amplicontraits_8-8-2025/amplicontraits/'

mkdir $wkdir/01_elegans/
mkdir $wkdir/01_elegans/00_usearch_out/


scp /Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/elegans_included_analysis/02_qiime_dada2/rep-seqs_fa/dna-sequences.fasta gcwoodruff@schooner.oscer.ou.edu:/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/amplicontraits_8-8-2025/amplicontraits/01_elegans/dna-sequences.fasta

mv $wkdir/traitdatabases/cellShape.fasta $wkdir/traitdatabases/cellshape.fasta


cd $wkdir/traitdatabases/

for i in *; do usearch -sinaps $wkdir/01_elegans/dna-sequences.fasta -db $i -attr ${i%.fasta} -tabbedout $wkdir/01_elegans/00_usearch_out/${i%.fasta}.txt -strand plus; done


scp -r gcwoodruff@schooner.oscer.ou.edu:/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/amplicontraits_8-8-2025/amplicontraits/01_elegans/00_usearch_out/ /Users/gavin/genome/microbiology_spectrum_revisions_16S_2025/amplicontraits/01_elegans/





