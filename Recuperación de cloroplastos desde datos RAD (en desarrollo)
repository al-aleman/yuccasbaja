# Homemade "fast" workflow to obtain potential plastid DNA from Next-Seq (or Rad-Seq) data

# First, we remove adapters and quality-trim sequences. We are used to utilize BBtools and force sequences to 100 bp because we use Stacks for denovo assembly, and any other tool to remove adapters is ok aswell. 

# Loop for clean-and-trim
for i in $(ls *.fastq.gz | sed 's/.fastq.gz//') 
do 
bbduk.sh in=$i\.fastq.gz out=$i\_cleaned.fastq.gz ktrim=r k=17 hdist=1 mink=8 ref=<path/to/adapters_or_contaminant_files> minlen=100 ow=t qtrim=rl trimq=10 forcetrimright=99
done

# Then, we download the reference genome that will serve as a template for the alignment of plastid-putative sequences (in our case and example, we used Yucca schidigera NCBI reference) and create a bowtie2 index database

bowtie2-build Yucshi.fa yucca

# result: 6 .bt2 database files
$ ls 
yucca.1.bt2
yucca.2.bt2
yucca.3.bt2
yucca.4.bt2
yucca.rev.1.bt2
yucca.rev.2.bt2

# Next, we align plastid-putative sequences from the cleaned fastq.gz files (which will result in a .sam file per sample)
for i in $(ls *.fastq.gz | sed 's/.fastq.gz//')
do 
bowtie2 -x <path/to/bowtie2/index/yucca> -U $i\.fastq.gz -S ./pathtoresults/$i.sam -p [threads]
done

# And transform the .sam files into binary
for i in $(ls *.sam | sed 's/.sam//')
do 
samtools view -S -b $i\.sam -@ [threads] | samtools sort -@ [threads] > $i.bam
done

# The next step is to obtain a fasta sequence from each sample via:
for i in 
$(ls *.bam | sed 's/.bam//') 
do 
samtools mpileup -uf <path/to/fasta/reference/yucca_schidigera.fa> $i.bam | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq | seqtk seq -a > $i.fa
done

# Although each file will respect its original name, the header inside each file will be something like >NC_032714.1, refering to the NCBI reference. We recommend to modify them by:
sed -i 's/>.*/>nameofthesample/g' nameofthefile.fa

# After that, we merge all the FASTA into one file

cat *.fa > <nameforthemultifasta>.fa

# And complete the tail of each sample to measure the same as our reference (the number may change of you use a different reference)

perl -e 'while(<>) {if (/^>/) {print;next}; chomp; print $_,"n" x (156158-length($_)),"\n"}'

# The alignment may be plot and visualized in R. Based on our results, we do not recommend the use of this information for intraspecific analysis (due to the levels of missing data you may experience), but if want to give it a try, one way could be removing the individuals with excess of missing data and then looking for variant positions, duties that may be done with tools like sortbyname.sh and snp-sites). 

# The we obtained the consensus sequence for the chloroplast of the species we are interested in
mothur 
> consensus.seqs(fasta=<nameforthemultifasta>.fa)

# The consensus is something like this (very short example):
# AGTCCGATCGGATCGATGCAGCTAGNNNNNNNNNNNNNTTTAGCTAGGATCGATCGGATCGAAATAGCATAGCTAGGATCGATTAGCTAYGCTAGCTAGGCTAGATCGAGCTAGCTAGSAGCTAGCTAGATCGWAGCTAGATCGATAKNNNNNNAGCTAGCTAGATCGGATCGATAGCTAGR
# Hopefully, each sequence between gaps ("N") measures more than 1000 bp. Then we extract each sequence in order that we obtain something like this:
>Consensus1
AGTCCGATCGGATCGATGCAGCTAG
>Consensus2
TTTAGCTAGGATCGATCGGATCGAAATAGCATAGCTAGGATCGATTAGCTAYGCTAGCTAGGCTAGATCGAGCTAGCTAGSAGCTAGCTAGATCGWAGCTAGATCGATAK
>Consensus3
AGCTAGCTAGATCGGATCGATAGCTAGR

# In our case, we obtained about 26 consensus sequences, each one of at least 1000 pb (we removed the ones with less than that). We checked for representation of each consensus on the multi-fasta alignment, ir order that we are sure it came from a significant number of samples.

# Next, we blasted the consensus sequences. We constructed a Blast database with a) the four species of Yucca on NCBI and Hesperoyucca whipplei and b) the consensus output of mothur via:
makeblastdb -in <referencesandonelineconsensusfile>.fasta -parse_seqids -dbtype nucl

# And blasted the consensus sequences
blastn -query <multilineconsensussequences>.fasta -subject <referencesandonelineconsensusfile>.fasta -outfmt "6 sseqid length sstart send sseq" -out results.tsv

# The format 6 of blast results outputs the aligned part of subject (reference) sequence (sseq option). Then we open the results.tsv file on any spreadsheet program. We verify that we obtain an alignment for each consensus*reference, and we removed the sequences that are reverse-aligned. Finally, we extract and concatenate the aligned parts for each reference. In our case, we obtained, approximately 100,000 pb per species. Then we use Clustalw to perfectly-align the concatenated sequences.
# We used the final output to construct a Maximum-Likelihood tree.
