Below is the simple variant annotation example using 2 exome datasets from SRA.

SAMPLES

https://www.ncbi.nlm.nih.gov/sra/SRX25324323[accn]
https://www.ncbi.nlm.nih.gov/sra/SRX25324326[accn]

5 min

DOWNLOAD

https://sra-pub-src-1.s3.amazonaws.com/SRR29825555/Preet_JME_1_R1.fastq.gz.1
https://sra-pub-src-1.s3.amazonaws.com/SRR29825555/Preet_JME_1_R2.fastq.gz.1
https://sra-pub-src-1.s3.amazonaws.com/SRR29825552/Kritika_4_R1.fastq.gz.1
https://sra-pub-src-1.s3.amazonaws.com/SRR29825552/Kritika_4_R2.fastq.gz.1

30 sec

UNCOMPRESSED

gunzip *gz

<1 sec

RENAMING

mv Preet_JME_1_R1.fastq.gz.1 sample1_R1.fq
mv Preet_JME_1_R2.fastq.gz.1 sample1_R1.fq
mv Kritika_4_R1.fastq sample2_R1.fq
mv Kritika_4_R2.fastq sample2_R2.fq

<30 sec

DATASTAT

(seqkit download:https://github.com/shenwei356/seqkit)

seqkit stat *fq
processed files:  4 / 4 [======================================] ETA: 0s. done
file           format  type  num_seqs    sum_len  min_len  avg_len  max_len
sample1_R1.fq  FASTQ   DNA     10,319  1,558,169      151      151      151
sample1_R2.fq  FASTQ   DNA     10,319  1,558,169      151      151      151
sample2_R1.fq  FASTQ   DNA     10,237  1,545,787      151      151      151
sample2_R2.fq  FASTQ   DNA     10,237  1,545,787      151      151      151

<1 sec

UCSC GENOME DOWNLOAD

https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/

~30min

genome.fasta

QC CHECK

fastqc download
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

fastqc *fq

<10 sec

ADAPTER REMOVAL & FILTERATION

Trim-Galore download:
https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

trim_galore --paired sample1_R1.fq sample1_R2.fq --fastqc --length 36
trim_galore --paired sample2_R1.fq sample2_R2.fq --fastqc --length 36

<30 sec

INDEXING REFERENCE USING BWA

(bwa download:https://github.com/lh3/bwa)
bwa index genome.fasta genome.fasta

~3 Hour

BWA MAPPING

bwa mem genome.fasta sample1_R1.fq sample1_R2.fq >out1_sample1_bwa_mapped.sam
bwa mem genome.fasta sample2_R1.fq sample2_R2.fq >out1_sample2_bwa_mapped.sam

<1 min

MAPPING STAT

(samtools download:http://www.htslib.org/download/)

SAMTOBAM

samtools sort -o out1_sample1_bwa_mapped.bam out1_sample1_bwa_mapped.sam
samtools sort -o out1_sample2_bwa_mapped.bam out1_sample2_bwa_mapped.sam

<5 sec

STAT

samtools flagstat out1_sample1_bwa_mapped.bam
samtools flagstat out1_sample2_bwa_mapped.bam

<5 sec

MARKDUPLICATES

(picard download:https://broadinstitute.github.io/picard/)
picard MarkDuplicates I=out1_sample1_bwa_mapped.bam O=out2_sample1_dedup.bam M=out2_sample1_dup_metrics.txt
picard MarkDuplicates I=out1_sample2_bwa_mapped.bam O=out2_sample2_dedup.bam M=out2_sample2_dup_metrics.txt

<30 sec

GENOME FASTA INDEX

samtools faidx genome.fasta

<30 min

.DICT FILE GENERATION

(GATK download:https://github.com/broadinstitute/gatk/releases)
gatk CreateSequenceDictionary -R genome.fasta -O genome.dict

<5 min

ADDING RG

picard AddOrReplaceReadGroups I=out2_sample1_dedup.bam O=out3_sample1_RG.bam SORT_ORDER=coordinate RGID=S1 RGLB=bar RGPL=illumina RGPU=Lane1 RGSM=Sample1 CREATE_INDEX=True
picard AddOrReplaceReadGroups I=out2_sample2_dedup.bam O=out3_sample2_RG.bam SORT_ORDER=coordinate RGID=S2 RGLB=bar RGPL=illumina RGPU=Lane1 RGSM=Sample2 CREATE_INDEX=True
<5 min

VARIANT CALLING

gatk Mutect2 -R genome.fasta -I out3_sample1_RG.bam -O out4_sample1.vcf.gz
gatk Mutect2 -R genome.fasta -I out3_sample2_RG.bam -O out4_sample2.vcf.gz

<30 min

SNP

gatk SelectVariants -V out4_sample1.vcf.gz --select-type-to-include SNP -O out4_sample1_snps.vcf
gatk SelectVariants -V out4_sample2.vcf.gz --select-type-to-include SNP -O out4_sample2_snps.vcf

<30 sec

INDEL

gatk SelectVariants -V out4_sample1.vcf.gz --select-type-to-include INDEL -O out4_sample1_indels.vcf
gatk SelectVariants -V out4_sample2.vcf.gz --select-type-to-include INDEL -O out4_sample2_indels.vcf

<30 sec

COMPARISON

sample1 and sample2 are having common SNP entries for positions 165990524, 166039668, 166046678, 166073704, 166127146, 166203939, 166226468, 166228614, 166233532, 166242648, 166243328, 166243736, 166259266, 166293183, 166305767, 166306533, 166306689, 166307153, 166311382, 166311407, 166311583, 166311954, 166321357, 166321478 and 166321645.
sample1 and sample2 are having common indel entries for positions 166204470, 166228569, 166228631, 166233203, 166321522, 166041479 and 166052612.

<5 min

ANNOTATION

java -jar snpEff.jar GRCh37.75 out4_sample1.vcf >out5_sample1_annotated.vcf
java -jar snpEff.jar GRCh37.75 out4_sample2.vcf >out5_sample2_annotated.vcf

<30 min

