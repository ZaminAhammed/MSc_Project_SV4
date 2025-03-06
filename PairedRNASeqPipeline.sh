trim_galore --illumina --dont_gzip --paired ES1_1.fastq.gz ES1_2.fastq.gz ES2_1.fastq.gz ES2_2.fastq.gz ES3_1.fastq.gz ES3_2.fastq.gz ET1_1.fastq.gz ET1_2.fastq.gz ET2_1.fastq.gz ET2_2.fastq.gz ET3_1.fastq.gz ET3_2.fastq.gz LS1_1.fastq.gz LS1_2.fastq.gz LS2_1.fastq.gz LS2_2.fastq.gz LS3_1.fastq.gz LS3_2.fastq.gz LT1_1.fastq.gz LT1_2.fastq.gz LT2_1.fastq.gz LT2_2.fastq.gz LT3_1.fastq.gz LT3_2.fastq.gz M1_1.fastq.gz M1_2.fastq.gz M2_1.fastq.gz M2_2.fastq.gz M3_1.fastq.gz M3_2.fastq.gz MT1_1.fastq.gz MT1_2.fastq.gz MT2_1.fastq.gz MT2_2.fastq.gz MT3_1.fastq.gz MT3_2.fastq.gz S1_1.fastq.gz S1_2.fastq.gz S2_1.fastq.gz S2_2.fastq.gz S3_1.fastq.gz S3_2.fastq.gz YR1_1.fastq.gz YR1_2.fastq.gz YR2_1.fastq.gz YR2_2.fastq.gz YR3_1.fastq.gz YR3_2.fastq.gz

STAR --runMode alignReads --runThreadN 4 --genomeDir ../PF3D7Genome/ --readFilesIn ES1_1_val_1.fq ES1_2_val_1.fq --outFileNamePrefix  ES1_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard

htseq-count ES1_Aligned.sortedByCoord.out.bam ../PF3D7_annotation.gtf --idattr gene_id -t mRNA -f bam > ES1_count.txt

grep "PF3D7" ES1_count.txt | awk '{print $1}' > genes.txt

INFILE=genes.txt
while read -r LINE
do
	gene="ID=$LINE;"
	grep $gene ../PF3D7_annotation.gff | awk '{print $5 - $4}' >> genelengths.txt
done < "$INFILE"

grep "PF3D7" ES1_count.txt | awk '{print $1 "," $2}' > ES1_count.csv

#Put ES1, ES2 and ES3 counts to one (ES_count.csv) 

count <- read.csv("ES_count.csv",header=FALSE,row.names=1)
gene_lengths <- read.csv("genelengths.txt", header=FALSE)
rownames(gene_lengths) <- rownames(count)
calculateTPM <- function(counts, lengths){
counts_per_base <- counts / lengths;
scaling_factors <- colSums(counts_per_base);
tpm <- t(t(counts_per_base) / scaling_factors) * 1e6;
return(round(tpm,4))}
tpm_matrix <- calculateTPM(count,gene_lengths)
write.csv(tpm_matrix, "ring_TPM.csv")

