#Installation: 
###Trim_Galore
#check for fastqc and cutadapt
fatsqc -v
cutadapt --version
#If cutadapt not there - sudo apt install cutadapt
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
#Edit ~/.bashrc to set alias trim_galore='Path/to/trim_galore' and source it
trim_galore --help

###AGAT
conda install bioconda::agat
agat_convert_sp_gff2gtf.pl -gff PF3D7_annotation.gff -o PF3D7_annotation.gtf

###STAR
#Installation:
conda install bioconda::star
STAR --help

###Qualimap
conda install bioconda::qualimap


#Generating Genome:
STAR --runMode genomeGenerate --genomeFastaFiles PF3D7_genome.fasta --genomeDir ./PfGenome

gzip ring_stage.fastq
trim_galore --illumina ring_stage.fastq.gz -o ring_stage
cd ring_stage/
gzip -d ring_stage_trimmed.fq.gz
STAR --runThreadN 4 --genomeDir ../PF3D7Genome/ --readFilesIn ring_stage_trimmed.fq --outFileNamePrefix ring_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard=
htseq-count ring_Aligned.sortedByCoord.out.bam ../PF3D7_annotation.gtf --idattr gene_id -t mRNA -f bam > ring_count.txt
grep "PF3D7" ring_count.txt | awk '{print $1}' > ring_genes.txt
INFILE=ring_genes.txt
while read -r LINE
do
	gene="ID=$LINE;"
	grep $gene ../PF3D7_annotation.gff | awk '{print $5 - $4}' >> ring_genelengths.txt
done < "$INFILE"
grep "PF3D7" ring_count.txt | awk '{print $1 "," $2}' > ring_count.csv
awk '{print $1}' ring_genelengths.txt > ring_genelengths.csv
#Ensure number of lines in both ring_count.csv and ring_genelengths.csv is the same

###Qualimap
qualimap rnaseq -bam ring_Aligned.sortedByCoord.out.bam -gtf ../PF3D7_annotation.gtf -outfile ring_qualimap_report.pdf

###R:
library('DESeq2')
count <- read.csv("ring_count.csv",header=FALSE,row.names=1)
gene_lengths <- read.csv("ring_genelengths.csv", header=FALSE)
rownames(gene_lengths) <- rownames(count)
calculateTPM <- function(counts, lengths){
counts_per_base <- counts / lengths;
scaling_factors <- colSums(counts_per_base);
tpm <- t(t(counts_per_base) / scaling_factors) * 1e6;
return(round(tpm,4))}
tpm_matrix <- calculateTPM(count,gene_lengths)
write.csv(tpm_matrix, "ring_TPM.csv")



