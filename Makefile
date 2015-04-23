# By default, run all steps. Triggers download if FastQ files are not found in 'input' dir.
all: check_bin trimmomatic ray_meta prodigal bowtie2-build bowtie2-run samtools-index bedtools-multicov

# Number of threads for assembly and mapping. Please adjust according to your system.
THREADS=32

###############
# Check progs
###############

.PHONY .SILENT: check_bin
check_bin:
	command -v wget >/dev/null 2>&1 || { echo >&2 "wget not found. Aborting."; exit 1; }
	command -v java >/dev/null 2>&1 || { echo >&2 "wget not found. Aborting."; exit 1; }
	command -v mpiexec >/dev/null 2>&1 || { echo >&2 "mpiexec not found. Aborting."; exit 1; }
	command -v Ray >/dev/null 2>&1 || { echo >&2 "Ray not found. Aborting."; exit 1; }
	command -v prodigal >/dev/null 2>&1 || { echo >&2 "prodigal not found. Aborting."; exit 1; }
	command -v bowtie2-build >/dev/null 2>&1 || { echo >&2 "bowtie2-build not found. Aborting."; exit 1; }
	command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "bowtie2 not found. Aborting."; exit 1; }
	command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools not found. Aborting."; exit 1; }
	command -v bedtools >/dev/null 2>&1 || { echo >&2 "bedtools not found. Aborting."; exit 1; }

###############
# Download data
###############

.PHONY: download
download: input/GAIIx_Lane6_R1.fastq.gz input/GAIIx_Lane6_R2.fastq.gz input/GAIIx_Lane7_R1.fastq.gz input/GAIIx_Lane7_R2.fastq.gz input/GAIIx_Lane8_R1.fastq.gz input/GAIIx_Lane8_R2.fastq.gz input/MiSeq_A1_R1.fastq.gz input/MiSeq_A1_R2.fastq.gz input/MiSeq_A2_R1.fastq.gz input/MiSeq_A2_R2.fastq.gz input/MiSeq_B1_R1.fastq.gz input/MiSeq_B1_R2.fastq.gz input/MiSeq_B2_R1.fastq.gz input/MiSeq_B2_R2.fastq.gz

input:
	mkdir -p input

input/GAIIx_Lane6_R1.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_6_1_sequence.txt.gz

input/GAIIx_Lane6_R2.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_6_2_sequence.txt.gz

input/GAIIx_Lane7_R1.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_7_1_sequence.txt.gz

input/GAIIx_Lane7_R2.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_7_2_sequence.txt.gz

input/GAIIx_Lane8_R1.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_8_1_sequence.txt.gz

input/GAIIx_Lane8_R2.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_8_2_sequence.txt.gz

input/MiSeq_A1_R1.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq_Biogas1_S1_L001_R1_001.fastq.gz

input/MiSeq_A1_R2.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq_Biogas1_S1_L001_R2_001.fastq.gz

input/MiSeq_A2_R1.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq_Biogas2_S2_L001_R1_001.fastq.gz

input/MiSeq_A2_R2.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq_Biogas2_S2_L001_R2_001.fastq.gz

input/MiSeq_B1_R1.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq2_Biogas1_S1_L001_R1_001.fastq.gz

input/MiSeq_B1_R2.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq2_Biogas1_S1_L001_R2_001.fastq.gz

input/MiSeq_B2_R1.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq2_Biogas2_S2_L001_R1_001.fastq.gz

input/MiSeq_B2_R2.fastq.gz: input
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq2_Biogas2_S2_L001_R2_001.fastq.gz

###############
# Read QC
###############

.PHONY: trimmomatic GAIIx_Lane6 GAIIx_Lane7 GAIIx_Lane8 MiSeq_A1 MiSeq_A2 MiSeq_B1 MiSeq_B2
# Trimmomatic plus adapter files must be in current working directory. Source: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
trimmomatic: trimmomatic-0.32.jar TruSeq2-PE.fa NexteraPE-PE.fa GAIIx_Lane6 GAIIx_Lane7 GAIIx_Lane8 MiSeq_A1 MiSeq_A2 MiSeq_B1 MiSeq_B2

GAIIx_Lane6: GAIIx_Lane6.trimmomatic_1P.fastq.gz GAIIx_Lane6.trimmomatic_1U.fastq.gz GAIIx_Lane6.trimmomatic_2P.fastq.gz GAIIx_Lane6.trimmomatic_2U.fastq.gz
GAIIx_Lane7: GAIIx_Lane7.trimmomatic_1P.fastq.gz GAIIx_Lane7.trimmomatic_1U.fastq.gz GAIIx_Lane7.trimmomatic_2P.fastq.gz GAIIx_Lane7.trimmomatic_2U.fastq.gz
GAIIx_Lane8: GAIIx_Lane8.trimmomatic_1P.fastq.gz GAIIx_Lane8.trimmomatic_1U.fastq.gz GAIIx_Lane8.trimmomatic_2P.fastq.gz GAIIx_Lane8.trimmomatic_2U.fastq.gz
MiSeq_A1: MiSeq_A1.trimmomatic_1P.fastq.gz MiSeq_A1.trimmomatic_1U.fastq.gz MiSeq_A1.trimmomatic_2P.fastq.gz MiSeq_A1.trimmomatic_2U.fastq.gz
MiSeq_A2: MiSeq_A2.trimmomatic_1P.fastq.gz MiSeq_A2.trimmomatic_1U.fastq.gz MiSeq_A2.trimmomatic_2P.fastq.gz MiSeq_A2.trimmomatic_2U.fastq.gz
MiSeq_B1: MiSeq_B1.trimmomatic_1P.fastq.gz MiSeq_B1.trimmomatic_1U.fastq.gz MiSeq_B1.trimmomatic_2P.fastq.gz MiSeq_B1.trimmomatic_2U.fastq.gz
MiSeq_B2: MiSeq_B2.trimmomatic_1P.fastq.gz MiSeq_B2.trimmomatic_1U.fastq.gz MiSeq_B2.trimmomatic_2P.fastq.gz MiSeq_B2.trimmomatic_2U.fastq.gz

GAIIx_%.trimmomatic_1P.fastq.gz GAIIx_%.trimmomatic_1U.fastq.gz GAIIx_%.trimmomatic_2P.fastq.gz GAIIx_%.trimmomatic_2U.fastq.gz: input/GAIIx_%_R1.fastq.gz input/GAIIx_%_R2.fastq.gz
	java -jar trimmomatic-0.32.jar PE -baseout GAIIx_$*.trimmomatic.fastq.gz $^ ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 TOPHRED33

MiSeq_%.trimmomatic_1P.fastq.gz MiSeq_%.trimmomatic_1U.fastq.gz MiSeq_%.trimmomatic_2P.fastq.gz MiSeq_%.trimmomatic_2U.fastq.gz: input/MiSeq_%_R1.fastq.gz input/MiSeq_%_R2.fastq.gz
	java -jar trimmomatic-0.32.jar PE -baseout MiSeq_$*.trimmomatic.fastq.gz $^ ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 TOPHRED33

###############
# Assembly (multithreaded)
###############

.PHONY: ray_meta
ray_meta: Contigs_gt1kb.fasta

RayMeta_k31/Contigs.fasta: GAIIx_Lane7 GAIIx_Lane8 MiSeq_A1 MiSeq_A2 MiSeq_B1 MiSeq_B2
	mpiexec -n $(THREADS) Ray -k 31 -p GAIIx_Lane7.trimmomatic_1P.fastq.gz GAIIx_Lane7.trimmomatic_2P.fastq.gz -s GAIIx_Lane7.trimmomatic_1U.fastq.gz -s GAIIx_Lane7.trimmomatic_2U.fastq.gz -p GAIIx_Lane8.trimmomatic_1P.fastq.gz GAIIx_Lane8.trimmomatic_2P.fastq.gz -s GAIIx_Lane8.trimmomatic_1U.fastq.gz -s GAIIx_Lane8.trimmomatic_2U.fastq.gz -p MiSeq_A1.trimmomatic_1P.fastq.gz MiSeq_A1.trimmomatic_2P.fastq.gz -s MiSeq_A1.trimmomatic_1U.fastq.gz -s MiSeq_A1.trimmomatic_2U.fastq.gz -p MiSeq_A2.trimmomatic_1P.fastq.gz MiSeq_A2.trimmomatic_2P.fastq.gz -s MiSeq_A2.trimmomatic_1U.fastq.gz -s MiSeq_A2.trimmomatic_2U.fastq.gz -p MiSeq_B1.trimmomatic_1P.fastq.gz MiSeq_B1.trimmomatic_2P.fastq.gz -s MiSeq_B1.trimmomatic_1U.fastq.gz -s MiSeq_B1.trimmomatic_2U.fastq.gz -p MiSeq_B2.trimmomatic_1P.fastq.gz MiSeq_B2.trimmomatic_2P.fastq.gz -s MiSeq_B2.trimmomatic_1U.fastq.gz -s MiSeq_B2.trimmomatic_2U.fastq.gz -o RayMeta_k31 -minimum-contig-length 1000

Contigs_gt1kb.fasta: RayMeta_k31/Contigs.fasta
	cp $^ $@

###############
# Gene prediction
###############

.PHONY: prodigal
prodigal: Contigs_gt1kb.prodigal.faa Contigs_gt1kb.prodigal.fna Contigs_gt1kb.prodigal.gff

%.prodigal.faa %.prodigal.fna %.prodigal.gff: Contigs_gt1kb.fasta
	prodigal -p meta -a Contigs_gt1kb.prodigal.faa -d Contigs_gt1kb.prodigal.fna -f gff -o Contigs_gt1kb.prodigal.gff -i $^

###############
# BLASTP against KEGG (fake)
###############

.PHONY .SILENT: kegg-blastp
kegg-blastp:
	echo 'We distributed BLASTP-jobs on our compute cluster. A simplified command, reproducing the results but running several weeks, looks like this:'
	echo 'blastp -max_target_seqs 1 -outfmt 6 -query Contigs_gt1kb.prodigal.faa -out Contigs_gt1kb.prodigal.faa.blastout.tab -db /path/to/kegg/blastdb'

###############
# Read mapping (multitreaded)
###############

.PHONY: bowtie2-build bowtie2-run samtools-index
bowtie2-build: Contigs_gt1kb.fasta.1.bt2 Contigs_gt1kb.fasta.2.bt2 Contigs_gt1kb.fasta.3.bt2 Contigs_gt1kb.fasta.4.bt2 Contigs_gt1kb.fasta.rev.1.bt2 Contigs_gt1kb.fasta.rev.2.bt2
bowtie2-run: GAIIx_Lane6.trimmomatic.bam GAIIx_Lane7.trimmomatic.bam GAIIx_Lane8.trimmomatic.bam MiSeq_A1.trimmomatic.bam MiSeq_A2.trimmomatic.bam MiSeq_B1.trimmomatic.bam MiSeq_B2.trimmomatic.bam
samtools-index: GAIIx_Lane6.trimmomatic.bam.bai GAIIx_Lane7.trimmomatic.bam.bai GAIIx_Lane8.trimmomatic.bam.bai MiSeq_A1.trimmomatic.bam.bai MiSeq_A2.trimmomatic.bam.bai MiSeq_B1.trimmomatic.bam.bai MiSeq_B2.trimmomatic.bam.bai

%.fasta.1.bt2 %.fasta.2.bt2 %.fasta.3.bt2 %.fasta.4.bt2 %.fasta.rev.1.bt2 %.fasta.rev.2.bt2: %.fasta
	bowtie2-build $^ $^

%.trimmomatic.bam: %.trimmomatic_1P.fastq.gz %.trimmomatic_2P.fastq.gz
	bowtie2 -X 1000 -p $(THREADS) --end-to-end --sensitive -x Contigs_gt1kb.fasta -1 $*.trimmomatic_1P.fastq.gz -2 $*.trimmomatic_2P.fastq.gz | samtools view -uT Contigs_gt1kb.fasta - | samtools sort - $*.trimmomatic

%.bam.bai: %.bam
	samtools index $^

###############
# Count reads in genes
###############

.PHONY: bedtools-multicov
bedtools-multicov: Contigs_gt1kb.prodigal.gff Contigs_gt1kb.prodigal.gff.bedtools.tsv samtools-index

Contigs_gt1kb.prodigal.gff.bedtools.tsv: GAIIx_Lane6.trimmomatic.bam GAIIx_Lane7.trimmomatic.bam GAIIx_Lane8.trimmomatic.bam MiSeq_A1.trimmomatic.bam MiSeq_A2.trimmomatic.bam MiSeq_B1.trimmomatic.bam MiSeq_B2.trimmomatic.bam
	bedtools multicov -bams $^ -bed Contigs_gt1kb.prodigal.gff > $@
