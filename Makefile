all: check_bin trimmomatic ray_meta prodigal bowtie2

###############
# Check progs
###############

.PHONY .SILENT: check_bin
check_bin:
	command -v mpiexec >/dev/null 2>&1 || { echo >&2 "I require mpiexec but it's not installed. Aborting."; exit 1; }
	command -v Ray >/dev/null 2>&1 || { echo >&2 "I require Ray but it's not installed. Aborting."; exit 1; }
	command -v prodigal >/dev/null 2>&1 || { echo >&2 "I require prodigal but it's not installed. Aborting."; exit 1; }
	command -v bowtie2-build >/dev/null 2>&1 || { echo >&2 "I require bowtie2-build but it's not installed. Aborting."; exit 1; }
	command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "I require bowtie2 but it's not installed. Aborting."; exit 1; }
	command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed. Aborting."; exit 1; }

###############
# Read QC
###############

.PHONY: trimmomatic GAIIx_Lane6 GAIIx_Lane7 GAIIx_Lane8 MiSeq_A1 MiSeq_A2 MiSeq_B1 MiSeq_B2
trimmomatic: trimmomatic-0.32.jar TruSeq2-PE.fa NexteraPE-PE.fa GAIIx_Lane6 GAIIx_Lane7 GAIIx_Lane8 MiSeq_A1 MiSeq_A2 MiSeq_B1 MiSeq_B2

GAIIx_Lane6: GAIIx_Lane6.trimmomatic_1P.fastq.gz GAIIx_Lane6.trimmomatic_1U.fastq.gz GAIIx_Lane6.trimmomatic_2P.fastq.gz GAIIx_Lane6.trimmomatic_2U.fastq.gz
GAIIx_Lane7: GAIIx_Lane7.trimmomatic_1P.fastq.gz GAIIx_Lane7.trimmomatic_1U.fastq.gz GAIIx_Lane7.trimmomatic_2P.fastq.gz GAIIx_Lane7.trimmomatic_2U.fastq.gz
GAIIx_Lane8: GAIIx_Lane8.trimmomatic_1P.fastq.gz GAIIx_Lane8.trimmomatic_1U.fastq.gz GAIIx_Lane8.trimmomatic_2P.fastq.gz GAIIx_Lane8.trimmomatic_2U.fastq.gz
MiSeq_A1: MiSeq_A1.trimmomatic_1P.fastq.gz MiSeq_A1.trimmomatic_1U.fastq.gz MiSeq_A1.trimmomatic_2P.fastq.gz MiSeq_A1.trimmomatic_2U.fastq.gz
MiSeq_A2: MiSeq_A2.trimmomatic_1P.fastq.gz MiSeq_A2.trimmomatic_1U.fastq.gz MiSeq_A2.trimmomatic_2P.fastq.gz MiSeq_A2.trimmomatic_2U.fastq.gz
MiSeq_B1: MiSeq_B1.trimmomatic_1P.fastq.gz MiSeq_B1.trimmomatic_1U.fastq.gz MiSeq_B1.trimmomatic_2P.fastq.gz MiSeq_B1.trimmomatic_2U.fastq.gz
MiSeq_B2: MiSeq_B2.trimmomatic_1P.fastq.gz MiSeq_B2.trimmomatic_1U.fastq.gz MiSeq_B2.trimmomatic_2P.fastq.gz MiSeq_B2.trimmomatic_2U.fastq.gz

GAIIx_%.trimmomatic_1P.fastq.gz GAIIx_%.trimmomatic_1U.fastq.gz GAIIx_%.trimmomatic_2P.fastq.gz GAIIx_%.trimmomatic_2U.fastq.gz: GAIIx_%_R1.fastq GAIIx_%_R2.fastq
	java -jar trimmomatic-0.32.jar PE -baseout GAIIx_$*.trimmomatic.fastq.gz $^ ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 TOPHRED33

MiSeq_%.trimmomatic_1P.fastq.gz MiSeq_%.trimmomatic_1U.fastq.gz MiSeq_%.trimmomatic_2P.fastq.gz MiSeq_%.trimmomatic_2U.fastq.gz: MiSeq_%_R1.fastq MiSeq_%_R2.fastq
	java -jar trimmomatic-0.32.jar PE -baseout MiSeq_$*.trimmomatic.fastq.gz $^ ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 TOPHRED33

###############
# Assembly
###############

.PHONY: ray_meta
ray_meta: Contigs_gt1kb.fasta

RayMeta_k31/Contigs.fasta: GAIIx_Lane7 GAIIx_Lane8 MiSeq_A1 MiSeq_A2 MiSeq_B1 MiSeq_B2
	mpiexec -n 40 Ray -k 31 -p GAIIx_Lane7.trimmomatic_1P.fastq.gz GAIIx_Lane7.trimmomatic_2P.fastq.gz -s GAIIx_Lane7.trimmomatic_1U.fastq.gz -s GAIIx_Lane7.trimmomatic_2U.fastq.gz -p GAIIx_Lane8.trimmomatic_1P.fastq.gz GAIIx_Lane8.trimmomatic_2P.fastq.gz -s GAIIx_Lane8.trimmomatic_1U.fastq.gz -s GAIIx_Lane8.trimmomatic_2U.fastq.gz -p MiSeq_A1.trimmomatic_1P.fastq.gz MiSeq_A1.trimmomatic_2P.fastq.gz -s MiSeq_A1.trimmomatic_1U.fastq.gz -s MiSeq_A1.trimmomatic_2U.fastq.gz -p MiSeq_A2.trimmomatic_1P.fastq.gz MiSeq_A2.trimmomatic_2P.fastq.gz -s MiSeq_A2.trimmomatic_1U.fastq.gz -s MiSeq_A2.trimmomatic_2U.fastq.gz -p MiSeq_B1.trimmomatic_1P.fastq.gz MiSeq_B1.trimmomatic_2P.fastq.gz -s MiSeq_B1.trimmomatic_1U.fastq.gz -s MiSeq_B1.trimmomatic_2U.fastq.gz -p MiSeq_B2.trimmomatic_1P.fastq.gz MiSeq_B2.trimmomatic_2P.fastq.gz -s MiSeq_B2.trimmomatic_1U.fastq.gz -s MiSeq_B2.trimmomatic_2U.fastq.gz -o RayMeta_k31 -minimum-contig-length 1000

Contigs_gt1kb.fasta: RayMeta_k31/Contigs.fasta
	cp $^ $@

###############
# Gene Prediction
###############

.PHONY: prodigal
prodigal: Contigs_gt1kb.prodigal.faa Contigs_gt1kb.prodigal.fna Contigs_gt1kb.prodigal.gff

%.prodigal.faa %.prodigal.fna %.prodigal.gff: Contigs_gt1kb.fasta
	prodigal -p meta -a Contigs_gt1kb.prodigal.faa -d Contigs_gt1kb.prodigal.fna -f gff -o Contigs_gt1kb.prodigal.gff -i $^

###############
# Read Mapping
###############

.PHONY: bowtie2 bowtie2-build
bowtie2: GAIIx_Lane6.trimmomatic.bam GAIIx_Lane7.trimmomatic.bam GAIIx_Lane8.trimmomatic.bam MiSeq_A1.trimmomatic.bam MiSeq_A2.trimmomatic.bam MiSeq_B1.trimmomatic.bam MiSeq_B2.trimmomatic.bam
bowtie2-build: Contigs_gt1kb.fasta.1.bt2 Contigs_gt1kb.fasta.2.bt2 Contigs_gt1kb.fasta.3.bt2 Contigs_gt1kb.fasta.4.bt2 Contigs_gt1kb.fasta.rev.1.bt2 Contigs_gt1kb.fasta.rev.2.bt2

%.trimmomatic.bam: %.trimmomatic_1P.fastq.gz %.trimmomatic_2P.fastq.gz
	bowtie2 -X 1000 -p 32 --end-to-end --sensitive -x Contigs_gt1kb.fasta -1 $*.trimmomatic_1P.fastq.gz -2 $*.trimmomatic_2P.fastq.gz | samtools view -uT Contigs_gt1kb.fasta - | samtools sort -m 32G - $*.trimmomatic

%.bam.bai: %.bam
	 samtools index $^

%.fasta.1.bt2 %.fasta.2.bt2 %.fasta.3.bt2 %.fasta.4.bt2 %.fasta.rev.1.bt2 %.fasta.rev.2.bt2: %.fasta
	bowtie2-build $^ $^

%.fasta.fai: %.fasta
	samtools faidx $^