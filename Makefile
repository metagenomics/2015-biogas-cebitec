AUTHOR  = Andreas Bremges
EMAIL   = abremges@cebitec.uni-bielefeld.de
VERSION = 1.1.0

# By default, this Makefile downloads all data and reproduces everything but the KEGG analysis
all: check_bin download_data trimmomatic ray_meta prodigal bowtie2-build bowtie2-run samtools-index bedtools-multicov

# Working directory
WORKING_DIR=data

# Number of threads, adjust accordingly. We used 1 machine with 48 cores.
# Ray is a bottleneck in the mackfile, and thus should use all cores
THREADS_RAY=48
# Using 'make -j', there will be 7 simultaneous Trimmomatic/Bowtie2 jobs
THREADS_MISC=8
# For non-parallel 'make', set both variables to e.g. the number of cores

###############
# Check dependencies
###############

.PHONY .SILENT: check_bin
check_bin:
	command -v wget >/dev/null 2>&1 || { echo >&2 "wget not found. Aborting."; exit 1; }
	command -v unzip >/dev/null 2>&1 || { echo >&2 "unzip not found. Aborting."; exit 1; }
	command -v java >/dev/null 2>&1 || { echo >&2 "java not found. Aborting."; exit 1; }
	command -v mpiexec >/dev/null 2>&1 || { echo >&2 "mpiexec not found. Aborting."; exit 1; }
	command -v Ray >/dev/null 2>&1 || { echo >&2 "Ray not found. Aborting."; exit 1; }
	command -v prodigal >/dev/null 2>&1 || { echo >&2 "prodigal not found. Aborting."; exit 1; }
	command -v bowtie2-build >/dev/null 2>&1 || { echo >&2 "bowtie2-build not found. Aborting."; exit 1; }
	command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "bowtie2 not found. Aborting."; exit 1; }
	command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools not found. Aborting."; exit 1; }
	command -v bedtools >/dev/null 2>&1 || { echo >&2 "bedtools not found. Aborting."; exit 1; }

###############
# Setup working directory
###############

.PHONY: check_dir
$(WORKING_DIR):
	mkdir -p $@

###############
# Download data
###############

.PHONY: download_data
download_data: $(WORKING_DIR)/GAIIx_Lane6_R1.fastq.gz $(WORKING_DIR)/GAIIx_Lane6_R2.fastq.gz $(WORKING_DIR)/GAIIx_Lane7_R1.fastq.gz $(WORKING_DIR)/GAIIx_Lane7_R2.fastq.gz $(WORKING_DIR)/GAIIx_Lane8_R1.fastq.gz $(WORKING_DIR)/GAIIx_Lane8_R2.fastq.gz $(WORKING_DIR)/MiSeq_A1_R1.fastq.gz $(WORKING_DIR)/MiSeq_A1_R2.fastq.gz $(WORKING_DIR)/MiSeq_A2_R1.fastq.gz $(WORKING_DIR)/MiSeq_A2_R2.fastq.gz $(WORKING_DIR)/MiSeq_B1_R1.fastq.gz $(WORKING_DIR)/MiSeq_B1_R2.fastq.gz $(WORKING_DIR)/MiSeq_B2_R1.fastq.gz $(WORKING_DIR)/MiSeq_B2_R2.fastq.gz

$(WORKING_DIR)/GAIIx_Lane6_R1.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_6_1_sequence.txt.gz

$(WORKING_DIR)/GAIIx_Lane6_R2.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_6_2_sequence.txt.gz

$(WORKING_DIR)/GAIIx_Lane7_R1.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_7_1_sequence.txt.gz

$(WORKING_DIR)/GAIIx_Lane7_R2.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_7_2_sequence.txt.gz

$(WORKING_DIR)/GAIIx_Lane8_R1.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_8_1_sequence.txt.gz

$(WORKING_DIR)/GAIIx_Lane8_R2.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/s_8_2_sequence.txt.gz

$(WORKING_DIR)/MiSeq_A1_R1.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq_Biogas1_S1_L001_R1_001.fastq.gz

$(WORKING_DIR)/MiSeq_A1_R2.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq_Biogas1_S1_L001_R2_001.fastq.gz

$(WORKING_DIR)/MiSeq_A2_R1.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq_Biogas2_S2_L001_R1_001.fastq.gz

$(WORKING_DIR)/MiSeq_A2_R2.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq_Biogas2_S2_L001_R2_001.fastq.gz

$(WORKING_DIR)/MiSeq_B1_R1.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq2_Biogas1_S1_L001_R1_001.fastq.gz

$(WORKING_DIR)/MiSeq_B1_R2.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq2_Biogas1_S1_L001_R2_001.fastq.gz

$(WORKING_DIR)/MiSeq_B2_R1.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq2_Biogas2_S2_L001_R1_001.fastq.gz

$(WORKING_DIR)/MiSeq_B2_R2.fastq.gz: $(WORKING_DIR)
	wget -O $@ ftp://ftp.sra.ebi.ac.uk/vol1/ERA427/ERA427694/fastq/MiSeq2_Biogas2_S2_L001_R2_001.fastq.gz

###############
# Read QC (multitreaded: THREADS_MISC)
###############

.PHONY: trimmomatic GAIIx_Lane6 GAIIx_Lane7 GAIIx_Lane8 MiSeq_A1 MiSeq_A2 MiSeq_B1 MiSeq_B2
trimmomatic: trimmomatic-0.32.jar TruSeq2-PE.fa NexteraPE-PE.fa GAIIx_Lane6 GAIIx_Lane7 GAIIx_Lane8 MiSeq_A1 MiSeq_A2 MiSeq_B1 MiSeq_B2

Trimmomatic-0.32.zip:
	wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip

trimmomatic-0.32.jar: Trimmomatic-0.32.zip
	unzip -p $^ Trimmomatic-0.32/trimmomatic-0.32.jar > $@

TruSeq2-PE.fa: Trimmomatic-0.32.zip
	unzip -p $^ Trimmomatic-0.32/adapters/TruSeq2-PE.fa > $@

NexteraPE-PE.fa: Trimmomatic-0.32.zip
	unzip -p $^ Trimmomatic-0.32/adapters/NexteraPE-PE.fa > $@

GAIIx_Lane6: $(WORKING_DIR)/GAIIx_Lane6.trimmomatic_1P.fastq.gz $(WORKING_DIR)/GAIIx_Lane6.trimmomatic_1U.fastq.gz $(WORKING_DIR)/GAIIx_Lane6.trimmomatic_2P.fastq.gz $(WORKING_DIR)/GAIIx_Lane6.trimmomatic_2U.fastq.gz
GAIIx_Lane7: $(WORKING_DIR)/GAIIx_Lane7.trimmomatic_1P.fastq.gz $(WORKING_DIR)/GAIIx_Lane7.trimmomatic_1U.fastq.gz $(WORKING_DIR)/GAIIx_Lane7.trimmomatic_2P.fastq.gz $(WORKING_DIR)/GAIIx_Lane7.trimmomatic_2U.fastq.gz
GAIIx_Lane8: $(WORKING_DIR)/GAIIx_Lane8.trimmomatic_1P.fastq.gz $(WORKING_DIR)/GAIIx_Lane8.trimmomatic_1U.fastq.gz $(WORKING_DIR)/GAIIx_Lane8.trimmomatic_2P.fastq.gz $(WORKING_DIR)/GAIIx_Lane8.trimmomatic_2U.fastq.gz
MiSeq_A1: $(WORKING_DIR)/MiSeq_A1.trimmomatic_1P.fastq.gz $(WORKING_DIR)/MiSeq_A1.trimmomatic_1U.fastq.gz $(WORKING_DIR)/MiSeq_A1.trimmomatic_2P.fastq.gz $(WORKING_DIR)/MiSeq_A1.trimmomatic_2U.fastq.gz
MiSeq_A2: $(WORKING_DIR)/MiSeq_A2.trimmomatic_1P.fastq.gz $(WORKING_DIR)/MiSeq_A2.trimmomatic_1U.fastq.gz $(WORKING_DIR)/MiSeq_A2.trimmomatic_2P.fastq.gz $(WORKING_DIR)/MiSeq_A2.trimmomatic_2U.fastq.gz
MiSeq_B1: $(WORKING_DIR)/MiSeq_B1.trimmomatic_1P.fastq.gz $(WORKING_DIR)/MiSeq_B1.trimmomatic_1U.fastq.gz $(WORKING_DIR)/MiSeq_B1.trimmomatic_2P.fastq.gz $(WORKING_DIR)/MiSeq_B1.trimmomatic_2U.fastq.gz
MiSeq_B2: $(WORKING_DIR)/MiSeq_B2.trimmomatic_1P.fastq.gz $(WORKING_DIR)/MiSeq_B2.trimmomatic_1U.fastq.gz $(WORKING_DIR)/MiSeq_B2.trimmomatic_2P.fastq.gz $(WORKING_DIR)/MiSeq_B2.trimmomatic_2U.fastq.gz

$(WORKING_DIR)/GAIIx_%.trimmomatic_1P.fastq.gz $(WORKING_DIR)/GAIIx_%.trimmomatic_1U.fastq.gz $(WORKING_DIR)/GAIIx_%.trimmomatic_2P.fastq.gz $(WORKING_DIR)/GAIIx_%.trimmomatic_2U.fastq.gz: $(WORKING_DIR)/GAIIx_%_R1.fastq.gz $(WORKING_DIR)/GAIIx_%_R2.fastq.gz
	java -jar trimmomatic-0.32.jar PE -threads $(THREADS_MISC) -baseout $(WORKING_DIR)/GAIIx_$*.trimmomatic.fastq.gz $^ ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 TOPHRED33

$(WORKING_DIR)/MiSeq_%.trimmomatic_1P.fastq.gz $(WORKING_DIR)/MiSeq_%.trimmomatic_1U.fastq.gz $(WORKING_DIR)/MiSeq_%.trimmomatic_2P.fastq.gz $(WORKING_DIR)/MiSeq_%.trimmomatic_2U.fastq.gz: $(WORKING_DIR)/MiSeq_%_R1.fastq.gz $(WORKING_DIR)/MiSeq_%_R2.fastq.gz
	java -jar trimmomatic-0.32.jar PE -threads $(THREADS_MISC) -baseout $(WORKING_DIR)/MiSeq_$*.trimmomatic.fastq.gz $^ ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 TOPHRED33

###############
# Assembly (multithreaded: THREADS_RAY)
###############

.PHONY: ray_meta
ray_meta: $(WORKING_DIR)/Contigs_gt1kb.fasta

$(WORKING_DIR)/RayMeta_k31/Contigs.fasta: GAIIx_Lane7 GAIIx_Lane8 MiSeq_A1 MiSeq_A2 MiSeq_B1 MiSeq_B2
	mpiexec -n $(THREADS_RAY) Ray -k 31 -p $(WORKING_DIR)/GAIIx_Lane7.trimmomatic_1P.fastq.gz $(WORKING_DIR)/GAIIx_Lane7.trimmomatic_2P.fastq.gz -s $(WORKING_DIR)/GAIIx_Lane7.trimmomatic_1U.fastq.gz -s $(WORKING_DIR)/GAIIx_Lane7.trimmomatic_2U.fastq.gz -p $(WORKING_DIR)/GAIIx_Lane8.trimmomatic_1P.fastq.gz $(WORKING_DIR)/GAIIx_Lane8.trimmomatic_2P.fastq.gz -s $(WORKING_DIR)/GAIIx_Lane8.trimmomatic_1U.fastq.gz -s $(WORKING_DIR)/GAIIx_Lane8.trimmomatic_2U.fastq.gz -p $(WORKING_DIR)/MiSeq_A1.trimmomatic_1P.fastq.gz $(WORKING_DIR)/MiSeq_A1.trimmomatic_2P.fastq.gz -s $(WORKING_DIR)/MiSeq_A1.trimmomatic_1U.fastq.gz -s $(WORKING_DIR)/MiSeq_A1.trimmomatic_2U.fastq.gz -p $(WORKING_DIR)/MiSeq_A2.trimmomatic_1P.fastq.gz $(WORKING_DIR)/MiSeq_A2.trimmomatic_2P.fastq.gz -s $(WORKING_DIR)/MiSeq_A2.trimmomatic_1U.fastq.gz -s $(WORKING_DIR)/MiSeq_A2.trimmomatic_2U.fastq.gz -p $(WORKING_DIR)/MiSeq_B1.trimmomatic_1P.fastq.gz $(WORKING_DIR)/MiSeq_B1.trimmomatic_2P.fastq.gz -s $(WORKING_DIR)/MiSeq_B1.trimmomatic_1U.fastq.gz -s $(WORKING_DIR)/MiSeq_B1.trimmomatic_2U.fastq.gz -p $(WORKING_DIR)/MiSeq_B2.trimmomatic_1P.fastq.gz $(WORKING_DIR)/MiSeq_B2.trimmomatic_2P.fastq.gz -s $(WORKING_DIR)/MiSeq_B2.trimmomatic_1U.fastq.gz -s $(WORKING_DIR)/MiSeq_B2.trimmomatic_2U.fastq.gz -o $(WORKING_DIR)/RayMeta_k31 -minimum-contig-length 1000

$(WORKING_DIR)/Contigs_gt1kb.fasta: $(WORKING_DIR)/RayMeta_k31/Contigs.fasta
	cp $^ $@

###############
# Gene prediction
###############

.PHONY: prodigal
prodigal: $(WORKING_DIR)/Contigs_gt1kb.prodigal.faa $(WORKING_DIR)/Contigs_gt1kb.prodigal.fna $(WORKING_DIR)/Contigs_gt1kb.prodigal.gff

$(WORKING_DIR)/%.prodigal.faa $(WORKING_DIR)/%.prodigal.fna $(WORKING_DIR)/%.prodigal.gff: $(WORKING_DIR)/Contigs_gt1kb.fasta
	prodigal -p meta -a $(WORKING_DIR)/Contigs_gt1kb.prodigal.faa -d $(WORKING_DIR)/Contigs_gt1kb.prodigal.fna -f gff -o $(WORKING_DIR)/Contigs_gt1kb.prodigal.gff -i $^

###############
# BLASTP against KEGG genes (only described)
###############

.PHONY: kegg-blastp
.SILENT kegg-blastp: $(WORKING_DIR)/Contigs_gt1kb.prodigal.faa.blastp.kegg.tsv

$(WORKING_DIR)/Contigs_gt1kb.prodigal.faa.blastp.kegg.tsv:
	echo 'We distributed BLASTP-jobs on our compute cluster. A simplified command, reproducing the results but running several weeks, looks like this:'
	echo 'blastp -max_target_seqs 1 -outfmt 6 -query Contigs_gt1kb.prodigal.faa -out Contigs_gt1kb.prodigal.faa.blastp.kegg.tsv -db /path/to/kegg/blastdb'

###############
# Read mapping (multithreaded: THREADS_MISC)
###############

.PHONY: bowtie2-build bowtie2-run samtools-index
bowtie2-build: $(WORKING_DIR)/Contigs_gt1kb.fasta.1.bt2 $(WORKING_DIR)/Contigs_gt1kb.fasta.2.bt2 $(WORKING_DIR)/Contigs_gt1kb.fasta.3.bt2 $(WORKING_DIR)/Contigs_gt1kb.fasta.4.bt2 $(WORKING_DIR)/Contigs_gt1kb.fasta.rev.1.bt2 $(WORKING_DIR)/Contigs_gt1kb.fasta.rev.2.bt2
bowtie2-run: $(WORKING_DIR)/GAIIx_Lane6.trimmomatic.bam $(WORKING_DIR)/GAIIx_Lane7.trimmomatic.bam $(WORKING_DIR)/GAIIx_Lane8.trimmomatic.bam $(WORKING_DIR)/MiSeq_A1.trimmomatic.bam $(WORKING_DIR)/MiSeq_A2.trimmomatic.bam $(WORKING_DIR)/MiSeq_B1.trimmomatic.bam $(WORKING_DIR)/MiSeq_B2.trimmomatic.bam
samtools-index: $(WORKING_DIR)/GAIIx_Lane6.trimmomatic.bam.bai $(WORKING_DIR)/GAIIx_Lane7.trimmomatic.bam.bai $(WORKING_DIR)/GAIIx_Lane8.trimmomatic.bam.bai $(WORKING_DIR)/MiSeq_A1.trimmomatic.bam.bai $(WORKING_DIR)/MiSeq_A2.trimmomatic.bam.bai $(WORKING_DIR)/MiSeq_B1.trimmomatic.bam.bai $(WORKING_DIR)/MiSeq_B2.trimmomatic.bam.bai

$(WORKING_DIR)/%.fasta.1.bt2 $(WORKING_DIR)/%.fasta.2.bt2 $(WORKING_DIR)/%.fasta.3.bt2 $(WORKING_DIR)/%.fasta.4.bt2 $(WORKING_DIR)/%.fasta.rev.1.bt2 $(WORKING_DIR)/%.fasta.rev.2.bt2: $(WORKING_DIR)/%.fasta
	bowtie2-build $^ $^

$(WORKING_DIR)/%.trimmomatic.bam: $(WORKING_DIR)/%.trimmomatic_1P.fastq.gz $(WORKING_DIR)/%.trimmomatic_2P.fastq.gz
	bowtie2 -X 1000 -p $(THREADS_MISC) --end-to-end --sensitive -x $(WORKING_DIR)/Contigs_gt1kb.fasta -1 $(WORKING_DIR)/$*.trimmomatic_1P.fastq.gz -2 $(WORKING_DIR)/$*.trimmomatic_2P.fastq.gz | samtools view -uT $(WORKING_DIR)/Contigs_gt1kb.fasta - | samtools sort - $(WORKING_DIR)/$*.trimmomatic

$(WORKING_DIR)/%.bam.bai: $(WORKING_DIR)/%.bam
	samtools index $^

###############
# Count reads in genes
###############

.PHONY: bedtools-multicov
bedtools-multicov: $(WORKING_DIR)/Contigs_gt1kb.prodigal.gff $(WORKING_DIR)/Contigs_gt1kb.prodigal.gff.bedtools.tsv samtools-index

$(WORKING_DIR)/Contigs_gt1kb.prodigal.gff.bedtools.tsv: $(WORKING_DIR)/GAIIx_Lane6.trimmomatic.bam $(WORKING_DIR)/GAIIx_Lane7.trimmomatic.bam $(WORKING_DIR)/GAIIx_Lane8.trimmomatic.bam $(WORKING_DIR)/MiSeq_A1.trimmomatic.bam $(WORKING_DIR)/MiSeq_A2.trimmomatic.bam $(WORKING_DIR)/MiSeq_B1.trimmomatic.bam $(WORKING_DIR)/MiSeq_B2.trimmomatic.bam
	bedtools multicov -bams $^ -bed $(WORKING_DIR)/Contigs_gt1kb.prodigal.gff > $@

###############
# Annotate results: KOs, modules, pathways
###############

.PHONY: kegg-annotate
kegg-annotate: $(WORKING_DIR)/Contigs_gt1kb.prodigal.faa.blastp.kegg.annotated.tsv

$(WORKING_DIR)/Contigs_gt1kb.prodigal.faa.blastp.kegg.annotated.tsv: annotate.pl genes_ko.list ko_module.list ko_pathway.list Contigs_gt1kb.prodigal.faa.blastp.kegg.tsv Contigs_gt1kb.prodigal.gff.bedtools.tsv
	perl annotate.pl > $@
