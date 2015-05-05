# 2015-biogas-cebitec

Work in progress, to be finished asap!

### Data for review

In the folder [raw_data](raw_data), you can find the following files:

File | Description | Analysis step
--- | --- | ---
[Contigs_gt1kb.fasta.gz](raw_data/Contigs_gt1kb.fasta.gz) | Assembled contigs | RayMeta assembly
[Contigs_gt1kb.prodigal.gff.gz](raw_data/Contigs_gt1kb.prodigal.gff.gz) | Genes in GFF format | Gene prediction
[Contigs_gt1kb.prodigal.fna.gz](raw_data/Contigs_gt1kb.prodigal.fna.gz) | Nucleotide sequences | Gene prediction
[Contigs_gt1kb.prodigal.faa.gz](raw_data/Contigs_gt1kb.prodigal.faa.gz) | Protein translations | Gene prediction
[Contigs_gt1kb.prodigal.faa.blastp.kegg.tsv.gz](raw_data/Contigs_gt1kb.prodigal.faa.blastp.kegg.tsv.gz) | Tabular BLAST output | BLASTP vs. KEGG
[Contigs_gt1kb.prodigal.faa.bedtools.tsv.gz](raw_data/Contigs_gt1kb.prodigal.faa.bedtools.tsv.gz) | Read counts per gene | BEDTools multicov
[Contigs_gt1kb.prodigal.faa.blastp.kegg.annotated.tsv.gz](raw_data/Contigs_gt1kb.prodigal.faa.blastp.kegg.annotated.tsv.gz) | Annotated results | Custom: [annotate.pl](annotate.pl)

These data will be submitted to GigaDB eventually.

### Availablility

Here, I will provide a description of the raw data, available from ENA.

g | h | i
--- | --- | ---
j | k | l

### Reproducibility

Please take a look at the [Makefile](Makefile). It downloads all data and re-runs all steps to reproduce my results (excluding the KEGG analyses, unfortunately). [@pbelmann](https://github.com/pbelmann) implemented and tested the accompanying Docker container.

#### How to run the docker container?

1. `docker pull metagenomics/2015-biogas-cebitec`
2. `docker run  -v /path/to/output/directory:/home/biogas/output 2015-biogas-cebitec`
     
Per default the container runs with 48 threads for Ray and 8 threads for Bowtie & Trimmomatic. You can change this by providing the following arguments after the docker name:
   
   * --threads-ray=NUMBER for Ray
   * --threads-misc=NUMBER for Bowtie and Trimmomatic

Example:

`docker run  -v /path/to/output/directory:/home/biogas/output 2015-biogas-cebitec --threads-ray=16 --threads-misc=16`

If you have any questions or run into problems, please file an [issue](https://github.com/abremges/2015-biogas-cebitec/issues)!

### The manuscript

In the folder [latex_src](latex_src), you can find the [manuscript's draft](latex_src/bremges_gigascience_2015.pdf) along with all LaTex source files and figures.
