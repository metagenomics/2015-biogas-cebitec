# Data availablility

### Metagenomic and metatranscriptomic sequencing

The datasets supporting the results of this article are available in the European Nucleotide Archive (ENA) under study accession [PRJEB8813](http://www.ebi.ac.uk/ena/data/view/PRJEB8813).

Please refer to Table 2 in [our manuscript](latex_src/bremges_gigascience_2015.pdf) for a detailed description of our metagenomic and metatranscriptomic sequencing efforts.

### Intermediate results for the review process

File | Description | Analysis step
--- | --- | ---
[Contigs_gt1kb.fasta.gz](raw_data/Contigs_gt1kb.fasta.gz) | Assembled contigs | RayMeta assembly
[Contigs_gt1kb.prodigal.gff.gz](raw_data/Contigs_gt1kb.prodigal.gff.gz) | Genes in GFF format | Gene prediction
[Contigs_gt1kb.prodigal.fna.gz](raw_data/Contigs_gt1kb.prodigal.fna.gz) | Nucleotide sequences | Gene prediction
[Contigs_gt1kb.prodigal.faa.gz](raw_data/Contigs_gt1kb.prodigal.faa.gz) | Protein translations | Gene prediction
[Contigs_gt1kb.prodigal.faa.blastp.kegg.tsv.gz](raw_data/Contigs_gt1kb.prodigal.faa.blastp.kegg.tsv.gz) | Tabular BLAST output | BLASTP vs. KEGG
[Contigs_gt1kb.prodigal.faa.bedtools.tsv.gz](raw_data/Contigs_gt1kb.prodigal.faa.bedtools.tsv.gz) | Read counts per gene | BEDTools multicov
[Contigs_gt1kb.prodigal.faa.blastp.kegg.annotated.tsv.gz](raw_data/Contigs_gt1kb.prodigal.faa.blastp.kegg.annotated.tsv.gz) | Annotated results | Custom: [annotate.pl](annotate.pl)

All data in [raw_data](raw_data) will be submitted to [GigaDB](http://gigadb.org/), the *GigaScience Database*, eventually.

# Reproducibility

Excluding the KEGG analysis, which relies on a commercial license of the KEGG database, all steps are performed using free and open-source software.

### Makefile

The complete workflow is organized in a single GNU [Makefile](Makefile). It downloads all data and re-runs all analysis steps (w/o the time-consuming BLASTP search against KEGG, for that please adjust the Makefile). All data and results can be reproduced by a simple invocation of `make`.

By default, the metagenome assembly (*Ray Meta*) will run with 48 threads. Read preprocessing (*Trimmomatic*) and mapping (*Bowtie2*) with 8 threads each. This suits e.g. a single-node machine with 48 cores and parallel execution of make with `make -j`. Please adjust the default values accordingly.

### Docker container

[@pbelmann](https://github.com/pbelmann) implemented and tested the accompanying Docker container.

1. `docker pull metagenomics/2015-biogas-cebitec`
2. `docker run  -v /path/to/output/directory:/home/biogas/output 2015-biogas-cebitec`
     
Per default the container runs with 48 threads for Ray and 8 threads for Bowtie & Trimmomatic. You can change this by providing the following arguments after the docker name:
   
   * --threads-ray=NUMBER for Ray
   * --threads-misc=NUMBER for Bowtie and Trimmomatic

Example:

`docker run  -v /path/to/output/directory:/home/biogas/output 2015-biogas-cebitec --threads-ray=16 --threads-misc=16`

--

**If you have any questions or run into problems, please file an [issue](https://github.com/abremges/2015-biogas-cebitec/issues)!**
