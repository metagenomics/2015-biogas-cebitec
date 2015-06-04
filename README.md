# Data availablility

### Metagenomic and metatranscriptomic sequencing

The datasets supporting the results of this article are available in the European Nucleotide Archive (ENA) under study accession [PRJEB8813](http://www.ebi.ac.uk/ena/data/view/PRJEB8813).

Please refer to Table 2 in [our manuscript](latex_src/bremges_gigascience_2015.pdf) for a detailed description of our metagenomic and metatranscriptomic sequencing efforts.

### Intermediate results for the review process

File | Description | Analysis step
--- | --- | ---
[Contigs.fasta.gz](raw_data/Contigs.fasta.gz) | Assembled contigs | RayMeta assembly
[Contigs.prodigal.gff.gz](raw_data/Contigs.prodigal.gff.gz) | Genes in GFF format | Gene prediction
[Contigs.prodigal.fna.gz](raw_data/Contigs.prodigal.fna.gz) | Nucleotide sequences | Gene prediction
[Contigs.prodigal.faa.gz](raw_data/Contigs.prodigal.faa.gz) | Protein translations | Gene prediction
[Contigs.prodigal.faa.blastp.kegg.tsv.gz](raw_data/Contigs.prodigal.faa.blastp.kegg.tsv.gz) | Tabular BLAST output | BLASTP vs. KEGG
[Contigs.prodigal.gff.bedtools.tsv.gz](raw_data/Contigs.prodigal.gff.bedtools.tsv.gz) | Read counts per gene | BEDTools multicov
[Contigs.prodigal.faa.blastp.kegg.annotated.tsv.gz](raw_data/Contigs.prodigal.faa.blastp.kegg.annotated.tsv.gz) | Annotated results | Custom: [annotate.pl](annotate.pl)

All data in [raw_data](raw_data) will be submitted to [GigaDB](http://gigadb.org/), the *GigaScience Database*, eventually.

# Reproducibility

Excluding the KEGG analysis, which relies on a commercial license of the KEGG database, all steps are performed using free and open-source software.

### Makefile

The complete workflow is organized in a single GNU [Makefile](Makefile). It downloads all data and re-runs all analysis steps (w/o the time-consuming BLASTP search against KEGG, for that please adjust the Makefile). All data and results can be reproduced by a simple invocation of `make`.

By default, the metagenome assembly (*Ray Meta*) will run with 48 threads. Read preprocessing (*Trimmomatic*) and mapping (*Bowtie2*) with 8 threads each. This suits e.g. a single-node machine with 48 cores and parallel execution of make with `make -j`. Please adjust the default values accordingly.

### Docker container

[@pbelmann](https://github.com/pbelmann) implemented and tested the [accompanying Docker container](https://registry.hub.docker.com/u/metagenomics/2015-biogas-cebitec).

1. `docker pull metagenomics/2015-biogas-cebitec`
2. `docker run  -v /path/to/wokspace/directory:/home/biogas/data metagenomics/2015-biogas-cebitec`

**Note:** The workspace directory, `/path/to/wokspace/directory`, mounted to the container should be on a volume with >83GB space.
After the container finished, all results can be found in here.

Per default the container runs with 8 threads (and a serial execution of make).
You can change this by specifying `--threads=NUMBER` after the name, e.g.
`docker run  -v /path/to/workspace/directory:/home/biogas/data metagenomics/2015-biogas-cebitec --threads=16`

### Docker on AWS

We tested the Docker container on an [r3.8xlarge](http://www.ec2instances.info/) instance with 32 Cores, 244GB RAM and a 320GB SSD volume. It should work on smaller instances, too: reproduction requires roughly **89GB memory** and **83GB storage**.

**Steps by step guide:**

1. Choose an instance with >83GB local volume size or mount an additional volume (>83GB) using the description provided by AWS: http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-add-volume-to-instance.html

2. Run `sudo apt-get update`

3. Install the newest docker version by using the description on docker.com:
https://docs.docker.com/installation/ubuntulinux/
(This image is tested with Docker version 1.6)

4. Start the container with
`sudo docker run  -v /path/to/workspace/directory:/home/biogas/data metagenomics/2015-biogas-cebitec`, where `/path/to/workspace/directory` is the path to a directory in your local storage volume or in a volume you mounted to your instance (see step 1).

--

**If you have any questions or run into problems, please file an [issue](https://github.com/abremges/2015-biogas-cebitec/issues)!**
