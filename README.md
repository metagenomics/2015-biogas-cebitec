# Data availablility

Raw sequencing data are available in the European Nucleotide Archive (ENA) under study accession [PRJEB8813](http://www.ebi.ac.uk/ena/data/view/PRJEB8813).
The datasets supporting the results of [our manuscript](LaTeX/bremges_gigascience_2015.pdf) are available in GigaDB:

Bremges et al. **(2015)**. *GigaScience Database*. [doi:10.5524/100151](http://dx.doi.org/10.5524/100151)

# Reproducibility

Excluding the KEGG analysis, which relies on a commercial license of the KEGG database, all steps are performed using free and open-source software.

### Makefile

The complete workflow is organized in a single GNU [Makefile](Makefile). It downloads all data and re-runs all analysis steps (w/o the time-consuming BLASTP search against KEGG, for that please adjust the Makefile). All data and results can be reproduced by a simple invocation of `make`.

By default, the metagenome assembly (*Ray Meta*) will run with 48 threads. Read preprocessing (*Trimmomatic*) and mapping (*Bowtie2*) with 8 threads each. This suits e.g. a single-node machine with 48 cores and parallel execution of make with `make -j`. Please adjust the default values accordingly.

**Note:** [*Ray Meta* is nondeterministic](http://sourceforge.net/p/denovoassembler/mailman/message/31936187/) on 2 or more cores, and thus the assembly results will slightly vary from run to run. Downstream analyses will be affected by this, but the results are of comparable quality and mostly consistent.

### Docker container

[@pbelmann](https://github.com/pbelmann) implemented and tested the [accompanying Docker container](https://registry.hub.docker.com/u/metagenomics/2015-biogas-cebitec).

1. `docker pull metagenomics/2015-biogas-cebitec`
2. `docker run  -v /path/to/workspace/directory:/home/biogas/data metagenomics/2015-biogas-cebitec`

**Note:** The workspace directory, `/path/to/workspace/directory`, mounted to the container should be on a volume with >83GB space.
After the container finished, all results can be found in here.

Per default the container runs with 8 threads (and a serial execution of make).
You can change this by specifying `--threads=NUMBER` after the name, e.g.
`docker run  -v /path/to/workspace/directory:/home/biogas/data metagenomics/2015-biogas-cebitec --threads=32`

### Docker on AWS

We tested the Docker container on an [r3.8xlarge](http://www.ec2instances.info/) instance with 32 Cores, 244GB RAM and a 320GB SSD volume. On such an instance, setting `--threads=32`, execution takes **less than 24 hours** to complete (21 hours and 16 minutes in our latest test). It should work on smaller instances, too: reproduction requires roughly **89GB memory** and **83GB storage**.

**Steps by step guide:**

1. Choose an instance with >83GB local volume size or mount an additional volume (>83GB) using the description provided by AWS: http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-add-volume-to-instance.html

2. Run `sudo apt-get update`

3. Install the newest docker version by using the description on docker.com:
https://docs.docker.com/installation/ubuntulinux/
(This image is tested with Docker version 1.6)

4. Start the container with
`sudo docker run  -v /path/to/workspace/directory:/home/biogas/data metagenomics/2015-biogas-cebitec`, where `/path/to/workspace/directory` is the path to a directory in your local storage volume or in a volume you mounted to your instance (see step 1). *Set the number of threads to the number of available cores to fully utilize your instance.*

--

**If you have any questions or run into problems, please file an [issue](https://github.com/abremges/2015-biogas-cebitec/issues)!**
