# 2015-biogas-cebitec

Work in progress, to be finished asap!

### Data for review

In the folder [raw_data](raw_data), you can find the following files:

a | b | c
--- | --- | ---
d | e | f

These data will be submitted to GigaDB eventually.

### Availablility

Here, I will provide a description of the raw data, available from ENA.

g | h | i
--- | --- | ---
j | k | l

### Reproducibility

Please take a look at the [Makefile](Makefile). It downloads all data and re-runs all steps to reproduce my results (excluding the KEGG analyses, unfortunately). [@pbelmann](https://github.com/pbelmann) implemented and tested the accompanying Docker container. I will add more information soon.

If you have any questions or run into problems, please file an [issue](https://github.com/abremges/2015-biogas-cebitec/issues)!

### The manuscript

In the folder [latex_src](latex_src), you can find the [manuscript's draft](latex_src/bremges_gigascience_2015.pdf) along with all LaTex source files and figures.

## How to run the docker container?

1. `docker pull 2015-biogas-cebitec`
2. `docker run  -v /path/to/output/directory:/home/biogas/output 2015-biogas-cebitec`

  where input_directory MUST contain the following files:
      
      1. GAIIx_Lane6_R1.fastq
      2. GAIIx_Lane6_R2.fastq
      3. GAIIx_Lane7_R1.fastq
      4. GAIIx_Lane7_R2.fastq
      5. GAIIx_Lane8_R1.fastq
      6. GAIIx_Lane8_R2.fastq
      7. MiSeq_A1_R1.fastq
      8. MiSeq_A1_R2.fastq
      9. MiSeq_A2_R1.fastq
      10. MiSeq_A2_R2.fastq
      11. MiSeq_B1_R1.fastq
      12. MiSeq_B1_R2.fastq
      13. MiSeq_B2_R1.fastq
      14. MiSeq_B2_R2.fastq
     
Per default the container runs with 32 threads. You can change this by providing a number after the docker name.

Example:

`docker run  -v /path/to/output/directory:/home/biogas/output 2015-biogas-cebitec 32`
