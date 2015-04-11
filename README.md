# 2015-biogas-cebitec

Work in progress, to be finished asap!

## How to run the docker container?


1. `docker pull 2015-biogas-cebitec`
2. `docker run -v /path/to/input_directory:/home/biogas/input -v /path/to/output/directory:/home/biogas/output 2015-biogas-cebitec`

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
     
Per default the container runs with 8 threads. You can change this by providing a number after the docker name.

Example:

`docker run -v /path/to/input_directory:/home/biogas/input -v /path/to/output/directory:/home/biogas/output 2015-biogas-cebitec 32`
     
     
