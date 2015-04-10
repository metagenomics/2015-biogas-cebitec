FROM ubuntu:latest
MAINTAINER Peter Belmann

#add data
RUN mkdir /home/biogas
ADD Makefile /home/biogas/Makefile
ADD run.sh /home/biogas/run.sh
ADD raw_data /home/biogas/raw_data
RUN gunzip /home/biogas/raw_data/*.gz && mv /home/biogas/raw_data/* /home/biogas

#install ray
RUN apt-get update
RUN apt-get install -y openssh-server openmpi-bin
ADD bin/Ray /opt/bin/Ray

#install prodigal
RUN wget -O /opt/bin/prodigal https://github.com/hyattpd/Prodigal/releases/download/v2.6.2/prodigal.linux
RUN chmod a+x /opt/bin/prodigal

#install bowtie
RUN apt-get install bowtie2

#install samtools
RUN apt-get install samtools

#install trimmomatic
RUN apt-get install -y openjdk-7-jre
RUN apt-get install unzip
RUN wget -O /opt/bin/trimmomatic.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
RUN unzip /opt/bin/trimmomatic.zip -d /opt/bin
RUN mv /opt/bin/Trimmomatic-0.32/trimmomatic-0.32.jar  /home/biogas/
RUN cp /opt/bin/Trimmomatic-0.32/adapters/* /home/biogas

#install bedtools
RUN apt-get install -y bedtools 

#install make
RUN apt-get install make

#set path
ENV PATH /opt/bin:$PATH

WORKDIR /home/biogas
ENTRYPOINT ["/bin/bash","/home/biogas/run.sh"]
CMD ["8"]
