#!/bin/bash

shopt -s extglob

#terminate after first line that fails
set -e 

OUTPUT=/home/biogas/output

if [ ! -d "$OUTPUT" ]; then
	echo "$OUTPUT does not exist or is empty."
	exit 1
fi

make THREADS_RAY=$1 THREADS_MISC=$2

mv !(run.sh|Makefile|trimmomatic-0.32.jar|Trimmomatic-0.32.zip|TruSeq*|NexteraPE-PE.fa|input|output) /home/biogas/output
