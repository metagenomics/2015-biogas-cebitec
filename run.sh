#!/bin/bash

shopt -s extglob

#terminate after first line that fails
set -e 

INPUT=/home/biogas/input
OUTPUT=/home/biogas/output

if [ ! -d "$INPUT" ] || [ ! "$(ls -A $INPUT)" ]; then
	echo " $INPUT does not exist or is empty."
	exit 1
fi

if [ ! -d "$OUTPUT" ]; then
	echo "$OUTPUT does not exist or is empty."
	exit 1
fi

make THREADS=$1

mv !(run.sh|Makefile|trimmomatic-0.32.jar|input|output) /home/biogas/output
