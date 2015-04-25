#!/bin/bash

shopt -s extglob

#terminate after first line that fails
set -e 

OUTPUT=/home/biogas/output

if [ ! -d "$OUTPUT" ]; then
	echo "$OUTPUT does not exist or is empty."
	exit 1
fi

make THREADS=$1

mv !(run.sh|Makefile|trimmomatic-0.32.jar|input|output) /home/biogas/output
