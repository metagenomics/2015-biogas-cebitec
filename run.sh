#!/bin/bash

shopt -s extglob

#terminate after first line that fails
set -e 

THREADS=8

for i in "$@"
do
case $i in
    -t=*|--threads=*)
    THREADS="${i#*=}"
    ;;
     *)
            # unknown option
    ;;
esac
done

OUTPUT=/home/biogas/output

if [ ! -d "$OUTPUT" ]; then
	echo "$OUTPUT does not exist or is empty."
	exit 1
fi

make THREADS_RAY=$THREADS THREADS_MISC=$THREADS

mv !(run.sh|Makefile|trimmomatic-0.32.jar|Trimmomatic-0.32.zip|TruSeq2-PE.fa|NexteraPE-PE.fa|input|output) /home/biogas/output
