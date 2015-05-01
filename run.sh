#!/bin/bash

shopt -s extglob

#terminate after first line that fails
set -e 

THREADS_RAY=48
THREADS_MISC=8

for i in "$@"
do
case $i in
    -tr=*|--threads-ray=*)
    THREADS_RAY="${i#*=}"
    ;;
    -tm=*|--threads-misc=*)
    THREADS_MISC="${i#*=}"
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

make THREADS_RAY=$THREADS_RAY THREADS_MISC=$THREADS_MISC

mv !(run.sh|Makefile|trimmomatic-0.32.jar|Trimmomatic-0.32.zip|TruSeq*|NexteraPE-PE.fa|input|output) /home/biogas/output
