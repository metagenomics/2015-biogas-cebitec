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

make THREADS_RAY=$THREADS THREADS_MISC=$THREADS WORKING_DIR="/home/biogas"
