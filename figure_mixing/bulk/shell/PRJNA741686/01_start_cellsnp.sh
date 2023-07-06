#!/bin/bash

for library in SRR*
do
	sbatch ./cellsnp-lite.sh $library
done

