#!/bin/bash

for exp in 1 5 10 50 100 500 1000
do
	sbatch ./extract_reads.sh Pool91-1_22/possorted_genome_bam.bam $exp/CLL_relapse1_1_barcodes.$exp $exp/CLL_relapse1_1.$exp
	sbatch ./extract_reads.sh Pool91-1_24/possorted_genome_bam.bam $exp/CLL_relapse3_1_barcodes.$exp $exp/CLL_relapse3_1.$exp
done

