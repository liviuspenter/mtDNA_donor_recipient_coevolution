#!/bin/bash
#SBATCH -c 4                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-1:59                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

module load gcc/6.2.0
module load samtools/1.15.1

for exp in 1 5 10 50 100 500 1000
do
	echo $exp
	cat $exp/CLL_relapse1_1_barcodes.$exp $exp/CLL_relapse3_1_barcodes.$exp | cut -d" " -f1 > $exp/barcodes.tsv
	samtools merge $exp/merge.bam $exp/CLL_relapse1_1.$exp/mixing.bam $exp/CLL_relapse3_1.$exp/mixing.bam -@ 4
	samtools index $exp/merge.bam
done
