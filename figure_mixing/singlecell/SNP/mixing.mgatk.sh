#!/bin/bash
#SBATCH -c 12                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-1:59                         # Runtime in D-HH:MM format
#SBATCH -p short                          # Partition to run in
#SBATCH --mem=64G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.mgatk.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.mgatk.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

module load gcc/6.2.0
module load python/3.6.0
module load java/jdk-1.8u112
module load R/3.6.1

ulimit -n 2000

# IMPORTANT:
# run before script
# python3 -m venv python3
# source python3/bin/activate
mgatk tenx -i $1/merge.bam -n mixing_$1 -o $1/mixing_${1}.mgatk -c 12 -bt CB -b $1/barcodes.tsv --mito-genome hg38 --nsamples 1900
