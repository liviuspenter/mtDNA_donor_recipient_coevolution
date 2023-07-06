#!/bin/bash
#SBATCH -c 12                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-11:59                         # Runtime in D-HH:MM format
#SBATCH -p short                          # Partition to run in
#SBATCH --mem=24G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.mgatk.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.mgatk.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

vireo -c $1/cellsnp/cellSNP.cells.vcf.gz -d donor_genotype/cellSNP.cells.vcf.gz -o $1/vireo/ -p 12 --randSeed 2

