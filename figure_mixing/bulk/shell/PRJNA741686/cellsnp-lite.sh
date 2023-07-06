#!/bin/bash
#SBATCH -c 12                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-00:20                         # Runtime in D-HH:MM format
#SBATCH -p short                          # Partition to run in
#SBATCH --mem=12G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%jcellsnp.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.cellsnp.err                 # File to which STDERR will be written, including job ID

cellsnp-lite -s $1/${1}.bam -I $1 -O $1/vcf/ -R /home/lp175/snp_reference/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz -p 20 --cellTAG None --UMItag None --gzip --genotype 
