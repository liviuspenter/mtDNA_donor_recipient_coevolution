#!/bin/bash
#SBATCH -c 12                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-15:59                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem=64G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.mgatk.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.mgatk.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

cellsnp-lite -s Pool91-1_22/possorted_genome_bam.bam,Pool91-1_24/possorted_genome_bam.bam -I CLL1,CLL3 -O donor_genotype -R /home/lp175/snp_reference/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz -p 20 --cellTAG None --UMItag None --gzip --genotype 

