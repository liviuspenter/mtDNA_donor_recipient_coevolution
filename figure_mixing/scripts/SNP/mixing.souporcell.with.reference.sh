#!/bin/bash
#SBATCH -c 12                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=24G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.mgatk.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.mgatk.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

#GRCH38_REF='/n/shared_db/GRCh38/uk/cellranger/5.0.0/5.0.0/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
#GRCH38_REF=/n/shared_db/GRCh38/uk/cellranger/3.0.0/3.0.0/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa
#GRCH38_REF=/home/lp175/GRCH38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
GRCH38_REF=/home/lp175/GRCH38/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa
#SNP_REF=/home/lp175/snp_reference/common_variants_grch38_with_chr.vcf
SNP_REF=/home/lp175/snp_reference/common_variants_grch38.vcf
SOUPORCELL_CONTAINER=/n/app/singularity/containers/lp175/souporcell.sif  

singularity exec $SOUPORCELL_CONTAINER souporcell_pipeline.py -i $1/merge.bam -b $1/barcodes.tsv -f $GRCH38_REF -t 12 \
	-o $1/mixing_${1}.souporcell.reference -k 2 --skip_remap SKIP_REMAP \
	--known_genotypes donor_genotype/cellSNP.cells.vcf --known_genotypes_sample_names CLL1 CLL3
