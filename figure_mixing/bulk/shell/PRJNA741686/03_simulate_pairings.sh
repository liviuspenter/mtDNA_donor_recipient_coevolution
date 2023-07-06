#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-11:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=3GB                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

for library in SRR*
do
	for library2 in SRR*
	do
		if [[ $library == $library2 ]]; then
			continue
		fi

		echo "processing $library $library2"

		bcftools isec $library/vcf/${library}.vcf.gz $library2/vcf/${library2}.vcf.gz -p test

		count_library=`cat test/0000.vcf | grep GT | wc -l`
		count_library2=`cat test/0001.vcf | grep GT | wc -l`

		echo "$library $library2 $count_library $count_library2" >> simulated_snps.csv
	done
done
