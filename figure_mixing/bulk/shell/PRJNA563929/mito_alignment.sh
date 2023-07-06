#!/bin/bash
#SBATCH -c 12                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-4:59                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=8GB                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

module load gcc/6.2.0
module load bowtie2/2.3.4.3

for library in SRR*
do
	echo $library

	if test -f $library/$library.chrM.bam; then
		echo "already processed"
		continue
	fi

	cd $library

	# alignment against hg38 chrM	
	bowtie2 -x /home/lp175/GRCH38/mito/chrM.hg38 --very-sensitive -p 12 -1 ${library}_1.fastq.gz -2 ${library}_2.fastq.gz -S tmp.sam

	# sort aligned reads
	samtools sort tmp.sam > tmp.sorted.sam

	# remove duplicates
	gatk MarkDuplicates -I tmp.sorted.sam -O ${library}.chrM.bam --REMOVE_DUPLICATES true --REMOVE_SEQUENCING_DUPLICATES true -M MarkDuplicates.metrics.txt

	# index file 
	samtools index ${library}.chrM.bam

	# cleanup
	rm tmp.sam
	rm tmp.sorted.sam

	cd ..
done
