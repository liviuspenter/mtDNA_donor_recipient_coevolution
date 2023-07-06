#!/bin/bash
for library in SRR*
do
	count=`gunzip -c $library/vcf/${library}.vcf.gz | grep GT | wc -l`
	echo "$library $count" >> SNP_count.csv
done
