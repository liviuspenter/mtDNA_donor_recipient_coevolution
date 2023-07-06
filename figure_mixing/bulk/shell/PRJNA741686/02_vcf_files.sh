#!/bin/bash

for library in SRR*
do
	echo $library
	bcftools sort $library/vcf/cellSNP.cells.vcf.gz -o $library/vcf/${library}.vcf.gz -Oz
	bcftools index $library/vcf/${library}.vcf.gz
done
