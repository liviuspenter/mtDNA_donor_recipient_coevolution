# Generation of in-silico mixing experiment data with CLL4 and CLL5

For historical reasons during the data analyses, the sample CLL4 is called CLL_relapse1_1 while the sample from CLL5 is called CLL_relapse3_1
in the scripts. 

## mtDNA-based deconvolution (mtscATAC-seq) 
To perform the mixing titration, first download the raw data (CLL4_1 and CLL5_1) from NCBI Geo (GSE163579) and process with cellranger-atac
as CLL_relapse1_1 and CLL_relapse3_1. Next, create an ArchR object (```00_create_ArchR.R```)

To generate the mixing steps, I first generated a list of cell barcodes from each sample for each titration step 
([```01_create_barcode_lists_atac.R```](R/01_create_barcode_lists_atac.R)).

Second, I used [sinto](https://github.com/timoast/sinto) to extract reads belonging to the cell barcodes 
([```01_processing.sh```](mtDNA/01_processing.sh)) into a new bam file for each sample and titration step.

Third, I merged the two bam files of each titration step using ([```02_merge_reads.sh```](mtDNA/02_merge_reads.sh)) 
and indexed the bam file.

Finally, I ran [mgatk](https://github.com/caleblareau/mgatk) on the merged bam files to perform discovery of 
mitochondrial DNA mutations ([```03_start_mgatk.sh```](mtDNA/03_start_mgatk.sh)).
