# Generation of in-silico mixing experiment data with CLL4 and CLL6

For historical reasons during the data analyses, the sample CLL4 is called CLL_relapse1_1 (mtDNA) or Pool91-1_22 (scRNA-seq),
while the sample from CLL6 is called CLL_relapse3_1 (mtDNA) or Pool91-1_24 (scRNA-seq) in the scripts. 

The mixing itself was performed using shell scripts on the HMS O2 cluster. 

## mtDNA-based deconvolution (mtscATAC-seq) 
To perform the mixing titration, first download the raw data (CLL4_1 and CLL6_1) from NCBI Geo (GSE163579) and process with cellranger-atac
as CLL_relapse1_1 and CLL_relapse3_1. Next, create an ArchR object ([```00_create_archr.R```](R/00_create_archr.R))

To generate the mixing steps, I first generated a list of cell barcodes from each sample for each titration step 
([```01_create_barcode_lists_atac.R```](R/01_create_barcode_lists_atac.R)).

Second, I used [sinto](https://github.com/timoast/sinto) to extract reads belonging to the cell barcodes 
([```01_processing.sh```](mtDNA/01_processing.sh)) into a new bam file for each sample and titration step.

Third, I merged the two bam files of each titration step using ([```02_merge_reads.sh```](mtDNA/02_merge_reads.sh)) 
and indexed the bam file.

Finally, I ran [mgatk](https://github.com/caleblareau/mgatk) on the merged bam files to perform discovery of 
mitochondrial DNA mutations ([```03_start_mgatk.sh```](mtDNA/03_start_mgatk.sh)).

## SNP-based deconvolution (scRNA-seq)
To perform the mixing titration, first download the raw data (CLL4_1 and CLL6_1) from NCBI Geo (GSE165087) and process with cellranger
as Pool91-1_22 and Pool91-1_24 for CLL4 and CLL6, respectively. Next, create a Seurat object 
([```10_create_seurat_object.R```](R/10_create_seurat_object.R)).

To generate the mixing steps, I first generated a list of cell barcodes from each sample for each titration step 
([```11_create_barcode_list_seurat.R```](R/11_create_barcode_list_seurat.R)).

Second, I used [sinto](https://github.com/timoast/sinto) to extract reads belonging to the cell barcodes 
([```01_processing.sh```](SNP/01_processing.sh)) into a new bam file for each sample and titration step.

Third, I merged the two bam files of each titration step using ([```02_merge_reads.sh```](mtDNA/02_merge_reads.sh)) 
and indexed the bam file.

Finally, I ran [souporcell](https://github.com/wheaton5/souporcell) and [vireo](https://github.com/single-cell-genetics/vireo)
for deconvolution with and without a germline reference. 

### souporcell without germline reference

### souporcell with germline reference
I was unable to make this work despite contacting 
