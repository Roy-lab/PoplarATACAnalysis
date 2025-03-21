# Poplar ATAC Analysis
Analysis of ATAC-Seq data from the stem of the Poplar (Populus trichocarpa) plant derived from the singleron-RD CeleScope platform using ArchR

> **Hardware requirements:**
1. Apple MacBook Pro or M3 laptop with 16 GB memory or equivalent decent configuration.
2. Unix OS with conda environment for ArchR. We used the latest Rocky CentOS and prebuilt conda environment for ArchR. For a new user, [ArchR website](https://www.archrproject.com/) is recommended for installing ArchR software.

> **Installing an ArchR conda environment:**

Currently, ArchR supports the R version [4.1 and 4.4](https://github.com/GreenleafLab/ArchR). Therefore, making a conda environment with ArchR with v4.4+ is recommended. However, some functionalities, like graph pdf making, require packages like Cairo, which need a lower version dependency, i.e., 4.3. Sometimes, downgrading the ArchR version may be necessary to get full functionalities. The current analysis installed the R version 4.4.2 but downgraded it to 4.3.1 for the latest functionalities. Some packages have the 4.4.2 dependencies, but this integrative dependency setting is recommended.
```
The following R session information is used for current analysis:
> sessionInfo()
R version 4.3.1 (2023-06-16)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Rocky Linux 8.10 (Green Obsidian)

Matrix products: default
BLAS/LAPACK: /mnt/dv/wid/projects6/Roy-singlecell3/bartholomew_lab/suvo_work/ArchR_analysis/scripts/suvo_ArchR/lib/libopenblasp-r0.3.29.so;  LAPACK version 3.12.0

Random number generation:
 RNG:     L'Ecuyer-CMRG 
 Normal:  Inversion 
 Sample:  Rejection 
 
locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Chicago
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] nabor_0.5.0                 lubridate_1.9.4            
 [3] forcats_1.0.0               dplyr_1.1.4                
 [5] purrr_1.0.4                 readr_2.1.5                
 [7] tidyr_1.3.1                 tibble_3.2.1               
 [9] tidyverse_2.0.0             AnnotationForge_1.48.0     
[11] rhdf5_2.50.2                SummarizedExperiment_1.36.0
[13] MatrixGenerics_1.18.1       Rcpp_1.0.14                
[15] Matrix_1.6-5                matrixStats_1.5.0          
[17] data.table_1.17.0           stringr_1.5.1              
[19] plyr_1.8.9                  magrittr_2.0.3             
[21] gtable_0.3.6                gtools_3.9.5               
[23] gridExtra_2.3               ArchR_1.0.2                
[25] ggplot2_3.4.4               Cairo_1.6-1                
[27] txdbmaker_1.2.1             GenomicFeatures_1.58.0     
[29] AnnotationDbi_1.68.0        Biobase_2.66.0             
[31] BSgenomeForge_1.6.0         BSgenome_1.74.0            
[33] rtracklayer_1.66.0          BiocIO_1.16.0              
[35] Biostrings_2.74.1           XVector_0.46.0             
[37] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
[39] IRanges_2.40.1              S4Vectors_0.44.0           
[41] BiocGenerics_0.52.0        

loaded via a namespace (and not attached):
 [1] DBI_1.2.3                bitops_1.0-9             httr2_1.1.1             
 [4] biomaRt_2.62.1           rlang_1.1.5              compiler_4.3.1          
 [7] RSQLite_2.3.9            png_0.1-8                vctrs_0.6.5             
[10] pkgconfig_2.0.3          crayon_1.5.3             fastmap_1.2.0           
[13] dbplyr_2.5.0             Rsamtools_2.22.0         tzdb_0.4.0              
[16] UCSC.utils_1.2.0         bit_4.6.0                zlibbioc_1.52.0         
[19] cachem_1.1.0             jsonlite_1.9.1           progress_1.2.3          
[22] blob_1.2.4               rhdf5filters_1.18.0      DelayedArray_0.32.0     
[25] Rhdf5lib_1.28.0          BiocParallel_1.40.0      parallel_4.3.1          
[28] prettyunits_1.2.0        R6_2.6.1                 stringi_1.7.12          
[31] timechange_0.3.0         tidyselect_1.2.1         abind_1.4-5             
[34] yaml_2.3.10              codetools_0.2-20         curl_6.2.1              
[37] lattice_0.22-6           withr_3.0.2              KEGGREST_1.46.0         
[40] BiocFileCache_2.14.0     xml2_1.3.5               pillar_1.10.1           
[43] filelock_1.0.3           generics_0.1.3           RCurl_1.98-1.16         
[46] hms_1.1.3                munsell_0.5.1            scales_1.3.0            
[49] glue_1.8.0               tools_4.3.1              GenomicAlignments_1.42.0
[52] XML_3.99-0.18            colorspace_2.1-1         GenomeInfoDbData_1.2.13 
[55] restfulr_0.0.15          cli_3.6.4                rappdirs_0.3.3          
[58] S4Arrays_1.6.0           digest_0.6.37            SparseArray_1.6.2       
[61] rjson_0.2.23             memoise_2.0.1            lifecycle_1.0.4         
[64] httr_1.4.7               bit64_4.5.2             
```

The making of the related conda environment can be followed in the following ways from UNIX terminal:
```
cd <your_conda_directory> e.g., <your_cconda_directory>=/mnt/dv/wid/projects6/Roy-singlecell3/bartholomew_lab/suvo_work/ArchR_analysis/scripts/
conda config --add channels conda-forge
conda config --set channel_priority strict
conda search r-base
conda create -n --prefix=<your_conda_environment_name> python=3.6 e.g. <your_conda_environment_name>=suvo_ArchR
# Output:
# To activate this environment, use
#     $ conda activate /mnt/dv/wid/projects6/Roy-singlecell3/bartholomew_lab/suvo_work/ArchR_analysis/scripts/suvo_ArchR
# To deactivate an active environment, use
#   $ conda deactivate
conda activate /mnt/dv/wid/projects6/Roy-singlecell3/bartholomew_lab/suvo_work/ArchR_analysis/scripts/suvo_ArchR
conda install -c conda-forge r-base=4.4.2 # The R base can be chosen to a lower version like 4.3.1 or 4.1.3, but it can downgraded later. The current analysis was later downgraded to 4.3.1 for Cairo functionality.
conda install -c conda-forge -c bioconda r-seurat=4 #Seurat tool install into conda environment
conda install -c conda-forge mamba #Mamba tool install into conda environment
conda install -c conda-forge -c bioconda bioconductor-chromvar #Chromvar tool install into conda environment
conda install -c conda-forge -c bioconda bioconductor-motifmatchr #Motifmachr tool install into conda environment
conda install -c conda-forge -c bioconda macs2 #Macs2 tool install into conda environment
conda install -c bioconda htslib #Tabix tool install into conda environment; this is required for proper index file generation against ATAC fragment.tsv.gz files.
conda install r-devtools #Devtools install into conda environment
##Next, open the R environment to install ArchR
R
 > if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
 > library(devtools)
 > devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
 > library(ArchR)
 > ArchR::installExtraPackages()
 > q()
## Reference webpage: http://biostars.org/p/498049/
```

> **ATAC files:**

Wendell sent us the three crucial ATAC files, which are as follows:

**1. ATAC-Seq fragment file:** fragments_corrected_dedup_count.tsv.gz
```
This is a fragment file with the following configuration:
> zcat fragments_corrected_dedup_count.tsv.gz | head -5
     Chr01	7504	7539	AACAGCGCTTGGCTTA	2
     Chr01	8106	8152	AACCGCCGCCAGTCTC	1
     Chr01	8107	8166	GGTAGCTGTAGGAGCG	1
     Chr01	8107	8173	CACTTGTTGTTCATGG	1
     Chr01	8114	8172	AATGTCACTCGGACCT	1
where,
Column 1: Sequence name (Chromosome or scaffold) information.
Column 2: Sequence (Chromosome or scaffold) start coordinate.
Column 3: Sequence (Chromosome or scaffold) end coordinate.
Column 4: ATAC fragment ID.
Column 5: Number of counts of the ATAC fragment.
```  
**2. Reference genome sequence file from which the fragments were referenced:** genome.fa
```
This is a fasta file of poplar genome which is used as reference in the sequencing:
> head -5 /mnt/dv/wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR/moved_files/genome.fa 
     >Chr01
     aaccctaaaccctaaaccctaaaccctaaaccctaaacctaaacccctaaaccctaaaccctaaaccctaaaccctaaac
     cctaaaccctaaaccctaaaccctaaaccctaaaccctaaaaccctaaaccctaaaaccctaaaccctaaaccctaaacc
     ctaaaccctaaaccaaaccctaaacctaacctaaaccctaaaccctaaaccctaaaccctagaaccctaaaccctaaccc
     taaccctaaaccctaaacctaacccctaaaccctaaaccctaaaccctaaacaaccctaaaccctaaaccctaaacccta
```  
**3. Reference General Feature Format (GTF) file from which the fragments were referenced:** gene.gtf
```
This is a genome information GTF file of poplar genome which is used as reference in the sequencing:
> head -5 gene.gtf 
     Chr01	phytozomev13	transcript	13697	17702	.	-	.	transcript_id "Potri.001G000400.5.v4.1"; gene_id "Potri.001G000400.v4.1"
     Chr01	phytozomev13	exon	13697	13947	.	-	.	transcript_id "Potri.001G000400.5.v4.1"; gene_id "Potri.001G000400.v4.1";
     Chr01	phytozomev13	exon	16186	16258	.	-	.	transcript_id "Potri.001G000400.5.v4.1"; gene_id "Potri.001G000400.v4.1";
     Chr01	phytozomev13	exon	16368	16463	.	-	.	transcript_id "Potri.001G000400.5.v4.1"; gene_id "Potri.001G000400.v4.1";
     Chr01	phytozomev13	exon	17511	17702	.	-	.	transcript_id "Potri.001G000400.5.v4.1"; gene_id "Potri.001G000400.v4.1";
where,
Column 1: Sequence (Chromosome or scaffold) information of where the feature is located.
Column 2: The source of annotation from Kirstlab. This indicates the annotation comes from phytozome v13 database.
Column 3: The type of feature. The value can be - "transcript", "exon" or "CDS" (coding sequence).
Column 4: The feature's start genomic coordinate (1-based index).
Column 5: The feature's end genomic coordinate (1-based index).
Column 6: The score of the feature. "." (dot) indicates no score is provided.
Column 7: The chromosomal strand information of the feature. The value can be - "+" or "-".
Column 8: The frame for the genomic feature. For "transcript" and "exon", the value is "." (dot) means the frame is not applicable because they are not CDS (coding sequence). CDS's  respective numerical frame score (0,1 and 2) is provided. 0 means the first base of the CDS aligns with the start codon., 1 means the first base of the CDS is shifted by one nucleotide (it would align differently), and 2 means two nucleotides shift the first base.
Column 9: The attributes of features that contain two key-value pairs: (i) transcript ID, a unique identifier for the transcript, and (ii) gene ID, a unique identifier for the gene associated with the transcript.
```
 




