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

Wendell sent us the three crucial ATAC files (see, attached file: [Kirstlab_Celescope_README.md](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/Kirstlab_Celescope_README.md)), which are as follows:

**1. ATAC-Seq fragment file:** fragments_corrected_dedup_count.tsv.gz [First 5000 lines](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/ExampleData/fragments_corrected_dedup_count.tsv.gz) #[fragments_corrected_dedup_count.tsv.gz](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/Kirstlab_Celescope_README.md)
```
This is a fragment file with the following configuration:
> zcat ExampleData/fragments_corrected_dedup_count.tsv.gz | head -5
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
**2. Reference genome sequence file from which the fragments were referenced:** genome.fa [First 5000 lines](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/ExampleData/genome.fa) #[Ptrichocarpa_533_v4.0.softmasked.fa](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/Kirstlab_Celescope_README.md)
```
This is a fasta file of poplar genome which is used as reference in the sequencing:
> head -5 ExampleData/genome.fa 
     >Chr01
     aaccctaaaccctaaaccctaaaccctaaaccctaaacctaaacccctaaaccctaaaccctaaaccctaaaccctaaac
     cctaaaccctaaaccctaaaccctaaaccctaaaccctaaaaccctaaaccctaaaaccctaaaccctaaaccctaaacc
     ctaaaccctaaaccaaaccctaaacctaacctaaaccctaaaccctaaaccctaaaccctagaaccctaaaccctaaccc
     taaccctaaaccctaaacctaacccctaaaccctaaaccctaaaccctaaacaaccctaaaccctaaaccctaaacccta
```  
**3. Reference General Feature Format (GTF) file from which the fragments were referenced:** gene.gtf [First 5000 lines](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/ExampleData/gene.gtf) #[Ptrichocarpa_533_v4.1.gtf](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/Kirstlab_Celescope_README.md)
```
This is a genome information GTF file of poplar genome which is used as reference in the sequencing:
> head -5 ExampleData/gene.gtf 
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

> **Critical assessment of the input ATAC files:**


**NOTE-1:**

There are some  discrepancies in the sequence names from three ATAC files of the Poplar ATAC data, which are as follows,<br/>
![image](https://github.com/user-attachments/assets/6df31770-defc-4c43-95c9-543c1a0f22c0)

As the scaffold information is mismatched in these three files, only chromosome 1 to 19 ATAC data is processed through this pipeline.

**NOTE-2:**

The sequence names for chromosomes in the poplar genome is represented as:\
Chr01, Chr02, Chr03, Chr04, Chr05, Chr06, Chr07, Chr08, Chr09, Chr10, Chr11, Chr12, Chr13, Chr14, Chr15, Chr16, Chr17, Chr18, Chr19.\
As [ArchR usage](https://www.archrproject.com/bookdown/getting-set-up.html) generally requires the prefix of annotation as "chr", the present study changes the sequence name in all three files to:<br/>
chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19

**NOTE-3:**

The fasta files describes that the length of chromosome is as following:
```
awk '/^>/ {if (seq) print chrom, length(seq); seq=""; chrom=$0; next} {seq = seq $0} END {if (seq) print chrom, length(seq)}' genome.fa
Chr01	49788581
Chr02	25242375
Chr03	21678634
Chr04	24140038
Chr05	24981103
Chr06	27516652
Chr07	15561420
Chr08	19195260
Chr09	12987399
Chr10	22799081
Chr11	19288771
Chr12	15589050
Chr13	15705617
Chr14	17801709
Chr15	15231745
Chr16	14619816
Chr17	15189755
Chr18	16264003
Chr19	15623655
scaffold_104	239390
scaffold_1549	75560
scaffold_1697	72059
scaffold_1865	68284
scaffold_1944	66818
scaffold_2151	63019
scaffold_2281	60697
scaffold_25	640640
scaffold_2555	56646
scaffold_2798	52914
scaffold_3255	45437
scaffold_3421	40476
scaffold_3435	40109
scaffold_3488	36978
scaffold_3518	35470
scaffold_3588	28433
scaffold_3607	24495
scaffold_3611	23661
scaffold_3618	21309
scaffold_3621	20187
scaffold_45	362030
```
However, some fragments have end sequence coordinate numbers that are higher than the chromosome length. like as the following example in Chromosome11, it is seen that,
```
zcat fragments_corrected_dedup_count.tsv.gz | grep Chr11 | awk '{diff = $3 - $2; print $0, diff}' | sort -k3,3nr -k6,6nr | grep 19288781

Chr11  19288645  19288781  GGCCTACTCCAGACTA  1 136
Chr11  19288749  19288781  TGTGCACCAGCCGTGA  2 32
```
The marked end coordinate of the fragments like "GGCCTACTCCAGACTA" or "TGTGCACCAGCCGTGA" (19288781) is higher than the length of chromosome 11 length, i.e., 19288771. Therefore, these creepy fragments with a dull end chromosome number must be cleaned from the fragment file.

<br/>
<br/>


> **Scripts usage for getting ATAC genescore matrix though ArchR protocol:**
<br/>

## Step-1: Updating the fragment file discrepancies:

There are three sub-steps to remove fragment file discrepancies,
**1A:** Remove the fragments with end sequence coordinates higher than the chromosomal sequence length.
**1B:** Rename the sequence name to chr1-chr19 from Chr01-Chr19.
**1C:** Making the updated tabix file for the fragment.

### Step-1A: Remove the fragments with end sequence coordinates higher than the chromosomal sequence length:

Script file: [scripts/filter_coordinates.sh](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/scripts/filter_coordinates.sh)
```
chmod 775 filter_coordinates.sh
./filter_coordinates.sh
wc -l ExampleData/fragments_corrected_dedup_count.tsv.gz
## 2069308
wc -l results/filtered_fragments_corrected_dedup_count.tsv.gz
## 2012859
```
**Output file:** [results/filtered_fragments_corrected_dedup_count.tsv.gz](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/results/filtered_fragments_corrected_dedup_count.tsv.gz) #First 5000 lines due to space issue in Github

### Step-1B: Rename the sequence name to chr1-chr19 from Chr01-Chr19:

Script file: [scripts/rename_fragment_file.sh](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/scripts/rename_fragment_file.sh)
```
chmod 775 rename_fragment_file.sh
./rename_fragment_file.sh
```
**Output file:** [results/renamed_filtered_fragments_corrected_dedup_count.tsv.gz](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/results/renamed_filtered_fragments_corrected_dedup_count.tsv.gz) #First 5000 lines due to space issue in Github

### Step-1C: Making the updated tabix file for the fragment:
```
#Here, we need the tabix function from the conda environment
for first time install tabix in conda environment use: conda activate your_environment; conda install -c bioconda htslib #for tabix
module load conda3-py311_23.11.0-2
conda activate /mnt/dv/wid/projects6/Roy-singlecell3/bartholomew_lab/suvo_work/ArchR_analysis/scripts/suvo_ArchR
gunzip renamed_filtered_fragments_corrected_dedup_count.tsv.gz # Decompress the fragment file
dos2unix renamed_filtered_fragments_corrected_dedup_count.tsv # Convert to Unix format (remove Windows-style line endings)
bgzip renamed_filtered_fragments_corrected_dedup_count.tsv # Recompress using bgzip (for Tabix compatibility)
tabix -s 1 -b 2 -e 3 renamed_filtered_fragments_corrected_dedup_count.tsv.gz # Index with Tabix
```
**Output file:** [results/renamed_filtered_fragments_corrected_dedup_count.tsv.gz.tbi](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/results/renamed_filtered_fragments_corrected_dedup_count.tsv.gz.tbi) #First 5000 lines due to space issue in Github



## Step-2: Updating the fasta file discrepancies:

There are two sub-steps to remove fragment file discrepancies,
**2A:** Rename the sequence name to chr1-chr19 from Chr01-Chr19.
**2B:** Converting .fa file .2bit file compatible for R environment.

### Step-2A: Rename the sequence name to chr1-chr19 from Chr01-Chr19:

Script file: [scripts/filter_fasta.sh](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/scripts/filter_fasta.sh)
```
chmod 775 filter_fasta.sh
./filter_fasta.sh
```
**Output file:** [results/filtered_genome.fa](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/results/filtered_genome.fa) #First 5000 lines due to space issue in Github

### Step-2B: Converting .fa file .2bit file compatible for R environment:

Script file: [scripts/convert_fasta_to_twobit.sh](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/scripts/convert_fasta_to_twobit.sh)
```
chmod 775 convert_fasta_to_twobit.sh
./convert_fasta_to_twobit.sh
```
**Output file:** [results/filtered_genome.2bit](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/results/filtered_genome.2bit) #First 5000 lines due to space issue in Github



## Step-3: Updating the GTF file discrepancies:

The sequence name is now renamed to chr1-chr19 from Chr01-Chr19 in the GTF file.

Script file: [scripts/python gtf_edit.py](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/scripts/gtf_edit.py)
```
python python gtf_edit.py
```
**Output file:** [results/renamed_gene.gtf](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/results/renamed_gene.gtf) #First 5000 lines due to space issue in Github

<br/>
<br/>

> Now the input files are prepared for ArchR input, which are as follows:

**Fragment file:** [results/renamed_filtered_fragments_corrected_dedup_count.tsv.gz](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/results/renamed_filtered_fragments_corrected_dedup_count.tsv.gz)

**Tabix index file for the fragment file:** [results/renamed_filtered_fragments_corrected_dedup_count.tsv.gz.tbi](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/results/renamed_filtered_fragments_corrected_dedup_count.tsv.gz.tbi)

**2bit fasta file:** [results/filtered_genome.2bit](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/results/filtered_genome.2bit)

**GTF file:** [results/renamed_gene.gtf](https://github.com/Roy-lab/PoplarATACAnalysis/blob/main/results/renamed_gene.gtf)

<br/>

## Step-4: Forging the BSGenome object for making a custom Poplar genome R package in the R environment:


Forging the BSGenome object for the Poplar genome in the R environment is essential as this custom Poplar genome R package will help to generate a custom **genome annotation** object in ArchR.
<br/>
**Requires seed file:** scripts/BSgenome-Ptrichocarpa-v4.1-seed
<br/>
**NOTE:** This seed file should not have an extension like .txt and can be made in Mac textEditor or other text editors.

```
cd <your_working_directory>
conda activate /mnt/dv/wid/projects6/Roy-singlecell3/bartholomew_lab/suvo_work/ArchR_analysis/scripts/suvo_ArchR
Next, open the R environment
R
 > library(BSgenomeForge)
 > seqs_srcdir <- "/mnt/dv/wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR/"
 > destdir <- "/mnt/dv/wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR"
 > forgeBSgenomeDataPkg("BSgenome-Ptrichocarpa-v4.1-seed", seqs_srcdir = seqs_srcdir, destdir = destdir, replace = TRUE)
 > q()
```

**Output forged R package for Poplar genome database:** BSgenome.pPopTri
<br/>
**NOTE:** This package folder could not be uploaded to Github due to its size; however, it is already [uploaded](https://drive.google.com/drive/folders/1r1JXM4sgSlrP5jL7wl6Y20Fc-wtKtX0e?usp=drive_link) in Roy Lab shared [Google Drive folder](https://drive.google.com/drive/folders/1VndDmLDSP_Xvv40U15cwJc7C1K4mz06_?usp=drive_link). Anyone must download the package into the working directory to make the custom "**genome annotation**" object in ArchR.

<br/>

## Step-5: Forging the Taxonomic (Txdb) and Organism (Org.db) database object for making a custom Poplar genome R package in R environment:


Now, making the Taxonomic database (Txdb) and Organism (Org.db) database for the Poplar genome in the R environment is essential as this custom Poplar genome R package will help to generate custom **gene annotation** object in ArchR.
<br/>
**Requires BSgenome object:** BSgenome.pPopTri
<br/>
```
cd <your_working_directory>
conda activate /mnt/dv/wid/projects6/Roy-singlecell3/bartholomew_lab/suvo_work/ArchR_analysis/scripts/suvo_ArchR
Next, open the R environment
R
 > library(txdbmaker)
 > txdb_ptri <- makeTxDbFromGFF("renamed_gene.gtf")
 > saveDb(txdb_ptri,file = "Ptrichocarpa_533_v4.1-genes.txdb")
 > library(AnnotationForge)
 > install.packages("BSgenome.pPopTri", repos=NULL, type="source", lib="/mnt/dv/wid/projects6/Roy-singlecell3/bartholomew_lab/suvo_work/ArchR_analysis/scripts/suvo_ArchR/lib/R/library")
 > library(BSgenome.pPopTri)
 > ptri_gtf_msymbols <- import("renamed_gene.gtf")
 > gene_info <- unique(mcols(subset(ptri_gtf_msymbols, type == "transcript"))[,c("gene_id","transcript_id")])
 > NAS <- is.na(gene_info$gene_id)
 > gene_info$gene_id[NAS] <- gene_info$gene_id[NAS]
 > fSym <- unique(gene_info[, c(1,2)])
 > gene_info$symbol <- gene_info$gene_id
 > fSym <- unique(gene_info[, c(1,3)])
 > colnames(fSym) <- c("GID", "SYMBOL")
 > ensembl_trans <- unique(gene_info[, c(1:2)])
 > colnames(ensembl_trans) <- c("GID", "ENSEMBLTRANS")
 > ensembl <- unique(gene_info[, c(1,1)])
 > colnames(ensembl) <- c("GID", "ENSEMBL")
 > library(AnnotationForge)
 > makeOrgPackage(gene_info = fSym,
  ensembl_trans = ensembl_trans,
  ensembl = ensembl,
  version = "0.1",
  maintainer = "Suvojit Hazra <shazra3@wisc.edu>",
  author = "shazra3@wisc.edu",
  outputDir = "/mnt/dv/wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR/",
  tax_id= "3694",
  genus = "Populus",
  species = "trichocarpa",
  goTable=NULL,
  verbose = T)
 > q()
```

**Output forged R package for Poplar Organism database:** org.Ptrichocarpa.eg.db
<br/>
**NOTE:** This package folder could not be uploaded in Github due to its size however, it is already [uploaded](https://drive.google.com/drive/folders/1WyguU6_rvE6Sq-2rybPj7uqXGuBeSpoZ?usp=drive_link) in Roy Lab shared [Google Drive folder](https://drive.google.com/drive/folders/1VndDmLDSP_Xvv40U15cwJc7C1K4mz06_?usp=drive_link). Anyone must download the package into the working directory to make the custom "**gene annotation**" object in ArchR.
<br/>
**Output forged R package for Poplar Taxonomic database:** Ptrichocarpa_533_v4.1-genes.txdb
<br/>
**NOTE:** This package folder could not be uploaded in Github due to its size however, it is already [uploaded](https://drive.google.com/file/d/1AAKP3uYh6yAoTULKWwcNMMA1OcV7velU/view?usp=drive_link) in Roy Lab shared [Google Drive folder](https://drive.google.com/drive/folders/1VndDmLDSP_Xvv40U15cwJc7C1K4mz06_?usp=drive_link). Anyone must download the package into the working directory to make the custom "**gene annotation**" object in ArchR.
<br/>


## Step-6: Running ArchR:

**Followed:** [latest ArchR manual](https://www.archrproject.com/bookdown/index.html) and [latest ArchR tutorial](https://www.archrproject.com/articles/Articles/tutorial.html)

As all the input files and packages are prepared, now the stepwise ArchR run is as follows:
```
#pwd=/mnt/dv/wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR/
module load gcc-13.2.0 
module load conda3-py311_23.11.0-2 
conda activate /mnt/dv/wid/projects6/Roy-singlecell3/bartholomew_lab/suvo_work/ArchR_analysis/scripts/suvo_ArchR
R
> library(ArchR)
> library(GenomicFeatures)
> getwd()
[1] "/mnt/dv/wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR"
> install.packages("org.Ptrichocarpa.eg.db", repos=NULL, type="source")
> library(org.Ptrichocarpa.eg.db) #adding oragnism database annotation library
> install.packages("BSgenome.pPopTri", repos=NULL, type="source")
> library(BSgenome.pPopTri) #adding genome annotation library
> txdb_obj <- loadDb("Ptrichocarpa_533_v4.1-genes.txdb") #adding taxonomy database annotation library
> library(tidyverse)
> geneAnnotation <- createGeneAnnotation(TxDb = txdb_obj, OrgDb = org.Ptrichocarpa.eg.db) #creating gene annotation object
> genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.pPopTri) #creating genome annotation object
> head_dir <- '/mnt/dv/wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR/'
> filenames = list('renamed_filtered_fragments_corrected_dedup_count.tsv.gz') #reading fragment file
> library(Cairo) #this helps in generating figure files (.pdf), for this R version needs downgradation from 4.4.2 to 4.3.1, otherwise, figures are not generating.
> library(ggplot2)
> filenames = list('renamed_filtered_fragments_corrected_dedup_count.tsv.gz') #checking fragment file is being read into R.
> filenames = as.character(filenames)
> input_files = paste0(head_dir, filenames)
> setnames = list('Poplar_stem')
> setnames = as.character(setnames)
> names(input_files) <- setnames
> file.exists(input_files)
[1] TRUE
> ArrowFiles <- createArrowFiles(inputFiles = input_files, geneAnnotation = geneAnnotation, genomeAnnotation = genomeAnnotation, minTSS = 1, minFrags = 1, addTileMat = TRUE, TileMatParams = list(tileSize = 1000), addGeneScoreMat = TRUE) #Creating arrowfiles.
# 'minTSS' refers to the minimum number of transcription start sites (TSS) required for a cell to be included in the analysis.
# 'minFrags' refers to the minimum number of fragments required for a cell to be included in the analysis.
# for human and mouse dataset, minTSS=4, minFrags=1000 filter application gives good results.
> library(nabor) #required package for synthetic doublet removal
> doubScores <- addDoubletScores(input = ArrowFiles,k = 10, knnMethod = "UMAP", LSIMethod = 1)
#'k' Refers to how many cells near a "pseudo-doublet" to count.
#'knnMethod' Refers to the embedding to use for nearest neighbor search with doublet projection.
> saveRDS(genomeAnnotation, "genomeAnnotation.rds") #saving genome annotation object for future use.
> saveRDS(geneAnnotation, "geneAnnotation.rds") #saving gene annotation object for future use.
> genomeAnnotation <- readRDS("genomeAnnotation.rds") #reading genome annotation object from R data file (.rds) which was just saved.
> geneAnnotation <- readRDS("geneAnnotation.rds") #reading gene annotation object from R data file (.rds) which was just saved.
> PoplarStem <- ArchR::ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "/mnt/dv/wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR/PoplarStem", geneAnnotation = geneAnnotation, genomeAnnotation = genomeAnnotation, copyArrows = TRUE) #Creating ArchR object 'PoplarStem' for future use.
> saveArchRProject(PoplarStem) #Saving ArchR object 'PoplarStem' for future use.
> summarize(as.data.frame(PoplarStem@cellColData), count = n(), .by=Sample) #Summarize the number of cells in ArchR object 'PoplarStem'.
##       Sample count
##1 Poplar_stem 12855
> p1 <- plotGroups(ArchRProj = PoplarStem, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
> p2 <- plotGroups(ArchRProj = PoplarStem, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
> p3 <- plotGroups(ArchRProj = PoplarStem, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "ridges")
> p4 <- plotGroups(ArchRProj = PoplarStem, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
> plotPDF(p1, p2, p3, p4, name = "QC-Sample-Statistics.pdf",ArchRProj = PoplarStem, width = 4, height = 4) #Generating Quality Control (QC) figure.
> getwd()
##[1] "/mnt/dv/wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR"
> save(file ='PoplarStem.Rdata', list = ls()) #Save the project into an R data file for future use.
> addGeneScoreMatrix(input = PoplarStem, geneModel = "exp(-abs(x)/5000) + exp(-1)", extendUpstream = c(1000, 1e+05), extendDownstream = c(1000, 1e+05), useTSS = TRUE, extendTSS = 1000, tileSize = 1000, force = TRUE) #Adding genescore matrix if previously 'addGeneScoreMat' set to "FALSE" during creation of ArrowFiles.
> save(file ='PoplarStem.Rdata', list = ls())
> Gene_score <- getMatrixFromProject(ArchRProj = PoplarStem, useMatrix = "GeneScoreMatrix")
> Gene_info<-as_tibble(rowData(Gene_score))
> Gene_info #Information of genes.
> write_tsv(Gene_info, 'gene_info.tsv')
> Gene_id<-as.data.frame(Gene_info$name)                                     
> Gene_id #Gene IDs.
> write.table(Gene_id, 'gene_ids.tsv', col.names = FALSE, row.names = FALSE)
> cell_info <- as_tibble(colData(Gene_score))
> cell_info #Information of cells.
> write_tsv(cell_info, 'cell_info.tsv')
> cell_ids <- colnames(Gene_score)                                           
> cell_ids
> cell_ids <- tibble(cell_ids)
> cell_ids <- cell_ids %>% rename('barcodes' = cell_ids)
> cell_ids #Cell IDs.
> write_tsv(cell_ids, 'cell_ids.tsv')
> scores <- as.matrix(assay(Gene_score, 1))                                  
#Warning message: In asMethod(object) : sparse->dense coercion: allocating vector of size 3.3 GiB
> scores #Genescore matrix.
> write.table(scores, 'gene_scores.tsv', sep = '\t', col.names =FALSE, row.names = FALSE)

Output folder: /mnt/dv/wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR/

Output sub-folders:

ArchR log files: ArchRLogs
QualityControl: QC stats
Temporary folder: tmp
Saved project: PoplarStem

Output files:

Genescore matrix (gene x cells; dimesion: 34488 x 12856): gene_scores.tsv 
Genes data from genescore matrix (Gene: 34488): gene_ids.tsv 
Cells data from genescore matrix (Cell: 12856): cell_ids.tsv
```

**NOTE:** This above-mentioned ArchR R script is demonstrated as it was run previously. The steps are replicative; make sure to change the filename and folders as per the analysis.




