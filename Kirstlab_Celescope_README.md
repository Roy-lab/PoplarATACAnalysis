# Executing and understanding the outputs of [CeleScope_ATAC](https://github.com/singleron-RD/CeleScope_ATAC)

## Installing CeleScope_ATAC

Start an interactive session on the HPC

```sh
srundev --time=08:00:00 --partition=hpg-dev --mem=96G --cpus-per-task=24
module load mamba
```

Creating the conda environment.

```sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install mamba -c conda-forge

cd CeleScope_ATAC
mamba create -n celescope_atac -y --file conda_pkgs.txt

mamba activate celescope_atac
pip install .
```

## Executing CeleScope_ATAC - Reference genome

```sh
mamba activate celescope_atac

cp /blue/kirst/share/Genomes/Ptrichocarpa_v4.1/Phytozome/PhytozomeV13/Ptrichocarpa/v4.1/assembly/Ptrichocarpa_533_v4.0.softmasked.fa.gz .
cp /blue/kirst/share/Genomes/Ptrichocarpa_v4.1/Phytozome/PhytozomeV13/Ptrichocarpa/v4.1/annotation/Ptrichocarpa_533_v4.1.gene_exons.gff3.gz .

gunzip -k Ptrichocarpa_533_v4.1.gene_exons.gff3.gz
gunzip -k Ptrichocarpa_533_v4.0.softmasked.fa.gz

ml gffread/0.12.7
gffread Ptrichocarpa_533_v4.1.gene_exons.gff3 -T -o Ptrichocarpa_533_v4.1.gtf

celescope atac mkref --fasta Ptrichocarpa_533_v4.0.softmasked.fa --gtf Ptrichocarpa_533_v4.1.gtf
```

The commands above will create a folder called `shell`. Within this folder, are the outputed scripts to execute the peak definition and counts generation.

## Executing - peak definition

The commands below are the same outputed in the `shell` folder. To try to understand what is happening in each step, I print the results of each command using `tree Poplar_tricocharpa_outputs`.


```sh
celescope atac sample --outdir ./Poplar_tricocharpa_outputs/00.sample --sample Poplar_tricocharpa_outputs --thread 24 --chemistry atac2 --wells 384 --fq2 fastq/24073FL-20-01-30_S2_L001_R2_001.fastq.gz 
```

```console
Poplar_tricocharpa_outputs/
├── 00.sample
│   └── stat.txt
└── Poplar_tricocharpa_outputs_report.html

1 directory, 2 files
```

```sh
celescope atac barcode --outdir ./Poplar_tricocharpa_outputs/01.barcode --sample Poplar_tricocharpa_outputs --thread 24 --chemistry atac2 --lowNum 2 --wells 384 --fq1 fastq/24073FL-20-01-30_S2_L001_R1_001.fastq.gz --fq2 fastq/24073FL-20-01-30_S2_L001_R2_001.fastq.gz --fq3 fastq/24073FL-20-01-30_S2_L001_R3_001.fastq.gz
```

```console
Poplar_tricocharpa_outputs/
├── 00.sample
│   └── stat.txt
├── 01.barcode
│   ├── Poplar_tricocharpa_outputs_S1_L001_R1_001.fastq
│   ├── Poplar_tricocharpa_outputs_S1_L001_R2_001.fastq
│   ├── Poplar_tricocharpa_outputs_S1_L001_R3_001.fastq
│   └── stat.txt
└── Poplar_tricocharpa_outputs_report.html
```

```sh
celescope atac atac --outdir ./Poplar_tricocharpa_outputs/02.atac --sample Poplar_tricocharpa_outputs --thread 24 --reference ./ --genomesize 402564050 --peak_cutoff auto --count_cutoff 500 --frip_cutoff 0.2 --cell_cutoff 1 --expected_target_cell_num 8000 --coef 0.1  --input_path ./Poplar_tricocharpa_outputs/01.barcode 
```

```console
Poplar_tricocharpa_outputs
├── 00.sample
│   └── stat.txt
├── 01.barcode
│   ├── Poplar_tricocharpa_outputs_S1_L001_R1_001.fastq
│   ├── Poplar_tricocharpa_outputs_S1_L001_R2_001.fastq
│   ├── Poplar_tricocharpa_outputs_S1_L001_R3_001.fastq
│   └── stat.txt
├── 02.atac
│   ├── cell_qc_metrics.tsv
│   ├── fragments_corrected_150bp.tsv
│   ├── fragments_corrected_count_sortedbybarcode.tsv
│   ├── fragments_corrected_dedup_count.tsv
│   ├── fragments_corrected_dedup_count.tsv.gz
│   ├── fragments_corrected_dedup_count.tsv.gz.tbi
│   ├── fragments_promoter_sortbybarcode.tsv
│   ├── fragments_promoter.tsv
│   ├── meta.csv
│   ├── peak
│   │   ├── Poplar_tricocharpa_outputs_150bp_control_lambda.bdg
│   │   ├── Poplar_tricocharpa_outputs_150bp_peaks.narrowPeak
│   │   ├── Poplar_tricocharpa_outputs_150bp_peaks.xls
│   │   ├── Poplar_tricocharpa_outputs_150bp_summits.bed
│   │   ├── Poplar_tricocharpa_outputs_150bp_treat_pileup.bdg
│   │   ├── Poplar_tricocharpa_outputs_control_lambda.bdg
│   │   ├── Poplar_tricocharpa_outputs_filtered_peak_count.h5
│   │   ├── Poplar_tricocharpa_outputs_final_peaks.bed
│   │   ├── Poplar_tricocharpa_outputs_peak_count.h5
│   │   ├── Poplar_tricocharpa_outputs_peaks.narrowPeak
│   │   ├── Poplar_tricocharpa_outputs_peaks.xls
│   │   ├── Poplar_tricocharpa_outputs_summits.bed
│   │   └── Poplar_tricocharpa_outputs_treat_pileup.bdg
│   ├── Poplar_tricocharpa_outputs_cluster.png
│   ├── Poplar_tricocharpa_outputs.rds
│   ├── singlecell_mapped_sortbybarcode.txt
│   ├── singlecell_mapped.txt
│   ├── singlecell_promoter_sortbybarcode.txt
│   ├── singlecell_promoter.txt
│   ├── singlecell.txt
│   ├── stat.txt
│   └── validcells.txt
└── Poplar_tricocharpa_outputs_report.html
```

```sh
celescope atac analysis --outdir ./Poplar_tricocharpa_outputs/03.analysis --sample Poplar_tricocharpa_outputs --thread 24  --analysis_dir ./Poplar_tricocharpa_outputs/02.atac --filtered_peak_count ./Poplar_tricocharpa_outputs/02.atac/peak/Poplar_tricocharpa_outputs_filtered_peak_count.h5 
```

```console
Poplar_tricocharpa_outputs/
├── 00.sample
│   └── stat.txt
├── 01.barcode
│   ├── Poplar_tricocharpa_outputs_S1_L001_R1_001.fastq
│   ├── Poplar_tricocharpa_outputs_S1_L001_R2_001.fastq
│   ├── Poplar_tricocharpa_outputs_S1_L001_R3_001.fastq
│   └── stat.txt
├── 02.atac
│   ├── cell_qc_metrics.tsv
│   ├── fragments_corrected_150bp.tsv
│   ├── fragments_corrected_count_sortedbybarcode.tsv
│   ├── fragments_corrected_dedup_count.tsv
│   ├── fragments_corrected_dedup_count.tsv.gz
│   ├── fragments_corrected_dedup_count.tsv.gz.tbi
│   ├── fragments_promoter_sortbybarcode.tsv
│   ├── fragments_promoter.tsv
│   ├── meta.csv
│   ├── peak
│   │   ├── Poplar_tricocharpa_outputs_150bp_control_lambda.bdg
│   │   ├── Poplar_tricocharpa_outputs_150bp_peaks.narrowPeak
│   │   ├── Poplar_tricocharpa_outputs_150bp_peaks.xls
│   │   ├── Poplar_tricocharpa_outputs_150bp_summits.bed
│   │   ├── Poplar_tricocharpa_outputs_150bp_treat_pileup.bdg
│   │   ├── Poplar_tricocharpa_outputs_control_lambda.bdg
│   │   ├── Poplar_tricocharpa_outputs_filtered_peak_count.h5
│   │   ├── Poplar_tricocharpa_outputs_final_peaks.bed
│   │   ├── Poplar_tricocharpa_outputs_peak_count.h5
│   │   ├── Poplar_tricocharpa_outputs_peaks.narrowPeak
│   │   ├── Poplar_tricocharpa_outputs_peaks.xls
│   │   ├── Poplar_tricocharpa_outputs_summits.bed
│   │   └── Poplar_tricocharpa_outputs_treat_pileup.bdg
│   ├── Poplar_tricocharpa_outputs_cluster.png
│   ├── Poplar_tricocharpa_outputs.rds
│   ├── singlecell_mapped_sortbybarcode.txt
│   ├── singlecell_mapped.txt
│   ├── singlecell_promoter_sortbybarcode.txt
│   ├── singlecell_promoter.txt
│   ├── singlecell.txt
│   ├── stat.txt
│   └── validcells.txt
├── 03.analysis
│   ├── cell_qc_metrics.tsv
│   ├── fragments_corrected_dedup_count.tsv.gz
│   ├── fragments_corrected_dedup_count.tsv.gz.tbi
│   ├── Poplar_tricocharpa_outputs_filtered_peak_count.h5
│   ├── Poplar_tricocharpa_outputs_final_peaks.bed
│   ├── Poplar_tricocharpa_outputs_peak_count.h5
│   ├── Poplar_tricocharpa_outputs.rds
│   └── stat.txt
├── outs
│   ├── cell_qc_metrics.tsv
│   ├── fragments_corrected_dedup_count.tsv.gz
│   ├── fragments_corrected_dedup_count.tsv.gz.tbi
│   └── Poplar_tricocharpa_outputs_filtered_peak_count.h5
└── Poplar_tricocharpa_outputs_report.html
```