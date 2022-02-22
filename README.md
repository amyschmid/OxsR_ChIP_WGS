# Data analysis for whole genome sequencing and ChIP-sequencing in Haloferax volcanii strains

## Corresponds to the paper "TrmB-like family transcription factor OxsR regulates the oxidative stress response in Haloferax volcanii" by Schmid group and Maupin-Furlow group.

### Dependencies and instructions in each analysis

#### 1) Whole genome sequencing
* fastqc v0.11.9
* bowtie2 v2.4.5
* breseq v0.35.1 <br/>
-Installation of packages, downloading raw sequencing data ([Accession #: PRJNA806939](https://www.ncbi.nlm.nih.gov/bioproject/806939)) and [reference genome of H. volcanii](https://www.ncbi.nlm.nih.gov/genome/?term=haloferax+volcanii) <br/>
-Execute a command: _breseq -p -r GENOME_FILE.gbff FORWARD_READ_FILE.fq REVERSE_READ_FILE.fq -o OUTPUT_ <br/>

#### 2) ChIP-seq
* The most recent versions of R, RStudio and R packages
  * [R](https://cran.r-project.org/) (R v4.1.0 in this study)
  * [RStudio](https://www.rstudio.com/products/rstudio/download/#download) (RStudio v1.3.1093 in this study)
  * R packages from Bioconductor, using `BiocManager::install()` function
```r
BiocManager::install("mosaics")
BiocManager::install("hexbin")
BiocManager::install("tidyverse")
BiocManager::install("ChIPQC")
BiocManager::install("DiffBind")
BiocManager::install("ChIPseeker")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicFeatures")
BiocManager::install("IRanges")
BiocManager::install("AnnotationHub")
```
-Installation of packages, downloading raw sequencing data ([GEO accession #GSE196894](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196894)) and [reference genome of H. volcanii](https://www.ncbi.nlm.nih.gov/genome/?term=haloferax+volcanii) <br/>
-Open the Rmd file contains command lines for the analysis <br/>
-Run analysis using raw data and [metadata](https://github.com/sungminhwang-duke/OxsR_ChIP_WGS/tree/master/Meta_data) files deposited <br/>

#### 2) arCOG functional enrichment analysis
in R or Rstudio, install packages with the following versions (or newer):
 * tidyverse 1.3.1
 * dplyr 1.0.7
 * openxlsx 4.2.4
 * ggplot2 3.3.5
 * bazar 1.0.11

The run "2022-02-07-arcog-OxsR.Rmd" with input files:
* oxsR-chipseq.xlsx
* arcogs-14-18.hvo.txt
* 20181113_hvol_GCF_000025685.1_ASM2568v1_genomic.gff.key.csv
