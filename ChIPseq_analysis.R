############# OxsR ChIP-seq analysis #############

### Unzip raw data
#   $ gunzip *.gz







### Trim adapter or low quality sequences in local laptop using Terminal followed by Quality check by Fastqc
#   $ trim_galore -o ./trimming *.fastq --paired --fastqc







### Index and Alignment in local laptop using Terminal
#   $ bowtie2-build HVO.fna HVO
#   $ bowtie2 -x ./Index/HVO -1 H26_1_IP_S74_S80_R1_val_1.fq -2 H26_1_IP_S74_S80_R2_val_2.fq -S ./sam/H26_1_IP_S74_S80.sam
#     = for f1 in *R1_val_1.fq; do for f2 in ${f1%%_R1_val_1.fq}"_R2_val_2.fq" ; do bowtie2 -x ./Index/HVO -1 $f1 -2 $f2 -S ./sam/$f1.sam ; done; done



### Bam file generation
# 1)  $ samtools view -bS S49_H26_1_glc_IP.sam > S49_H26_1_glc_IP.bam
#     = for file in ./*.sam; do samtools view -bS $file -o ./bam/$file.bam ; done


# 2)  $ samtools sort S49_H26_1_glc_IP.bam > S49_H26_1_glc_IP_sorted.bam
#     = for file in ./*.bam; do samtools sort $file -o ./bam_sort/$file.sort.bam ; done
#     *Change the file name from ~.sam.bam.sort.bam to(->) ~_sort.bam


# 3)  $ samtools index S49_H26_1_glc_IP_sorted.bam
#     = for file in ./*.bam; do samtools index $file ; done











### Peak calling from pair-end data by MACS2
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html
# https://github.com/taoliu/MACS


#macs2 callpeak -t H26_1_IP_S74_S80_sort.bam -c H26_1_WCE_S65_S71_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n H26_1_no 2> macs2/H26_1_no-macs2.log
#macs2 callpeak -t H26_2_IP_S91_S97_sort.bam -c H26_2_WCE_S95_S43_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n H26_2_no 2> macs2/H26_2_no-macs2.log
#macs2 callpeak -t oxsRHA_1_IP_S40_S46_sort.bam -c oxsRHA_1_WCE_S90_S96_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n oxsRHA_1_no 2> macs2/oxsRHA_1_no-macs2.log
#macs2 callpeak -t oxsRHA_2_IP_S57_S63_sort.bam -c oxsRHA_2_WCE_S70_S76_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n oxsRHA_2_no 2> macs2/oxsRHA_2_no-macs2.log
#macs2 callpeak -t oxsRHA_3_IP_S59_S65_sort.bam -c oxsRHA_3_WCE_S96_S44_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n oxsRHA_3_no 2> macs2/oxsRHA_3_no-macs2.log
#macs2 callpeak -t oxsRHA_4_IP_S82_S88_sort.bam -c oxsRHA_4_WCE_S80_S86_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n oxsRHA_4_no 2> macs2/oxsRHA_4_no-macs2.log

#macs2 callpeak -t H26_1_NaOCl_IP_S81_S87_sort.bam -c H26_1_NaOCl_WCE_S86_S92_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n H26_1_NaOCl 2> macs2/H26_1_NaOCl-macs2.log
#macs2 callpeak -t H26_2_NaOCl_IP_S60_S66_sort.bam -c H26_2_NaOCl_WCE_S51_S57_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n H26_2_NaOCl 2> macs2/H26_2_NaOCl-macs2.log
#macs2 callpeak -t oxsRHA_1_NaOCl_IP_S83_S89_sort.bam -c oxsRHA_1_NaOCl_WCE_S49_S55_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n oxsRHA_1_NaOCl 2> macs2/oxsRHA_1_NaOCl-macs2.log
#macs2 callpeak -t oxsRHA_2_NaOCl_IP_S89_S95_sort.bam -c oxsRHA_2_NaOCl_WCE_S42_S48_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n oxsRHA_2_NaOCl 2> macs2/oxsRHA_2_NaOCl-macs2.log
#macs2 callpeak -t oxsRHA_3_NaOCl_IP_S77_S83_sort.bam -c oxsRHA_3_NaOCl_WCE_S39_S45_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n oxsRHA_3_NaOCl 2> macs2/oxsRHA_3_NaOCl-macs2.log
#macs2 callpeak -t oxsRHA_4_NaOCl_IP_S87_S93_sort.bam -c oxsRHA_4_NaOCl_WCE_S84_S90_sort.bam -f BAMPE -g 4.1e6 -q 0.05 --outdir macs2 -n oxsRHA_4_NaOCl 2> macs2/oxsRHA_4_NaOCl-macs2.log


     







### Peakcalling quality control I by ChIP_QC 
#BiocManager::install(c("GenomicRanges","rtracklayer", "ChIPQC"))

# To fix a bug in QCmetrics(ChIPQCsample) which will cause the error: "names' attribute [9] must be the same length as the vector [7]" (https://github.com/shengqh/ChIPQC)
library(devtools)
install_github("shengqh/ChIPQC")

library(tidyverse)
library(ChIPQC)

# Creating a custom annotation track from .gff file
txdb <- GenomicFeatures::makeTxDbFromGFF("GCF_000025685.1_ASM2568v1_genomic.gff", format = "gff")

#reduce(unique(unlist(GenomicFeatures::cdsBy(txdb, "tx"))))
txn <- GenomicFeatures::transcripts(txdb)
gene <- unlist(GenomicFeatures::cdsBy(txdb, "tx"))
pro500 <- GenomicFeatures::promoters(txdb, upstream = 500, downstream = 0)
pro250 <- GenomicFeatures::promoters(txdb, upstream = 250, downstream = 0)

hvo <- list(version="",
            gff.features = txn,
            genes = gene,
            promoters250 = pro250,
            promoters500 = pro500)

# load sample file
samples <- as.data.frame(read_csv("Meta_chipqc_NaOCl_1.csv"))

# create and execute the experiment!
#register(SerialParam()) # prevents BiocParallel error (Run this if you use Windows!)
exp <- ChIPQC(experiment = samples, annotation = hvo)
exp
ChIPQCreport(exp, reportFolder ="ChIPQC report", reportName="ChIPQC_gc")

mean(head(as.data.frame(QCmetrics(exp))$FragL, n=8))
mean(tail(as.data.frame(QCmetrics(exp))$FragL, n=8))

# Use fragment lengths of 150 for IP samples and 220 for WCE sampels in final peak calling. 
plotCC(exp, facetBy = c("Condition", "Factor"))

#FragmentLengthCrossCoverage(exp)[1:8]/FragmentLengthCrossCoverage(exp)[9:16]
#ReadLengthCrossCoverage(exp)[1:8]/ReadLengthCrossCoverage(exp)[9:16]
#RelativeCrossCoverage(exp)
abs(RelativeCrossCoverage(exp)[1:8]/RelativeCrossCoverage(exp)[9:16])
abs(RelativeCrossCoverage(exp2)[1:8]/RelativeCrossCoverage(exp2)[9:16])

QCmetadata(exp)
#for (i in 1:5) {
#  plot(crosscoverage(QCsample(exp,i)),type='l',
#       ylab="Cross-coverage",
#       xlab="Fragment length", main = rownames(QCmetadata(exp))[i])}

# Single sample assessment 
tbsPHA_1_glc_IP = ChIPQCsample("data/tbsPHA_1_glc_IP_S79_S85_sorted.bam", 
                               peaks="macs2/tbsPHA_1_glc_peaks.narrowPeak", annotation = hvo)
tbsPHA_1_glc_WCE = ChIPQCsample("data/tbsPHA_1_glc_WCE_S48_S54_sorted.bam", 
                               peaks="macs2/tbsPHA_1_glc_peaks.narrowPeak", annotation = hvo)
tbsPHA_1_glc_IP
tbsPHA_1_glc_WCE

QCmetrics(tbsPHA_1_glc_IP)
QCmetrics(tbsPHA_1_glc_WCE)







### Peak calling quality control II by IDR
# https://github.com/nboley/idr
# Required r files: 1) batch-consistency-analysis.r, 2) batch-consistency-plot.r, 3) functions-all-clayton-12-13.r
# Required genome_table for the strain (i.e., Haloferax volcanii)
# Run the code below in terminal

#---- no oxidative stress
#  $ Rscript batch-consistency-analysis.r H26_1_no_peaks.narrowPeak H26_2_no_peaks.narrowPeak -1 ./IDR_analysis/H26_1vs_H26_2 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_1_no_peaks.narrowPeak oxsRHA_2_no_peaks.narrowPeak -1 ./IDR_analysis/oxsR_1_vs_oxsR_2 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_1_no_peaks.narrowPeak oxsRHA_3_no_peaks.narrowPeak -1 ./IDR_analysis/oxsR_1_vs_oxsR_3 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_1_no_peaks.narrowPeak oxsRHA_4_no_peaks.narrowPeak -1 ./IDR_analysis/oxsR_1_vs_oxsR_4 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_2_no_peaks.narrowPeak oxsRHA_3_no_peaks.narrowPeak -1 ./IDR_analysis/oxsR_2_vs_oxsR_3 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_2_no_peaks.narrowPeak oxsRHA_4_no_peaks.narrowPeak -1 ./IDR_analysis/oxsR_2_vs_oxsR_4 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_3_no_peaks.narrowPeak oxsRHA_4_no_peaks.narrowPeak -1 ./IDR_analysis/oxsR_3_vs_oxsR_4 0 F p.value

#  $ Rscript batch-consistency-plot.r 6 ./IDR_output/tbsP ./IDR_analysis/oxsR_1_vs_oxsR_2 ./IDR_analysis/oxsR_1_vs_oxsR_3 ./IDR_analysis/oxsR_1_vs_oxsR_4 ./IDR_analysis/oxsR_2_vs_oxsR_3 ./IDR_analysis/oxsR_2_vs_oxsR_4 ./IDR_analysis/oxsR_3_vs_oxsR_4

#---- NaOCl present
#  $ Rscript batch-consistency-analysis.r oxsRHA_1_NaOCl_peaks.narrowPeak oxsRHA_2_NaOCl_peaks.narrowPeak -1 ./IDR_analysis/oxsR_1_vs_oxsR_2 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_1_NaOCl_peaks.narrowPeak oxsRHA_3_NaOCl_peaks.narrowPeak -1 ./IDR_analysis/oxsR_1_vs_oxsR_3 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_1_NaOCl_peaks.narrowPeak oxsRHA_4_NaOCl_peaks.narrowPeak -1 ./IDR_analysis/oxsR_1_vs_oxsR_4 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_2_NaOCl_peaks.narrowPeak oxsRHA_3_NaOCl_peaks.narrowPeak -1 ./IDR_analysis/oxsR_2_vs_oxsR_3 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_2_NaOCl_peaks.narrowPeak oxsRHA_4_NaOCl_peaks.narrowPeak -1 ./IDR_analysis/oxsR_2_vs_oxsR_4 0 F p.value
#  $ Rscript batch-consistency-analysis.r oxsRHA_3_NaOCl_peaks.narrowPeak oxsRHA_4_NaOCl_peaks.narrowPeak -1 ./IDR_analysis/oxsR_3_vs_oxsR_4 0 F p.value

#  $ Rscript batch-consistency-plot.r 6 ./IDR_output/tbsP ./IDR_analysis/oxsR_1_vs_oxsR_2 ./IDR_analysis/oxsR_1_vs_oxsR_3 ./IDR_analysis/oxsR_1_vs_oxsR_4 ./IDR_analysis/oxsR_2_vs_oxsR_3 ./IDR_analysis/oxsR_2_vs_oxsR_4 ./IDR_analysis/oxsR_3_vs_oxsR_4




### Annotation and identification of peak on genome
# http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
#BiocManager::install("ChIPseeker")
#BiocManager::install("clusterProfiler")

library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(GenomicRanges)
library(IRanges)

# Load data (*Remove "the header" in the bed file or message like as Error in .Call2("solve_user_SEW0", start, end, width, PACKAGE = "IRanges") : In range 1: at least two out of 'start', 'end', and 'width', must be supplied.)
samplefiles <- list.files("./1_macs2", pattern= ".narrowPeak", full.names=TRUE) #Input file I
#samplefiles <- list.files("./macs2", pattern= ".bed", full.names=TRUE)        #Input file II
samplefiles <- as.list(samplefiles)
#names(samplefiles) <- c("H26_1", "H26_2", "oxsRHA_1", "oxsRHA_2", "oxsRHA_3", "oxsRHA_4")
names(samplefiles) <- c("oxsRHA_1", "oxsRHA_2", "oxsRHA_3", "oxsRHA_4")

# Database construction for the index genome from gff file deposited in NCBI
library(GenomicFeatures)
hvo_TxDb <- makeTxDbFromGFF("GCF_000025685.1_ASM2568v1_genomic.gff", 
                            organism="Haloferax volcanii", format = "gff")

# Analyze data one-by-one -> Recommend to go below for the group analysis
peakAnnoList <- lapply(samplefiles[[3]], annotatePeak, TxDb=hvo_TxDb, 
                       tssRegion=c(-500, 500), verbose=FALSE)

peak <- readPeakFile(samplefiles[[3]])
peak

covplot(peak, weightCol="V5")

promoter <- getPromoters(TxDb=hvo_TxDb, upstream=500, downstream=500)
tagMatrix <- getTagMatrix(peak, windows=promoter)

tagHeatmap(tagMatrix, xlim=c(-500, 500), color="red")

peakHeatmap(samplefiles[[3]], TxDb=hvo_TxDb, upstream=500, downstream=500, color="red")
plotAvgProf2(samplefiles[[3]], TxDb=hvo_TxDb, upstream=500, downstream=500,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagMatrix, xlim=c(-500, 500), conf = 0.95, resample = 500)

peakAnno <- annotatePeak(samplefiles[[3]], tssRegion=c(-500, 500), TxDb=hvo_TxDb, 
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                         addFlankGeneInfo = TRUE, verbose=FALSE)

gene.list <- as.data.frame(peakAnno)
write.csv(gene.list, "ChIP_seq_genes.csv")  # Results

plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

# Group analysis
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=hvo_TxDb, 
                       tssRegion=c(-500, 500), verbose=FALSE)

promoter <- getPromoters(TxDb=hvo_TxDb, upstream=500, downstream=500)
tagMatrixList <- lapply(samplefiles, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList, xlim=c(-500, 500))

plotAvgProf(tagMatrixList, xlim=c(-500, 500), conf=0.95,resample=500, facet="row")

tagHeatmap(tagMatrixList, xlim=c(-500, 500), color=NULL)

plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

gene.list <- as.data.frame(peakAnnoList[[4]])
write.csv(gene.list, "Peak_anno_oxsRHA_4.csv")  # Results






### Merge the bed files and identify the whole peaks profile - Codes from Rylee
# for two replicate comparisons:
#     $ bedtools intersect -f 0.5 -r -a replicate1.bed -b replicate2.bed > OUTFILE.bed

# for 3+ replicate comparisons, we need to approach this combinatorically:
#     $ bedtools intersect -f 0.5 -r -wb -names B C D -a oxsRHA_1_NaOCl_summits.bed -b oxsRHA_2_NaOCl_summits.bed oxsRHA_3_NaOCl_summits.bed oxsRHA_4_NaOCl_summits.bed > intersectABCD_bed.bed
#     $ bedtools intersect -f 0.5 -r -wb -names B C D -a oxsRHA_1_NaOCl_peaks.narrowPeak -b oxsRHA_2_NaOCl_peaks.narrowPeak oxsRHA_3_NaOCl_peaks.narrowPeak oxsRHA_4_NaOCl_peaks.narrowPeak > intersectABCD_narrowPeak.bed
# Code above gives peaks that are in ABCD, ABC, ABD, ACD, AB, AC, AD

#     $ bedtools intersect -f 0.5 -r -wb -names C D A -a oxsRHA_2_NaOCl_summits.bed -b oxsRHA_3_NaOCl_summits.bed oxsRHA_4_NaOCl_summits.bed oxsRHA_1_NaOCl_summits.bed > intersectBCDA_bed.bed
#     $ bedtools intersect -f 0.5 -r -wb -names C D A -a oxsRHA_2_NaOCl_peaks.narrowPeak -b oxsRHA_3_NaOCl_peaks.narrowPeak oxsRHA_4_NaOCl_peaks.narrowPeak oxsRHA_1_NaOCl_peaks.narrowPeak > intersectBCDA_narrowPeak.bed
# Code above gives peaks that are in BCD, BC, BD

# This is not essential if you are only interested in peaks that are conserved in 3+ replicates
#     $ bedtools intersect -f 0.5 -r -wa -a rep_C.bed -b rep_D.bed > intersectCD.bed

#I combine the peak files by hand. Doing so, I update the peak range to the narrowest range shared across replicates, 
#take the highest score, and update the peak name to reflect the replicates that shared the peak.







### Analysis of the merged bed file (at least three bio reps)
# Load data 
samplefiles <- list.files("./temp/", pattern= ".bed", full.names=TRUE)
samplefiles <- as.list(samplefiles)

# Database construction for the index genome from gff file deposited in NCBI
library(GenomicFeatures)
hvo_TxDb <- makeTxDbFromGFF("GCF_000025685.1_ASM2568v1_genomic.gff", 
                            organism="Haloferax volcanii", format = "gff")

# Analyze data one-by-one -> Recommend to go below for the group analysis
peakAnnoList <- lapply(samplefiles[[1]], annotatePeak, TxDb=hvo_TxDb, 
                       tssRegion=c(-500, 500), verbose=FALSE)

peak <- readPeakFile(samplefiles[[1]])
peak

covplot(peak, weightCol="V5")

promoter <- getPromoters(TxDb=hvo_TxDb, upstream=500, downstream=500)
tagMatrix <- getTagMatrix(peak, windows=promoter)

peakHeatmap(samplefiles[[1]], TxDb=hvo_TxDb, upstream=500, downstream=500, color="red")
plotAvgProf2(samplefiles[[1]], TxDb=hvo_TxDb, upstream=500, downstream=500,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagMatrix, xlim=c(-500, 500), conf = 0.95, resample = 500)

peakAnno <- annotatePeak(samplefiles[[1]], tssRegion=c(-500, 500), TxDb=hvo_TxDb, 
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                         addFlankGeneInfo = TRUE, verbose=FALSE)

gene.list <- as.data.frame(peakAnno)
write.csv(gene.list, "Significant_peaks_from3biorep.csv")  # Results


# Matching for a gene name replcement and a significant gene note
# you need to add column names as "Old_locus", "Uniprot" and "Protein_name" in the .csv file
raw.reads <- read.csv("Significant_peaks_from3biorep_noNaOCl_manualcheck.csv")
reference <- read.csv("DS2_uniprot_2.csv")

#indicate a significant gene by "Yes" in the raw.reads file
# a: gene name in the signif.gene file
# b: significant gene in the signif.gene file
# c: gene name in the raw.reads file
# d: Yes in the raw.reads file
#    d    <-                 b  [match( c , a)]

raw.reads$Old_locus <- reference$Old_locus[match(raw.reads$transcriptId, reference$Locus)]
raw.reads$Uniprot <- reference$Entry[match(raw.reads$Old_locus, reference$Gene_uniprot)]
raw.reads$Protein_name <- reference$Protein_name[match(raw.reads$Uniprot, reference$Entry)]

write.csv(raw.reads, "Significant_peaks_from3biorep_anno.csv") 







### BigWig file generation (normalized by WCE)
#------no oxidative stress
#   $ bamCompare -b1 ./0_data/H26_1_IP_S74_S80_sort.bam -b2 ./0_data/H26_1_WCE_S65_S71_sort.bam -o ./6_bigWig/H26_1_no_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/H26_1_no_bamCompare.log
#   $ bamCompare -b1 ./0_data/H26_2_IP_S91_S97_sort.bam -b2 ./0_data/H26_2_WCE_S95_S43_sort.bam -o ./6_bigWig/H26_2_no_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/H26_2_no_bamCompare.log
#   $ bamCompare -b1 ./0_data/oxsRHA_1_IP_S40_S46_sort.bam -b2 ./0_data/oxsRHA_1_WCE_S90_S96_sort.bam -o ./6_bigWig/oxsRHA_1_no_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/oxsRHA_1_no_bamCompare.log
#   $ bamCompare -b1 ./0_data/oxsRHA_2_IP_S57_S63_sort.bam -b2 ./0_data/oxsRHA_2_WCE_S70_S76_sort.bam -o ./6_bigWig/oxsRHA_2_no_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/oxsRHA_2_no_bamCompare.log
#   $ bamCompare -b1 ./0_data/oxsRHA_3_IP_S59_S65_sort.bam -b2 ./0_data/oxsRHA_3_WCE_S96_S44_sort.bam -o ./6_bigWig/oxsRHA_3_no_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/oxsRHA_3_no_bamCompare.log
#   $ bamCompare -b1 ./0_data/oxsRHA_4_IP_S82_S88_sort.bam -b2 ./0_data/oxsRHA_4_WCE_S80_S86_sort.bam -o ./6_bigWig/oxsRHA_4_no_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/oxsRHA_4_no_bamCompare.log

#------with NaOCl
#   $ bamCompare -b1 ./0_data/H26_1_NaOCl_IP_S81_S87_sort.bam -b2 ./0_data/H26_1_NaOCl_WCE_S86_S92_sort.bam -o ./6_bigWig/H26_1_NaOCl_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/H26_1_NaOCl_bamCompare.log
#   $ bamCompare -b1 ./0_data/H26_2_NaOCl_IP_S60_S66_sort.bam -b2 ./0_data/H26_2_NaOCl_WCE_S51_S57_sort.bam -o ./6_bigWig/H26_2_NaOCl_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/H26_2_NaOCl_bamCompare.log
#   $ bamCompare -b1 ./0_data/oxsRHA_1_NaOCl_IP_S83_S89_sort.bam -b2 ./0_data/oxsRHA_1_NaOCl_WCE_S49_S55_sort.bam -o ./6_bigWig/oxsRHA_1_NaOCl_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/oxsRHA_1_NaOCl_bamCompare.log
#   $ bamCompare -b1 ./0_data/oxsRHA_2_NaOCl_IP_S89_S95_sort.bam -b2 ./0_data/oxsRHA_2_NaOCl_WCE_S42_S48_sort.bam -o ./6_bigWig/oxsRHA_2_NaOCl_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/oxsRHA_2_NaOCl_bamCompare.log
#   $ bamCompare -b1 ./0_data/oxsRHA_3_NaOCl_IP_S77_S83_sort.bam -b2 ./0_data/oxsRHA_3_NaOCl_WCE_S39_S45_sort.bam -o ./6_bigWig/oxsRHA_3_NaOCl_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/oxsRHA_3_NaOCl_bamCompare.log
#   $ bamCompare -b1 ./0_data/oxsRHA_4_NaOCl_IP_S87_S93_sort.bam -b2 ./0_data/oxsRHA_4_NaOCl_WCE_S84_S90_sort.bam -o ./6_bigWig/oxsRHA_4_NaOCl_Norm.bw --binSize 20 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 60 --extendReads 150 --centerReads -p 6 2> ./6_bigWig/oxsRHA_4_NaOCl_bamCompare.log






### Identify the conserved motif by MEME or DREME (http://meme-suite.org)
# Extract the first three columns for 
#     $ cut -f 1,2,3 intersect_3biorep_narrowPeak.bed  > intersect_3biorep_narrowPeak-simple.bed

#     $ bedtools getfasta -fi HVO.fa -bed intersect_3biorep_narrowPeak-simple.bed -fo intersect_3biorep-simple-dreme.fasta







### A differential enrichment by DiffBind
# http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
# https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/08_diffbind_differential_peaks.html

library(DiffBind)
library(tidyverse)

samples <- read.csv("Meta_diffbind.csv")

dbObj <- dba(sampleSheet = samples)

#Take the alignment files and compute count information for each of the peaks/regions in the consensus set
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

#PCA plot and correlation heatmap using all the consensus sites
dba.plotPCA(dbObj,  attributes=DBA_CONDITION, label=DBA_ID)
plot(dbObj)

#Establishing a contrast to show which samples we want to compare to one another
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION)

#Performing the differential enrichment analysis
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
dba.show(dbObj, bContrasts=T)
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_CONDITION, label=DBA_ID)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

#MA plots - the effect of normalization on data, as well as seeing which of the data points are being identified as differentially bound.
dba.plotMA(dbObj, method=DBA_DESEQ2) #red points: the sites identified by DESeq2 as differentially bound (FDR < 0.05)
dba.plotMA(dbObj, bXY=TRUE)

#How the reads are distributed amongst the different classes of differentially bound sites and sample groups
pvals <- dba.plotBox(dbObj)

#Extract the full results from DESeq2
res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
out <- as.data.frame(res_deseq)
write.csv(out, "DiffBind_result.csv")
write.table(out, file="DiffBind_result.bed", sep="\t", quote=F, row.names=F, col.names=F)
# - Conc: mean read concentration over all the samples (the default calculation uses log2 normalized ChIP read counts with control read counts subtracted)
# - Conc_no_NaOCl: mean concentration over the first group
# - Fold: shows the difference in mean concentrations between the two groups, 
#   with a positive value indicating increased binding affinity in the first group
#   with a negative value indicating increased binding affinity in the second group.

DiffBind_sig_enrich <- out %>% 
  filter(FDR < 0.05 & abs(Fold) > 1)
write.csv(DiffBind_sig_enrich, "DiffBind_result_sig.csv")

# Create bed files for each keeping only significant peaks (FDR < 0.05, Fold >1)
DiffBind_sig_enrich <- out %>% 
  filter(FDR < 0.05 & abs(Fold) > 1) %>% 
  select(seqnames, start, end)
write.table(DiffBind_sig_enrich, file="DiffBind_result_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

### Followed by the annotation of data from the DiffBind
library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(GenomicRanges)
library(IRanges)

# Load data 
samplefiles <- list.files("./8_DiffBind/", pattern= ".bed", full.names=TRUE)
samplefiles <- as.list(samplefiles)

# Database construction for the index genome from gff file deposited in NCBI
library(GenomicFeatures)
hvo_TxDb <- makeTxDbFromGFF("GCF_000025685.1_ASM2568v1_genomic.gff", 
                            organism="Haloferax volcanii", format = "gff")

# Analyze data one-by-one -> Recommend to go below for the group analysis
peakAnnoList <- lapply(samplefiles[[2]], annotatePeak, TxDb=hvo_TxDb, 
                       tssRegion=c(-500, 500), verbose=FALSE)

peak <- readPeakFile(samplefiles[[2]])
peak

covplot(peak, weightCol="V6")

promoter <- getPromoters(TxDb=hvo_TxDb, upstream=500, downstream=500)
tagMatrix <- getTagMatrix(peak, windows=promoter)

peakHeatmap(samplefiles[[2]], TxDb=hvo_TxDb, upstream=500, downstream=500, color="red")

plotAvgProf2(samplefiles[[2]], TxDb=hvo_TxDb, upstream=500, downstream=500,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

peakAnno <- annotatePeak(samplefiles[[2]], tssRegion=c(-500, 500), TxDb=hvo_TxDb, 
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                         addFlankGeneInfo = TRUE, verbose=FALSE)

gene.list <- as.data.frame(peakAnno)
write.csv(gene.list, "DiffBind_peaks.csv")  # Results


# Matching for a gene name replcement and a significant gene note
rm(list = ls()) 

# you need to add column names as "Old_locus", "Uniprot" and "Protein_name" in the .csv file
raw.reads <- read.csv("DiffBind_peaks.csv")
reference <- read.csv("DS2_uniprot_2.csv")

#indicate a significant gene by "Yes" in the raw.reads file
# a: gene name in the signif.gene file
# b: significant gene in the signif.gene file
# c: gene name in the raw.reads file
# d: Yes in the raw.reads file
#    d    <-                 b  [match( c , a)]

raw.reads$Old_locus <- reference$Old_locus[match(raw.reads$transcriptId, reference$Locus)]
raw.reads$Uniprot <- reference$Entry[match(raw.reads$Old_locus, reference$Gene_uniprot)]
raw.reads$Protein_name <- reference$Protein_name[match(raw.reads$Uniprot, reference$Entry)]

write.csv(raw.reads, "DiffBind_peaks_anno.csv") 


#Input the peak location from the manual investigation into the DiffBind result
raw.reads <- read.csv("DiffBind_peaks_anno.csv")
reference1 <- read.csv("Significant_peaks_noNaOCl_manualcheck.csv")
reference2 <- read.csv("Significant_peaks_NaOCl_manualcheck.csv")

raw.reads$annotation_noNaOCl <- reference1$annotation[match(raw.reads$geneId, reference1$geneId)]
raw.reads$annotation_NaOCl <- reference2$annotation[match(raw.reads$geneId, reference2$geneId)]

write.csv(raw.reads, "DiffBind_peaks_anno_2.csv") 






####################################### Other information ####################################### 

### Peak calling by MOSAiCS R package
# Cynthia L. Darnell, Rylee K. Hackley, and Amy K. Schmid
# https://github.com/amyschmid/rosr-chip-utils/tree/master/RosR_ChIP-seq_pipeline
# https://bioconductor.org/packages/release/bioc/html/mosaics.html

# Package installation and libarary loading
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("mosaics")

library(mosaics)
library(hexbin)
library(tidyverse)

sample_file <- read_csv("Meta_mosaics_1.csv", col_names = F)

# For single-end ChIP-seq,,,
# Construct bins for IP, you can adjust fragLen and binSize depending on the ChIPQC later. 
IP_files <- unique(sample_file$X1)

for (i in 1:length(IP_files)){
  constructBins(infile=paste("data/", IP_files[i], sep=""),
                fileFormat="bam",
                outfileLoc="mosaics/bins/",
                byChr=FALSE,
                fragLen=200,
                binSize=200,
                capping=0,
                PET=FALSE)
}

# Construct bins for WCE
WCE_files <- sample_file$X2

for (i in 1:length(WCE_files)){
  constructBins(infile=paste("data/", WCE_files[i], sep=""),
                fileFormat="bam",
                outfileLoc="mosaics/bins/",
                byChr=FALSE,
                fragLen=200,
                binSize=200,
                capping=0,
                PET=FALSE)
}

for (i in 1:nrow(sample_file)) {
  sample_name <- paste("mosaics/bins/", sample_file[i,1], sep = "")
  sample_name <- str_replace(string = sample_name, pattern = ".bam", replacement = ".bam_fragL200_bin200.txt")
  ref_name <- paste("mosaics/bins/", sample_file[i,2], sep = "")
  ref_name <- str_replace(string = ref_name, pattern = ".bam", replacement = ".bam_fragL200_bin200.txt")
  
  print(paste("analyzing", sample_name, "against", ref_name))
  
  binTest <- readBins(type=c("chip", "input"), fileName= c(sample_name, ref_name))
  count_data <- hexbin (binTest@input, binTest@tagCount, xbins=100)
  control <- plot(count_data, trans=log, inv=exp, colramp=rainbow, xlab="WCE", ylab="ChIP", lcex=0.9)
  hexVP.abline(control$plot.vp, a=0, b=sum(binTest@tagCount)/sum(binTest@input), lwd=0.2)
  
  dev.copy(png, paste("mosaics/plots/", sample_file$X3[i], "_counts.png", sep=""))
  dev.off()
  
  fitTest <- mosaicsFit(binTest, analysisType="IO", bgEst="rMOM")
  plot(fitTest)
  
  dev.copy(png, paste("mosaics/plots/", sample_file$X3[i], "_fit.png", sep=""))
  dev.off()
  
  peakTest <- mosaicsPeak(fitTest, signalModel="2S", FDR=0.01)
  
  export(peakTest, type="bed", filename=paste("mosaics/", sample_file$X3[i], ".bed", sep=""))
  export(peakTest, type="txt", filename=paste("mosaics/", sample_file$X3[i],  ".txt", sep=""))
  
  sample_name_narrowP <- paste("data/", sample_file[i,1], sep = "")
  ref_name_narrowP <- paste("data/", sample_file[i,2], sep = "")
  
  peakNarrow <- extractReads(peakTest, 
                             chipFile=sample_name_narrowP, chipFileFormat="bam", chipPET=FALSE, chipFragLen=200, 
                             controlFile=ref_name_narrowP, controlFileFormat="bam", controlPET=FALSE, controlFragLen=200)
  
  peakNarrow <- findSummit(peakNarrow)
  
  export(peakNarrow, type="narrowPeak", filename=paste("mosaics/", sample_file$X3[i],  ".narrowPeak", sep=""))
}


#Initial nearby peaks are merged if the distance (in bp) between them is less than 'maxgap'. 
#Some initial peaks are removed if their lengths are shorter than 'minsize'.

#If you use a bin size shorter than the average fragment length in the experiment, 
#we recommend to set 'maxgap' to the average fragment length and 'minsize' to the bin size. 
#This setting removes peaks that are too narrow (e.g., singletons). 
#If you set the bin size to the average fragment length (or maybe bin size is larger than the average fragment length), 
#we recommend setting 'minsize' to a value smaller than the average fragment length while leaving 'maxgap' the same as the average fragment length. 
#This is to prevent filtering using 'minsize' because initial peaks would already be at a reasonable width. 


# For pair-end ChIP-seq,,, (You have to use "sam" files, or error meaasage is shown!)

sample_file <- read_csv("Meta_mosaics_2.csv", col_names = F)

IP_files <- unique(sample_file$X1)
for (i in 1:length(IP_files)){
  constructBins(infile=paste("data/", IP_files[i], sep=""),
                fileFormat="bam",
                outfileLoc="mosaics/bins/",
                byChr=FALSE,
                binSize=200,
                capping=0,
                PET=TRUE)
}

WCE_files <- sample_file$X2
for (i in 1:length(WCE_files)){
  constructBins(infile=paste("data/", WCE_files[i], sep=""),
                fileFormat="bam",
                outfileLoc="mosaics/bins/",
                byChr=FALSE,
                binSize=200,
                capping=0,
                PET=TRUE)
}

for (i in 1:nrow(sample_file)) {
  sample_name <- paste("mosaics/bins/", sample_file[i,1], sep = "")
  sample_name <- str_replace(string = sample_name, pattern = ".bam", replacement = ".sam_bin200.txt")
  ref_name <- paste("mosaics/bins/", sample_file[i,2], sep = "")
  ref_name <- str_replace(string = ref_name, pattern = ".bam", replacement = ".sam_bin200.txt")
  
  print(paste("analyzing", sample_name, "against", ref_name))
  
  binTest <- readBins(type=c("chip", "input"), fileName= c(sample_name, ref_name))
  count_data <- hexbin (binTest@input, binTest@tagCount, xbins=100)
  control <- plot(count_data, trans=log, inv=exp, colramp=rainbow, xlab="WCE", ylab="ChIP", lcex=0.9)
  hexVP.abline(control$plot.vp, a=0, b=sum(binTest@tagCount)/sum(binTest@input), lwd=0.2)
  
  dev.copy(png, paste("mosaics/plots/", sample_file$X3[i], "_counts.png", sep=""))
  dev.off()
  
  fitTest <- mosaicsFit(binTest, analysisType="IO", bgEst="rMOM")
  plot(fitTest)
  
  dev.copy(png, paste("mosaics/plots/", sample_file$X3[i], "_fit.png", sep=""))
  dev.off()
  
  peakTest <- mosaicsPeak(fitTest, signalModel="2S", FDR=0.01)
  
  export(peakTest, type="bed", filename=paste("mosaics/", sample_file$X3[i], ".bed", sep=""))
  export(peakTest, type="txt", filename=paste("mosaics/", sample_file$X3[i],  ".txt", sep=""))
  #export(peakNarrow, type="narrowPeak", filename=paste("mosaics/", sample_file$X3[i],  ".narrowPeak", sep=""))
  
  #sample_name_narrowP <- paste("data/", sample_file[i,1], sep = "")
  #ref_name_narrowP <- paste("data/", sample_file[i,2], sep = "")
  
  #peakNarrow <- extractReads(peakTest, 
  #                           chipFile=sample_name_narrowP, chipFileFormat="sam", chipPET=TRUE,
  #                           controlFile=ref_name_narrowP, controlFileFormat="sam", controlPET=TRUE)
  
  #peakNarrow <- findSummit(peakNarrow)
  
  #export(peakNarrow, type="narrowPeak", filename=paste("mosaics/", sample_file$X3[i],  ".narrowPeak", sep=""))
}



