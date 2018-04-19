
biol_reps_data<-read.table("C:\\Users\\Shamita\\Desktop\\PSU\\biol_reps.txt",header = T)
dim(biol_reps_data)
## [1] 17564    21 

biol_reps_data[1:5,1:5]
##             CDIPT_2 CDIPT_3 CDIPT_4 KCTD13_1 KCTD13_2
## FBgn0000003       0       1       0        0        0
## FBgn0000008    1739    2176    1445     1628     1741
## FBgn0000014       3       3       0        2        2
## FBgn0000015       1       0       0        1        1
## FBgn0000017    9976    8824    6422    11300    12067


logbiolrepsdata=log2(biol_reps_data+.25)
par(mfrow=c(3,3))
for (i in 1:21) hist(logbiolrepsdata[,i],main=colnames(biol_reps_data)[i])

condition= factor(paste(c(rep("CCDIPT",3),rep("KCTD13",3),rep("MAPK3",3),rep("CORO1A",3),rep("C16orf53",3),rep("DOC2A",3),rep("WT",3))))

sample<-c(1:21)
sample.data.sep<-data.frame(sample,condition)

#Convert to matrix and convert characters to numeric type/round to integers

biol_reps_numeric<-as.matrix(biol_reps_data)

class(biol_reps_numeric)<-"numeric"

biol_reps_numeric<-round(biol_reps_numeric)

require(DESeq2)

## Loading required package: DESeq2
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, xtabs
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     grep, grepl, intersect, is.unsorted, lapply, lengths, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unlist, unsplit
## Loading required package: IRanges
## Loading required package: GenomicRanges
## Loading required package: GenomeInfoDb
## Loading required package: SummarizedExperiment
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## Loading required package: Rcpp
## Warning: package 'Rcpp' was built under R version 3.2.4
## Loading required package: RcppArmadillo
## Warning: package 'RcppArmadillo' was built under R version 3.2.4

dds_sep<-DESeqDataSetFromMatrix(countData=biol_reps_numeric, colData=sample.data.sep, design= ~ condition)

## converting counts to integer mode
#Pre-filter to remove rows with <2 total counts
nrow(dds_sep)
## [1] 17564

dds<-dds_sep[rowSums(counts(dds_sep))>1,]
nrow(dds)
## [1] 15526

rld <- rlog(dds)
head(assay(rld), 3)

##                       1           2           3          4           5
## FBgn0000003 -1.28225052 -1.25451921 -1.28041958 -1.2820187 -1.28230945
## FBgn0000008 10.48612503 10.92848804 10.80617384 10.4889437 10.46966470
## FBgn0000014  0.04604333  0.05830042 -0.01788328  0.0265776  0.02308518
##                        6          7           8            9          10
## FBgn0000003 -1.281871409 -1.2816845 -1.23132141 -1.254496763 -1.28169183
## FBgn0000008 10.637145865 10.4405693 10.40379646 10.354309725 10.35974571
## FBgn0000014  0.003791175  0.0554062 -0.02118674  0.006456863  0.03050138
##                      11          12           13          14         15
## FBgn0000003 -1.28174391 -1.25667947 -1.234331197 -1.28164771 -1.2816318
## FBgn0000008 10.20684033 10.26915141 10.242151140 10.30698594 10.3732368
## FBgn0000014 -0.02104935 -0.02135919  0.002247937 -0.02081928 -0.0207813
##                      16          17           18           19          20
## FBgn0000003 -1.28182949 -1.28154941 -1.256541679 -1.255509254 -1.28123944
## FBgn0000008 10.26733547 10.29305552 10.311360060 10.305537076 10.13166720
## FBgn0000014 -0.02125402  0.03221021  0.003945482  0.005213501 -0.01984304
##                      21
## FBgn0000003 -1.28128016
## FBgn0000008 10.27462316
## FBgn0000014 -0.01994039

par( mfrow = c( 1, 2 ) )
ddse <- estimateSizeFactors(dds)
plot(log2(counts(ddse, normalized=TRUE)[,1:2] + 1), main="log2 transform of normalized counts", pch=16, cex=0.3)
plot(assay(rld)[,1:2], pch=16, main="rlog", cex=0.3)

sampleDists <- dist( t( assay(rld) ) )
sampleDists

##           1        2        3        4        5        6        7        8
## 2  27.61928                                                               
## 3  29.71277 14.13129                                                      
## 4  38.71147 44.35692 43.47135                                             
## 5  36.01839 43.23140 42.16990 12.35890                                    
## 6  39.57496 46.06544 44.90782 15.51598 12.35682                           
## 7  39.13591 44.19701 43.73252 29.32311 28.18615 29.33120                  
## 8  39.56841 45.31146 44.52298 25.57149 25.05698 27.11296 17.23592         
## 9  38.49753 43.30060 42.14930 31.09093 30.21101 31.78310 18.22619 17.51730
## 10 40.29812 43.77126 43.90130 38.06318 36.36784 38.61217 37.75845 37.14569
## 11 43.82767 47.71844 46.71589 21.75449 24.27907 25.85080 29.45420 24.29681
## 12 42.88324 47.37239 46.35366 22.68804 24.12074 25.54474 28.00613 22.28100
## 13 37.05900 44.47096 43.39274 22.08711 21.51948 24.44932 26.21171 21.47871
## 14 37.59886 43.40854 42.19397 26.01151 23.73220 24.73259 22.82306 24.10289
## 15 40.61163 47.39895 46.22397 30.21816 27.25807 26.27493 24.62963 26.67354
## 16 34.46304 42.72667 41.75042 25.79098 22.96202 24.84519 23.72619 23.39078
## 17 37.96813 46.09766 44.98880 24.18575 22.68200 24.02584 24.58418 22.16374
## 18 35.68775 45.23060 44.08956 27.29327 24.20856 26.10042 25.89722 24.99480
## 19 39.64298 48.92177 47.44405 30.31105 26.68750 27.09675 27.74845 27.07709
## 20 37.56177 46.74597 45.47031 33.67993 30.23726 32.27017 28.29235 27.75797
## 21 35.79480 44.13561 43.37537 30.52379 28.50698 30.38235 25.69463 26.33825
##           9       10       11       12       13       14       15       16
                                                                       
## 10 37.65367                                                               
## 11 29.48071 37.21525                                                      
## 12 26.17893 36.79145 13.44307                                             
## 13 26.32025 34.32230 19.95640 19.28156                                    
## 14 26.10195 33.43626 25.20266 24.51575 17.50477                           
## 15 28.67607 36.77174 31.37915 29.41554 23.33906 17.43432                  
## 16 25.85212 33.31942 24.92464 24.55664 18.22633 17.40139 20.52888         
## 17 25.43985 34.89178 23.50874 21.74588 16.45450 18.50893 20.57653 13.75189
## 18 27.24482 34.19906 28.23018 26.84328 18.10737 18.35586 21.08028 14.39219
## 19 29.62240 37.87361 31.61721 29.17606 23.34753 21.81618 21.12387 19.51533
## 20 27.82479 38.90865 33.87441 30.69253 24.93756 23.66498 23.61043 20.11916
## 21 26.93473 36.29862 31.65360 29.19062 22.60023 20.69373 21.26757 18.37784
##          17       18       19       20
                                  
## 18 13.10717                           
## 19 18.03295 17.21076                  
## 20 19.00424 18.84864 18.32563         
## 21 17.12967 17.64551 19.25563 14.74276

library("pheatmap")
## Warning: package 'pheatmap' was built under R version 3.2.4
library("RColorBrewer")

sampleDistMatrix <- as.matrix(sampleDists)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

dds_analysis<-DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
ddss <- estimateSizeFactors(dds)
#Normalize samples

par(mfrow=c(1,1))
rldb <- rlogTransformation(ddss, blind = TRUE)
print(plotPCA(rld, intgroup = c("condition")))

#Sample QC: Transform data to log space and visualize samples
ddss <- estimateDispersions(ddss)
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
plotDispEsts(ddss)


#Estimate biological variance and visualize
ddss <- nbinomWaldTest(ddss)
#Determine differential expression

# Get results and write to output files
res1 <- results(ddss)
res1 <- res1[order(res1$padj), ]
write.table(cbind('#ID'=rownames(as.data.frame(res1)), as.data.frame(res1)),
file='deseq2_results.txt', quote=FALSE, sep='\t', row.names=FALSE) 

# Display the top few genes
head(res1)
## log2 fold change (MAP): condition WT vs C16orf53 
## Wald test p-value: condition WT vs C16orf53 
## DataFrame with 6 rows and 6 columns
##                baseMean log2FoldChange      lfcSE       stat       pvalue
##               <numeric>      <numeric>  <numeric>  <numeric>    <numeric>
## FBgn0262983   113.36936       2.730910 0.20243457  13.490334 1.782917e-41
## FBgn0085224    93.49030      -2.432127 0.18134799 -13.411378 5.186595e-41
## FBgn0032549    53.59738      -3.113570 0.25960942 -11.993287 3.853058e-33
## FBgn0003996 14843.21813      -1.502959 0.12798596 -11.743157 7.657355e-32
## FBgn0050295   576.42352      -1.094101 0.09717824 -11.258702 2.098275e-29
## FBgn0037672   197.01456       1.632855 0.16533191   9.876225 5.278421e-23
##                     padj
##                <numeric>
## FBgn0262983 1.894884e-37
## FBgn0085224 2.756157e-37
## FBgn0032549 1.365010e-29
## FBgn0003996 2.034559e-28
## FBgn0050295 4.460093e-26
## FBgn0037672 9.349842e-20
#Visualize with MA plot
plotMA(ddss,ylim=c(-2,2),main="DESeq2")

res <- results(dds_analysis)

results<-results(dds_analysis, alpha=0.05) 

#P-value is set to 0.05, default in DESeq is 0.1
table(results$pv < 0.05)
## 
## FALSE  TRUE 
## 13066  2300

table(results$padj < 0.05)
## 
## FALSE  TRUE 
## 10003   921

resultnames= resultsNames(dds_analysis)
resultnew= removeResults(dds_analysis)

summary(results)
## 
## out of 15526 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)     : 434, 2.8% 
## LFC < 0 (down)   : 487, 3.1% 
## outliers [1]     : 160, 1% 
## low counts [2]   : 4442, 29% 
## (mean count < 6)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
#results(dds_analysis, contrast=c("condition","WT"))

resLFC1 <- results(dds_analysis, lfcThreshold=1)
table(resLFC1$padj < 0.1)
## 
## FALSE  TRUE 
## 15362     4
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")


plotMA(results, main="Result of Differential expression analysis")


#MA Plot comparing mean expression to log-fold change (quality check)

plotMA(resLFC1, ylim=c(-5,5))
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

require('genefilter')
## Loading required package: genefilter
## 
## Attaching package: 'genefilter'
## The following object is masked from 'package:base':
## 
##     anyNA
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)

mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld))
pheatmap(mat, annotation_col=df)


#resOrdered <- res[order(res$padj),]
#head(resOrdered)
#Get an ordered report of results
#library("ReportingTools")
#resOrderedDF <- as.data.frame(resOrdered)[1:100,]
#write.csv(resOrderedDF, file="results.csv")

#htmlRep <- HTMLReport(shortName="report", title="My report", reportDirectory="./report")
#publish(resOrderedDF, htmlRep)
#url <- finish(htmlRep)
#browseURL(url) - to open the html file


#Export results to text files
#write.csv(as.data.frame(res), file="results_deseq2.csv")

sig=subset(results, results$padj<0.05)
siggenes=rownames(sig)
#write.csv(siggenes, file="significant gene names.csv")
#Export FlyBase IDs of significant genes for analysis in DAVID 

#Export normalized counts to files
#write.csv(counts(dds_sep, normalized=TRUE), file="counts_norm_deseq2.csv")
library("AnnotationDbi")
library('org.Dm.eg.db')
## Loading required package: DBI
## 
#Genome wide annotation for Fly

res$symbol <- mapIds(org.Dm.eg.db, keys=row.names(res), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

res$entrez <- mapIds(org.Dm.eg.db, keys=row.names(res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")

cols <- c("SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID")
               
resannotated= select(org.Dm.eg.db, keys=row.names(res), columns=cols, keytype="ENSEMBL")
## 'select()' returned 1:many mapping between keys and columns

##                      ENSEMBL                        SYMBOL
## 1                FBgn0000003                          <NA>
## 2                FBgn0000008                             a
## 3                FBgn0000014                         abd-A
## 4                FBgn0000015                         Abd-B
## 5                FBgn0000017                           Abl
## 6                FBgn0000018                           abo
## 7                FBgn0000024                           Ace
## 8                FBgn0000028                          acj6
## 9                FBgn0000032                        Acph-1
## 10               FBgn0000036                   nAChRalpha1



write.csv(resannotated, file="annotated genes from result.csv")
GOanalysisBP= read.csv("C:\\Users\\Shamita\\Desktop\\PSU\\spring 2016\\stat 555\\DEseq results\\GO chart_5DD20948D1C21462043470521.csv", header = T)
head(GOanalysisBP)
##        Category                                    Term Count        X.   PValue
## 1 GOTERM_BP_FAT   GO:0000022~mitotic spindle elongation    39  4.659498   7.02e-25
## 2 GOTERM_BP_FAT           GO:0051231~spindle elongation    39  4.659498   1.27e-24
## 3 GOTERM_BP_FAT                  GO:0006412~translation    87 10.394265   2.69e-16
## 4 GOTERM_BP_FAT GO:0007052~mitotic spindle organization    41  4.898447   4.87e-11
## 5 GOTERM_BP_FAT         GO:0007051~spindle organization    44  5.256870   1.20e-10
## 6 GOTERM_BP_FAT    GO:0007010~cytoskeleton organization    68  8.124253   4.18e-10
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            FBGN0033699, FBGN0001942, FBGN0032518, FBGN0031980, FBGN0034968, FBGN0032987, FBGN0029897, FBGN0017579, FBGN0003517, FBGN0011272, FBGN0034138, FBGN0027948, FBGN0039857, FBGN0036213, FBGN0034743, FBGN0010409, FBGN0038834, FBGN0002590, FBGN0010408, FBGN0014026, FBGN0029785, FBGN0002626, FBGN0002607, FBGN0086710, FBGN0036825, FBGN0011284, FBGN0035422, FBGN0010078, FBGN0010265, FBGN0020910, FBGN0005593, FBGN0035753, FBGN0026250, FBGN0017545, FBGN0025286, FBGN0010411, FBGN0039359, FBGN0015756, FBGN0010198
## 3 FBGN0032518, FBGN0001942, FBGN0033699, FBGN0037686, FBGN0032987, FBGN0029897, FBGN0086472, FBGN0017579, FBGN0010803, FBGN0000100, FBGN0011272, FBGN0031869, FBGN0004867, FBGN0034654, FBGN0003279, FBGN0039300, FBGN0036213, FBGN0034743, FBGN0039406, FBGN0002590, FBGN0037351, FBGN0038834, FBGN0010409, FBGN0010408, FBGN0003274, FBGN0025582, FBGN0014026, FBGN0002593, FBGN0015834, FBGN0002622, FBGN0029785, FBGN0002626, FBGN0033480, FBGN0086710, FBGN0028697, FBGN0064225, FBGN0010078, FBGN0005593, FBGN0035753, FBGN0003941, FBGN0003942, FBGN0025286, FBGN0030616, FBGN0023519, FBGN0019936, FBGN0001995, FBGN0039359, FBGN0260441, FBGN0031980, FBGN0036135, FBGN0034968, FBGN0034967, FBGN0002579, FBGN0037135, FBGN0037899, FBGN0003517, FBGN0034138, FBGN0039857, FBGN0039757, FBGN0003149, FBGN0026372, FBGN0039713, FBGN0033912, FBGN0002607, FBGN0005533, FBGN0011284, FBGN0036825, FBGN0035422, FBGN0020910, FBGN0010265, FBGN0034915, FBGN0030136, FBGN0026250, FBGN0017545, FBGN0042712, FBGN0037328, FBGN0028737, FBGN0034087, FBGN0010412, FBGN0035423, FBGN0030802, FBGN0010411, FBGN0016726, FBGN0004404, FBGN0015756, FBGN0004403, FBGN0010198, FBGN0015521
                                                                                                                                                                                                                                                                     FBGN0033699, FBGN0001942, FBGN0032518, FBGN0040235, FBGN0000044, FBGN0029897, FBGN0000042, FBGN0032987, FBGN0038294, FBGN0017579, FBGN0011272, FBGN0027948, FBGN0004009, FBGN0036213, FBGN0034743, FBGN0002590, FBGN0038834, FBGN0010409, FBGN0010408, FBGN0000273, FBGN0014026, FBGN0029785, FBGN0053556, FBGN0002626, FBGN0014020, FBGN0086710, FBGN0259108, FBGN0010078, FBGN0005593, FBGN0035753, FBGN0000723, FBGN0086906, FBGN0025286, FBGN0039359, FBGN0031980, FBGN0034968, FBGN0020440, FBGN0002736, FBGN0023081, FBGN0014141, FBGN0027492, FBGN0003517, FBGN0034138, FBGN0039857, FBGN0000578, FBGN0000463, FBGN0086359, FBGN0003149, FBGN0032409, FBGN0000152, FBGN0002607, FBGN0011284, FBGN0036825, FBGN0035422, FBGN0010265, FBGN0020910, FBGN0003447, FBGN0010314, FBGN0013733, FBGN0026250, FBGN0017545, FBGN0010411, FBGN0085447, FBGN0000163, FBGN0001404, FBGN0015756, FBGN0000547, FBGN0010198
##   List.Total Pop.Hits Pop.Total Fold.Enrichment Bonferroni Benjamini  FDR
## 1        523       78      7937        7.587954   1.28e-21  1.28e-21   1.19e-21
## 2        523       79      7937        7.491904   2.31e-21  1.16e-21   2.16e-21
## 3        523      520      7937        2.539046   4.05e-13  1.35e-13   3.77e-13
## 4        523      194      7937        3.207280   8.88e-08  2.22e-08   8.28e-08
## 5        523      225      7937        2.967733   2.19e-07  4.38e-08   2.04e-07
## 6        523      465      7937        2.219273   7.62e-07  1.27e-07   7.12e-07


#GO analysis of all significant genes showing biological process GO

neuroGOanalysisBP= read.csv("C:\\Users\\Shamita\\Desktop\\PSU\\DEseq results\\neurodevelopmental subset chart_6706354EB0811462045395343.csv", header = T)
neuroGOanalysisBP[1:30,]
##         Category
## 1  GOTERM_BP_FAT
## 2  GOTERM_BP_FAT
## 3  GOTERM_BP_FAT
## 4  GOTERM_BP_FAT
## 5  GOTERM_BP_FAT
## 6  GOTERM_BP_FAT
## 7  GOTERM_BP_FAT
## 8  GOTERM_BP_FAT
## 9  GOTERM_BP_FAT
## 10 GOTERM_BP_FAT
## 11 GOTERM_BP_FAT
## 12 GOTERM_BP_FAT
## 13 GOTERM_BP_FAT
## 14 GOTERM_BP_FAT
## 15 GOTERM_BP_FAT
## 16 GOTERM_BP_FAT
## 17 GOTERM_BP_FAT
## 18 GOTERM_BP_FAT
## 19 GOTERM_BP_FAT
## 20 GOTERM_BP_FAT
## 21 GOTERM_BP_FAT
## 22 GOTERM_BP_FAT
## 23 GOTERM_BP_FAT
## 24 GOTERM_BP_FAT
## 25 GOTERM_BP_FAT
## 26 GOTERM_BP_FAT
## 27 GOTERM_BP_FAT
## 28 GOTERM_BP_FAT
## 29 GOTERM_BP_FAT
## 30 GOTERM_BP_FAT
##                                                                Term Count
## 1                                 GO:0030182~neuron differentiation    36
## 2                                     GO:0048666~neuron development    31
## 3                        GO:0048812~neuron projection morphogenesis    27
## 4                          GO:0031175~neuron projection development    27
## 5                           GO:0030030~cell projection organization    29
## 6                       GO:0032989~cellular component morphogenesis    32
## 7         GO:0000904~cell morphogenesis involved in differentiation    27
## 8                                            GO:0006928~cell motion    27
## 9                                     GO:0000902~cell morphogenesis    30
## 10                         GO:0048858~cell projection morphogenesis    27
## 11 GO:0048667~cell morphogenesis involved in neuron differentiation    26
## 12                               GO:0032990~cell part morphogenesis    27
## 13                                          GO:0007409~axonogenesis    22
## 14                                         GO:0007411~axon guidance    19
## 15                        GO:0060284~regulation of cell development    13
## 16                             GO:0007423~sensory organ development    19
## 17              GO:0051960~regulation of nervous system development    11
## 18                                        GO:0016477~cell migration    13
## 19                                         GO:0048870~cell motility    13
## 20                                  GO:0051674~localization of cell    13
## 21                                          GO:0035282~segmentation    14
## 22                                           GO:0042063~gliogenesis     8
## 23                         GO:0007389~pattern specification process    17
## 24                                  GO:0016358~dendrite development    10
## 25                                GO:0048813~dendrite morphogenesis    10
## 26                    GO:0046530~photoreceptor cell differentiation    10
## 27                                       GO:0003002~regionalization    16
## 28                            GO:0050767~regulation of neurogenesis     8
## 29                          GO:0030029~actin filament-based process    10
## 30                              GO:0048749~compound eye development    13
##          X.   PValue
## 1  65.45455 1.47e-32
## 2  56.36364 2.53e-27
## 3  49.09091 7.08e-24
## 4  49.09091 7.75e-24
## 5  52.72727 7.84e-24
## 6  58.18182 1.92e-23
## 7  49.09091 3.19e-23
## 8  49.09091 3.79e-23
## 9  54.54545 7.67e-23
## 10 49.09091 1.32e-22
## 11 47.27273 2.21e-22
## 12 49.09091 3.16e-22
## 13 40.00000 1.83e-20
## 14 34.54545 2.87e-19
## 15 23.63636 4.68e-11
## 16 34.54545 8.24e-11
## 17 20.00000 3.36e-10
## 18 23.63636 1.39e-09
## 19 23.63636 2.97e-09
## 20 23.63636 5.08e-09
## 21 25.45455 1.24e-08
## 22 14.54545 2.68e-08
## 23 30.90909 6.20e-08
## 24 18.18182 7.35e-08
## 25 18.18182 7.35e-08
## 26 18.18182 1.60e-07
## 27 29.09091 2.00e-07
## 28 14.54545 2.87e-07
## 29 18.18182 4.40e-07
## 30 23.63636 7.46e-07
##                                                                                                                                                                                                                                                                                                                                 Genes
## 1  1990532, 1992610, 1977118, 1979999, 1974045, 1984317, 1991109, 1994460, 1980150, 1985146, 1991973, 1985129, 1977625, 1973359, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1989410, 1979036, 1983532, 1978563, 1974058, 1975492, 1980536, 1988685, 1984964, 1988909, 1974221, 1974346, 1984480, 1991416, 1984865, 1994597
## 2                                               1990532, 1979999, 1974045, 1984317, 1991109, 1994460, 1980150, 1985146, 1991973, 1985129, 1977625, 1973359, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1979036, 1983532, 1978563, 1974058, 1975492, 1988685, 1980536, 1984964, 1988909, 1974221, 1974346, 1991416, 1984865
## 3                                                                                   1990532, 1979999, 1991109, 1984317, 1994460, 1980150, 1991973, 1973359, 1977625, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1979036, 1983532, 1978563, 1974058, 1975492, 1988685, 1980536, 1984964, 1988909, 1974346, 1991416, 1984865
## 4                                                                                   1990532, 1979999, 1991109, 1984317, 1994460, 1980150, 1991973, 1973359, 1977625, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1979036, 1983532, 1978563, 1974058, 1975492, 1988685, 1980536, 1984964, 1988909, 1974346, 1991416, 1984865
## 5                                                                 1990532, 1979999, 1988551, 1984317, 1991109, 1994460, 1980150, 1991973, 1973359, 1977625, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1979036, 1983532, 1978563, 1974058, 1975492, 1988685, 1980536, 1984964, 1988909, 1974346, 1981813, 1991416, 1984865
## 6                                      1990532, 1992610, 1979999, 1984317, 1991109, 1994460, 1980150, 1991973, 1977625, 1973359, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1990016, 1979036, 1983532, 1978563, 1974058, 1975492, 1974209, 1988685, 1980536, 1984964, 1988909, 1974346, 1984480, 1981813, 1991416, 1984865
## 7                                                                                   1992610, 1979999, 1984317, 1991109, 1994460, 1980150, 1991973, 1973359, 1977625, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1979036, 1983532, 1978563, 1974058, 1975492, 1988685, 1980536, 1984964, 1988909, 1974346, 1991416, 1984865
## 8                                                                                   1990532, 1992610, 1979999, 1984317, 1991109, 1994460, 1980150, 1991973, 1973359, 1977625, 1985876, 1990629, 1980219, 1979036, 1983532, 1975168, 1978563, 1975492, 1988685, 1980536, 1984964, 1988909, 1984480, 1974346, 1991416, 1984865, 1994597
## 9                                                        1990532, 1992610, 1979999, 1984317, 1991109, 1994460, 1980150, 1991973, 1973359, 1977625, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1979036, 1983532, 1978563, 1974058, 1975492, 1974209, 1988685, 1980536, 1984964, 1988909, 1974346, 1981813, 1991416, 1984865
## 10                                                                                  1990532, 1979999, 1991109, 1984317, 1994460, 1980150, 1991973, 1973359, 1977625, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1979036, 1983532, 1978563, 1974058, 1975492, 1988685, 1980536, 1984964, 1988909, 1974346, 1991416, 1984865
## 11                                                                                           1979999, 1991109, 1984317, 1994460, 1980150, 1991973, 1973359, 1977625, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1979036, 1983532, 1978563, 1974058, 1975492, 1988685, 1980536, 1984964, 1988909, 1974346, 1991416, 1984865
## 12                                                                                  1990532, 1979999, 1991109, 1984317, 1994460, 1980150, 1991973, 1973359, 1977625, 1985876, 1984423, 1990629, 1993294, 1981250, 1980219, 1979036, 1983532, 1978563, 1974058, 1975492, 1988685, 1980536, 1984964, 1988909, 1974346, 1991416, 1984865
## 13                                                                                                                               1979999, 1979036, 1983532, 1978563, 1974058, 1984317, 1991109, 1975492, 1994460, 1980150, 1988685, 1980536, 1984964, 1988909, 1977625, 1974346, 1985876, 1984423, 1991416, 1984865, 1990629, 1980219
## 14                                                                                                                                                          1979999, 1979036, 1983532, 1984317, 1991109, 1975492, 1994460, 1980150, 1988685, 1980536, 1984964, 1988909, 1977625, 1974346, 1985876, 1991416, 1990629, 1984865, 1980219
## 15                                                                                                                                                                                                                1992610, 1974863, 1979999, 1986373, 1975168, 1974058, 1977286, 1984964, 1973655, 1991973, 1981813, 1991416, 1990016
## 16                                                                                                                                                          1990532, 1992610, 1977118, 1979999, 1978563, 1974058, 1974045, 1988551, 1987321, 1984964, 1973655, 1974221, 1985146, 1977625, 1984480, 1981813, 1994597, 1989410, 1990016
## 17                                                                                                                                                                                                                                  1984964, 1992610, 1991973, 1974863, 1979999, 1986373, 1975168, 1974058, 1977286, 1991416, 1990016
## 18                                                                                                                                                                                                                1984964, 1990532, 1988909, 1992610, 1991973, 1973359, 1977625, 1984480, 1979036, 1975168, 1978563, 1991109, 1994597
## 19                                                                                                                                                                                                                1984964, 1990532, 1988909, 1992610, 1991973, 1973359, 1977625, 1984480, 1979036, 1975168, 1978563, 1991109, 1994597
## 20                                                                                                                                                                                                                1984964, 1990532, 1988909, 1992610, 1991973, 1973359, 1977625, 1984480, 1979036, 1975168, 1978563, 1991109, 1994597
## 21                                                                                                                                                                                                       1992610, 1994875, 1974863, 1977588, 1978563, 1977286, 1987321, 1973655, 1991973, 1988402, 1977625, 1973359, 1984480, 1981813
## 22                                                                                                                                                                                                                                                             1984964, 1992610, 1985058, 1977625, 1978563, 1979428, 1994597, 1974209
## 23                                                                                                                                                                            1990532, 1992610, 1994875, 1974863, 1977588, 1978563, 1977286, 1987321, 1984964, 1973655, 1991973, 1988402, 1977625, 1973359, 1984480, 1985876, 1981813
## 24                                                                                                                                                                                                                                           1984964, 1991973, 1973359, 1979999, 1979036, 1974058, 1975492, 1993294, 1981250, 1980536
## 25                                                                                                                                                                                                                                           1984964, 1991973, 1973359, 1979999, 1979036, 1974058, 1975492, 1993294, 1981250, 1980536
## 26                                                                                                                                                                                                                                           1974221, 1985146, 1992610, 1977118, 1984480, 1974346, 1978563, 1974045, 1994597, 1989410
## 27                                                                                                                                                                                     1990532, 1992610, 1994875, 1974863, 1977588, 1978563, 1977286, 1987321, 1973655, 1991973, 1988402, 1977625, 1973359, 1984480, 1985876, 1981813
## 28                                                                                                                                                                                                                                                             1984964, 1992610, 1974863, 1979999, 1986373, 1975168, 1977286, 1990016
## 29                                                                                                                                                                                                                                           1984964, 1990532, 1992610, 1979999, 1974346, 1979036, 1986373, 1988551, 1981813, 1990016
## 30                                                                                                                                                                                                                1984964, 1973655, 1974221, 1985146, 1992610, 1977118, 1984480, 1978563, 1974058, 1974045, 1981813, 1994597, 1989410
##    List.Total Pop.Hits Pop.Total Fold.Enrichment Bonferroni Benjamini
## 1          55      409      7937       12.702023   1.36e-29  1.36e-29
## 2          55      347      7937       12.892167   2.34e-24  1.17e-24
## 3          55      286      7937       13.623586   6.54e-21  2.18e-21
## 4          55      287      7937       13.576117   7.16e-21  1.79e-21
## 5          55      365      7937       11.465654   7.24e-21  1.45e-21
## 6          55      518      7937        8.914847   1.77e-20  2.95e-21
## 7          55      303      7937       12.859226   2.95e-20  4.21e-21
## 8          55      305      7937       12.774903   3.50e-20  4.38e-21
## 9          55      442      7937        9.794735   7.08e-20  7.87e-21
## 10         55      320      7937       12.176080   1.22e-19  1.22e-20
## 11         55      288      7937       13.027904   2.04e-19  1.85e-20
## 12         55      331      7937       11.771436   2.92e-19  2.43e-20
## 13         55      198      7937       16.034343   1.69e-17  1.30e-18
## 14         55      136      7937       20.160829   2.66e-16  1.90e-17
## 15         55      131      7937       14.320749   4.32e-08  2.88e-09
## 16         55      410      7937        6.687494   7.62e-08  4.76e-09
## 17         55       90      7937       17.637778   3.10e-07  1.82e-08
## 18         55      175      7937       10.720104   1.28e-06  7.12e-08
## 19         55      187      7937       10.032183   2.75e-06  1.45e-07
## 20         55      196      7937        9.571521   4.69e-06  2.35e-07
## 21         55      260      7937        7.770490   1.14e-05  5.45e-07
## 22         55       48      7937       24.051515   2.48e-05  1.13e-06
## 23         55      480      7937        5.110947   5.73e-05  2.49e-06
## 24         55      117      7937       12.334110   6.79e-05  2.83e-06
## 25         55      117      7937       12.334110   6.79e-05  2.83e-06
## 26         55      128      7937       11.274148   1.48e-04  5.93e-06
## 27         55      454      7937        5.085783   1.85e-04  7.10e-06
## 28         55       67      7937       17.230936   2.65e-04  9.82e-06
## 29         55      144      7937       10.021465   4.07e-04  1.45e-05
## 30         55      308      7937        6.090968   6.89e-04  2.38e-05
##             FDR
## 1  2.300000e-29
## 2  3.960000e-24
## 3  1.110000e-20
## 4  1.210000e-20
## 5  1.230000e-20
## 6  3.000000e-20
## 7  4.990000e-20
## 8  5.930000e-20
## 9  1.200000e-19
## 10 2.060000e-19
## 11 3.450000e-19
## 12 4.940000e-19
## 13 2.860000e-17
## 14 4.490000e-16
## 15 7.320000e-08
## 16 1.290000e-07
## 17 5.250000e-07
## 18 2.170000e-06
## 19 4.650000e-06
## 20 7.950000e-06
## 21 1.940000e-05
## 22 4.200000e-05
## 23 9.690000e-05
## 24 1.150000e-04
## 25 1.150000e-04
## 26 2.510000e-04
## 27 3.130000e-04
## 28 4.490000e-04
## 29 6.890000e-04
## 30 1.167229e-03

#subset of above document showing genes essential for neurodevelopmental function

functional= read.csv("C://Users//Shamita//Desktop//PSU//DEseq results//PIR rt_6706354EB0811462045902018.csv")
head(functional)
##   X.        Category                             Term     Kappa
## 1  1 SP_PIR_KEYWORDS                     neurogenesis 1.0000000
## 2  2 SP_PIR_KEYWORDS                  differentiation 0.6994536
## 3  3  UP_SEQ_FEATURE topological domain:Extracellular 0.6496815
## 4  4  UP_SEQ_FEATURE   topological domain:Cytoplasmic 0.6496815
## 5  5  UP_SEQ_FEATURE                      domain:Sema 0.6496815
## 6  6  UP_SEQ_FEATURE                   disulfide bond 0.6373626
#Functional categories information

library(png)
img = readPNG("C:\\Users\\Shamita\\Desktop\\PSU\\DEseq results\\KEGG dme04310.png")
par(mfrow=c(1,1))
grid::grid.raster(img)


#KEGG pathway results for subset
sessionInfo()
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 10586)
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] png_0.1-7                  org.Dm.eg.db_3.2.3        
##  [3] RSQLite_1.0.0              DBI_0.3.1                 
##  [5] AnnotationDbi_1.32.3       genefilter_1.52.1         
##  [7] RColorBrewer_1.1-2         pheatmap_1.0.8            
##  [9] DESeq2_1.10.1              RcppArmadillo_0.6.700.3.0 
## [11] Rcpp_0.12.4                SummarizedExperiment_1.0.2
## [13] Biobase_2.30.0             GenomicRanges_1.22.4      
## [15] GenomeInfoDb_1.6.3         IRanges_2.4.8             
## [17] S4Vectors_0.8.11           BiocGenerics_0.16.1       
## 
## loaded via a namespace (and not attached):
##  [1] locfit_1.5-9.1       splines_3.2.3        lattice_0.20-33     
##  [4] colorspace_1.2-6     htmltools_0.3.5      yaml_2.1.13         
##  [7] survival_2.39-2      XML_3.98-1.4         foreign_0.8-66      
## [10] BiocParallel_1.4.3   lambda.r_1.1.7       plyr_1.8.3          
## [13] stringr_1.0.0        zlibbioc_1.16.0      munsell_0.4.3       
## [16] gtable_0.2.0         futile.logger_1.4.1  evaluate_0.8.3      
## [19] labeling_0.3         latticeExtra_0.6-28  knitr_1.12.3        
## [22] geneplotter_1.48.0   acepack_1.3-3.3      xtable_1.8-2        
## [25] scales_0.4.0         formatR_1.3          Hmisc_3.17-3        
## [28] annotate_1.48.0      XVector_0.10.0       gridExtra_2.2.1     
## [31] ggplot2_2.1.0        digest_0.6.9         stringi_1.0-1       
## [34] grid_3.2.3           tools_3.2.3          magrittr_1.5        
## [37] Formula_1.2-1        cluster_2.0.4        futile.options_1.0.0
## [40] Matrix_1.2-5         rmarkdown_0.9.5      rpart_4.1-10        
## [43] nnet_7.3-11
