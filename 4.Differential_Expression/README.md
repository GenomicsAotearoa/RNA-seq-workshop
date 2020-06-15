# RNA-seq : recap

- In the last section we worked through the process of quality assessment and alignment 
for RNA-seq data
- This typically takes place on the command line, but can also be done from within R.
- The end result was the generation of count data (counts of reads aligned to each gene, per sample) using the FeatureCounts command from Subread/Rsubread.
- Now that we've got count data in R, we can begin our differental expression analysis.

## Data set reminder

 - Data obtained from yeast RNA-seq experiment, <a href="https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000299">Lee et al 2008 </a>
 - Wild-type versus RNA degradation mutants 
 - Six samples (3 WT / 3 MT)
 - We are working with data from chromosome 1 to keep the files sizes relatively small.

---

## Create an RStudio project

One of the first benefits we will take advantage of in RStudio is something called an RStudio Project. An RStudio project allows you to more easily:

- Save data, files, variables, packages, etc. related to a specific analysis project
- Restart work where you left off
- Collaborate, especially if you are using version control such as git.

To create a project,
- Open RStudio and go to the File menu, and click New Project.
- In the window that opens select Existing Project, and browse to the RNA-seq folder we have created on our Desktop.
- Finally click Create Project


(save source from untitled to yeast_data.R and continue saving reguralry as we work) 

---

## Count data

- Note: I have now aligned the data for ALL CHROMOSOMES and generated counts, so we are working with data from all 7127 genes.

*Let's look at our dataset and manipulate it is as we prepare for differential expression.*

```R

> getwd()
[1] "/Users/ngonifaya/Desktop/RNA_seq"

> library(dplyr)

> fcData = read.table('yeast_counts_all_chr.txt', sep='\t', header=TRUE)

> fcData %>% head()
     Geneid Chr Start   End Strand Length ...STAR.SRR014335.Aligned.out.sam
1   YDL248W  IV  1802  2953      +   1152                                52
2 YDL247W-A  IV  3762  3836      +     75                                 0
3   YDL247W  IV  5985  7814      +   1830                                 2
4   YDL246C  IV  8683  9756      -   1074                                 0
5   YDL245C  IV 11657 13360      -   1704                                 0
6   YDL244W  IV 16204 17226      +   1023                                 6
  ...STAR.SRR014336.Aligned.out.sam ...STAR.SRR014337.Aligned.out.sam
1                                46                                36
2                                 0                                 0
3                                 4                                 2
4                                 0                                 1
5                                 3                                 0
6                                 6                                 5
  ...STAR.SRR014339.Aligned.out.sam ...STAR.SRR014340.Aligned.out.sam
1                                65                                70
.
.
.

```

Further checking our dataset

```R

> dim(fcData)
[1] 7127   12

> names(fcData)
 [1] "Geneid"                                 "Chr"                                   
 [3] "Start"                                  "End"                                   
 [5] "Strand"                                 "Length"                                
 [7] "...STAR.SRR014335.chr1.Aligned.out.sam" "...STAR.SRR014336.chr1.Aligned.out.sam"
 [9] "...STAR.SRR014337.chr1.Aligned.out.sam" "...STAR.SRR014339.chr1.Aligned.out.sam"
[11] "...STAR.SRR014340.chr1.Aligned.out.sam" "...STAR.SRR014341.chr1.Aligned.out.sam"

```
## Rename data columns

```R

> names(fcData)[7:12] = c("SRR014335", "SRR014336", "SRR014337", 
                        "SRR014339", "SRR014340", "SRR014341")
 
> fcData %>% head()
     Geneid Chr Start   End Strand Length SRR014335 SRR014336 SRR014337 SRR014339 SRR014340 SRR014341
1   YDL248W  IV  1802  2953      +   1152        52        46        36        65        70        78
2 YDL247W-A  IV  3762  3836      +     75         0         0         0         0         1         0
3   YDL247W  IV  5985  7814      +   1830         2         4         2         6         8         5
4   YDL246C  IV  8683  9756      -   1074         0         0         1         1         2         0
5   YDL245C  IV 11657 13360      -   1704         0         3         0         5         7         4
6   YDL244W  IV 16204 17226      +   1023         6         6         5        20        30        19

```

## Extract count data

 - Remove annotation columns
 - Add row names

```R

> counts = fcData[, 7:12]

> rownames(counts) = fcData$Geneid

> counts %>% head()
          SRR014335 SRR014336 SRR014337 SRR014339 SRR014340 SRR014341
YDL248W          52        46        36        65        70        78
YDL247W-A         0         0         0         0         1         0
YDL247W           2         4         2         6         8         5
YDL246C           0         0         1         1         2         0
YDL245C           0         3         0         5         7         4
YDL244W           6         6         5        20        30        19


```

## Read counts per sample

 - Normalisation process (slightly different for each analysis method) takes 
 "library size" (number of reads generated for each sample) into account.

```R

> colSums(counts)
SRR014335 SRR014336 SRR014337 SRR014339 SRR014340 SRR014341 
  4915975   4892227   4778158   4618409   4719413   4554283 

```

## Visualise via bar plot

```R

> colSums(counts) %>% barplot(., las=3, ylab="Reads mapped per sample")
```
![Alt text](https://github.com/foreal17/RNA-seq-workshop/blob/master/Prep_Files/Images/Rplot_bar1.png)


*Now we are ready for differential expression analysis*

---

# Detecting differential expression: 

We are going to identify genes that are differential expressed using 3 different packages (time allowing) and compare the results. 

---
# DESeq2

 - The DESeq2 package uses the *Negative Binomial* distribution to model the count data from each sample.
 - A statistical test based on the Negative Binomial distribution (via a generalized linear model, GLM) can be used to assess differential expression for each gene.
 - Use of the Negative Binomial distribution attempts to accurately capture the variation 
 that is observed for count data.

More information about DESeq2: <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8">article by Love et al, 2014</a> 


```R

> library(DESeq2)

# Specify "conditions" (groups: WT and MT)
> conds <- c("WT","WT","WT","MT","MT","MT")

# Create object of class CountDataSet derived from eSet class
> dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = data.frame(conds=factor(conds)), design = formula(~conds))

# CountDataSet has similar accessor methods as eSet class.
> knitr::kable(counts(dds)[1:4, ]) 
|          | SRR014335| SRR014336| SRR014337| SRR014339| SRR014340| SRR014341|
|:---------|---------:|---------:|---------:|---------:|---------:|---------:|
|YDL248W   |        52|        46|        36|        65|        70|        78|
|YDL247W-A |         0|         0|         0|         0|         1|         0|
|YDL247W   |         2|         4|         2|         6|         8|         5|
|YDL246C   |         0|         0|         1|         1|         2|         0|

```

#### Fit DESeq model to identify DE transcripts
```R

> dds <- DESeq(dds)
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing

> res <- DESeq2::results(dds)

> knitr::kable(res[1:6,])
|          |   baseMean| log2FoldChange|     lfcSE|       stat|    pvalue|      padj|
|:---------|----------:|--------------:|---------:|----------:|---------:|---------:|
|YDL248W   | 56.2230316|     -0.2339470| 0.2278417| -1.0267961| 0.3045165| 0.3663639|
|YDL247W-A |  0.1397714|     -0.5255769| 4.0804729| -0.1288030| 0.8975136| 0.9173926|
|YDL247W   |  4.2428211|     -0.8100552| 0.8332543| -0.9721584| 0.3309718| 0.3948536|
|YDL246C   |  0.6182409|     -1.0326364| 2.1560899| -0.4789394| 0.6319817| 0.6915148|
|YDL245C   |  2.8486580|     -1.9781751| 1.1758845| -1.6822869| 0.0925132| 0.1230028|
|YDL244W   | 13.0883354|     -1.5823354| 0.5128788| -3.0852037| 0.0020341| 0.0032860|

```

#### Note: p-value adjustment

 - The "padj" column of the DESeq2 results (`res`) contains adjusted p-values (FDR).
 - Can use the `p.adjust` function to manually adjust the DESeq2 p-values if needed (e.g., to use Holm correction)

```R

# Remove rows with NAs
> res = na.omit(res)

# Get the rows of "res" with significant adjusted p-values
> resPadj<-res[res$padj <= 0.05 , ]

# Get dimensions
> dim(resPadj)
[1] 4811    6

```

```R

# Number of adjusted p-values less than 0.05
> sum(res$padj <= 0.05)
[1] 4811

# Check that this is the same using p.adjust with FDR correction
> sum(p.adjust(res$pvalue, method="fdr") <= 0.05)
[1] 4811

# Number of Holm adjusted p-values less than 0.05
> sum(p.adjust(res$pvalue, method="holm") <= 0.01)
[1] 3429

```
#### Volcano plot

```R

#reset par
> par(mfrow=c(1,1))

# Make a basic volcano plot
> with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
> with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
> with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
```

---

![Alt text](https://github.com/foreal17/RNA-seq-workshop/blob/master/Prep_Files/Images/Volcano%20plot.png)

# edgeR

 - The edgeR package also uses the negative binomial distribution to model the RNA-seq count data.
 - Takes an empirical Bayes approach to the statistical analysis, in a similar way to how the `limma` package handles microarray data.
 
*Identify DEGs with edgeR’s Exact Method*
 
 ```R
 
> library(edgeR)

# Construct DGEList object
> y <- DGEList(counts=counts, group=conds)

# Calculate library size (counts per sample)
> y <- calcNormFactors(y)

# Estimate common dispersion (overall variability)
> y <- estimateCommonDisp(y)

# Estimate tagwise dispersion (per gene variability)
> y <- estimateTagwiseDisp(y)

```

*edgeR analysis and output

```R

# Compute exact test for the negative binomial distribution.
> et <- exactTest(y) 

> knitr::kable(topTags(et, n=4)$table)
|        |     logFC|    logCPM| PValue| FDR|
|:-------|---------:|---------:|------:|---:|
|YCR077C | 13.098599|  6.629306|      0|   0|
|snR71   | -9.870811|  8.294459|      0|   0|
|snR59   | -8.752555| 10.105453|      0|   0|
|snR53   | -8.503101|  7.687908|      0|   0|

```
*adjusted p-values*

```R

> edge <- as.data.frame(topTags(et, n=nrow(counts))) 

> sum(edge$FDR <= 0.05)
[1] 5242

> sum(p.adjust(edge$PValue, method="fdr") <= 0.05)
[1] 5242

## Get the rows of "edge" with significant adjusted p-values
> edgePadj <- edge[edge$FDR <= 0.05, ]

```
---

### DESeq2 vs edgeR

 - Generate Venn diagram to compare DESeq2 and edgeR results.

```R

> library(systemPipeR)

> setlist <- list(edgeRexact=rownames(edgePadj), DESeq2=rownames(resPadj))

> vennset <- overLapper(setlist=setlist[1:2], type="vennsets")

> vennPlot(vennset, mymain="DEG Comparison")

```
![Alt text](https://github.com/foreal17/RNA-seq-workshop/blob/master/Prep_Files/Images/DEG_vs_edgeR.png)

---

# Limma

<!-- NEED TO EXPLAIN "TREND=TRUE" -->

 - Limma can be used for analysis (log-scale normality-based assumption rather than Negative Binomial for count data)
 - Use data transformation and log to satisfy normality assumptions (CPM = Counts per Million).

```R
> design <- model.matrix(~conds)

> dge <- DGEList(counts=counts)

> dge <- calcNormFactors(dge)

> logCPM <- cpm(dge, log=TRUE, prior.count=3)

> options(width=100)

> head(logCPM, 3)
           SRR014335  SRR014336  SRR014337  SRR014339  SRR014340  SRR014341
YDL248W    3.7199528  3.5561232  3.2538405  3.6446399  3.7156488  3.9155366
YDL247W-A -0.6765789 -0.6765789 -0.6765789 -0.6765789 -0.3140297 -0.6765789
YDL247W    0.1484688  0.6727144  0.1645731  0.7843936  1.0395626  0.6349276

```

### Limma: voom

 - The "voom" function estimates relationship between the mean and the variance of the logCPM data, normalises the data, and 
 creates "precision weights" for each observation that are incorporated into the limma analysis.

```R

> v <- voom(dge, design, plot=TRUE)

```
![Alt text](https://github.com/foreal17/RNA-seq-workshop/blob/master/Prep_Files/Images/Voom_Mean_Variance.png)


*Limma: voom* (impact on first three samples)

```R

> par(mfrow=c(1,3))

> for(i in 1:3){plot(logCPM[,i], v$E[,i], xlab="LogCPM", ylab="Voom", main=colnames(logCPM)[i])abline(0,1)}

```
![Alt text](https://github.com/foreal17/RNA-seq-workshop/blob/master/Prep_Files/Images/voom_3_samples.png)

```R

> fit <- lmFit(v, design)

> fit <- eBayes(fit)

> tt <- topTable(fit, coef=ncol(design), n=nrow(counts))

head(tt)
           logFC  AveExpr        t      P.Value    adj.P.Val        B
YAL038W 2.313985 10.80214 319.5950 3.725452e-13 8.850433e-10 21.08951
YOR161C 2.568389 10.80811 321.9628 3.574510e-13 8.850433e-10 21.08187
YML128C 1.640664 11.40819 286.6167 6.857932e-13 9.775297e-10 20.84520
YMR105C 2.772539  9.65092 331.8249 3.018547e-13 8.850433e-10 20.16815
YHL021C 2.034496 10.17510 269.4034 9.702963e-13 1.152550e-09 20.07857
YDR516C 2.085424 10.05426 260.8061 1.163655e-12 1.184767e-09 19.87217

```


*limma: adjusted p-values*

```R

> sum(tt$adj.P.Val <= 0.05)
[1] 5140

> sum(p.adjust(tt$P.Value, method="fdr") <= 0.05)
[1] 5140

## Get the rows of top table with significant adjusted p-values
> limmaPadj <- tt[tt$adj.P.Val <= 0.05, ]

```

---

## Limma vs edgeR vs DESeq2

```R

> setlist <- list(edgeRexact=rownames(edgePadj), DESeq2=rownames(resPadj), LimmaVoom=rownames(limmaPadj))

> vennset <- overLapper(setlist=setlist[1:3], type="vennsets")

vennPlot(vennset, mymain="DEG Comparison")

```
![Alt text](https://github.com/foreal17/RNA-seq-workshop/blob/master/Prep_Files/Images/DEG_Comparison.png)

---
