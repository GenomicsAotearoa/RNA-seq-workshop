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

```

> getwd()
[1] "/Users/ngonifaya/Desktop/RNA_seq"

> library(dplyr)

> fcData = read.table('yeast_chr1_counts.txt', sep='\t', header=TRUE)

> fcData %>% head()
     Geneid Chr Start   End Strand Length ...STAR.SRR014335.chr1.Aligned.out.sam
1   YDL248W  IV  1802  2953      +   1152                                      0
2 YDL247W-A  IV  3762  3836      +     75                                      0
3   YDL247W  IV  5985  7814      +   1830                                      0
4   YDL246C  IV  8683  9756      -   1074                                      0
5   YDL245C  IV 11657 13360      -   1704                                      0
6   YDL244W  IV 16204 17226      +   1023                                      0
  ...STAR.SRR014336.chr1.Aligned.out.sam ...STAR.SRR014337.chr1.Aligned.out.sam
1                                      0                                      0
2                                      0                                      0
3                                      0                                      0
4                                      0                                      0
5                                      0                                      0
6                                      0                                      0
  ...STAR.SRR014339.chr1.Aligned.out.sam ...STAR.SRR014340.chr1.Aligned.out.sam
  .
  .
  .

```

Further checking our dataset

```{r}

> dim(fcData)

> names(fcData)
 [1] "Geneid"                                 "Chr"                                   
 [3] "Start"                                  "End"                                   
 [5] "Strand"                                 "Length"                                
 [7] "...STAR.SRR014335.chr1.Aligned.out.sam" "...STAR.SRR014336.chr1.Aligned.out.sam"
 [9] "...STAR.SRR014337.chr1.Aligned.out.sam" "...STAR.SRR014339.chr1.Aligned.out.sam"
[11] "...STAR.SRR014340.chr1.Aligned.out.sam" "...STAR.SRR014341.chr1.Aligned.out.sam"

```
## Rename data columns

```{r}

> names(fcData)[7:12] = c("SRR014335", "SRR014336", "SRR014337", 
                        "SRR014339", "SRR014340", "SRR014341")
 
> fcData %>% head()
     Geneid Chr Start   End Strand Length SRR014335 SRR014336 SRR014337 SRR014339
1   YDL248W  IV  1802  2953      +   1152         0         0         0         0
2 YDL247W-A  IV  3762  3836      +     75         0         0         0         0
3   YDL247W  IV  5985  7814      +   1830         0         0         0         0
4   YDL246C  IV  8683  9756      -   1074         0         0         0         0
5   YDL245C  IV 11657 13360      -   1704         0         0         0         0
6   YDL244W  IV 16204 17226      +   1023         0         0         0         0
  SRR014340 SRR014341
1         0         0
2         0         0
3         0         0
4         0         0
5         0         0
6         0         0

```

## Extract count data

 - Remove annotation columns
 - Add row names

```{r}

> counts = fcData[, 7:12]

> rownames(counts) = fcData$Geneid

> counts %>% head()
          SRR014335 SRR014336 SRR014337 SRR014339 SRR014340 SRR014341
YDL248W           0         0         0         0         0         0
YDL247W-A         0         0         0         0         0         0
YDL247W           0         0         0         0         0         0
YDL246C           0         0         0         0         0         0
YDL245C           0         0         0         0         0         0
YDL244W           0         0         0         0         0         0


```

## Read counts per sample

 - Normalisation process (slightly different for each analysis method) takes 
 "library size" (number of reads generated for each sample) into account.

```{r}

> colSums(counts)
SRR014335 SRR014336 SRR014337 SRR014339 SRR014340 SRR014341 
    66112     66607     75730     52490     52915     51802 

```

## Visualise via bar plot

```{r}

> colSums(counts) %>% barplot(., las=3, ylab="Reads mapped per sample")
```
![Alt text](https://github.com/foreal17/RNA-seq-workshop/blob/master/Prep_Files/Images/Rplot_bar.png)


*Now we are ready for differential expression analysis*

---

# Detecting differential expression: DESeq2

 - The DESeq2 package uses the *Negative Binomial* distribution to model the count data from each sample.
 - A statistical test based on the Negative Binomial distribution (via a generalized linear model, GLM) can be used to assess differential expression for each gene.
 - Use of the Negative Binomial distribution attempts to accurately capture the variation 
 that is observed for count data.

More information about DESeq2: <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8">article by Love et al, 2014</a> 


```{r, comment=FALSE, message=FALSE}

> library(DESeq2)

# Specify "conditions" (groups: WT and MT)
> conds <- c("WT","WT","WT","MT","MT","MT")

# Create object of class CountDataSet derived from eSet class
> dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = data.frame(conds=factor(conds)), design = formula(~conds))

# CountDataSet has similar accessor methods as eSet class.
> knitr::kable(counts(dds)[1:4, ]) 
```

```{r, message=FALSE, warning=FALSE}

# Fit DESeq model to identify DE transcripts
> dds <- DESeq(dds)
> res <- DESeq2::results(dds)
> knitr::kable(res[1:6,])

```

#### Note: p-value adjustment

 - The "padj" column of the DESeq2 results (`res`) contains adjusted p-values (FDR).
 - Can use the `p.adjust` function to manually adjust the DESeq2 p-values if needed (e.g., to use Holm correction)

```{r}

# Remove rows with NAs
> res = na.omit(res)

# Get the rows of "res" with significant adjusted p-values
> resPadj<-res[res$padj <= 0.05 , ]

# Get dimensions
> dim(resPadj)

```

```{r, warning=FALSE}
# Number of adjusted p-values less than 0.05
> sum(res$padj <= 0.05)

# Check that this is the same using p.adjust with FDR correction
> sum(p.adjust(res$pvalue, method="fdr") <= 0.05)

# Number of Holm adjusted p-values less than 0.05
> sum(p.adjust(res$pvalue, method="holm") <= 0.01)

```

---

## Detecting differential expression: edgeR

 - The edgeR package also uses the negative binomial distribution to model the RNA-seq count data.
 - Takes an empirical Bayes approach to the statistical analysis, in a similar way to how the `limma` package handles microarray data.

```{r}

# Compute exact test for the negative binomial distribution.
> et <- exactTest(y) 

> knitr::kable(topTags(et, n=4)$table)

```
#### edgeR: adjusted p-values

```{r}

> edge <- as.data.frame(topTags(et, n=nrow(counts))) 

> sum(edge$FDR <= 0.05)

> sum(p.adjust(edge$PValue, method="fdr") <= 0.05)

## Get the rows of "edge" with significant adjusted p-values
> edgePadj <- edge[edge$FDR <= 0.05, ]

```

### DESeq2 vs edgeR

 - Generate Venn diagram to compare DESeq2 and edgeR results.

```{r}

> library(systemPipeR)

> setlist <- list(edgeRexact=rownames(edgePadj), DESeq2=rownames(resPadj))

> vennset <- overLapper(setlist=setlist[1:2], type="vennsets")

> vennPlot(vennset, mymain="DEG Comparison")

```
