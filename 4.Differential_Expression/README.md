## RNA-seq : recap

- In the last section we worked through the process of quality assessment and alignment 
for RNA-seq data
- This typically takes place on the command line, but can also be done from within R.
- The end result was the generation of count data (counts of reads aligned to each gene, per sample) using the FeatureCounts command from Subread/Rsubread.
- Now that we've got count data in R, we can begin our differental expression analysis.

---

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


(save source from untitled to yeast_data.R) 
---

## Count data

- Note: I have now aligned the data for ALL CHROMOSOMES and generated counts, so we are working with data from all 7127 genes.

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

---

## Extract count data

 - Remove annotation columns
 - Add row names

```{r}
counts = fcData[, 7:12]
rownames(counts) = fcData$Geneid
counts %>% head()
```

---

## Read counts per sample

 - Normalisation process (slightly different for each analysis method) takes 
 "library size" (number of reads generated for each sample) into account.

```{r}
colSums(counts)
```

---

## Visualise via bar plot

