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

## Count data

- Note: I have now aligned the data for ALL CHROMOSOMES and generated counts, so we are working with data from all 7127 genes.

```{r, message=FALSE, warning=FALSE}
library(dplyr)
fcData = read.table('YeastData/Counts/yeast_counts_all_chr.txt', sep='\t', header=TRUE)
fcData %>% head()
```

---

## Count data

```{r}
dim(fcData)
names(fcData)
```

---
