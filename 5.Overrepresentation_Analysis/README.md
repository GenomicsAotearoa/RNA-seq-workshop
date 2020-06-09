# Overrepresentation analysis (Gene Ontology)

Over Representation Analysis (<a href="https://academic.oup.com/bioinformatics/article/20/18/3710/202612">Boyle et al. 2004</a>) is a widely used approach to determine whether known biological functions or processes are over-represented (= enriched) in an experimentally-derived gene list, e.g. a list of differentially expressed genes (DEGs).

 - Can perform overrepresentation analysis (e.g., GATHER, GeneSetDB, PantherDB) in R.
 - The basic principles are to:
    - identify a collection of differentially expressed genes
    - test to see if genes that are members of specific gene sets (e.g.Reactome pathways, Gene Ontologies categories) are differentially expressed more often than would
    be expected by chance.

*Some caveats for RNA-seq data*

 - The gene-set analysis methods are applicable to transcriptomic data from both microarrays and RNA-seq.
 - One caveat, however, is that the results need to take _gene length_ into account.
    - RNA-seq tends to produce higher expression levels (i.e., greater counts) for longer genes: a longer transcript implies more aligned fragments, and thus higher counts. This also gives these genes a great chance of being statistically differentially expressed.
    - Some gene sets (pathways, GO terms) tend to involve families of long genes: if long genes have a great chance of being detected as differentially expressed, then _gene sets_ consisting of long genes will have a great chance of appeared to be enriched in the analysis. 

---

## GOseq

 - The GOseq methodology (<a href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-2-r14">Young et al., 2010</a>) overcomes this issue by allowing the over-representation analysis to be adjusted for gene length.
     - modification to hypergeometric sampling probability 
     - exact method (resampling) and an approximation-based method
 - More recent publications have also applied this gene-length correction to GSEA-based methods.
 - Still not widely understood to be an issue when performing RNA-seq pathway analysis, but REALLY important to take into account.

*Proportion DE by gene length and reads* 

![Alt text](https://github.com/foreal17/RNA-seq-workshop/blob/master/Prep_Files/Images/goseq-1.png)


*Gene set ranks by standard analysis*

![Alt text](https://github.com/foreal17/RNA-seq-workshop/blob/master/Prep_Files/Images/goseq-2.png)


*Gene set ranks by GOseq*

![Alt text](https://github.com/foreal17/RNA-seq-workshop/blob/master/Prep_Files/Images/goseq-3.png)

---

## GOseq analysis

 - Need to figure out if our organism is supported... (code is "sacCer")

```{r}

> library(goseq)

> supportedOrganisms()

```


*Define differentially expressed genes*

 - Create a vector of 0's and 1's to denot whether or not genes are differentially expressed (limma analysis: topTable).
 - Add gene names to the vector so that GOSeq knows which gene each data point relates to.

```{r}

> genes <- ifelse(tt$adj.P.Val < 0.05, 1, 0)

> names(genes) <- rownames(tt)

> head(genes)

> table(genes)

```

#### Calculate gene weights
  
 - Put genes into length-based "bins", and plot length vs proportion differentially expressed
 - Likely restricts to only those genes with GO annotation
  
```{r, message=FALSE, warning=FALSE, fig.height=5}
pwf=nullp(genes,"sacCer1","ensGene")
```

```{r, echo=FALSE, eval=FALSE}

> library(GenomicFeatures)

> txdb <- makeTxDbFromGFF(file="Saccharomyces_cerevisiae.R64-1-1.99.gtf", format = "gtf", dataSource = "SGD", organism = "Saccharomyces cerevisiae")

> txsByGene=transcriptsBy(txdb,"gene")

> lengthData=median(width(txsByGene))

```


#### Inspect output

 - Report length (bias) and weight data per gene.

```{r}

> head(pwf)

```

#### Gene lengths and weights

```{r, fig.height=5.5, fig.width=10}

> par(mfrow=c(1,2))

> hist(pwf$bias.data,30)

>hist(pwf$pwf,30)

```


#### Gene length vs average expression

```{r, fig.height=5.5, fig.width=9, warning=FALSE, message=FALSE}

> library(ggplot2)

> data.frame(logGeneLength = log2(pwf$bias.data), avgExpr = tt$AveExpr) %>% 
  ggplot(., aes(x=logGeneLength, y=avgExpr)) + geom_point(size=0.2) + 
  geom_smooth(method='lm')
  
```

#### Length correction in GOSeq

 - Uses "Wallenius approximation" to perform correction.
 - Essentially it is performing a weighted Fisher's Exact Test, but each gene in the 
 2x2 data does not contribute equally to the per cell count - it instead contributes its weight (based on its length).
 - This means that a gene set (e.g., GO term) containing lots of significant short genes will 
 be considered more likely to be enriched that a gene set with a similar proportion of long genes
 that are differentially expressed.


#### Run GOSeq with gene length correction

```{r, cache=TRUE}

> GO.wall=goseq(pwf, "sacCer1", "ensGene")

```


#### Output: Wallenius method

```{r}

> head(GO.wall)

```

#### P-value adjustment

```{r}

> GO.wall.padj <- p.adjust(GO.wall$over_represented_pvalue, method="fdr")

> sum(GO.wall.padj < 0.05)

> GO.wall.sig <- GO.wall$category[GO.wall.padj < 0.05]

> length(GO.wall.sig)

> head(GO.wall.sig)

```

---

## GO terms

 - Can use the `GO.db` package to get more information about the significant gene sets.

```{r}

> library(GO.db)

> GOTERM[[GO.wall.sig[1]]]

```

---

#### Run GOSeq without gene length correction

```{r, cache=TRUE}

> GO.nobias=goseq(pwf, "sacCer1", "ensGene", method="Hypergeometric")

```

*Output: Hypergeomtric (Fisher) method*

```{r}

> head(GO.nobias)

```


#### P-value adjustment

```{r}

> GO.nobias.padj <- p.adjust(GO.nobias$over_represented_pvalue, method="fdr")

> sum(GO.nobias.padj < 0.05)

> GO.nobias.sig <- GO.nobias$category[GO.nobias.padj < 0.05]

> length(GO.nobias.sig)

> head(GO.nobias.sig)

```


## Compare with and without adjustment

```{r}

> venn(list(GO.wall=GO.wall.sig, GO.nobias=GO.nobias.sig))

```

 - Extract out the different parts of the Venn diagram (yes, there are definitely better ways to do this).
 

```{r}

## Only significant in Hypergeomtric analysis
> onlySig.nobias <- setdiff(GO.nobias.sig, GO.wall.sig)

## Only significant in Wallenius analysis
> onlySig.wall <- setdiff(GO.wall.sig, GO.nobias.sig)

## Significant in both
> sig.wall.nobias <- intersect(GO.wall.sig, GO.nobias.sig)

```


#### Gene lengths and GO term membership

 - Can also extract gene length and GO membership information.

```{r}

> len=getlength(names(genes),"sacCer1","ensGene")

> head(len)

> go = getgo(names(genes),"sacCer1","ensGene")

> head(go)

```

---

#### Getting fancy...

 - Figure out which genes are in the significant GO groups, and then gets their lengths.

```{r}

> lengths.onlySig.nobias <- list()
> for(i in 1:length(onlySig.nobias)){
  inGo <- lapply(go, function(x)  onlySig.nobias[i] %in% x) %>% unlist()
  lengths.onlySig.nobias[[i]] <- len[inGo]
}

> lengths.onlySig.wall <- list()
> for(i in 1:length(onlySig.wall)){
  inGo <- lapply(go, function(x)  onlySig.wall[i] %in% x) %>% unlist()
  lengths.onlySig.wall[[i]] <- len[inGo]
}

```

#### Significant: Hypergeometric vs Wallenius 

 - Only Hypergeometric (pink) vs only Wallenius (blue)
 - Hypergeometric method is findings GO terms containing longer genes.

```{r, fig.width=9, fig.height=5}

> cols <- rep(c("lightpink", "lightblue"), c(10,7))

> boxplot(c(lengths.onlySig.nobias, lengths.onlySig.wall), col=cols)

```

#### All significant GO terms

```{r, echo=FALSE}

> lengths.sig.wall.nobias <- list()

> for(i in 1:length(sig.wall.nobias)){
  inGo <- lapply(go, function(x)  sig.wall.nobias[i] %in% x) %>% unlist()
  lengths.sig.wall.nobias[[i]] <- len[inGo]
}

```

```{r, echo=FALSE}

> cols <- rep(c("lightpink", grey(0.7), "lightblue"), c(10,37,7))

> avgLength <- lapply(c(lengths.onlySig.nobias, lengths.sig.wall.nobias, lengths.onlySig.wall),
                    median) %>% unlist()

> oo <- order(avgLength, decreasing=TRUE)

```

```{r, echo=FALSE, fig.width=12, fig.height=7}

> boxplot(c(lengths.onlySig.nobias, lengths.sig.wall.nobias, lengths.onlySig.wall)[oo],
        col=cols[oo], ylab="Gene Length", xlab = "GO term")
```


#### Gene length versus P-value

```{r, echo=FALSE}

> avgLength.wall <- lapply(c(lengths.onlySig.wall, lengths.sig.wall.nobias), median)

> avgLength.nobias <- lapply(c(lengths.onlySig.nobias, lengths.sig.wall.nobias), median)

> cols <- rep(c("blue", "lightblue", "red","lightpink"),
            c(length(lengths.onlySig.wall), length(lengths.sig.wall.nobias),
              length(lengths.onlySig.nobias), length(lengths.sig.wall.nobias)))

> plot(c(avgLength.wall, avgLength.nobias), 
     -log(c(GO.nobias.padj[GO.nobias.padj < 0.05], GO.wall.padj[GO.wall.padj < 0.05])), 
     col=cols, pch=16, xlab="Median Gene Length", ylab ="-log(FDR adj-pval)")
legend('topright', c("Only sig in NoBias", "Sig in both (nobias adjp)", 
                     "Sig in both (wal adjp)", "Only sig in Wall"), 
       fill=c("red", "pink", "lightblue", "blue"))
       
```

---

## Summary

 - Once we've generated count data, there are a number of ways to perform a differential expression analysis.
   - DESeq2 and edgeR model the count data, and assume a Negative Binomial distribution
   - Limma transforms (and logs) the data and assumes normality
 - Here we've seen that these three approches give quite similar results.
 - For Gene Set analysis, gene length needs to be accounted for, since longer transcripts are more likely to be found to be differentially expressed.
   - GOSeq adjusts for transcript length to take this into account.
   - It is also possible to use GOSeq with other types of annotation (e.g., Reactome or KEGG pathways).
