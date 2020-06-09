# Workshop Preparation

Things to check:
- Are you able to login into NeSI (If not please get in touch with Dini - dinindu.senanayake@nesi.org.nz)
- Do you have RStudio installed on your laptop

Once you login to NeSI, you will see a folder RNA_seq in your home directory:

```
$ cd RNA_seq

$ ls -F
Genome/ Raw/

```

If you see the RNA_seq directory with 2 sub-directories: Genome/ and Raw/, yayyy you are good to go.

If not, no pressure, we can still copy the RNA_seq directory to your home directory.

```
$ cp -r /nesi/nobackup/nesi02659/RNA_seq/ /home/Your_Username

$ ls
RNA_seq

```

---


## Differential expression prep
These are the packages needed for the differential expression analysis

```
# install DESeq2, edgeR and limma

# Older versions of R (older than 3.5)
> source("https://bioconductor.org/biocLite.R")

> biocLite("DESeq2")


# Newer versions of R
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

> BiocManager::install("DESeq2")

# Installing edgeR will also install limma as it is a dependency package
> BiocManager::install("edgeR")

```


```
# install systemPipeR: R package for building and running automated end-to-end analysis workflows for a wide range of next generation sequence (NGS) applications such as RNA-Seq, ChIP-Seq, VAR-Seq and Ribo-Seq. 

> BiocManager::install("systemPipeR")

```





