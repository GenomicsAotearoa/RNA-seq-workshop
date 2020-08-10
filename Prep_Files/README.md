# Workshop Preparation

---
## HPC and dataset check
Things to check:
- Are you able to login into NeSI (If not please get in touch with Dini - dinindu.senanayake@nesi.org.nz)
- Do you have RStudio installed on your laptop

Once you login to NeSI, you will see a folder RNA_seq in your home directory:

```bash

$ ssh -y [your_username]@lander.nesi.org.nz

#logging to the training platform
$ ssh ga-vl01

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

#### install `DESeq2`, `edgeR`, `limma`, `systemPipeR`

```

# Older versions of R (older than 3.5)
> source("https://bioconductor.org/biocLite.R")

> biocLite("DESeq2")

# Newer versions of R
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

> BiocManager::install("DESeq2")

# Installing edgeR will also install limma as it is a dependency package
> BiocManager::install("edgeR")

# systemPipeR: R package for building and running automated end-to-end analysis workflows for a wide range of next generation sequence (NGS) applications such as RNA-Seq, ChIP-Seq, VAR-Seq and Ribo-Seq.

> BiocManager::install("systemPipeR")

```

---

## Overrepresentation analysis prep

#### install `goseq`, `ggplot2`, `org.Sc.sgd.db`, `venn`

```
# goseq: R package for detecting Gene Ontology and/or other user defined categories which are over/under represented in RNA-seq data
> BiocManager::install("goseq")

# ggplot: is a data visualization package
> install.packages("ggplot2")

# org.Sc.sgd.db: Genome wide annotation for Yeast, primarily based on mapping using ORF identifiers from SGD.
> BiocManager::install("org.Sc.sgd.db")

# venn: Draws and displays Venn diagrams up to 7 sets, and any Boolean union of set intersections.

> install.packages("venn")

```











