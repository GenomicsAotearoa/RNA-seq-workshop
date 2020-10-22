# Introduction to R for Genomics

### Why use R for genomics?
- R, with its statistical analysis heritage, plotting features, and rich user-contributed packages is one of the best languages for the task of analyzing genomic data.
- High-dimensional genomics datasets are usually suitable to be analyzed with core R packages and functions.
- Bioconductor and CRAN have an array of specialized tools for doing genomics-specific analysis.

---
### Installing Packages
- CRAN packages can be installed using install.packages()

```{r}
# install package called `tidyverse`
> install.packages("tidyverse")
```
- You can install bioconductor packages with a specific installer script.

```{r}
# get the installer package if you don't have
> install.packages("BiocManager")

# install bioconductor package "DESeq2"
> BiocManager::install("DESeq2")
```

- You can also install packages from Github, or from the source by downloading the source file

### Updating CRAN and Bioconductor packages.

```{r}
# updating CRAN packages
update.packages()

# updating bioconductor packages
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install()
```

### Getting help on functions and packages
You can get help on functions by using `help()` and ?package_name.

```
> ?hist

> help("hist")
```

Above all, google search can be your best friend to quickly find answers.

---

### R as a calculator
- One of the most basic ways to use R is as a calculator.
Take a minute tp play around with R to perform some calculations.

### Data structures
- R has multiple data structures
- R deals with named data structures

#### Vectors
- A list of elements of the same type e.g numeric, character, logical
- Operation on vectors propagates to all the elements

```{r}
# create vectors named x and y with 5 elements

> x <- c(1,2,3,4,5)
[1] 1 2 3 4 5

> y <- 6:10
[1]  6  7  8  9 10

# scalar addition and multiplication
> y + 2
[1]  8  9 10 11 12

> x * 2
[1]  2  4  6  8 10

# create a vector of 1s, length 3
> r1 <- rep(1,3)

> length(r1)
[1] 3

> class(r1)
[1] "numeric"

```

#### Matrices

