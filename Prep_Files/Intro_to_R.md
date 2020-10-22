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
- A numeric array of columns and rows (tabular data structures)

```{r}
# creating matrices from vectors
> x <- c(1,2,3,4,5)
> y <- c(7,7,8,9,10)
> m1 <- cbind(x,y)
> m1
     x  y
[1,] 1  6
[2,] 2  7
[3,] 3  8
[4,] 4  9
[5,] 5 10

> dim(m1)
[1] 5 2

#transpose m1
> m2 <- t(m1)
> m2
  [,1] [,2] [,3] [,4] [,5]
x    1    2    3    4    5
y    6    7    8    9   10
> dim(m2)
[1] 2 5


# Creating matrices by directly adding elements
> m3<-matrix(c(1,3,2,5,-1,2,2,3,9),nrow=3)
> m3
     [,1] [,2] [,3]
[1,]    1    5    2
[2,]    3   -1    3
[3,]    2    2    9

```

#### Data frames
- A tabular structure that can have different modes ( numeric, character etc)

```{r}
# creating data frame
> chr <- c("chr1", "chr1", "chr2", "chr2")
> strand <- c("-","-","+","+")
> start<- c(200,4000,100,400)
> end<-c(250,410,200,450)
> mydata <- data.frame(chr,start,end,strand)

#change column names
names(mydata) <- c("chr","start","end","strand")
> mydata 
   chr start end strand
1 chr1   200 250      -
2 chr1  4000 410      -
3 chr2   100 200      +
4 chr2   400 450      +
```

- Accessing elements in a data frame
```
#looking at columns 2 to 4
> mydata[,2:4]
  start end strand
1   200 250      -
2  4000 410      -
3   100 200      +
4   400 450      +

# columns chr and start from data frame
> mydata[,c("chr","start")]
   chr start
1 chr1   200
2 chr1  4000
3 chr2   100
4 chr2   400

# variable start in the data frame
> mydata$start
[1]  200 4000  100  400

# get 1st and 3rd rows
> mydata[c(1,3),]
   chr start end strand
1 chr1   200 250      -
3 chr2   100 200      +

# get all rows where start>400
> mydata[mydata$start>400,]
   chr start end strand
2 chr1  4000 410      -
```

#### Lists
- An ordered collection of objects.
- In a list, you can store a variety of (possibly unrelated) objects under one name. 
- Each object or element in a list has a numbered position and can have names.

```{r}
# creating a list with 4 objects
> mylist
$name
[1] "Bob"

$mynumbers
[1] 1 2 3

$mymatrix
     [,1] [,2]
[1,]    1    3
[2,]    2    4

$age
[1] 5.3


# Objects in a list can be extracted using double square-bracket `[[]]` with either position or name in the brackets
> mylist[[2]]
[1] 1 2 3

> mylist[["mymatrix"]]
     [,1] [,2]
[1,]    1    3
[2,]    2    4


# This also works
> mylist$age
[1] 5.3

```



