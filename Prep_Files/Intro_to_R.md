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

#### Factors
- This data structure is used to store categorial data.
- Very important in statistical models as a categorical data is treated differently from continuous variables.

```{r}
# creating factors
> features=c("promoter","exon","intron")
> class(features)
[1] "character"
> f.feat=factor(features)
> class(f.feat)
[1] "factor"
```
- It is important to note that when you are reading a data frame with read.table() or creating a data frame with data.frame() function, the character columns are stored as factors by default.

---

### Data types
- There are four common data types in R, they are numeric, logical, character and integer. All these data types can be used to create vectors natively.

#### Numeric data type

```{r}
#create a numeric vector x with 5 components
> x<-c(1,3,2,10,5)
> x
[1]  1  3  2 10  5
```
#### Logical data type
```
#create a logical vector x
> x<-c(TRUE,FALSE,TRUE)
> x
[1]  TRUE FALSE  TRUE
```

#### Character data type

```
# create a character vector
> x<-c("sds","sd","as")
> x
[1] "sds" "sd"  "as"
```

#### Integer data type
```
# create an integer vector
> x<-c(1L,2L,3L)
> x
[1] 1 2 3
```

---
### Reading and writing data
- Most of the genomics datasets come in table format e.g. gene expression count data or BED format.
- In R you can easily read tabular data with the read.table() and will be converted to a dataframe structure

*We are going to read the count data file located in RNA-seq-workshop/Prep_Files/Files*

```
> cd1 <- read.table("test_table.csv", sep=",", header=TRUE)
> head(cd1)
          X WT1 WT2 WT3 MT1 MT2 MT3
1   YDL248W  52  46  36  65  70  78
2 YDL247W-A   0   0   0   0   1   0
3   YDL247W   2   4   2   6   8   5
4   YDL246C   0   0   1   1   2   0
5   YDL245C   0   3   0   5   7   4
6   YDL244W   6   6   5  20  30  19

> dim(cd1)
[1] 7127    7

# check column names
> names(cd1)
[1] "X"   "WT1" "WT2" "WT3" "MT1" "MT2" "MT3"

# change column name
> names(cd1)[1] = c("Gene_Names")
> names(cd1)
[1] "Gene_Names" "WT1"        "WT2"        "WT3"        "MT1"        "MT2"        "MT3"  
```

*We are also going to write out our data*

- A data frame or matrix can be saved using write.table()
```
> write.table(cd1, file="new_test_data.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
```

### Plotting in R with base graphics

- sample 50 values from normal distribution and store them in vector x
-

#### Histogram

```
> x <- rnorm(50)
> hist(x)
```
![Alt_text](https://github.com/GenomicsAotearoa/RNA-seq-workshop/blob/master/Prep_Files/Files/Hist.png)

- using help (?hist), we can see how to manipulate the figure
```
> hist(x,main="Hello histogram!!!",col="red")
```
![Alt_text](https://github.com/GenomicsAotearoa/RNA-seq-workshop/blob/master/Prep_Files/Files/Hist_2.png)

#### Scatter plot
- Scatter plots are one of the most common plots you will meet in data analysis
- They are useful to visualise the relationships between two variables

*We are going to create another vector y and compare with x*

```

> y <- rnorm(50)

# plot scatter plot
> plot(x,y,main="scatterplot of random samples", ylab="y values",xlab="x values")

```
![](https://github.com/GenomicsAotearoa/RNA-seq-workshop/blob/master/Prep_Files/Files/Scatter_plot.png)


#### Box plot
- These depict groups of numerical data through their quartiles.
- The edges of the box denote the 1st and 3rd quartiles, and the line that crosses the box is the median. 
- The distance between the 1st and the 3rd quartiles is called interquartile range. 
- Outliers can be depicted as dots

*Using the same data from vectors x and y, we can plot a bocplot*

```
> boxplot(x,y,main="boxplots of random samples")
```
![](https://github.com/GenomicsAotearoa/RNA-seq-workshop/blob/master/Prep_Files/Files/boxplot.png)


#### Bar plot
- We are going to plot lengths of 4 genes

```
# Creating vectors with the variables
> gene_len <- c(250, 400, 320, 100)
> gene_names <- c("lipase1","lipase2", "protease1", "protease2")

# plotting the barplot
> barplot(gene_len, names.arg=gene_names, ylab = "Gene Length", main = "Lengths of experimental genes", col = c("red", "red", "blue", "blue"))

# creating the legend 
> legend("topright",legend=c("test","control"), fill=c("red","blue"))
```
![](https://github.com/GenomicsAotearoa/RNA-seq-workshop/blob/master/Prep_Files/Files/Barplot.png)






        

















