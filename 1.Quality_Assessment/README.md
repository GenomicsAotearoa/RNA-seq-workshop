# Quality control of the sequencing data.

Several tools available to do so. For this workshop, we will use fastqc.

First, it is always good to verify where we are:

```
$ pwd
/home/[your_username]
# good I am ready to work

```

Checking to make sure we have the Raw files for the workshop.

```
$ ls
modload.sh RNA_seq ...
```
Creating a directory where to store the QC data:

```
$ cd RNA_seq

$ ls
Raw

$ mkdir - QC

```

Since we are working on the NeSI HPC, we need to search and load the package before we start using it.

Search
```
$ module spider fastqc
```

and then load 

```
$ module load FastQC
```

Now we can start the quality control:

```
$ fastqc -o QC/ /Raw/*

```
You will see an automatically updating output message telling you the progress of the analysis. It will start like this:

```
Started analysis of SRR014335-chr1.fastq
Approx 5% complete for SRR014335-chr1.fastq
Approx 10% complete for SRR014335-chr1.fastq
Approx 15% complete for SRR014335-chr1.fastq
Approx 20% complete for SRR014335-chr1.fastq
Approx 25% complete for SRR014335-chr1.fastq
Approx 30% complete for SRR014335-chr1.fastq
Approx 35% complete for SRR014335-chr1.fastq

```

The FastQC program has created several new files within our RNA_seq/QC/ directory.

```
$ ls QC
SRR014335-chr1_fastqc.html  SRR014336-chr1_fastqc.zip   SRR014339-chr1_fastqc.html  SRR014340-chr1_fastqc.zip
SRR014335-chr1_fastqc.zip   SRR014337-chr1_fastqc.html  SRR014339-chr1_fastqc.zip   SRR014341-chr1_fastqc.html
SRR014336-chr1_fastqc.html  SRR014337-chr1_fastqc.zip   SRR014340-chr1_fastqc.html  SRR014341-chr1_fastqc.zip

```

## Viewing the FastQC results

If we were working on our local computers, we’d be able to look at each of these HTML files by opening them in a web browser.

However, these files are currently sitting on our remote NeSI HPC, where our local computer can’t see them. And, since we are only logging into NeSI via the command line - it doesn’t have any web browser setup to display these files either.

So the easiest way to look at these webpage summary reports will be to transfer them to our local computers (i.e. your laptop).

To transfer a file from a remote server to our own machines, we will use scp.

First we will make a new directory on our computer to store the HTML files we’re transferring. Let’s put it on our desktop for now. Open a new tab in your terminal program (you can use the pull down menu at the top of your screen or the Cmd+t keyboard shortcut) and type:

```

mkdir -p ~/Desktop/fastqc_html 

scp -r fayfa80p@login.mahuika.nesi.org.nz:/home/fayfa80p/RNA_seq/QC/ ~/Desktop/fastqc_html

```
