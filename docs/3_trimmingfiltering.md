# Cleaning Reads

!!! info "Objectives"
    
    - Understand what adapters are
    - Understand that trimming can be done based on Phred scores or sequence

<center>
![image](./theme_images/cleaning_reads.png){width="300"}
</center>

In the previous section, we took a high-level look at the quality of each of our samples using FastQC. We visualized per-base quality graphs showing the distribution of read quality at each base across all reads in a sample and extracted information about which samples fail which quality checks. Some of our samples failed quite a few quality metrics used by FastQC. This doesn’t mean, though, that our samples should be thrown out! It’s very common to have some quality metrics fail, and this may or may not be a problem for your downstream application. 

## Adapter removal

!!! danger ""

    - "Adapters" are short DNA sequences that are added to each read as part of the sequencing process (we won't get into "why" here).
    - These are removed as part of the data generation steps that occur during the sequencing run, but sometimes there is still a non-trivial amount of adapter sequence present in the FASTQ files.
    - Since the sequence is not part of the target genome (i.e., the genome if the species from which teh samples were derived) then we need to remove it to prevent it affecting the downstream analysis.
    - The FastQC application get detection adapter contamination in samples.

We will use a program called CutAdapt to filter poor quality reads and trim poor quality bases from our samples.


## How to act on fastq after QC.

We can do several trimming:

  * on quality using Phred score. What will be the Phred score?
  * on the sequences, if they contain adaptor sequences.

To do so, we can use on tools: The cutadapt application is often used to remove adapter sequence
from FASTQ files.
- The following syntax will remove the adapter sequence AACCGGTT from the file SRR014335-chr1.fastq, create a new file called SRR014335-chr1_trimmed.fastq, and write a summary to the log file SRR014335-chr1.log:

```bash
pwd
```
```
    /home/$USER/RNA_seq
```
```bash
mkdir Trimmed
module load cutadapt/4.1-gimkl-2022a-Python-3.10.5
cutadapt -q 20 -a AACCGGTT -o Trimmed/SRR014335-chr1_cutadapt.fastq Raw/SRR014335-chr1.fastq > Trimmed/SRR014335-chr1.log
```
!!! quote ""
    * **-q** (`--quality-cutoff`)  parameter can be used to trim low-quality ends from reads. If you specify a single cutoff value, the 3’ end of each read is trimmed.

```bash
less Trimmed/SRR014335-chr1.log
```

Now we should trim all samples.

```bash 
cd Raw
ls
```
```
SRR014335-chr1.fastq  SRR014336-chr1.fastq  SRR014337-chr1.fastq  SRR014339-chr1.fastq  SRR014340-chr1.fastq  SRR014341-chr1.fastq
```
```bash
for filename in *.fastq
do 
base=$(basename ${filename} .fastq)
cutadapt -q 20 -a AACCGGTT -o ../Trimmed/${base}.trimmed.fastq ${filename} > ../Trimmed/${base}.log
done
```

---

## MultiQC: `cutadapt` log files

 - If the log files from `cutadapt` are added to the directory containing the FastQC output, this information will also be incorporated into the MultiQC report the next time it is run.
 
```bash
cd ../MultiQC
cp ../Trimmed/*log .
multiqc .
```
![image](./Images/MQC2.png)

- - - 


<p align="center"><b><a class="btn" href="https://genomicsaotearoa.github.io/RNA-seq-workshop/" style="background: var(--bs-dark);font-weight:bold">Back to homepage</a></b></p>
