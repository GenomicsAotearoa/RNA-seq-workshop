## How to act on fastq after QC.

We can do several trimming:

  * on quality using Phred score: we want an accuracy of 99%. What will be the Phred score?
  * on the sequences, if they contain adaptor sequences.

To do so, we can use on tools: cutadapt.

```
cutadapt -q 20 -a file:/mnt/RNAseq_Workshop_Data/QC/adaptors.fa -A file:/mnt/RNAseq_Workshop_Data/QC/adaptors.fa -o trimmed/WT1_R1_trimmed.fastq -p trimmed/WT1_R2_trimmed.fastq /mnt/RNAseq_Workshop_Data/sequencing/WT1_R1.fastq.gz /mnt/RNAseq_Workshop_Data/sequencing/WT1_R2.fastq.gz

```

Now we should trim all samples.

```
for i in `seq 1 3`; 
do
# WT first:
cutadapt -q 20 -a file:/mnt/RNAseq_Workshop_Data/QC/adaptors.fa -A file:/mnt/RNAseq_Workshop_Data/QC/adaptors.fa -o trimmed/WT$i\_R1_trimmed.fastq -p trimmed/WT$i\_R2_trimmed.fastq /mnt/RNAseq_Workshop_Data/sequencing/WT$i\_R1.fastq.gz /mnt/RNAseq_Workshop_Data/sequencing/WT$i\_R2.fastq.gz

# then comes the KOs:
cutadapt -q 20 -a file:/mnt/RNAseq_Workshop_Data/QC/adaptors.fa -A file:/mnt/RNAseq_Workshop_Data/QC/adaptors.fa -o trimmed/KO$i\_R1_trimmed.fastq -p trimmed/KO$i\_R2_trimmed.fastq /mnt/RNAseq_Workshop_Data/sequencing/KO$i\_R1.fastq.gz /mnt/RNAseq_Workshop_Data/sequencing/KO$i\_R2.fastq.gz

done
```
