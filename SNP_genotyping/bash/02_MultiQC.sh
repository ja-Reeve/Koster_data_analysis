#!/bin/bash

### Generating read quality reports for 99 Littorina saxatilis samples that
### were sequenced by BGI with 5X coverage. Reports are generated for each
### sample using MultiQC. This job was run on the server 'Rackham', part of
### the National Academic Infrastructure for Supercomputing in Sweden.
### James Reeve - Göteborgs Universitet
### 30/06/2021

### Job parameters
#SBATCH -A snic2022-5-266
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 3:00:00
#SBATCH -J KT_MultiQC

### Required software
module load bioinfo-tools
module load FastQC/0.11.9
module load MultiQC/1.10

### Command
fastqc --dir $TMPDIR -t 20 /file/path/trimmed_read/*.fastq.gz --outdir /file/path/fastqc
multiqc /file/path/fastqc -o /file/path/MultiQC
