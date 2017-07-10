## RACE-SEQ lite

This is a custom R script for the downstream analysis of RACE-seq data.

## Dependencies

The following programs and packages must be installed on your computer.

Package            | Version
--------------     | --------------
Cutadapt (optional)| >=1.12
Bowtie             | 1.0.0
Samtools           | >=1.3.1
Bedtools           | >=2.17.0
Tmap (optional)    | >=3.4.1

## Scripts

There are two stable scripts for RACE-seq analysis

 - *stable.r* is a simple R script that uses Cutadapt, Bowtie, Samtools and Bedtools and the options need to be edited manually in the file.

 - *combined_stable.r* is an optimized version of the stable.r script that parses command line arguments and incorporates the Tmap aligner. 

## Configuration 

To use the *stable.r* open the file using a text editor of your choice and edit the script variables according to your preferences.

- Set the start (str) and finish (end) variables to the specific nucleotide positions of your reference in between which the graphs will be plotted. 

- Input the specific adapter (RACE_adapter) nucleotide sequence which needs to be trimmed from the reads before alignment.

- Set the file name (prefix) of the output and graph that will be created.

- Set the mismatch tolerence (mismatch) during Bowtie alignment, set preferably between "0" and "3"

## Usage 
You will need to have in the same directory the following files:
- The RACE-SEQ lite *stable.r* script downloded and edited as stated above or the *combined_stable.r* as it is.
- Only one sequencing data file in .fastq format.
- Only one reference sequence file in .fasta format.
  
To run the *stable.r* you can either use Rscript from the command line or run it from within Rstudio.
To run the *combined_stable.r* you will need to run it using Rscript from the command line
