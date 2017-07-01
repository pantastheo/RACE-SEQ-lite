## RACE-SEQ lite

This is a custom R script for the downstream analysis of RACE-seq data.

## Dependencies

The following programs and packages must be installed on your computer.

Package       | Written with version
------------- | --------------------
Cutadapt      | >=1.12
Bowtie        | 1.0.0
Samtools      | >=1.3.1
Bedtools      | >=2.17.0

## Configuration 
Open the *stable.r* using a standart text editor of your choice.
Edit the script variables acording to your preferences and based on the comments provided.

- Set the start and finish nucleotide positions in between which the graphs will be plotted. 
	- To do that edit the "str" and "end" variables and set them to the specific nucleotide positions in the reference sequence you provide.

- Input the specific adapter nucleotide sequence which needs to be trimmed from the reads before alignment.
	- To do that edit the "RACE_adapter" variable and paste the specific nucleotide sequence of your adapter.

- Set the name of the output files and graphs that will be created.
	- To do that edit the "prefix" variable according to your prefered file name.

- Set the mismatch tolerence that you want the BOWTIE aligner to perform alignment with.
	- To do that edit the "mismatch" variable and set it to a numeric value, preferably between "0" and "3"

## Usage 
You will need to have in the same directory the following files:
- The RACE-SEQ lite *stable.r* script downloded and edited as stated above
- Only one sequencing data file in .fastq format
- Only one reference sequence file in .fasta format
  
To run the script you can either use Rscript from the command line or run it from within Rstudio.
