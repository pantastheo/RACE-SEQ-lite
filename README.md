# rnai_mismatch_tolerance

This is a custom R script for the downstream analysis of RACE-seq data.

UNIX/LINUX package dependencies:
	- CUTADAPT
	- BOWTIE
	- SAMTOOLS
	- BEDTOOLS

HOW TO CONFIGURE:
- Open the script and edit the variables acording to your RNA binding region, adapter sequence and other preferences.

- Set the start and finish nucleotide positions in between which the graphs will be plotted. 
	- To do that edit the "str" and "end" variables and set them to the specific nucleotide positions in the reference sequence you provide.

- Input the specific adapter nucleotide sequence which needs to be trimmed from the reads before alignment.
	- To do that edit the "RACE_adapter" variable and paste the specific nucleotide sequence of your adapter.

- Set the name of the output files and graphs that will be created.
	- To do that edit the "prefix" variable according to your prefered file name.

- Set the mismatch tolerence that you want the BOWTIE aligner to perform alignment with.
	- To do that edit the "mismatch" variable and set it to a numeric value, preferably between "0" and "3"

HOW TO RUN:
- To run the analysis you will need to have in the same folder these specific files:
	- The R script downloded and edited according to your options
	- Input dataset in .FASTQ format
	- Reference sequence in .FASTA format
  
- To run the script you can use either run it using the Rscript command from the command line/terminal or open it in Rstudio and run it from there.
