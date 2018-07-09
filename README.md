## RACE-SEQ lite

This is a custom R script for the downstream analysis of RACE-seq data.

## Dependencies

The following command line toold must be installed on your UNIX machine.

Package            | Version
--------------     | --------------
Cutadapt (optional)| >=1.12
Bowtie             | 1.0.0
Samtools           | >=1.3.1
Bedtools           | >=2.17.0
Tmap (optional)    | >=3.4.1

## Downloads

#### Cutadapt 

Cutdapt is a tool writen in Python that is used for the trimmming of specific nucleotides sequences and adapters from high-throughput sequencing reads.
Cutadapt helps with these trimming tasks by finding the adapter or primer sequences in an error-tolerant way. It can also modify and filter reads in various ways. 

The easiest way to download and install Cutadapt system-wide is using the <b>pip</b> command.
To install <b>pip</b> :
```
sudo apt-get install python-pip
```
Then to install Cutadapt :
```
sudo pip install cutadapt
```
For more information and troubleshooting please read Cutadapt documentation: (https://cutadapt.readthedocs.io/)

#### Bowtie 

Bowtie is an ultrafast, memory-efficient short read aligner. It aligns short DNA sequences to the human genome at a rate of over 25 million 35-bp reads per hour.

To download and install Bowtie system-wide:
```
sudo apt-get install bowtie
```
For specific versions and installation instructions please visit: (http://bowtie-bio.sourceforge.net/index.shtml)

#### Samtools

Samtools is a set of tools for high-throughput sequencing data manipulation. It also contains the libraries for several sequence and variant calling data formats.
To download and install Samtolls sustem-wide:
```
sudo apt-get install samtools
```
For documentation please visit: (http://samtools.github.io/)

#### Bedtools

Bedtools is a powerful toolset for genome arithmetic and contains many utilities for a wide-range of genomics analysis tasks.

To download and install bedtools system-wide :
```
sudo apt-get install bedtools
```
For documetation please visit: (http://bedtools.readthedocs.io/en/latest/)

#### Tmap

TMAP is a fast and accurate alignment software for short and long nucleotide sequences produced by next-generation sequencing technologies such as Ion Torrent PGM and Proton sequencers from ThermoFisher Scientific.
The latest TMAP is unsupported. To use a supported version, please see the TMAP version associated with a Torrent Suite.
*The download and use of the Tmap aligner is optional.*

To download and install system wide please follow the instructions on their official repository:
(https://github.com/iontorrent/TS/tree/master/Analysis/TMAP)

#### R
To run the RACE-SEQ-lite script you will need to have an R version ( >=3.4) installed on your system. All the necessary R packages will get installed and loaded automatically.

## Scripts

There is one main script for the RACE-seq analysis

 - *RACEseqMM.r* is an optimized R script that combines and configures all the necessary command line packages in a complete pipeline solution. The script will get all the configuration options as command line parameters parsed with Rscript.

In the *legacy_scripts* folder you can find the vanilla RACE-SEQ-lite script *stable.r*.
 - *stable.r* is a simple R script that uses the Cutadapt, Bowtie, Samtools and Bedtools packages and the configuration options need to be edited manually in the file. 

## Usage 

To use the **stable.r** you will need to have in the same directory the following files:
- The RACE-SEQ-lite *stable.r* script downloded and edited.
- Only one sequencing data file in .fastq format.
- Only one reference sequence file in .fasta format.

Open the *stable.r* using a text editor of your choice or Rstudio and edit the script variables according to your preferences:
- Set the start (str) and finish (end) variables to the specific nucleotide positions of your reference in between which the graphs will be plotted. 
- Input the specific adapter (RACE_adapter) nucleotide sequence which needs to be trimmed from the reads before alignment.
- Set the file name (prefix) of the output and graph that will be created.
- Set the mismatch tolerence (mismatch) during Bowtie alignment, set preferably between "0" and "3"

To run the *stable.r* you can either use Rscript from the command line or run it from within Rstudio.

## Usage 

To use the **RACEseqMM.r** you will need to run it using Rscript from the command line
```
Rscript RACEseqMM.r
```
To get information and help about the script use option *--help* or *-h*
```
Rscript RACEseqMM.r --help
```

Options           | Description
----------------- | --------------
-s, --start       | "Input the start nucleotide position"
-e, --end         | "Input the end nucleotide position"
-a, --adapter     | "Input the RACE adapter sequence [default: NO_trimming]"
-m, --mismatch    | "Input number of mismatches during alignement [default: %default]"
-p, --plot        | "Print output graph plot [default: NO]"
-t, --tmap        | "Use tmap aligner instead of bowtie [default: NO]"
--nocsv           | "Do not print output CSV file [default: NO]"
-i, --iterate     | "Create the alternative references to cover all the possible SNPs of the reference between the genomics locations specified by -s and -e and then perform global alignment and generate graph.[default: NO]"
-h, --help        | "Print this help page"              




