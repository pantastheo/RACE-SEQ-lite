# RACE-SEQ-lite example run

In this folder you will find all the necessary files to perform an example run of the **RACEseqMM.r** script.
To run the example you will need to download or clone the *RACE-SEQ-lite* repository on your working directory.
The example dataset is a fastq file containing short sequencing reads from a RACEseq experiment that have already been trimmed of the RACEseq adapter using the Cutadapt adapter trimming software.  

## Example dataset

The only files you will need to perform an example run are the:
- RACEseqMM.r
- sample_data.fastq
- sample_reference.fasta

## Download 

Clone the Github repository which includes the example dataset on your system

```
cd /path/to/your/home/directory
```
```
git clone https://github.com/pantastheo/RACE-SEQ-lite.git 
```
```
cd src/example_data
```
## Do a test run

First run the Rscript command to check that you can access the file and use R.
```
Rscript RACEseqMM.r --help
```
The help screen will be printed on your terminal.
Now perform a test run using the example dataset provided in the *example_data* folder by using the minimum options.
The file **sampe_data_mm0.tsv** file will be generated.

```
Rscript -s 9478 -e 9498
```
This file is a tab delimited table that gives you the cleavage incidents that occurred in the genomic locations specified between *-s 9478* and *-e 9498* in actual numbers, in percentage and in a logarithmic scale as well.

## Run with options
Now this time run the script but with the *-p* option and the *--notsv* option as well.
Using the *-p* option the script will to generate and print a graph of the cleavage incidents between the genomic locations.
Using the *--notsv* option the script will not output a table file with the cleavage incidents.
The file **sample_data_bowtie_mm0.pdf** file will be generated.

```
Rscript -s 9478 -e 9498 -p --notsv
```
Now run the script but give the option to perform alignment using Bowtie and allow for a 1 mismatch during alignment. You can use this option by allow for one or more mismatches during alignment. 
The file **sample_data_bowtie_mm1.pdf** file will be generated.

```
Rscript -s 9478 -e 9498 -p --notsv -m 1
```

## Run the iterate mode
By invoking the *-i* option will generate multiple reference sequences using all the possible nucleotide combinations found between the genomic locations provided and then try to align the dataset against all these new references. Using this option the default output will be two graphs, one for the seed region and one around the cleavage region of RNAi.   
The files **sample_multiplot_cleavage_mm0.pdf** and **sample_multiplot_seed_mm0.pdf** will be generated.

```
Rscript -s 9478 -e 9498 -i

```
## All outputs
The iterate *-i* option can be used along with generating an output table (without using *--notsv*) and a plot (*-p*). You can also allow for one or more mismatches during this process (-m 1). The command that will output all the possible files is:
```
Rscript -s 9478 -e 9498 -p --notsv -m 1 -p -i 
```






