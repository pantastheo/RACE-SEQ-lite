#Function to download, install and load the required libraries only when needed
packages <- function(x) {
  x <- as.character(match.call()[[2]])
  if (!require(x, character.only = TRUE)) {
    install.packages(pkgs = x, repos = "http://cran.r-project.org")
    require(x, character.only = TRUE)
  }
}
#load optparser package first so as to print msg and read vars
suppressMessages(packages(optparse))


#List of pasring options loaded
option_list<- list(
  make_option(c("-s", "--start"), type="integer", action = "store", default = NA,
              help="Input the first nucleotide reference genomic location."),

  make_option(c("-e", "--end"), type="integer", action = "store" ,default = NA,
              help="Input the last nucleotide reference genomic location."),

  make_option(c("-a", "--adapter"), type="character", action = "store" ,default = NA,
              help="Input the RACE adapter nucleotide sequence. [default: NO]"),

  make_option(c("-m", "--mismatch"), type="integer", default = 0,
              help="Input number of mismatches that are allowed during alignement. [default: %default]"),

  make_option(c("-p", "--plot"), action="store_true", default = FALSE,
              help="Print output graph between the specified genomic locations. [default: NO]"),

  make_option(c("-t", "--tmap"), action="store_true", default = FALSE,
              help="Use the Tmap aligner instead of Bowtie. [default: NO]"),
  
  make_option(c("--notsv"), action="store_true", default = FALSE ,
              help="Do not write output tsv file. [default: NO]"),
  
  make_option(c("-i", "--iterate"), action="store_true", default= FALSE,
              help="Create the alternative references to cover all the possible SNPs \n
              of the reference between the genomics locations specified by -s and -e\n
              and then perform global alignment and generate graph.[default: NO] ")
  
)

opt = parse_args(OptionParser(description = "RACE-SEQ-lite\n
This is a custom R script for the downstream analysis of RACE-seq data.\n
The pipeline uses common bioinformatics command line packages such as BOWTIE, SAMTOOLS and BEDTOOLS that should be installed system-wide.\n
The script reads the necessary input files from your workind directory and outputs a graph or a tsv file.\n
One fasta and one fastq file can only be in the working directory",
                              usage = "Rscript %prog [options] -s <integer> -e <integer> -m <integer> \n",
                              option_list = option_list,
                              add_help_option = TRUE,
                              epilogue = "Thank you for using RACE-SEQ lite.
                              \nFor documentation visit: https://github.com/pantastheo/RACE-SEQ-lite.
                              \nAuth: Pantazis Theotokis 2018
                              \nContact: p.theotokis@imperial.ac.uk
                              \n"))




str<- opt$s
end<- opt$e
mismatch<- opt$mismatch
RACE_adapter<- opt$a

#example values to be used
#str<- 9478
#end<- 9498
#RACE_adapter<- "GGACACTGACATGGACTGAAGGAGTAGAAA"

#Check if the necessary fasta and fastq files are located in the working directory.
#If TRUE load read the files.
#If FALSE exit with a message.
if(!is.na(opt$s) & !is.na(opt$e)) {
  #input the reference sequence in .fasta format
  reference<- list.files(".", pattern ="fasta", all.files = F, full.names = F)
  print(paste0("Reading reference file in fasta format."))
  print(paste0("File " ,reference, " found in working directory."))
  
    if ((length(reference))==0) {
      stop("No input .fasta reference file available in working directory.")} else if ((length(reference))>=2) {
      stop("More than one .fasta reference file in working directory.")
    }

  #input the data in .fastq or .fastq.gz format
  input_data<- list.files(".", pattern="fastq", all.files = F, full.names = F)
  print(paste0("Reading data input file in fastq format."))
  print(paste0("File " ,input_data, " found in working directory."))
    if ((length(input_data))==0) {
      stop("No input .fastq file available in working directory.")
    } else if ((length(input_data))>=2) {
      stop("More than one .fastq file in working directory. \nFor paired end reads please concatenate and run again")
    } 
} else 
  stop("Please input Start and End nucleotide reference genomic locations \nOr type [option] -h for help")

print(paste0("Loading R packages"))



#List of required libraries to be loaded
suppressMessages(source("https://bioconductor.org/biocLite.R"))
suppressMessages(if (!"ShortRead" %in% installed.packages()) biocLite(ShortRead))
suppressMessages(if (!"Biostrings" %in% installed.packages()) biocLite(Biostrings))
#List of required libraries to be loaded
suppressMessages(packages(tools))
suppressMessages(packages(ShortRead))
suppressMessages(packages(Biostrings))
suppressMessages(packages(stringi))
suppressMessages(packages(ggplot2))
suppressMessages(packages(cowplot))


#If iterate option is TRUE run the script that will generate all the alternative references and will perform alignment.
if (opt$i==FALSE){
  
#reading and transforming reference sequence
nt_reference <-strsplit((toString(readBStringSet(reference))), NULL , fixed = T)
nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)

#set output names
input_name<- file_path_sans_ext(input_data)
filename <- paste("mm", mismatch, sep = "")
out_name <- paste("read_count_", filename, sep="")

#If tmap TRUE perform alignmment using TMAP
if (opt$t==TRUE){
  #check if the tmap aligner is installed
  if (system("which tmap")==0) {
  print(paste0("Generating tmap index files"))
  prefix<-"tmap"
  #build the index
  CMD_tmapindex<- paste("tmap index -f", reference , sep=" ")
  system(CMD_tmapindex)

  #perform alignment with tmap and read count using bedtools
  if (is.na(opt$a)){
    
    #no adapter trimming
    print(paste0("Performing alignment with ", mismatch, " mismatch using tmap"))
    CMD_tmap<- paste("tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", reference," -r ", input_data, " | samtools view -bt ", reference," - | genomeCoverageBed -d -5 -ibam stdin > ",out_name, sep="")
    system(CMD_tmap)
  } else {
    #Check if cutadapt installed and perform adapter trimming
    if (system("which cutadapt")==0){
    #adapter trimming using cutadapt
    print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using tmap"))
    CMD_tmap<- paste("cutadapt -g ", RACE_adapter, " -e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data," |tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", reference," -i fastq | samtools view -bt ", reference," - | genomeCoverageBed -d -5 -ibam stdin > ",out_name, sep="")
    system(CMD_tmap) 
    } else {
      stop("Cutadapt software is not installed or not in $PATH. Please see documentation for installation.")}
    }
  } else {
    stop("Tmap software is not installed or not in $PATH. Please see documentation for installation.")}
} 
else {
  #If tmap FALSE perform alignment using bowtie, which is the default option.
  #check if bowtie aligner is installed
  if (system("which bowtie")==0) {
  print(paste0("Generating bowtie index files"))
  prefix<-"bowtie"
  #build the index 
  CMD_bowindex<- paste("bowtie-build -q -f", reference, "index", sep=" ")
  system(CMD_bowindex)

  #perform alignment with bowtie and read count using bedtools
  if (is.na(opt$a)){
    #no adapter trimming
    print(paste0("Performing alignment with ", mismatch, " mismatch using bowtie"))
    CMD_bow<- paste("bowtie -p 4 -S -k 1 -v", mismatch, "index", input_data," | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
    system(CMD_bow)

  } else {
    if (system("which cutadapt")==0) {
    #adapter trimming using cutadapt
    print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using bowtie"))
    CMD_bow<- paste("cutadapt -g", RACE_adapter, "-e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data,"|bowtie -p 4 -S -k 1 -v", mismatch, "index - | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
    system(CMD_bow)
    }else {
      stop("Cutadapt software is not installed or not in $PATH. Please see documentation for installation.")}
    }
  } else {
    stop("Bowtie software is not installed or not in $PATH. Please see documentation for installation.")}}

#read and merge ref and reads
reads<- read.delim(out_name, header = F )
dataframe<- data.frame(reads, nt_reference , stringsAsFactors = F)
    
#calculating the % and log10 columns
dataframe[,5] <- (dataframe[,3]/sum(dataframe[,3])*100)
dataframe[,6] <- (log10(dataframe[,3]))
dataframe[dataframe== -Inf] <-0

#focusing on target region can be ajusted acording to experiment
binding_region <- dataframe[str:end,]

#function to delete files created
del_files<- function(pattern){
  fl_rm<-list.files(".", pattern = pattern, all.files = F, full.names = F)
  for(i in fl_rm){
    i<-paste("rm", i , sep = " ")
    system(i)
  }
}
#delete generated files using fuunction
del_files("read_count")
del_files("fasta.tmap.")
del_files("aligned.bam")
del_files("out.sam")
del_files("index")

#print the wildtype alignment in tsv format table
if (opt$notsv==FALSE){
print(paste0("Writing results to output tsv"))
write.table(binding_region, file = paste0(input_name,"_",prefix, "_", filename, ".tsv") , sep = "\t",
            col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ),
            row.names = F )
}

if (opt$p==TRUE){
#create wildtype linear & log scale graph
pdf(paste0(input_name,"_",prefix, "_", filename, ".pdf"), width=20)
print(paste0("Generating graph with ",mismatch," mismatches."))
#in 100% linear scale
mp <- barplot(binding_region[,5],
              xlab="Binding site", 
              names.arg=(binding_region[,4]), 
              las=1, 
              cex.names = 2.2, 
              col="darkgrey" , 
              main="Novel 5' Ends in linear",
              cex.main=2.3,
              cex.lab=1.5,
              ylim=c(0,100))
title(ylab="Novel 5\' Ends (%)", line=2, cex.lab=1.5)
text(mp,binding_region[,5]+5 ,cex = 1.3, adj = 0 ,labels=binding_region[,3] ,srt=90)
    
#in log10  logarithmic scale
mp <- barplot(binding_region[,6],
              xlab="Binding site", 
              names.arg=(binding_region[,4]), 
              las=1, 
              cex.names = 2.2, 
              col="darkgrey", 
              main="Novel 5' Ends in logarithmic",
              cex.main=2.3,
              cex.lab=1.5,
              ylim=c(0,10))
title(ylab=expression("Novel 5\' Ends (log"[10]*")"), line=2, cex.lab=1.5)
text(mp,binding_region[,6]+0.5 ,cex = 1.3, adj = 0 ,labels=binding_region[,3] ,srt=90)

dev.off()
}

} else {
  
  #Use the combined script to iterate the all the possible SNPs in the reference sequence between the specified genomic locations
  #This is a massive code repetition that needs tidying up.
    
  #reading and transforming reference sequence
  nt_reference <-strsplit((toString(readBStringSet(reference))), NULL , fixed = T)
  nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)
    
  #set output names
  input_name<- file_path_sans_ext(input_data)
  filename <- paste("mm", mismatch, sep = "")
  out_name <- paste("read_count_", filename, sep="")  
    
  #read the original wildtype reference
  replicon_str <- (toString(readBStringSet(reference)))
  
  #read and transform the reference
  ref_replace <- function(str, end, reference) {
    nt_reference <-strsplit((toString(readBStringSet(reference))), NULL , fixed = T)
    nt <-data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)
    
    
    a = 0
    count_names <- data.frame(NA)
    nt_sub<- data.frame(NA)
    for (i in nt_reference[str:end]) {
      if (nt[str + a, 1] == "A") {
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "C"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "C", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"C"
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "T"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "T", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"T"
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "G"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "G", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"G"
        a = a + 1
      }
      if (nt[str + a, 1] == "C") {
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "A"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "A", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"A"
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "T"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "T", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"T"
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "G"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "G", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"G"
        a = a + 1
      }
      if (nt[str + a, 1] == "T") {
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "A"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "A", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"A"
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "C"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "C", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"C"
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "G"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "G", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"G"
        a = a + 1
      }
      if (nt[str + a, 1] == "G") {
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "A"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "A", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"A"
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "T"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "T", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"T"
        nt[, (ncol(nt) + 1)] <- nt[, 1]
        nt[str + a, ncol(nt)] <- "C"
        count_names[nrow(count_names) + 1,] <-
          paste((str + a), nt[(str + a), 1], "to", "C", sep = "_")
        nt_sub[nrow(nt_sub) + 1,] <-"C"
        a = a + 1
      }
    }
    
    
    names_list <-
      data.frame(count_names[((((end - str) + 1) * 3) + 1):1,], stringsAsFactors = F)
    names_list[(nrow(names_list)),] <- "wildtype_siRNA"
    
    nt_sub_list <-
      data.frame(nt_sub[((((end - str) + 1) * 3) + 1):1,], stringsAsFactors = F)
    #nt_sub_list[(nrow(nt_sub_list)),] <- "substitution"
    
    ref_list <- nt[str:end, 1:((((end - str) + 1) * 3) + 1)]
    trans_list <-
      as.data.frame(t(ref_list[, ncol(ref_list):1]), stringsAsFactors = F)
    
    for (i in (c(1:nrow(names_list)))) {
      names_list[i, 2] <-
        paste(trans_list[i, 1:(ncol(trans_list))],sep = "", collapse = "")
    }
    names_list<- cbind(names_list, nt_sub_list)
    colnames(names_list) <- c("name", "target", "nucleotide")
    names_list <-
      rbind((names_list[(nrow(names_list)),]), names_list)
    names_list <- names_list[-c(nrow(names_list)),]
    return(names_list)
  }
  
  #call the ref_replace function
  target <- ref_replace(str, end, reference)
  
  #transform the reference and create the output tables
  dataframe_counts <- data.frame(nt_reference, stringsAsFactors = F)
  dataframe_log10 <- data.frame(nt_reference , stringsAsFactors = F)
  dataframe_linear <- data.frame(nt_reference , stringsAsFactors = F)
  
  #extract the wildtype siRNA sequence
  siRNA_ref <- subseq((replicon_str), start = str, end = end)
  
  #check if packages are installed on system
  #If not exit with a message
  if (system("which bowtie")==0) {
    print("Bowtie aligner installed and on $PATH")} 
  else {stop("Bowtie software is not installed or not on $PATH. Please see documentation for installation.")}
  if (opt$t==TRUE ){
    if (system("which tmap")==0) {
      print("Tmap aligner installed and on $PATH")}
    else if (system("which tmap")==1) {
      print(paste0("Tmap software is not installed or not on $PATH."))
      stop("Please see documentation for installation.")}
    }
  if (!is.na(opt$a)){
    if (system("which cutadapt")==0) {
      print("Cutadapter trimmer installed and on $PATH")}
    else if (system("which Cutadapt")==1){
      print(paste0("Cutadapt trimmer is not installed or not on $PATH."))
      stop("Please see documentation for installation.")}
    }
  
  
  
  #set the counter
  counter <- 100
  #run the script for mismatch references in loop
  for (i in target$target) {
    print(paste("working on alignment",counter - 99 ,"of",length(target$target)))
    
    #reference and output prefix name
    prefix <- i
    
    #read the mismatch and substitute the original sequence with the mismach sequence
    sub_ref <-stri_replace_all_fixed(replicon_str, pattern = siRNA_ref, replacement = i)
    new_fasta_ref <- paste(prefix, "_ref.fasta", sep = "")
    
    #write the new reference in .fasta format
    writeFasta (DNAStringSet(sub_ref), new_fasta_ref, mode = "w")
    
    #input the reference sequence in .fasta format
    mm_ref <-list.files(".",pattern = prefix,all.files = F,full.names = F)
    
    #reading and transforming reference sequence
    ref_str <-strsplit((toString(readBStringSet(mm_ref))), NULL , fixed = T)
    ref_str <- data.frame(lapply(ref_str, function(x) toupper(x)), stringsAsFactors = F)
    
      if (opt$t==TRUE){
      
        prefix<-"tmap"
        #build the index
        CMD_tmapindex<- paste("tmap index -f", mm_ref , sep=" ")
        system(CMD_tmapindex)
        
        #perform alignment with tmap and read count using bedtools
        if (is.na(opt$a)){
          #no adapter trimming
          print(paste0("Performing alignment with ", mismatch, " mismatch using tmap"))
          CMD_tmap<- paste("tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", mm_ref," -r ", input_data, " | samtools view -bt ", mm_ref," - | genomeCoverageBed -d -5 -ibam stdin > ",out_name, sep="")
          system(CMD_tmap)
        } else {
            #adapter trimming using cutadapt
            print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using tmap"))
            CMD_tmap<- paste("cutadapt -g ", RACE_adapter, " -e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data," |tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", mm_ref," -i fastq | samtools view -bt ", mm_ref," - | genomeCoverageBed -d -5 -ibam stdin > ",out_name, sep="")
            system(CMD_tmap)}
        } else {
        
        prefix<-"bowtie"
        #build the index 
        CMD_bowindex<- paste("bowtie-build -q -f", mm_ref, "index", sep=" ")
        system(CMD_bowindex)
        
        #perform alignment with bowtie and read count using bedtools
        if (is.na(opt$a)){
          #no adapter trimming
          print(paste0("Performing alignment with ", mismatch, " mismatch using bowtie"))
          CMD_bow<- paste("bowtie -p 2 -S -k 1 -v", mismatch, "index", input_data," | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
          system(CMD_bow)
           } else {
            #adapter trimming using cutadapt
            print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using bowtie"))
            CMD_bow<- paste("cutadapt -g", RACE_adapter, "-e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data,"|bowtie -p 8 -S -k 1 -v", mismatch, "index - | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
            system(CMD_bow)}
      } 
    #read and merge ref and reads
    reads <- read.delim(out_name, header = F, stringsAsFactors = F)
    dataframe_counts <-data.frame(dataframe_counts, reads[, 3], stringsAsFactors = F)
    dataframe_log10 <-data.frame(dataframe_log10, (log10(reads[, 3])), stringsAsFactors = F)
    dataframe_linear <-data.frame(dataframe_linear, (reads[, 3] / sum(reads[, 3]) * 100), stringsAsFactors = F)
    
    #function to delete read_count and index file created
    del_files<- function(pattern){
      fl_rm<-list.files(".", pattern = pattern, all.files = F, full.names = F)
      for(i in fl_rm){
        i<-paste("rm", i , sep = " ")
        system(i)
      }
    }
    #delete files generated during alignment
    del_files("read_count")
    del_files("index")
    del_files(new_fasta_ref)
    
    #remove R environment variables
    rm(reads)
    rm(sub_ref)
    rm(mm_ref)
    rm(ref_str)
    rm(refs)
    counter <- counter + 1
  }
  
  #create CSV files#
  nt_reference <- data.frame(nt_reference, check.rows = T)
  dataframe_counts[1] <- NULL
  dataframe_log10[1] <- NULL
  dataframe_linear[1] <- NULL
  
  #name the mismatches columns acording to position and nt tranformation
  colnames(nt_reference) <- "nucleotide"
  colnames(dataframe_counts) <- target$name
  colnames(dataframe_log10) <- target$name
  colnames(dataframe_linear) <- target$name
  
  #merge and write csv tables in log and linear
  CSV_log <-data.frame(nt_reference,dataframe_log10,check.names = T,check.rows = T)
  CSV_linear <-data.frame(nt_reference,dataframe_linear,check.names = T,check.rows = T)
  CSV_log[CSV_log == -Inf] <- 0
  
  #focusing on target region can be ajusted acording to experiment
  N <- CSV_log[str:end,]
  
  #create and write wildtype table in .csv
  binding_region <-data.frame(nt_reference[str:end,],
                   dataframe_counts[str:end, 1],
                   dataframe_linear[str:end, 1],
                   dataframe_log10[str:end, 1],
                   stringsAsFactors = F)
  binding_region[binding_region == -Inf] <- 0
  colnames(binding_region) <- c("nucleotide", "counts", "linear", "log10")
  
  #write log10 siRNA region in .tsv
  if (opt$notsv==FALSE){
    write.table(binding_region, file = paste0(input_name,"_", filename, ".tsv"),sep = "\t",quote = F,row.names = F)
    #binding_region <- read.table("siRNA22_1_mm0.tsv", sep = "\t", quote = F, row.names = F, header = T)
  }
  
  if (opt$p==TRUE){
    #create wildtype linear & log scale graph
    pdf(paste0(input_name,"_",prefix, "_", filename, ".pdf"), width=20)
    print(paste0("Generating wildtype graph with ",mismatch," mismatches."))
    #in 100% linear scale
    mp <- barplot(binding_region[,3],
                  xlab="Binding site", 
                  names.arg=(binding_region[,1]), 
                  las=1, 
                  cex.names = 2.2, 
                  col="darkgrey" , 
                  main="Novel 5' Ends in linear",
                  cex.main=2.3,
                  cex.lab=1.5,
                  ylim=c(0,100))
    title(ylab="Novel 5\' Ends (%)", line=2, cex.lab=1.5)
    text(mp,binding_region[,3]+5 ,cex = 1.3, adj = 0 ,labels=binding_region[,2] ,srt=90)
    
    #in log10  logarithmic scale
    mp <- barplot(binding_region[,4],
                  xlab="Binding site", 
                  names.arg=(binding_region[,1]), 
                  las=1, 
                  cex.names = 2.2, 
                  col="darkgrey", 
                  main="Novel 5' Ends in logarithmic",
                  cex.main=2.3,
                  cex.lab=1.5,
                  ylim=c(0,10))
    title(ylab=expression("Novel 5\' Ends (log"[10]*")"), line=2, cex.lab=1.5)
    text(mp,binding_region[,4]+0.5 ,cex = 1.3, adj = 0 ,labels=binding_region[,2] ,srt=90)
    
    dev.off()
  }
  
  
  #ggplot2 graph function
  print(paste0("Generating multiplot with ",mismatch," mismatches."))
  
  plot_RACEseq<- function(N){
    
    N[,1]<- as.character(N[,1])
    values = c("C"="blue", "G"="black", 
               "A"="green", "T"= "red", "black"="black")
    
    ggraph4 <- function(N, a, b) {
      
      gg<- ggplot(data = N, aes(x = seq_along(N[,1]))) +
        geom_col(aes(y = N[,2]), size=1.4,width = 0.4, fill= "lightgrey", colour="black") +
        geom_line(aes(y = N[,a], colour = (strsplit((colnames(N)), "_")[[a]][[4]])), size=1.4) +
        geom_line(aes(y = N[,(a+1)], colour = (strsplit((colnames(N)), "_")[[(a+1)]][[4]])), size=1.4) +
        geom_line(aes(y = N[,(a+2)], colour = (strsplit((colnames(N)), "_")[[(a+2)]][[4]])), size=1.4) +
        scale_colour_manual("Wildtype",
                            breaks = c("C", "G", "A", "T", "black"),
                            values = values)+
        xlab("Binding Site") +
        scale_x_discrete() +
        coord_cartesian(xlim = c(1, (nrow(N))) ) +
        scale_y_continuous(expression("Novel 5\' Ends (log"[10] * ")"), limits = c(-0.3,8) ) + 
        geom_text(data = NULL,x = c(1:(nrow(N))),y = -0.4,label = N[, 1], size=5) +
        geom_text(data = NULL,x = (((nrow(N)) + 1) - b),y = -0.4,label = N[(((nrow(N)) + 1) - b), 1],size = 5, colour = "orange2") +
        theme(legend.position = c(.95, .95),
              legend.justification = c("right", "top"),
              
              legend.margin = margin(1, 1, 1, 1))+
        labs(title=paste0("Nucleotide position ", b))
      
      return(gg)
    }
    
    si01 <- ggraph4(N, 3, 1)
    si02 <- ggraph4(N, 6, 2)
    si03 <- ggraph4(N, 9, 3)
    si04 <- ggraph4(N, 12, 4)
    si05 <- ggraph4(N, 15, 5)
    si06 <- ggraph4(N, 18, 6)
    si07 <- ggraph4(N, 21, 7)
    si08 <- ggraph4(N, 24, 8)
    si09 <- ggraph4(N, 27, 9)
    si10 <- ggraph4(N, 30, 10)
    si11 <- ggraph4(N, 33, 11)
    si12 <- ggraph4(N, 36, 12)
    si13 <- ggraph4(N, 39, 13)
    si14 <- ggraph4(N, 43, 14)
    si15 <- ggraph4(N, 45, 15)
    
    plot2_8<- plot_grid(si02, si03, si04, si05, si06, si07, labels = "AUTO", label_size = 14 , hjust = 0, ncol = 2)
    ggsave(filename=paste0(input_name,"_multiplot_seed_mm",mismatch,".pdf"), plot = plot2_8, scale = 1.8 )
    
    plot9_12<- plot_grid(si09, si10, si11, si12, labels = "AUTO", label_size = 14 , hjust = 0)
    ggsave(filename=paste0(input_name,"_multiplot_cleavage_mm",mismatch,".pdf"), plot = plot9_12, scale = 0.8, width= 20.4 ,height = 7.8 )
  }
  plot_RACEseq(N)
}

