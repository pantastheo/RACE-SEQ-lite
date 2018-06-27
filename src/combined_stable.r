
#Function to download, install and load the required libraries only when needed
packages <- function(x) {
  x <- as.character(match.call()[[2]])
  if (!require(x, character.only = TRUE)) {
    install.packages(pkgs = x, repos = "http://cran.r-project.org")
    require(x, character.only = TRUE)
  }
}
#List of required libraries to be loaded
suppressMessages(packages(optparse))
suppressMessages(packages(Biostrings))

option_list<- list(
  make_option(c("-s", "--start"), type="integer", action = "store", default = NA,
              help="Input the first nucleotide reference genomic location."),

  make_option(c("-e", "--end"), type="integer", action = "store" ,default = NA,
              help="Input the last nucleotide reference genomic location."),

  make_option(c("-a", "--adapter"), type="character", action = "store" ,default = NA,
              help="Input the RACE adapter nucleotide sequence. [default: no_trimming]"),

  make_option(c("-m", "--mismatch"), type="integer", default = 0,
              help="Input number of mismatches that are allowed during alignement. [default: %default]"),

  make_option(c("-p", "--plot"), action="store_true", default = FALSE,
              help="Print output graph between the specified genomic locations. [default: NO]"),

  make_option(c("-t", "--tmap"), action="store_true", default = FALSE,
              help="Use the Tmap aligner instead of Bowtie. [default: BOWTIE]"),
  
  make_option(c("--no_csv"), action="store_true", default = FALSE ,
              help="Do not print output CSV file. [default: <filename.csv>]")
  
)

opt = parse_args(OptionParser(description = "RACE-SEQ-lite 
                              \nThis is a custom R script for the downstream analysis of RACE-seq data. 
                              \nThe pipeline uses common bioinformatics command line packages such as BOWTIE, SAMTOOLS and BEDTOOLS that should be installed system-wide.
                              \nThe script reads the necessary input files from your workind directory and outputs a graph or a csv file. 
                              \nOne fasta and one fastq file can only be in the working directory",
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

if(!is.na(opt$s) & !is.na(opt$e)) {
  #input the reference sequence in .fasta format
  refname<- list.files(".", pattern ="fasta", all.files = F, full.names = F)
    if ((length(refname))==0) {
      stop("No input .fasta reference file available in working directory.")
    } else if ((length(refname))>=2) {
      stop("More than one .fasta reference file in working directory.")
    } else if ((length(refname))==1) {
      reference<- as.character(refname)
    } 

  #input the data in .fastq or .fastq.gz format
  data_fastq<- list.files(".", pattern="fastq", all.files = F, full.names = F)
    if ((length(data_fastq))==0) {
      stop("No input .fastq file available in working directory.")
    } else if ((length(data_fastq))>=2) {
      stop("More than one .fastq file in working directory. For paired end reads please concatenate and run again")
    } else if ((length(data_fastq))==1) {
      input_data<- data_fastq
    } 
}else {stop("Please input Start and End nucleotide reference genomic locations \nOr type [option] -h for help")}

#reading and transforming reference sequence
nt_reference <-strsplit((toString(readBStringSet(reference))), NULL , fixed = T)
nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)

#set output names
filename <- paste("mm", mismatch, sep = "")
out_name <- paste("read_count_", filename, sep="")

if (opt$t==TRUE){
  if (system("which tmap",show.output.on.console = F  )==0) {
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
    if (system("which cutadapt",show.output.on.console = F )==0){
    #adapter trimming using cutadapt
    print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using tmap"))
    CMD_tmap<- paste("cutadapt -g ", RACE_adapter, " -e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data," |tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", reference," -i fastq | samtools view -bt ", reference," - | genomeCoverageBed -d -5 -ibam stdin > ",out_name, sep="")
    system(CMD_tmap) 
    } else {
      stop("Cutadapt software is not installed or not in $PATH. Please see documentation for installation.")
    }
  }
  
  } else {
    stop("Tmap software is not installed or not in $PATH. Please see documentation for installation.")
  }
  
} else {
  
  if (system("which bowtie",show.output.on.console = F )==0) {
  prefix<-"bowtie"
  #build the index 
  CMD_bowindex<- paste("bowtie-build -q -f", reference, "index", sep=" ")
  system(CMD_bowindex)

  #perform alignment with bowtie and read count using bedtools
  if (is.na(opt$a)){
    #no adapter trimming
    print(paste0("Performing alignment with ", mismatch, " mismatch using bowtie"))
    CMD_bow<- paste("bowtie -p 2 -S -k 1 -v", mismatch, "index", input_data," | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
    system(CMD_bow)
  } else {
    if (system("which cutadapt",show.output.on.console = F )==0) {
    #adapter trimming using cutadapt
    print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using bowtie"))
    CMD_bow<- paste("cutadapt -g", RACE_adapter, "-e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data,"|bowtie -p 8 -S -k 1 -v", mismatch, "index - | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
    system(CMD_bow)
    }else {
      stop("Cutadapt software is not installed or not in $PATH. Please see documentation for installation.")
    }
  }
  } else {
    stop("Bowtie software is not installed or not in $PATH. Please see documentation for installation.")
  }
}

#read and merge ref and reads
reads<- read.delim(out_name, header = F )
dataframe<- data.frame(reads, nt_reference , stringsAsFactors = F)
    
#calculating the % and log10 columns
dataframe[,5] <- (dataframe[,3]/sum(dataframe[,3])*100)
dataframe[,6] <- (log10(dataframe[,3]))
dataframe[dataframe== -Inf] <-0

#focusing on target region can be ajusted acording to experiment
binding_region <- dataframe[str:end,]

#remove read_count and index file created

del_files<- function(pattern){
  fl_rm<-list.files(".", pattern = pattern, all.files = F, full.names = F)
  for(i in fl_rm){
    i<-paste("rm", i , sep = " ")
    system(i)
  }
}

del_files("read_counts")
del_files("fasta.tmap.")
del_files("aligned.bam")
del_files("out.sam")
del_files("index")

if (opt$no_csv==FALSE){
write.table(binding_region, file = paste0(prefix, "_", filename, ".csv") , sep = "\t",
            col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ),
            row.names = F )
}


if (opt$p==TRUE){
#create wildtype linear & log scale graph
pdf(paste0(prefix, "_", filename, ".pdf"), width=15)
    
#in 100% linear scale
mp <- barplot(binding_region[,5],
              ylab="Novel 5\' Ends (%)",
              xlab="Binding site", 
              names.arg=(binding_region[,4]), 
              las=1, 
              cex.names = 2.2, 
              col="darkgrey" , 
              main="Novel 5' Ends in linear",
              cex.main=2.3,
              cex.lab=1.3,
              ylim=c(0,100))
text(mp,binding_region[,5]+5 ,cex = 1.3, adj = 0 ,labels=binding_region[,3] ,srt=90)
    
#in log10  logarithmic scale
mp <- barplot(binding_region[,6],
              ylab=expression("Novel 5\' Ends (log"[10]*")"),
              xlab="Binding site", 
              names.arg=(binding_region[,4]), 
              las=1, 
              cex.names = 2.2, 
              col="darkgrey", 
              main="Novel 5' Ends in logarithmic",
              cex.main=2.3,
              cex.lab=1.3,
              ylim=c(0,10))
text(mp,binding_region[,6]+0.5 ,cex = 1.3, adj = 0 ,labels=binding_region[,3] ,srt=90)
dev.off()
}

