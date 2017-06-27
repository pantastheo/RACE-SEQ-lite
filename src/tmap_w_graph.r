
#List of required libraries to be loaded
suppressMessages(library(Biostrings))
suppressMessages(library(optparse))

option_list<- list(
  make_option(c("-s", "--start"), type="integer", default = NA, 
              help="Input the start nucleotide position"),
  
  make_option(c("-e", "--end"), type="integer", default = NA, 
              help="Input the end nucleotide position"),
  
  make_option(c("-m", "--mismatch"), type="integer", default = 0, 
              help="Input number of mismatches during alignement [default: %default]"),
  
  make_option(c("-a", "--adapter"), action="store", default=NA, type="character", 
              help="Input RACE adapter sequence [default: %default]") ,
  
  make_option(c("-r", "--reference"), action="store", default= NA, type="character", 
              help="Input .fasta reference file [default: .fasta]") ,
  
  make_option(c("-i", "--input"), action="store", default= NA , type="character", 
              help="Input .fastaq file [default: .fastq]") ,
  
  make_option(c("-o", "--output"), action="store", default="substitutions.txt", type="character", 
              help="Output file [default: %default]"),
  
  make_option(c("-g", "--graph"), action="store_true", default=TRUE,
              help="Print output graph [default: %default]"),
  
  make_option(c("-t", "--tmap"), action="store_false", default=FALSE,
              help="Use tmap aligner instead of bowtie [default: bowtie]")
  
)

opt = parse_args(OptionParser(usage = "Usage: %prog [options] -r <fasta> -i <fastq> \n", option_list = option_list,
                              add_help_option = TRUE, 
                              description = "This is a custom R script for the downstream analysis of RACE-seq data.", 
                              epilogue = "Thank you for using RACE-SEQ lite"))

print(opt$s)
print(opt$e)
print(opt$m)
print(opt$a)
print(opt$r)
print(opt$i)
print(opt$o)
print(opt$g)

str<- opt$s
end<- opt$e

if(!is.na(opt$s) & !is.na(opt$e)) {
  #input the reference sequence in .fasta format
  refname<- list.files(".", pattern ="fasta", all.files = F, full.names = F)
  if (length(refname)==0) {
    stop("No input .fasta reference file available")
  } else if (is.na(opt$r)) {
    replicon_ref<- opt$r
  } else {
    replicon_ref<- refname
  }
  #input the data in .fastq or .fastq.gz format
  data_fastq<- list.files(".", pattern="fastq", all.files = F, full.names = F)
  if (length(data_fastq)==0) {
    stop("No input .fastq file available")
  } else if (is.na(opt$i)) {
    replicon_ref<- opt$i
  } else {
    input_data<- data_fastq
  }
  #set output name prefix
  output_name<- "your_data"
  if (opt$o!=output_name) {
    prefix<- opt$o
  } else {
    prefix<-output_name
  }
  
}else {stop("Please input Start and End nucleotide positions \n Or type [option] -h for help")}

#reading and transforming reference sequence
nt_reference <-strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)

#set output names
filename <- paste("mm", mismatch, sep="")
out_name <- paste("read_count_", filename, ".txt", sep="")
csv_name <- paste(prefix, "_", filename, ".csv", sep="")
csv_br_name <- paste(prefix, "_", filename, "_br.csv", sep="")


if (opt$t==TRUE){
  #build the index and perform alignment with tmap and read count using bedtools
  
  CMD_tmapindex<- paste("tmap index -f", mm_ref , sep=" ")
  system(CMD_tmapindex)
  
  CMD_tmap<- paste("tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", replicon_ref," -r ", input_data, " | samtools view -bt replicon.fasta - | genomeCoverageBed -d -5 -ibam stdin > ",out_name, sep="")
  system(CMD_tmap)
} else {
  
  #build the index and perform alignment with bowtie and read count using bedtools
  CMD_bowindex<- paste("bowtie-build -q -f", replicon_ref, "index", sep=" ")
  system(CMD_bowindex)
  
  CMD_bow<- paste("bowtie -p 2 -S -k 1 -v", mismatch, "index", input_data," | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
  system(CMD_bow)
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
write.table(binding_region, file = csv_br_name , sep = "\t",col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ),row.names = F )

#remove read_count and index file created
counts<- list.files(".", pattern="read_count", all.files = F, full.names = F)
for(i in counts ){
  i<- paste("rm", i , sep=" ")
  system(i)
}
index_files<- list.files(".", pattern=".fasta.tmap.", all.files = F, full.names = F)
for(i in index_files ){
  i<- paste("rm", i , sep=" ")
  system(i)
}
index_files<- list.files(".", pattern="aligned.bam", all.files = F, full.names = F)
for(i in index_files ){
  i<- paste("rm", i , sep=" ")
  system(i)
}
index_files<- list.files(".", pattern="out.sam", all.files = F, full.names = F)
for(i in index_files ){
  i<- paste("rm", i , sep=" ")
  system(i)
}
index_files<- list.files(".", pattern="index", all.files = F, full.names = F)
for(i in index_files ){
  i<- paste("rm", i , sep=" ")
  system(i)
}

if (opt$g==TRUE){
#create wildtype linear & log scale graph
pdf(paste0(filename, "_graph.pdf"), width=15)
    
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
