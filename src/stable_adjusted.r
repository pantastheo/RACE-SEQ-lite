#Function to download, install and load the required libraries only when needed
packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}
#List of required libraries to be loaded
suppressMessages(packages(Biostrings))

#set output name prefix
prefix<- "your_file_name_here"

#select start and end positions
str<- 9478
end<- 9498

#set the mismatch tolerance for the bowtie alligner
mismatch<- 0

#read the original wildtype reference
replicon_ref<- list.files(".", pattern =".fasta", all.files = F, full.names = F)

#read nucleotide reference and convert to character string
nt <- data.frame(lapply(replicon_ref, function(v) {
  if (is.character(v)) return(strsplit(toupper((toString(readBStringSet(replicon_ref)))), NULL , fixed = T))
  else return(v)
}), stringsAsFactors = F)

#input the data in .fastq or .fastq.gz format
input_data<- list.files(".", pattern="fastq", all.files = F, full.names = F)

#build the index
CMD_index<- paste("bowtie-build -q -f", replicon_ref, "index", sep=" ")
system(CMD_index)

#set output names
filename <- paste(mismatch, "mm", sep="")
out_name <- paste("read_count_", filename, ".txt", sep="")
csv_name <- paste("CSV_",prefix, "_", filename, ".csv", sep="")
csv_br_name <- paste("CSV_",prefix, "_", filename, "_region.csv", sep="")

#perform alignment with bowtie and read count using bedtools
CMD_bow<- paste("bowtie -p 2 -S -k 1 -v", mismatch, "index", input_data," | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
system(CMD_bow)

#read and merge ref and reads
reads<- read.delim(out_name, header = F )
dataframe<- data.frame(reads, nt_reference , stringsAsFactors = F)

#calculating the % and log10 columns
dataframe[,5] <- (dataframe[,3]/sum(dataframe[,3])*100)
dataframe[,6] <- (log10(dataframe[,3]))
dataframe[dataframe== -Inf] <-0

write.table(dataframe, file = csv_name , sep = "\t",col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ),row.names = F )
write.table(dataframe[str:end,], file = csv_br_name , sep = "\t",col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ),row.names = F )

#remove read_count and index file created
counts<- list.files(".", pattern="read_count", all.files = F, full.names = F)
for(i in counts ){
  i<- paste("rm", i , sep=" ")
  system(i)
}
index_files<- list.files(".", pattern="index", all.files = F, full.names = F)
for(i in index_files ){
  i<- paste("rm", i , sep=" ")
  system(i)
}


