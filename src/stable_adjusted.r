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
suppressMessages(packages(tools))
suppressMessages(packages(stringi))

#set output name prefix
prefix<- "your_file_name_here"

#select start and end positions
str<- 9478
end<- 9498

#set the mismatch tolerance for the bowtie alligner
mismatch<- 0

#read the original wildtype reference
replicon_ref<- list.files(".", pattern =".fasta", all.files = F, full.names = F)

nt_reference<- strsplit((toString(readBStringSet(replicon_ref))), NULL ,fixed = T )
nt_reference<- data.frame(lapply(nt_reference, function(v) {
  if (is.character(v)) return(toupper(v))
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
dataframe<- read.delim(out_name, header = F, stringsAsFactors = F )
dataframe[,4] <- (dataframe[,3])
dataframe[,3] <- (nt_reference)

#calculating the % and log10 columns
dataframe[,5] <- (dataframe[,4]/sum(dataframe[,4])*100)
dataframe[,6] <- (log10(dataframe[,4]))

dataframe[,1]<- NULL
dataframe[dataframe== -Inf] <-0
colnames(dataframe) <- c("position","nucleotide","counts", "linear", "log10")

write.table(dataframe, csv_name, sep = "\t", quote = F, row.names = F)
write.table(dataframe[str:end,], csv_br_name, sep = "\t", quote = F, row.names = F)

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


