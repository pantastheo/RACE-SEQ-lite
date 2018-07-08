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
prefix<- "bowtie"

#select start and end positions
str<- 9478
end<- 9498

#set the mismatch tolerance for the bowtie alligner
mismatch<- 0

#read .fasta reference and .fastq data
replicon_ref<- list.files(".", pattern ="fasta", all.files = F, full.names = F)
input_data<- list.files(".", pattern="fastq", all.files = F, full.names = F)

#reading and transforming reference sequence
nt_reference <-strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)

#set output name
filename <- paste("mm", mismatch, sep = "")
out_name <- paste("read_count_", filename, sep="")

#build the index and perform alignment with bowtie and read count using bedtools
  
print(paste0("Performing alignment with ", mismatch, " mismatch using bowtie"))

CMD_bowindex<- paste("bowtie-build -q -f", replicon_ref, "index", sep=" ")
system(CMD_bowindex)

CMD_bow<- paste("bowtie -p 2 -S -k 1 -v", mismatch, "index", input_data," | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
system(CMD_bow)

#read and merge ref and reads
reads<- read.delim(out_name, header = F )
dataframe<- data.frame(reads, nt , stringsAsFactors = F)

#calculating the % and log10 columns
dataframe[,5] <- (dataframe[,3]/sum(dataframe[,3])*100)
dataframe[,6] <- (log10(dataframe[,3]))
dataframe[dataframe== -Inf] <-0

#focusing on target region can be ajusted acording to experiment
binding_region <- dataframe[str:end,]

#write.table(dataframe, file = paste0(prefix, "_", filename, "_whole.csv") , sep = "\t",col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ),row.names = F )
write.table(binding_region, file = paste0(prefix, "_", filename, ".csv") , sep = "\t",
            col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ), row.names = F )

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


