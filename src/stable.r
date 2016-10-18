suppressPackageStartupMessages(suppressWarnings(suppressMessages(library("Biostrings"))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library("tools"))))

#set output name prefix
prefix<- "ENTER_NAME_PREFIX_HERE"
#set the mismatch tolerance for alignment 
mismatch<- 0
#select start and end positions of the graph
str<- 9095
end<- 9117
#input the data in .fastq or .fastq.gz format
input_data<- list.files(".", pattern="fastq", all.files = F, full.names = F)
#input the reference sequence in .fasta format
mm_ref<- list.files(".", pattern ="fasta", all.files = F, full.names = F)
#reading and transforming reference sequence
ref_str<- strsplit((toString(readBStringSet(mm_ref))), NULL ,fixed = T )
#build the index
CMD_index<- paste("bowtie-build -q -f", mm_ref, "index", sep=" ")
system(CMD_index)
#set output names
filename <- paste("mm", mismatch, sep="")
out_name <- paste("read_count_", filename, ".txt", sep="")
pdf_name <- paste(prefix, "_", filename, ".pdf", sep="")
csv_name <- paste(prefix, "_", filename, ".csv", sep="")
csv_br_name <- paste(prefix, "_", filename, "_br.csv", sep="")
#perform adapter trimming with cutadapt, alignment with bowtie and read count using bedtools
CMD_bow<- paste("cutadapt -gGGACACTGACATGGACTGAAGGAGTAGAAA -e0 --no-indels -m10 --discard-untrimmed", input_data,"|bowtie -p 2 -S -k 1 -v", mismatch, "index - | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
system(CMD_bow)
#read and merge ref and reads
reads<- read.delim(out_name, header = F )
dataframe<- data.frame(reads, ref_str )
#calculating the % and log10 columns
dataframe[,5] <- (dataframe[,3]/sum(dataframe[,3])*100)
dataframe[,6] <- (log10(dataframe[,3]))
#write CSV table
write.table(dataframe, file = csv_name , sep = "\t",col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ),row.names = F )
#focusing on target region can be ajusted acording to experiment
binding_region <- dataframe[str:end,]
names <- as.vector(binding_region[,4])
write.table(binding_region, file = csv_br_name , sep = "\t",col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ),row.names = F )
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
#ploting the tables
pdf(pdf_name, width=20)
#in 100%
mp <- barplot(binding_region[,5],ylab="Novel 5\' Ends (%)" ,xlab="Binding site", names.arg=names , las=2 , cex.names = 2, col="darkgrey" , main=filename ,cex.main=2.7,cex.lab=1.2,ylim=c(0,100))
text(mp,binding_region[,5]+5,cex = 1.3, adj = 0 ,labels=binding_region[,3] ,srt=90)
#in log10
mp <- barplot(binding_region[,6],ylab=("Novel 5\' Ends (log10)") ,xlab="Binding site", names.arg=names , las=2 , cex.names = 2 , col="darkgrey" , main=filename ,cex.main=2.7,cex.lab=1.2,ylim=c(0,10))
text(mp,binding_region[,6]+0.5 ,cex = 1.3, adj = 0 ,labels=binding_region[,3] ,srt=90)
dev.off()


