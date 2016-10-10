suppressPackageStartupMessages(suppressWarnings(suppressMessages(library("Biostrings"))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library("tools"))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library("ShortRead"))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library("stringi"))))

#select start and end positions
str<- 9478
end<- 9498
#set the mismatch tolerance for the bowtie alligner
mismatch<- 0
#read the mismatch list created in .txt
ref_list_01_12<- read.delim("ref_list_01_12.txt", header = T, stringsAsFactors = F )
#read the original wildtype reference
replicon_ref<- list.files(".", pattern ="replicon", all.files = F, full.names = F)
replicon_str<- (toString(readBStringSet(replicon_ref)))
siRNA_ref<- subseq((replicon_str), start = str, end = end )
#input the data in .fastq or .fastq.gz format
input_data<- list.files(".", pattern="fastq", all.files = F, full.names = F)
counter<- 100
#substitute the wildtype with the mismatch reference 
#and run the script for the specific reference in a loop
for (i in ref_list_01_12$target) {
  #set reference and output prefix name
  prefix <-i
  #read the mismatch and substitute the original sequence with the mismach sequence
  sub_ref<- stri_replace_all_fixed(replicon_str, pattern = siRNA_ref, replacement = i )
  new_fasta_ref<- paste(prefix,"_ref.fasta", sep = "")
  #write the new reference in .fasta format
  writeFasta (DNAStringSet(sub_ref), new_fasta_ref, mode="w" )

  #input the reference sequence in .fasta format
  mm_ref<- list.files(".", pattern =prefix, all.files = F, full.names = F)
  #reading and transforming reference sequence
  ref_str<- strsplit((toString(readBStringSet(mm_ref))), NULL ,fixed = T )
  #build the index
  CMD_index<- paste("bowtie-build -q -f", mm_ref, "index", sep=" ")
  system(CMD_index)
  
  #set output names
  filename <- paste("mm", mismatch, sep="")  
  out_name <- paste("read_count_", filename, ".txt", sep="")
  csv_name <- paste(counter, "_", filename, ".csv", sep="")
  #perform alignment with bowtie and read count using bedtools
  CMD_bow<- paste("bowtie -p 2 -S -k 1 -v", mismatch, "index", input_data," | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
  system(CMD_bow)
  #read and merge ref and reads
  reads<- read.delim(out_name, header = F, stringsAsFactors = F )
  dataframe<- data.frame(reads, ref_str , stringsAsFactors = F )
  #calculating the % and log10 columns
  dataframe[,5] <- (dataframe[,3]/sum(dataframe[,3])*100)
  dataframe[,6] <- (log10(dataframe[,3]))
  #write CSV table
  write.table(dataframe, file = csv_name , sep = "\t",col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ),row.names = F )
  counter<- counter+1
  #remove read_count and index file created
  counts<- list.files(".", pattern="read_count", all.files = F, full.names = F)
  for(h in counts ){
    h<- paste("rm", h , sep=" ")
    system(h)
  }
  refs<- list.files(".", pattern=new_fasta_ref, all.files = F, full.names = F)
  for(h in refs ){
    h<- paste("rm", h , sep=" ")
    system(h)
  }
  index_files<- list.files(".", pattern="index", all.files = F, full.names = F)
  for(h in index_files ){
    h<- paste("rm", h , sep=" ")
    system(h)
  }
  rm(reads)
  rm(sub_ref)
  rm(mm_ref)
  rm(ref_str)
  rm(refs)
}
#read all .csv files
csv_list<- list.files(".", pattern = ".csv", all.files = F, full.names = F)
csv_N_name <- paste("100_", filename, ".csv", sep="")
N<- read.csv(csv_N_name, header = T, sep = "\t", quote = "" )

#create csv data frame with all the log values
CSV_df_log <- N[6]
a=2
for (i in csv_list ) {
  log_col<- read.csv(i, header = T, sep = "\t", quote = "" )
  CSV_df_log[,a]<- data.frame(log_col[6])
  a=a+1
}
#create csv data frame with all the linear values
CSV_df_linear <- N[5]
a=2
for (i in csv_list ) {
  linear_col<- read.csv(i, header = T, sep = "\t", quote = "" )
  CSV_df_linear[,a]<- data.frame(linear_col[5])
  a=a+1
}
CSV_df_linear[1]<- NULL
CSV_df_log[1]<- NULL
nucleotides<- data.frame(N[,4])
#name the mismatches columns acording to position and nt tranformation
colnames(nucleotides) <- "nucleotide"
colnames(CSV_df_log) <- c(ref_list_01_12$name)
colnames(CSV_df_linear) <- c(ref_list_01_12$name)
#merge and write csv tables in log
CSV_log<- data.frame(nucleotides, CSV_df_log)
write.csv(CSV_log, "CSV_log10.csv" )
CSV_br_log <- CSV_log[str:end,]
write.csv(CSV_br_log, file = "CSV_log10_region.csv")
#merge and write csv tables in linear
CSV_linear<- data.frame(nucleotides, CSV_df_linear)
write.csv(CSV_linear, "CSV_linear.csv" )
CSV_br_linear <- CSV_linear[str:end,]
write.csv(CSV_br_linear, "CSV_linear_region.csv")
#delete csv files
csv_patern<- paste(filename, ".csv", sep="")
csv_files<- list.files(".", pattern=csv_patern, all.files = F, full.names = F)
for(h in csv_files ){
  h<- paste("rm", h , sep=" ")
  system(h)
}

