
time_start<- Sys.time()
print("operation started on")
print(time_start)

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
suppressMessages(packages(ShortRead))
suppressMessages(packages(stringi))
suppressMessages(packages(ggplot2))
suppressMessages(packages(gridExtra))
suppressMessages(packages(grid))


#####set parameters and read input files####

#select start and end positions
str<- 9478
end<- 9498
#set the mismatch tolerance for the bowtie alligner
mismatch<- 0
#read the mismatch list created in .txt
reference_list<- list.files(".", pattern ="ref_list", all.files = F, full.names = F)
ref_list_01_12<- read.table(reference_list, sep = "\t", quote = "\"", header = T, stringsAsFactors = F )
#read the original wildtype reference
replicon_ref<- list.files(".", pattern ="replicon", all.files = F, full.names = F)
replicon_str<- (toString(readBStringSet(replicon_ref)))

nt_reference<- strsplit((toString(readBStringSet(replicon_ref))), NULL ,fixed = T )
nt_reference<- data.frame(lapply(nt_reference, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}), stringsAsFactors = F)

dataframe_counts<- data.frame(nt_reference, stringsAsFactors = F)
dataframe_log10<- data.frame(nt_reference , stringsAsFactors = F)
dataframe_linear<- data.frame(nt_reference , stringsAsFactors = F)

siRNA_ref<- subseq((replicon_str), start = str, end = end )
#input the data in .fastq or .fastq.gz format
input_data<- list.files(".", pattern="fastq", all.files = F, full.names = F)
counter<- 100

#####run the stable.r script for mismatch references in loop#####

for (i in ref_list_01_12$target) {
  
  print(paste("working on file", counter-99 , "of", length(ref_list_01_12$target)))
  #reference and output prefix name
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
  ref_str<- data.frame(lapply(ref_str, function(v) {
   if (is.character(v)) return(toupper(v))
   else return(v)
  }))
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
  
  dataframe_counts<- data.frame(dataframe_counts, reads[,3], stringsAsFactors = F)
  dataframe_log10<- data.frame(dataframe_log10, (log10(reads[,3])), stringsAsFactors = F)
  dataframe_linear<- data.frame(dataframe_linear, (reads[,3]/sum(reads[,3])*100), stringsAsFactors = F)
  
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
  counter<- counter+1
  
}

#####create CSV files####
nt_reference<- data.frame(nt_reference, check.rows = T)
dataframe_counts[1]<- NULL
dataframe_log10[1]<- NULL
dataframe_linear[1]<- NULL

#name the mismatches columns acording to position and nt tranformation
colnames(nt_reference) <- "nucleotide"
colnames(dataframe_counts)<- c(ref_list_01_12$name)
colnames(dataframe_log10) <- c(ref_list_01_12$name)
colnames(dataframe_linear) <- c(ref_list_01_12$name)

#merge and write csv tables in log and linear
CSV_log<- data.frame(nt_reference, dataframe_log10, check.names = T, check.rows = T)
CSV_linear<- data.frame(nt_reference, dataframe_linear, check.names = T, check.rows = T)
#set -inf values to 0.00 so there is no problem plotting
CSV_log[CSV_log== -Inf] <- 0
#write.table(CSV_log, "CSV_log10.csv", sep = "\t", quote = F, row.names = F)
#write.table(CSV_linear, "CSV_linear.csv" , sep = "\t", quote = F, row.names = F)
write.table(CSV_log[str:end,], "CSV_log10_region.csv", sep = "\t", quote = F, row.names = F)
write.table(CSV_linear[str:end,], "CSV_linear_region.csv", sep = "\t", quote = F, row.names = F)

mm0<- data.frame(nt_reference[str:end,], dataframe_counts[str:end,1], dataframe_linear[str:end,1], dataframe_log10[str:end,1] , stringsAsFactors = F)
mm0[mm0== -Inf] <-0
colnames(mm0) <- c("nucleotide","counts", "linear", "log10")
write.table(mm0, "mm0_wildtype.csv", sep = "\t", quote = F, row.names = F)

##### create wildtype linear & log scale graph #####
pdf("mm0_graph.pdf", width=15)
#in 100% linear scale
mp <- barplot(mm0[,3],
              ylab="Novel 5\' Ends (%)",
              xlab="Binding site", 
              names.arg=mm0[,1], 
              las=1, 
              cex.names = 2, 
              col="darkgrey" , 
              main="Novel 5' Ends in linear",
              cex.main=2.7,
              cex.lab=1.2,
              ylim=c(0,100))
text(mp,mm0[,3]+5 ,cex = 1.3, adj = 0 ,labels=mm0[,2] ,srt=90)
#in log10  logarithmic scale
mp <- barplot(mm0[,4],
              ylab=expression("Novel 5\' Ends (log"[10]*")"),
              xlab="Binding site", 
              names.arg=mm0[,1], 
              las=1, 
              cex.names = 2, 
              col="darkgrey", 
              main="Novel 5' Ends in logarithmic",
              cex.main=2.7,
              cex.lab=1.2,
              ylim=c(0,10))
text(mp,mm0[,4]+0.5 ,cex = 1.3, adj = 0 ,labels=mm0[,2] ,srt=90)
dev.off()
##end of wildtype graph##

##### ggplot2 graph#####


#ref_list_01_12$index<- c(0,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12)

# "G"<- "black"
# "A"<- "green"
# "T"<- "red"
# "C"<- "blue"

N<- CSV_log[str:end,]

N$nt2<- N[,1]
N[20,36]<- NA
N$nt3<- N[,1]
N[19,37]<- NA
N$nt4<- N[,1]
N[18,38]<- NA
N$nt5<- N[,1]
N[17,39]<- NA
N$nt6<- N[,1]
N[16,40]<- NA
N$nt7<- N[,1]
N[15,41]<- NA
N$nt8<- N[,1]
N[14,42]<- NA
N$nt9<- N[,1]
N[13,43]<- NA
N$nt10<- N[,1]
N[12,44]<- NA
N$nt11<- N[,1]
N[11,45]<- NA
N$nt12<- N[,1]
N[10,46]<- NA

N$legend_mix<- as.factor(N[,1]) 
levels(N$legend_mix)<- c(levels(N$legend_mix), "wildtype")

write.table(N, file = "graph_log_table.csv", sep = "\t", quote = F, row.names = F)

setEPS()
postscript("seed_region.eps", width=10)

si02<- ggplot(N, aes(siRNA_00_N, x = seq_along(siRNA_00_N)))+
  ylim(-0.4,8)+
  xlab("Binding site")+
  ylab(expression("Novel 5\' Ends (log"[10]*")"))+
  geom_bar(stat = "identity", colour="black",fill="white", width = 0.4) +
  geom_line(aes(y=N$siRNA_02_g), colour="black", size= 1) + 
  geom_line(aes(y=N$siRNA_02_c), colour="blue", size=1) + 
  geom_line(aes(y=N$siRNA_02_t), colour="red", size=1)+
  scale_x_discrete()+
  geom_text(data = NULL, x=c(1:21), y=-0.4, label= N$nt2)+
  coord_cartesian(xlim = c(1,21))+
  geom_text(data = NULL, x=20, y=-0.4, label= N[20,1], size=5)


si03<- ggplot(N, aes(siRNA_00_N, x = seq_along(siRNA_00_N)))+
  ylim(-0.4,8)+
  xlab("Binding site")+
  ylab(expression("Novel 5\' Ends (log"[10]*")"))+
  geom_bar(stat = "identity", colour="black",fill="white", width = 0.4) + 
  geom_line(aes(y=N$siRNA_03_g), colour="black", size= 1) + 
  geom_line(aes(y=N$siRNA_03_c), colour="blue", size=1) + 
  geom_line(aes(y=N$siRNA_03_t), colour="red", size=1)+
  scale_x_discrete()+
  geom_text(data = NULL, x=c(1:21), y=-0.4, label= N$nt3)+
  coord_cartesian(xlim = c(1,21))+
  geom_text(data = NULL, x=19, y=-0.4, label= N[19,1], size=5)

si04<- ggplot(N, aes(siRNA_00_N, x = seq_along(siRNA_00_N)))+
  ylim(-0.4,8)+
  xlab("Binding site")+
  ylab(expression("Novel 5\' Ends (log"[10]*")"))+
  geom_bar(stat = "identity", colour="black",fill="white", width = 0.4) + 
  geom_line(aes(y=N$siRNA_04_a), colour="green", size= 1) + 
  geom_line(aes(y=N$siRNA_04_g), colour="black", size=1)+
  geom_line(aes(y=N$siRNA_04_t), colour="red", size=1) + 
  scale_x_discrete()+
  geom_text(data = NULL, x=c(1:21), y=-0.4, label= N$nt4)+
  coord_cartesian(xlim = c(1,21))+
  geom_text(data = NULL, x=18, y=-0.4, label= N[18,1], size=5)

si05<- ggplot(N, aes(siRNA_00_N, x = seq_along(siRNA_00_N)))+
  ylim(-0.4,8)+
  xlab("Binding site")+
  ylab(expression("Novel 5\' Ends (log"[10]*")"))+
  geom_bar(stat = "identity", colour="black",fill="white", width = 0.4) + 
  geom_line(aes(y=N$siRNA_05_a), colour="green", size= 1) + 
  geom_line(aes(y=N$siRNA_05_g), colour="black", size=1)+
  geom_line(aes(y=N$siRNA_05_t), colour="red", size=1) + 
  scale_x_discrete() +
  geom_text(data = NULL, x=c(1:21), y=-0.4, label= N$nt5)+
  coord_cartesian(xlim = c(1,21))+
  geom_text(data = NULL, x=17, y=-0.4, label= N[17,1], size=5)

si06<- ggplot(N, aes(siRNA_00_N, x = seq_along(siRNA_00_N)))+
  ylim(-0.4,8)+
  xlab("Binding site")+
  ylab(expression("Novel 5\' Ends (log"[10]*")"))+
  geom_bar(stat = "identity", colour="black",fill="white", width = 0.4) + 
  geom_line(aes(y=N$siRNA_06_a), colour="green", size=1)+
  geom_line(aes(y=N$siRNA_06_g), colour="black", size= 1) + 
  geom_line(aes(y=N$siRNA_06_c), colour="blue", size=1) + 
  scale_x_discrete() +
  geom_text(data = NULL, x=c(1:21), y=-0.4, label= N$nt6)+
  coord_cartesian(xlim = c(1,21))+
  geom_text(data = NULL, x=16, y=-0.4, label= N[16,1], size=5)

si07<- ggplot(N, aes(siRNA_00_N, x = seq_along(siRNA_00_N)))+
  ylim(-0.4,8)+
  xlab("Binding site")+
  ylab(expression("Novel 5\' Ends (log"[10]*")"))+
  geom_bar(stat = "identity", colour="black",fill="white", width = 0.4) + 
  geom_line(aes(y=N$siRNA_07_a), colour="green", size=1) + 
  geom_line(aes(y=N$siRNA_07_g), colour="black", size= 1) + 
  geom_line(aes(y=N$siRNA_07_t), colour="red", size=1)+
  scale_x_discrete() +
  geom_text(data = NULL, x=c(1:21), y=-0.4, label= N$nt7)+
  coord_cartesian(xlim = c(1,21))+
  geom_text(data = NULL, x=15, y=-0.4, label= N[15,1], size=5)

grid.arrange(si02, si03, si04, si05, si06, si07)

dev.off()

setEPS()
postscript("cleavage_region.eps", width=10)

si09<- ggplot(N, aes(siRNA_00_N, x = seq_along(siRNA_00_N)))+
  ylim(-0.4,8)+
  xlab("Binding site")+
  ylab(expression("Novel 5\' Ends (log"[10]*")"))+
  geom_bar(stat = "identity", colour="black",fill="white", width = 0.4) +
  geom_line(aes(y=N$siRNA_09_a), colour="green", size=1)+
  geom_line(aes(y=N$siRNA_09_g), colour="black", size=1) + 
  geom_line(aes(y=N$siRNA_09_t), colour="red", size= 1) + 
  scale_x_discrete() +
  geom_text(data = NULL, x=c(1:21), y=-0.4, label= N$nt9)+
  coord_cartesian(xlim = c(1,21))+
  geom_text(data = NULL, x=13, y=-0.4, label= N[13,1], size=5)

si10<- ggplot(N, aes(siRNA_00_N, x = seq_along(siRNA_00_N)))+
  ylim(-0.4,8)+
  xlab("Binding site")+
  ylab(expression("Novel 5\' Ends (log"[10]*")"))+
  geom_bar(stat = "identity", colour="black",fill="white", width = 0.4) +
  geom_line(aes(y=N$siRNA_10_a), colour="green", size= 1) + 
  geom_line(aes(y=N$siRNA_10_g), colour="black", size=1) + 
  geom_line(aes(y=N$siRNA_10_c), colour="blue", size=1)+
  scale_x_discrete() +
  geom_text(data = NULL, x=c(1:21), y=-0.4, label= N$nt10)+
  coord_cartesian(xlim = c(1,21))+
  geom_text(data = NULL, x=12, y=-0.4, label= N[12,1], size=5)

si11<- ggplot(N, aes(siRNA_00_N, x = seq_along(siRNA_00_N)))+
  ylim(-0.4,8)+
  xlab("Binding site")+
  ylab(expression("Novel 5\' Ends (log"[10]*")"))+
  geom_bar(stat = "identity", colour="black",fill="white", width = 0.4) +
  geom_line(aes(y=N$siRNA_11_a), colour="green", size=1) + 
  geom_line(aes(y=N$siRNA_11_g), colour="black", size= 1) + 
  geom_line(aes(y=N$siRNA_11_t), colour="red", size=1)+
  scale_x_discrete() +
  geom_text(data = NULL, x=c(1:21), y=-0.4, label= N$nt11)+
  coord_cartesian(xlim = c(1,21))+
  geom_text(data = NULL, x=11, y=-0.4, label= N[11,1], size=5)

si12<- ggplot(N, aes(siRNA_00_N, x = seq_along(siRNA_00_N)))+
  ylim(-0.4,8)+
  xlab("Binding site")+
  ylab(expression("Novel 5\' Ends (log"[10]*")"))+
  geom_bar(stat = "identity", colour="black",fill="white", width = 0.4) +
  geom_line(aes(y=N$siRNA_12_g), colour="black", size=1) + 
  geom_line(aes(y=N$siRNA_12_c), colour="blue", size= 1) + 
  geom_line(aes(y=N$siRNA_12_t), colour="red", size=1)+
  scale_x_discrete() +
  geom_text(data = NULL, x=c(1:21), y=-0.4, label= N$nt12)+
  coord_cartesian(xlim = c(1,21))+
  geom_text(data = NULL, x=10, y=-0.4, label= N[10,1], size=5)

grid.arrange(si09, si10, si11, si12)

dev.off()

time_finish<- Sys.time()
print("operation finished on")
print(time_start)
time_duration<- time_finish - time_start 
print(time_duration)





