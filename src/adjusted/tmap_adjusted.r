
#Function to download, install and load the required libraries only when needed
packages <- function(x) {
  x <- as.character(match.call()[[2]])
  if (!require(x, character.only = TRUE)) {
    install.packages(pkgs = x, repos = "http://cran.r-project.org")
    require(x, character.only = TRUE)
  }
}
#List of required libraries to be loaded
suppressMessages(packages(Biostrings))

replicon_ref<- list.files(".", pattern ="fasta", all.files = F, full.names = F)
input_data<- list.files(".", pattern="fastq", all.files = F, full.names = F)

#reading and transforming reference sequence
nt_reference <-strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)


#set output name
filename <- paste("mm", mismatch, sep = "")
out_name <- paste("read_count_", filename, sep="")
  
#build the index and perform alignment with tmap and read count using bedtools
  
print(paste0("Performing alignment with ", mismatch, " mismatch using tmap"))

#build tmap index
CMD_tmapindex<- paste("tmap index -f", replicon_ref , sep=" ")
system(CMD_tmapindex)

#align using tmap with 
# map1 :mapping procedure #1 (bwa-short variant)
# -a 0 :unique best hits
# -g 3 :no softclipping
# --maz-mismatches :<mismatch>

CMD_tmap<- paste("tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", replicon_ref," -r ", input_data, " | samtools view -bt replicon.fasta - | genomeCoverageBed -d -5 -ibam stdin > ",out_name, sep="")
system(CMD_tmap)

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

write.table(binding_region, file = paste0(prefix, "_", filename, ".csv") , sep = "\t",
            col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ),
            row.names = F )

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
