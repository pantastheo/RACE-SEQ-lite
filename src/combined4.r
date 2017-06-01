
time_start <- Sys.time()
print("operation started on")
print(time_start)

#Function to download, install and load the required libraries only when needed
packages <- function(x) {
  x <- as.character(match.call()[[2]])
  if (!require(x, character.only = TRUE)) {
    install.packages(pkgs = x, repos = "http://cran.r-project.org")
    require(x, character.only = TRUE)
  }
}
#List of required libraries to be loaded
source("https://bioconductor.org/biocLite.R")
suppressMessages(biocLite("Biostrings"))
suppressMessages(biocLite("ShortRead"))

suppressMessages(packages(Biostrings))
suppressMessages(packages(tools))
suppressMessages(packages(ShortRead))
suppressMessages(packages(stringi))
suppressMessages(packages(ggplot2))
suppressMessages(packages(gridExtra))
suppressMessages(packages(grid))


#set parameters and read input files#

#select start and end positions
str <- 9478
end <- 9498
#set the mismatch tolerance for the bowtie alligner
mismatch <- 0

#read the original wildtype reference
replicon_ref <- list.files(".",pattern = "replicon",all.files = F,full.names = F)
replicon_str <- (toString(readBStringSet(replicon_ref)))

#read and transform the reference

ref_replace <- function(str, end, replicon_ref) {
  nt_reference <-
    strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
  nt <- data.frame(lapply(nt_reference, function(v) {
    if (is.character(v))
      return(toupper(v))
    else
      return(v)
  }), stringsAsFactors = F)
  
  
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



target <- ref_replace(str, end, replicon_ref)

#transform the reference and create the output tables#


nt_reference <-strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
nt_reference <- data.frame(lapply(nt_reference, function(v) {
  if (is.character(v))
    return(toupper(v))
  else
    return(v)
}), stringsAsFactors = F)

dataframe_counts <- data.frame(nt_reference, stringsAsFactors = F)
dataframe_log10 <- data.frame(nt_reference , stringsAsFactors = F)
dataframe_linear <- data.frame(nt_reference , stringsAsFactors = F)

siRNA_ref <- subseq((replicon_str), start = str, end = end)
#input the data in .fastq or .fastq.gz format
input_data <-list.files(".",pattern = "fastq",all.files = F,full.names = F)
counter <- 100

#run the stable.r script for mismatch references in loop

for (i in target$target) {
  print(paste("working on alignment",counter - 99 ,"of",length(target$target)
  ))
  #reference and output prefix name
  prefix <- i
  #read the mismatch and substitute the original sequence with the mismach sequence
  sub_ref <-
    stri_replace_all_fixed(replicon_str, pattern = siRNA_ref, replacement = i)
  new_fasta_ref <- paste(prefix, "_ref.fasta", sep = "")
  #write the new reference in .fasta format
  writeFasta (DNAStringSet(sub_ref), new_fasta_ref, mode = "w")
  
  #input the reference sequence in .fasta format
  mm_ref <-list.files(".",pattern = prefix,all.files = F,full.names = F)
  #reading and transforming reference sequence
  ref_str <-strsplit((toString(readBStringSet(mm_ref))), NULL , fixed = T)
  ref_str <- data.frame(lapply(ref_str, function(v) {
    if (is.character(v))
      return(toupper(v))
    else
      return(v)
  }))
  #build the index
  CMD_index <-
    paste("bowtie-build -q -f", mm_ref, "index", sep = " ")
  system(CMD_index)
  
  #set output names
  filename <- paste("mm", mismatch, sep = "")
  out_name <- paste("read_count_", filename, ".txt", sep = "")
  csv_name <- paste(counter, "_", filename, ".csv", sep = "")
  #perform alignment with bowtie and read count using bedtools
  CMD_bow <-
    paste(
      "bowtie -p 4 -S -k 1 -v",
      mismatch,
      "index",
      input_data,
      " | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >",
      out_name,
      sep = " "
    )
  system(CMD_bow)
  #read and merge ref and reads
  reads <- read.delim(out_name, header = F, stringsAsFactors = F)
  
  dataframe_counts <-data.frame(dataframe_counts, reads[, 3], stringsAsFactors = F)
  dataframe_log10 <-data.frame(dataframe_log10, (log10(reads[, 3])), stringsAsFactors = F)
  dataframe_linear <-data.frame(dataframe_linear, (reads[, 3] / sum(reads[, 3]) * 100), stringsAsFactors = F)
  
  counts <-list.files(".",pattern = "read_count",all.files = F,full.names = F)
  for (h in counts) {
    h <- paste("rm", h , sep = " ")
    system(h)
  }
  refs <-list.files(".",pattern = new_fasta_ref,all.files = F,full.names = F)
  for (h in refs) {
    h <- paste("rm", h , sep = " ")
    system(h)
  }
  index_files <-
    list.files(".",pattern = "index",all.files = F,full.names = F)
  for (h in index_files) {
    h <- paste("rm", h , sep = " ")
    system(h)
  }
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

#write log10 siRNA region in .csv
write.table(CSV_log[str:end,],"CSV_log10_region.csv",sep = "\t",quote = F,row.names = F)

#create wildtype table
mm0 <-data.frame(nt_reference[str:end,],
                 dataframe_counts[str:end, 1],
                 dataframe_linear[str:end, 1],
                 dataframe_log10[str:end, 1],
                 stringsAsFactors = F)
mm0[mm0 == -Inf] <- 0
colnames(mm0) <- c("nucleotide", "counts", "linear", "log10")

#write wildtype table in .csv
write.table(mm0,"mm0_wildtype.csv",sep = "\t",quote = F,row.names = F)

#mm0 <- read.csv("mm0_wildtype.csv", sep = "\t")
#log10_region <- read.csv("CSV_log10_region.csv", sep = "\t")

#create wildtype linear & log scale graph
pdf <- paste("mm0", "graph.pdf", sep = "_")

graph <- function(mm0, pdf_name) {
  pdf(pdf_name, width = 15)
  
  #in 100% linear scale
  mp <- barplot(mm0[, 3],
    ylab = "Novel 5\' Ends (%)",
    xlab = "Binding site",
    names.arg = mm0[, 1],
    las = 1,
    cex.names = 2,
    col = "darkgrey" ,
    main = "Novel 5' Ends in linear",
    cex.main = 2.7,
    cex.lab = 1.2,
    ylim = c(0, 100))
  text(mp,mm0[, 3] + 5 ,cex = 1.3,adj = 0 ,labels = mm0[, 2] ,srt = 90)
  
  #in log10  logarithmic scale
  mp <- barplot(mm0[, 4],
    ylab = expression("Novel 5\' Ends (log"[10] * ")"),
    xlab = "Binding site",
    names.arg = mm0[, 1],
    las = 1,
    cex.names = 2,
    col = "darkgrey",
    main = "Novel 5' Ends in logarithmic",
    cex.main = 2.7,
    cex.lab = 1.2,
    ylim = c(0, 10))
  text(mp,mm0[, 4] + 0.5 ,cex = 1.3,adj = 0 ,labels = mm0[, 2] ,srt = 90)
  
  dev.off()
}

graph(mm0, pdf)

#ggplot2 graph function


N <- CSV_log[str:end,]

write.table(N,file = "graph_log_table.csv",sep = "\t",quote = F,row.names = F)
#N<-read.csv("graph_log_table.csv", sep = "\t")

ggraph2 <- function(N, a, b) {
  si02 <- ggplot(N, aes(N[, 2], x = seq_along(N[, 2]))) +
    ylim(-0.4, 8) +
    xlab("Binding site") +
    ylab(expression("Novel 5\' Ends (log"[10] * ")")) +
    geom_bar(stat = "identity",colour = "black",fill = "white",width = 0.4) +
    geom_line(aes(y = N[, a]), colour = "lightgrey", size = 1) +
    geom_line(aes(y = N[, (a + 1)]), colour = "darkgrey", size = 1) +
    geom_line(aes(y = N[, (a + 2)]), colour = "black", size = 1) +
    scale_x_discrete() +
    geom_text(data = NULL,x = c(1:(nrow(N))),y = -0.4,label = N[, 1]) +
    coord_cartesian(xlim = c(1, (nrow(N)))) +
    geom_text(data = NULL,x = (((nrow(N)) + 1) - b),y = -0.4,label = N[(((nrow(N)) + 1) - b), 1],size = 5)
  return(si02)
}
#ggraph2(N,3,1)


gtable <- data.frame(1)
col <- 3
row <- 1
for (i in (c(1:(nrow(N))))) {
  gtable[i,] <- paste("ggraph2(N,", col, ",", row, ")", sep = "")
  col = col + 3
  row = row + 1
  
}

si01 <- ggraph2(N, 3, 1)
si02 <- ggraph2(N, 6, 2)
si03 <- ggraph2(N, 9, 3)
si04 <- ggraph2(N, 12, 4)
si05 <- ggraph2(N, 15, 5)
si06 <- ggraph2(N, 18, 6)
si07 <- ggraph2(N, 21, 7)
si08 <- ggraph2(N, 24, 8)
si09 <- ggraph2(N, 27, 9)
si10 <- ggraph2(N, 30, 10)
si11 <- ggraph2(N, 33, 11)
si12 <- ggraph2(N, 36, 12)
si13 <- ggraph2(N, 39, 13)
si14 <- ggraph2(N, 43, 14)
si15 <- ggraph2(N, 45, 15)

pdf("seed_reg.pdf", width = 13)

grid.arrange(si01, si02, si03, si04, si05, si06, si07, si08, si09)
dev.off()


pdf("cleavage_reg.pdf", width = 13)
grid.arrange(si10, si11, si12, si13, si14, si15)
dev.off()

# setEPS()
# postscript("cleavage_region.eps", width=10)

# grid.arrange(si09, si10, si11, si12)
#
# dev.off()


time_finish <- Sys.time()
print("operation finished on")
print(time_start)
time_duration <- time_finish - time_start
print(time_duration)
