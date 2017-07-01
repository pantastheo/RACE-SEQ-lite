args <- commandArgs(trailingOnly = TRUE)

str = as.integer(args[1])
end = as.integer(args[2])

#start and end nucleotide positions on reference genome
# str<- 9478
# end<- 9498

# replicon_ref<- list.files(".", pattern ="replicon", all.files = F, full.names = F)


ref_replace <- function( str, end, replicon_ref ) {
  #Function to download, install and load the required libraries only when needed
  packages<-function(x){
    x<-as.character(match.call()[[2]])
    if (!require(x,character.only=TRUE)){
      install.packages(pkgs=x,repos="https://bioconductor.org/biocLite.R")
      require(x,character.only=TRUE)
    }
  }
  #List of required libraries to be loaded
  suppressMessages(packages(Biostrings))
  
  #read nucleotide reference and convert to character string
  nt <-strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
  nt<- data.frame(lapply(nt, function(x) toupper(x)), stringsAsFactors = F)
  
  
  #loop the reference for start and end location and create nucleotide substitutions
  a = 0
  rlength<- list(1:nrow(nt))
  count_names <- data.frame(NA)
  nt_sub<- data.frame(NA)
  for (i in rlength[str:end]) {
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
  
  #create dataframe columns
  names_list <-
    data.frame(count_names[((((end - str) + 1) * 3) + 1):1,], stringsAsFactors = F)
  names_list[(nrow(names_list)),] <- "wildtype_siRNA"
  
  nt_sub_list <-
    data.frame(nt_sub[((((end - str) + 1) * 3) + 1):1,], stringsAsFactors = F)
  #nt_sub_list[(nrow(nt_sub_list)),] <- "substitution"
  
  #transformand write table
  ref_list <- nt[str:end, 1:((((end - str) + 1) * 3) + 1)]
  trans_list <-
    as.data.frame(t(ref_list[, ncol(ref_list):1]), stringsAsFactors = F)
  
  for (i in (c(1:nrow(names_list)))) {
    names_list[i, 2] <-
      paste(trans_list[i, 1:(ncol(trans_list))],sep = "", collapse = "")
  }
  #write column names
  names_list<- cbind(names_list, nt_sub_list)
  colnames(names_list) <- c("name", "target", "nucleotide")
  names_list <-
    rbind((names_list[(nrow(names_list)),]), names_list)
  names_list <- names_list[-c(nrow(names_list)),]
  
  return(names_list)
}

#target<- ref_replace(str,end, replicon_ref)
#write.table(target, "ref_list.csv", sep = "\t", quote = F, row.names = F)

main <- function() {
  
  if (length(args)==0) {
    stop("At least the first 2 arguments must be provided.
         START nucleotide position.
         END nucleotide position.
         Input .fasta file name (default = *.fasta)
         Output file name (default = substitutions.txt)", call.=FALSE)
  }
  else if (length(args)==2) {
    # stdin input read

    refname<- list.files(".", pattern ="fasta", all.files = F, full.names = F)
    if (length(refname)==0) {
      stop("No .fasta file available for input")
    }
    else if (length(refname)==1) {
      target<- ref_replace(str,end, refname)
    }
    args[4] = "substitutions.txt"
  }
  else if (length(args)==3) {
    # default output file
    args[4] = "substitutions.txt"
    target<- ref_replace(str,end, args[3])
  }
  else if (length(args)==4) {
    filenames = as.character(args[3])
    target<- ref_replace(str,end, filenames)
  }
  else if (length(args)>4) {
    stop("Too many arguments")
  }
  #cat(args, sep = "\n")
  write.table(target, file = args[4] , sep = "\t", quote = F, row.names = F)
}

main()
