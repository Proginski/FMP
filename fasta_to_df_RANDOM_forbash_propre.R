#! /usr/bin/Rscript

# This script takes as input a fasta file and output a csv file containing 
# the codon counts needed for further prediction


# In order to pass arguments with bash
args <- commandArgs(trailingOnly = TRUE)



main <- function(fasta_path, 
                 features = None,
                 log = ""){
  

  library(stringr)
  
  # Return the codon frequency table
  get_codon_count <- function(fasta_path){
    
    library(coRdon)
    set <- readSet(file = fasta_path)
    codon_Table <- codonTable(set)
    codon_Counts <- codonCounts(codon_Table)
    row.names(codon_Counts) <- codon_Table@ID
    
    # Remove the stop codons before scaling
    exit_codons <- c("TAA","TAG","TGA")
    codon_Counts[ ,exit_codons] <- 0
    
    return(codon_Counts)
  }
  
  
  # Return the mean frequency of a given reference set of nucleotidic sequences
  getSeqFreq <- function(codon_count, by){
    
    codon_Counts_ref <- codon_count
    
    
    if(by == "same_position"){
      # Transform codon_Counts to a table so that for position 1 of the codon it displays the proportion of A T C and G.
      
      library(stringr)
      meanBasesFractions_ref <- data.frame(matrix(nrow = nrow(codon_count), ncol = 12))
      colnames <- c()
      for(pos in c(1,2,3)){
        for(base in c("A","C","G","T")){
          colnames <- c(colnames, paste(pos,base,sep=""))
        }
      }
      
      colnames(meanBasesFractions_ref) <- colnames
      
      for(pos in c(1,2,3)){
        for(base in c("A","C","G","T")){
          
          codon_givenBase_givenPosition_ref <- substring(colnames(codon_Counts_ref),pos,pos) == base
          
          
          sum_givenBase_givenPosition_ref <- rowSums(t(t(codon_Counts_ref)*codon_givenBase_givenPosition_ref))
          # relFreqBasePos_ref <- sum_givenBase_givenPosition_ref/rowSums(codon_Counts_ref)
          relFreqBasePos_ref <- sum_givenBase_givenPosition_ref
          
          meanBasesFractions_ref[,paste(pos, base, sep="")] <- relFreqBasePos_ref
          
          
        }
      }
    }
    
    
    
    
    if(by == "same_base"){
      # Transform codon_Counts to a table so that for the A base it displays its fraction in firstn second and third position of the codon.
      library(stringr)
      # Mean fraction of each base among the 3 codon position
      meanBasesFractions_ref <- data.frame(A1=rep(as.double(NA),dim( codon_Counts_ref)[1]),A2=rep(as.double(NA),dim( codon_Counts_ref)[1]),A3=rep(as.double(NA),dim( codon_Counts_ref)[1]),
                                           C1=rep(as.double(NA),dim( codon_Counts_ref)[1]),C2=rep(as.double(NA),dim( codon_Counts_ref)[1]),C3=rep(as.double(NA),dim( codon_Counts_ref)[1]),
                                           G1=rep(as.double(NA),dim( codon_Counts_ref)[1]),G2=rep(as.double(NA),dim( codon_Counts_ref)[1]),G3=rep(as.double(NA),dim( codon_Counts_ref)[1]),
                                           T1=rep(as.double(NA),dim( codon_Counts_ref)[1]),T2=rep(as.double(NA),dim( codon_Counts_ref)[1]),T3=rep(as.double(NA),dim( codon_Counts_ref)[1]))
      
      
      
      for (base in c("A","C","G","T")){
        for (pos in c(1,2,3)){
          codons_givenBase_ref <- str_count(colnames(codon_Counts_ref), base)
          
          codon_givenBase_givenPosition_ref <- substring(colnames(codon_Counts_ref),pos,pos) == base
          
          
          sum_givenBase_givenPosition_ref <- rowSums(t(t(codon_Counts_ref)*codon_givenBase_givenPosition_ref))
          
          all_givenBase_ref <- rowSums(t(t(codon_Counts_ref)*codons_givenBase_ref))
          
          # relFreqBasePos_ref <- sum_givenBase_givenPosition_ref/all_givenBase_ref
          relFreqBasePos_ref <- sum_givenBase_givenPosition_ref
          
          meanBasesFractions_ref[,paste(base,pos, sep="")] <- relFreqBasePos_ref
          
        }
      }
    }
    
    return(meanBasesFractions_ref)
    
  }
  
  
  
  # Return the mean frequency of a given reference set of nucleotidic sequences
  get_NNX_Freq <- function(codon_count){
    
    # Retrieve the first two letters of codon_count colnames
    NNNs <- colnames(codon_count)
    NNX_colnames <- substr(NNNs ,start = 1,stop = 2)
    
    
    # Generate all possible NNX codons first two letters
    NNXs <- c()
    for (base1 in c("T","C","A","G")) {
      for (base2 in c("T","C","A","G")){
        NNX <- paste(base1, base2, sep="")
        NNXs <- c(NNXs, NNX)
      }
    }
    
    # Create an empty dataframe to return
    NNX_count <- data.frame(matrix(ncol = 16, nrow = nrow(codon_count)))
    colnames(NNX_count) <- NNXs
    
    
    # Compute codon_count row sums
    rowsums <- rowSums(codon_count)
    
    for (NNX in NNXs){
      col_tosum <- grep(NNX, NNX_colnames)
      NNX_count[ ,NNX] <- rowSums(codon_count[ ,col_tosum])#/rowsums
      
    }
    
    colnames(NNX_count) <- paste(colnames(NNX_count),"X", sep="")
    
    return(NNX_count)
  }
  
  
  
  # Return the mean frequency of a given reference set of nucleotidic sequences
  getAAfreq <- function(codon_count){
    
    
    #### Generate the AAcodons table ###
    codons <- c()
    for (base1 in c('T','C','A','G')) {
      for (base2 in c('T','C','A','G')) {
        for (base3 in c('T','C','A','G')) {
          codons <- c(codons, paste(base1, base2, base3, sep=""))
        }
      }
    }
    
    letterCode1 <- c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", "X", "X", "C", "C", "X", "W", "L", "L", "L", "L", "P", "P", "P", "P", "H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I", "M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", "G", "G")
    
    letterCOde3 <- c("Phe", "Phe", "Leu", "Leu", "Ser", "Ser", "Ser", "Ser", "Tyr", "Tyr", "Xxx", "Xxx", "Cys", "Cys", "Xxx", "Trp", "Leu", "Leu", "Leu", "Leu", "Asp", "Asp", "Pro", "Pro", "Pro", "Pro", "His", "His", "Gln", "Gln", "Arg", "Arg", "Arg", "Arg", "Ile", "Ile", "Ile", "Met", "Val", "Thr", "Thr", "Thr", "Thr", "Asn", "Asn", "Lys", "Lys", "Ser", "Ser", "Arg", "Arg", "Val", "Val", "Val", "Ala", "Ala", "Ala", "Ala", "Glu", "Glu", "Gly", "Gly", "Gly", "Gly")
    
    AAcodons <- data.frame(codon = codons, letterCode1, letterCOde3)
    
    ####################################
    
    
    # AAcodons <- readRDS("AAcodons2.rds")
    AAcodons$letterCOde3 <- as.factor(AAcodons$letterCOde3)
    AAs <- levels(AAcodons$letterCOde3)
    
    # The dataframe that will receive the frequency of each AA
    AAfreqDF <- data.frame(matrix(ncol = length(AAs),
                                  nrow = nrow(codon_count)))
    # rownames(AAfreqDF) <- rownames(codon_count)
    colnames(AAfreqDF) <- AAs
    
    for (AA in AAs){
      
      # print(AA)
      # The codons that code for this AA
      assoc_codons <- AAcodons$codon[AAcodons$letterCOde3 == AA]
      
      # print(assoc_codons)
      
      if (length(assoc_codons) > 1){
        AAfreqDF[ ,AA] <- rowSums(codon_count[ ,assoc_codons])#/rowSums(codon_count)
      } 
      else{
        AAfreqDF[ ,AA] <- codon_count[ ,assoc_codons]#/rowSums(codon_count)
      }
    }
    
    # Add the cation info
    AAfreqDF$LysArg <- AAfreqDF$Lys+AAfreqDF$Arg
    
    return(AAfreqDF)

  }
  
  
  
  # Return the base content of a given set
  getBaseContent <- function(codon_count){
    
    library(stringr)
    bases <- c("A","C","G","T")
    
    df <- data.frame(matrix(nrow = nrow(codon_count),
                            ncol = length(bases)))
    colnames(df) <- paste(bases,"content", sep = "_")
    colsums <- colSums(codon_count)
    
    for (base in bases){
      base_nb_by_codon <- c()
      for (codon in colnames(codon_count)){
        base_nb_by_codon <- c(base_nb_by_codon, str_count(codon, base))
      }
      df[ ,paste(base,"content",sep="_")] <- rowSums(t(apply(codon_count,1,function(x) x*base_nb_by_codon)))
    }
    
    # Get relative frequency
    df <- as.data.frame(t(scale(t(df), center = FALSE, scale = colSums(t(df)))))
    df$AT_content <- df$A_content+df$T_content
    df$GC_content <- df$G_content+df$C_content
    # df$GC_skew <- (df$G_content - df$C_content) / df$GC_content
    
    return(df)
  }
         
  
  
  
  codon_count <- get_codon_count(fasta_path)

  # Get sequences length and codons relative frequency 
  seq_length <- rowSums(codon_count)
  codon_count <- codon_count/seq_length


  
  # The dataframe that will receive the features
  # First Generate the anwser vector
  category <- grepl("random", row.names(codon_count), fixed = TRUE)
  category <- as.numeric(category) # 0/1 version

  
  df <- data.frame(Name = row.names(codon_count), category, seq_length_nucl = seq_length)

  
  # If we need the codon count
  if("codon_count" %in% features){
    print("Codon count...")
    df <- cbind(df, codon_count)
    
  }
  
  
  # If we need the same base count
  if("same_base" %in% features){ 
    print("Same base...")
    df <- cbind(df, getSeqFreq(codon_count = codon_count, by = "same_base")) }
  
  
  # If we need the same position count
  if("same_position" %in% features){ 
    print("Same position...")
    df <- cbind(df, getSeqFreq(codon_count = codon_count, by = "same_position")) }
    
  
  # If we need the NNX count
  if("NNX" %in% features){ 
    print("NNXs...")
    df <- cbind(df, get_NNX_Freq(codon_count = codon_count))}
  
  # If we need the AA count
  if("AAs" %in% features){ 
    print("AAs...")
    df <- cbind(df, getAAfreq(codon_count = codon_count))}
  
  # If we need the Base content count
  if("BaseContent" %in% features){ 
    print("Base content...")
    df <- cbind(df, getBaseContent(codon_count = codon_count))}
  

  
  # Generate a csv file from the input file name.
  write.csv(df, file = paste(str_split(fasta_path,".nfasta")[[1]][1],".csv", sep=""), row.names = FALSE)

  print('DONE.')
  
  return(df)

}

# For now we always need to get all the features :
output <- main(fasta_path = args[1], features=c("codon_count","same_base","same_position","NNX","BaseContent","AAs"))
