library(dplyr)
library(biogram)
library(seqinr)
library(pbapply)

source("./functions/extract_strains.R")

raw_dat <- read.csv("./data/dump_17-07-23.csv") 

rna_seqs <- read.fasta("./data/rRNA_.txt")

write.csv(extract_strains(raw_dat, rna_seqs), 
          file = "./results/doubtful_strains.csv", row.names = FALSE)

mcra_seqs <- c(read.fasta("./raw_data/McrA_1.txt"),
               read.fasta("./raw_data/McrA_2.txt"),
               read.fasta("./raw_data/McrA_3.txt")) %>% 
  lapply(function(i) {
    res <- i[i != "-"]
    mostattributes(res) <- attributes(i)
    res
  })

write.csv(extract_strains(raw_dat, mcra_seqs), 
          file = "./results/doubtful_strains_mcra.csv", row.names = FALSE)
