library(dplyr)
library(biogram)
library(seqinr)
library(pbapply)

raw_dat <- read.csv("./data/dump_17-07-23.csv") 

rna_seqs <- read.fasta("./data/rRNA_.txt")

validated_seq_names <- lapply(rna_seqs, function(i) 
  strsplit(attr(i, "name"), split = "|", fixed = TRUE)[[1]][2]) %>% 
  unlist %>% 
  unname %>% 
  sub("_", " ", .) %>% 
  sub("_", "|", .) %>% 
  strsplit(split = "|", fixed = TRUE) %>% 
  lapply(function(i) if(length(i) == 1) {
    c(i, "") 
  } else {
    i
  }) %>% 
  do.call(rbind, .) %>% 
  data.frame(stringsAsFactors = FALSE) %>% 
  mutate(X2 = sub("_", "-", X2)) %>% 
  left_join(read.csv("./results/doubtful_strains_modif.csv",
                     stringsAsFactors = FALSE)[c("Name", "X2", "is.reference")],
            by = c("X1" = "Name", "X2" = "X2"))  %>% 
  mutate(is.reference = ifelse(is.na(is.reference), TRUE, is.reference)) 
  
filter(validated_seq_names, is.reference) %>% 
  select(X1) %>% 
  unlist %>% 
  unique %>% 
  length
# 141 species

# remove following nucleotides: c("r", "n", "b", "s", "m", "d", "w", "y", "k", "v")
validated_seqs <- list2matrix(rna_seqs[validated_seq_names[["is.reference"]]])

rownames(validated_seqs) <- strsplit(rownames(validated_seqs), split = "|", fixed = TRUE) %>% 
  sapply(function(i) i[2]) %>% 
  strsplit(split = "_") %>% 
  sapply(function(i) paste0(i[1L:2], collapse = " "))

ngram_matrix <- count_ngrams(validated_seqs, 4, u = c("a", "c", "g", "t"), scale = TRUE)

which(validated_seqs == "r", arr.ind = TRUE)
validated_seqs[140, ] 
