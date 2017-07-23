library(dplyr)
library(biogram)
library(seqinr)
library(pbapply)

raw_dat <- read.csv("./data/dump_17-07-23.csv") 

rna_seqs <- read.fasta("./data/rRNA_.txt")

# extracting strain info
strain_df <- lapply(rna_seqs, function(i) 
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
  data.frame

# species who have more than one strain
more_then_one_strain <- strain_df %>% 
  group_by(X1, X2) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  select(X1) %>% 
  table %>% 
  data.frame %>% 
  filter(Freq > 1) %>% 
  select(-Freq) %>% 
  unlist %>% 
  as.character 

# save strains and their strain ID from sequences
filter(raw_dat, Name %in% more_then_one_strain) %>% 
  select(Name, Type.strain, DSM.strain.number) %>% 
  inner_join(strain_df, by = c("Name" = "X1")) %>% 
  filter(!duplicated(.)) %>% 
  mutate(X2 = sub("_", "-", X2)) %>% 
  filter(X2 != "") %>% 
  group_by(Name, X2) %>% 
  mutate(found = grepl(X2, Type.strain) + grepl(X2, DSM.strain.number)) %>% 
  filter(found == 0) %>% 
  write.csv(file = "./results/doubtful_strains.csv", row.names = FALSE)


