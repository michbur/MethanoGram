library(dplyr)

read.csv("./results/ngram_benchmark_multi.csv") %>% 
  mutate(mean_error = sqrt(mean_error),
         sd_error = sqrt(sd_error)) %>% 
  left_join(read.csv("./data/full_names.csv")) %>% 
  droplevels() %>% 
  mutate(mtry_nice = round(mtry/((4^ngram_length)*feature_prop*ifelse(seq_type == "both", 2, 1))*12, 0)/12,
         seq_type = factor(seq_type, labels = c("Both", "16 rRNA", "mcrA"))) %>% 
  filter(rna_type == "RNA2", mcra_type == "McrA3", feature_prop == 0.25) %>% 
  group_by(nice) %>% 
  filter(mean_error == min(mean_error)) %>% 
  ungroup %>% 
  select(-mean_error, -sd_error, -mtry_nice)
