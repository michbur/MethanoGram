library(dplyr)
library(ggplot2)
library(xtable)

res <- read.csv("./results/ngram_benchmark_multi.csv") %>% 
  mutate(mean_error = sqrt(mean_error),
         sd_error = sqrt(sd_error)) %>% 
  left_join(read.csv("./data/full_names.csv"))

group_by(res, nice) %>% 
  filter(mean_error == min(mean_error)) %>% 
  mutate(mean_error = formatC(mean_error, 4)) %>% 
  select(nice, mean_error)

group_by(res, nice) %>% 
  filter(mean_error == min(mean_error)) %>% 
  mutate(train_seq = ifelse(seq_type == "rna_seqs", 
                            as.character(rna_type), 
                            as.character(mcra_type))) %>% 
  select(-seq_type, -rna_type, -mcra_type) 

dat_type <- mutate(res, train_seq = ifelse(seq_type == "rna_seqs", 
                            as.character(rna_type), 
                            as.character(mcra_type))) %>% 
  select(-seq_type, -rna_type, -mcra_type) %>% 
  group_by(nice, train_seq) %>% 
  filter(mean_error == min(mean_error)) %>% 
  ungroup %>% 
  arrange(task.id, mean_error) %>% 
  ungroup %>% 
  filter(!duplicated(.)) 
  
ggplot(dat_type, aes(x = task.id, y = mean_error, 
                     color = train_seq, shape = factor(feature_prop))) + 
  geom_point(position = "jitter") +
  facet_wrap(~ nice, scales = "free") +
  theme_bw()

# select McrA3 i RNA2, feature_prop 0.25

filter(res, rna_type == "RNA2", mcra_type == "McrA3", feature_prop == 0.25) %>% 
  group_by(nice) %>% 
  filter(mean_error == min(mean_error)) %>% 
  mutate(mean_error = formatC(mean_error, digits = 2, format = "g")) %>% 
  select(Name = nice, `Mean error` = mean_error) %>% 
  ungroup() %>% 
  slice(c(2, 4, 3, 1, 5)) %>% 
  xtable(caption = "Results of nested cross-validation of MethanoGram.",
         label = "tab:nested_cv") %>% 
  print(include.rownames = FALSE)

filter(res, rna_type == "RNA2", mcra_type == "McrA3", feature_prop == 0.25) %>% 
  group_by(nice) %>% 
  filter(mean_error == min(mean_error)) %>% 
  mutate(mean_error = formatC(mean_error, digits = 2, format = "g")) %>% 
  select(Name = nice, `Mean error` = mean_error) %>% 
  ungroup() %>% 
  slice(c(2, 4, 3, 1, 5)) %>% 
  write.csv(file = "nested_cv.csv", row.names = FALSE)

