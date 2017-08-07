library(dplyr)
library(ggplot2)

res <- read.csv("./results/ngram_benchmark_multi.csv") %>% 
  mutate(error = sqrt(mean),
         error_r = round(error, 2)) %>% 
  inner_join(read.csv("./data/full_names.csv"), by = c("task" = "task.id"))

group_by(res, nice) %>% 
  filter(error == min(error)) %>% 
  select(nice, error_r)

group_by(res, seq_type, nice, ngram_length) %>% 
  filter(error == min(error)) %>% 
  ggplot(aes(x = factor(ngram_length), y = error, color = seq_type))  +
  geom_point(size = 3, position = position_dodge(0.2)) +
  #geom_errorbar() +
  scale_x_discrete("n-gram length") +
  facet_wrap(~ nice, scales = "free_y") +
  theme_bw()
