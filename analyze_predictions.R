library(dplyr)
library(ggplot2)

benchmark_ngram_length <- read.csv("./results/ngram_benchmark.csv")

mutate(benchmark_ngram_length, error = sqrt(mse)) %>%
  group_by(task.id, ngram_length, seq_type) %>%
  summarise(mean_error = mean(error),
            sd_error = sd(error)) %>%
  mutate(upper = mean_error + sd_error,
         lower = mean_error - sd_error) %>% 
  ggplot(aes(x = factor(ngram_length), y = mean_error, color = seq_type,
             ymax = upper, ymin = lower)) +
  geom_point(position = position_dodge(0.2)) +
  geom_errorbar(position = position_dodge(0.2)) +
  facet_wrap(~ task.id, scales = "free_y") +
  theme_bw()
