library(dplyr)
library(ggplot2)

benchmark_ngram_length <- rbind(read.csv("./results/ngram_benchmark_opt.csv") %>% 
                                  mutate(error = sqrt(mse),
                                         opt = TRUE) %>% 
                                  inner_join(read.csv("./data/full_names.csv")),
                                read.csv("./results/ngram_benchmark.csv") %>% 
                                  mutate(error = sqrt(mse),
                                         opt = FALSE) %>% 
                                  inner_join(read.csv("./data/full_names.csv"))
)



benchmark_ngram_length %>% 
  group_by(nice, ngram_length, seq_type, opt) %>%
  summarise(mean_error = mean(error),
            sd_error = sd(error)) %>%
  mutate(upper = mean_error + sd_error,
         lower = mean_error - sd_error) %>% 
  ggplot(aes(x = factor(ngram_length), y = mean_error, color = opt, shape = seq_type,
             ymax = upper, ymin = lower)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  #geom_errorbar() +
  scale_x_discrete("n-gram length") +
  facet_wrap(~ nice, scales = "free_y") +
  theme_bw()


lapply(list("growth_doubl",
            "growth_rate",
            c("mean_gt", "mean_ogt"),
            c("mean_gn", "mean_ogn"),
            c("mean_gp", "mean_ogp")),
       function(i)
         filter(benchmark_ngram_length, task.id %in% i) %>% 
         group_by(nice, ngram_length, seq_type) %>%
         summarise(mean_error = mean(error),
                   sd_error = sd(error)) %>%
         mutate(upper = mean_error + sd_error,
                lower = mean_error - sd_error) %>% 
         ggplot(aes(x = factor(ngram_length), y = mean_error, color = seq_type,
                    ymax = upper, ymin = lower)) +
         geom_point(position = position_dodge(0.2)) + 
         scale_x_discrete("n-gram length") +
         geom_errorbar(position = position_dodge(0.2)) +
         facet_wrap(~ nice, nrow = 1) +
         theme_bw()) 
  

rbind(lapply(list("growth_doubl",
                  "growth_rate",
                  c("mean_gt", "mean_ogt"),
                  c("mean_gn", "mean_ogn"),
                  c("mean_gp", "mean_ogp")),
             function(i)
               filter(benchmark_ngram_length_opt, task.id %in% i) %>% 
               group_by(nice, ngram_length, seq_type) %>%
               summarise(mean_error = mean(error),
                         sd_error = sd(error))) %>% 
        do.call(rbind, .) %>% 
        mutate(only_optim = TRUE),
      lapply(list("growth_doubl",
                  "growth_rate",
                  c("mean_gt", "mean_ogt"),
                  c("mean_gn", "mean_ogn"),
                  c("mean_gp", "mean_ogp")),
             function(i)
               filter(benchmark_ngram_length, task.id %in% i) %>% 
               group_by(nice, ngram_length, seq_type) %>%
               summarise(mean_error = mean(error),
                         sd_error = sd(error))) %>% 
        do.call(rbind, .) %>% 
        mutate(only_optim = FALSE)) %>% 
  write.csv(file = "./results/all_res.csv", row.names = FALSE)

library(xtable)
read.csv("./results/ngram_benchmark_opt.csv") %>% 
  mutate(error = sqrt(mse),
         opt = TRUE) %>% 
  inner_join(read.csv("./data/full_names.csv")) %>% 
  group_by(nice) %>% 
  summarise(`Error (mean)` = mean(error),
            `Error (SD)` = sd(error)) %>% 
  xtable %>% 
  print(include.rownames = FALSE)
