library(dplyr)
library(biogram)
library(seqinr)
library(pbapply)
library(mlr)

source("./functions/validate_seqs.R")

raw_dat <- read.csv("./data/dump_17-07-23.csv") 

rna_seqs <- read.fasta("./data/rRNA_.txt") %>% 
  validate_seqs()


mcra_seqs <- c(read.fasta("./raw_data/McrA_1.txt"),
               read.fasta("./raw_data/McrA_2.txt"),
               read.fasta("./raw_data/McrA_3.txt")) %>% 
  lapply(function(i) {
    res <- i[i != "-"]
    mostattributes(res) <- attributes(i)
    res
  }) %>% 
  validate_seqs()

seq_data <- list(mcra_seqs = mcra_seqs,
                 rna_seqs = rna_seqs)

# raw_dat[raw_dat[["Name"]] %in% normalized_ngrams[["source"]], 
#         c("Growth.doubling.time..h.", "Growth.rate", 
#           "Min..growth.temp.", "Max..growth.temp.", 
#           "Min..optimal.growth.temp.", "Max..optimal.growth.temp.",  
#           "Min..growth.NaCl", "Max..growth.NaCl", 
#           "Min..optimal.growth.NaCl", "Max..optimal.growth.NaCl", 
#           "Min..growth.pH", "Max..growth.pH", 
#           "Min..optimal.growth.pH", "Max..optimal.growth.pH")]

conditions_dat <- raw_dat[c("Name", 
                            "Growth.doubling.time..h.", "Growth.rate", 
                            "Min..growth.temp.", "Max..growth.temp.", 
                            "Min..optimal.growth.temp.", "Max..optimal.growth.temp.",  
                            "Min..growth.NaCl", "Max..growth.NaCl", 
                            "Min..optimal.growth.NaCl", "Max..optimal.growth.NaCl", 
                            "Min..growth.pH", "Max..growth.pH", 
                            "Min..optimal.growth.pH", "Max..optimal.growth.pH")] %>% 
  rename(growth_doubl = Growth.doubling.time..h.,
         growth_rate = Growth.rate,
         min_gt = Min..growth.temp.,
         max_gt = Max..growth.temp.,
         min_ogt = Min..optimal.growth.temp.,
         max_ogt = Max..optimal.growth.temp.,
         min_gn = Min..growth.NaCl,
         max_gn = Max..growth.NaCl,
         min_ogn = Min..optimal.growth.NaCl,
         max_ogn = Max..optimal.growth.NaCl,
         min_gp = Min..growth.pH,
         max_gp = Max..growth.pH,
         min_ogp = Min..optimal.growth.pH,
         max_ogp = Max..optimal.growth.pH) %>% 
  mutate(mean_gt = (min_gt + max_gt)/2,
         mean_ogt = (min_ogt + max_ogt)/2,
         mean_gn = (min_gn + max_gn)/2,
         mean_ogn = (min_ogn + max_ogn)/2,
         mean_gp = (min_gp + max_gp)/2,
         mean_ogp = (min_ogp + max_ogp)/2) %>% 
  na.omit


both_mcra_rna <- intersect(unique(rownames(rna_seqs)), unique(rownames(mcra_seqs)))
both_mcra_conditions <- intersect(as.character(conditions_dat[["Name"]]), unique(rownames(mcra_seqs)))
both_rna_conditions <- intersect(as.character(conditions_dat[["Name"]]), unique(rownames(rna_seqs)))
all_three <- intersect(as.character(conditions_dat[["Name"]]), both_mcra_rna)


# library(VennDiagram)
# 
# venn.plot <- draw.triple.venn(length(unique(rownames(rna_seqs))), 
#                               length(unique(rownames(mcra_seqs))), 
#                               length(as.character(conditions_dat[["Name"]])),
#                               length(both_mcra_rna), 
#                               length(both_mcra_conditions), 
#                               length(both_rna_conditions), 
#                               length(all_three), 
#                               c("Known RNA", "Known mcrA", "Known conditions"),
#                               fill = c("#f1a340", "#998ec3", "#f7f7f7"))
# grid.draw(venn.plot)
# grid.newpage()





benchmark_ngram_length <- lapply(1L:2, function(ngram_length) {
  
  lapply(c("mcra_seqs", "rna_seqs"), function(ith_seqs) {
    ith_seqs_data <- seq_data[[ith_seqs]] 
    ith_seqs_data <- ith_seqs_data[rownames(ith_seqs_data) %in% all_three, ]

    ngram_matrix <- count_ngrams(ith_seqs_data, ngram_length, u = c("a", "c", "g", "t"), scale = TRUE)
    
    normalized_ngrams <- ngram_matrix %>% 
      as.matrix %>% 
      data.frame %>% 
      mutate(source = rownames(ith_seqs_data)) %>% 
      group_by(source) %>% 
      summarise_all(mean) 
    
    
    bench_res <- lapply(c("growth_doubl", "growth_rate", "mean_gt", "mean_ogt", 
                          "mean_gn", "mean_ogn", "mean_gp", "mean_ogp"),
                        function(ith_condition) {
                          dat <- conditions_dat[, c("Name", ith_condition)] %>% 
                            inner_join(normalized_ngrams, by = c("Name" = "source")) %>% 
                            select(-Name)
                          
                          predict_par <- makeRegrTask(id = ith_condition, 
                                                      data = dat, 
                                                      target = ith_condition)
                          
                          learnerRF <- makeLearner("regr.ranger")
                          
                          set.seed(1410)
                          benchmark(learnerRF, predict_par, makeResampleDesc("CV", iters = 5L))
                        })
    
    
    lapply(bench_res, function(i) 
      data.frame(i)) %>% 
      do.call(rbind, .) %>% 
      mutate(ngram_length = ngram_length,
             seq_type = ith_seqs)  
  }) %>% 
    do.call(rbind, .)
}) %>%
  do.call(rbind, .)

write.csv(benchmark_ngram_length, file = "./results/ngram_benchmark.csv")

# mutate(benchmark_ngram_length, error = sqrt(mse)) %>%
#   group_by(task.id, ngram_length, seq_type) %>%
#   summarise(mean_error = mean(error)) %>% 
#   ggplot(aes(x = factor(ngram_length), y = mean_error, color = seq_type)) +
#   geom_point(position = position_jitter(width = 0.2)) +
#   facet_wrap(~ task.id, scales = "free_y") +
#   theme_bw()
 
# group_by(task.id) %>% 
#   summarise(mse = mean(mse)) %>% 
#   mutate(error = sqrt(mse),
#          ngram_length = ngram_length)

