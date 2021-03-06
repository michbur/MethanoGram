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
                            "Min..optimal.growth.temp.", "Max..optimal.growth.temp.",
                            "Min..optimal.growth.NaCl", "Max..optimal.growth.NaCl", 
                            "Min..optimal.growth.pH", "Max..optimal.growth.pH")] %>% 
  rename(growth_doubl = Growth.doubling.time..h.,
         growth_rate = Growth.rate,
         min_ogt = Min..optimal.growth.temp.,
         max_ogt = Max..optimal.growth.temp.,
         min_ogn = Min..optimal.growth.NaCl,
         max_ogn = Max..optimal.growth.NaCl,
         min_ogp = Min..optimal.growth.pH,
         max_ogp = Max..optimal.growth.pH) %>% 
  mutate(mean_ogt = (min_ogt + max_ogt)/2,
         mean_ogn = (min_ogn + max_ogn)/2,
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





benchmark_ngram_length <- lapply(1L:7, function(ngram_length) {
  
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
    
    
    bench_res <- lapply(c("growth_doubl", "growth_rate", "mean_ogt", 
                          "mean_ogn", "mean_ogp"),
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

write.csv(benchmark_ngram_length, file = "./results/ngram_benchmark_opt.csv", row.names = FALSE)

 
# group_by(task.id) %>% 
#   summarise(mse = mean(mse)) %>% 
#   mutate(error = sqrt(mse),
#          ngram_length = ngram_length)

