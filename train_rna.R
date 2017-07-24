library(dplyr)
library(biogram)
library(seqinr)
library(pbapply)
library(mlr)

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
only_standard_nucleotides <- sapply(rna_seqs, function(i)
  !any(c("r", "n", "b", "s", "m", "d", "w", "y", "k", "v") %in% i))

validated_seqs <- list2matrix(rna_seqs[validated_seq_names[["is.reference"]] & only_standard_nucleotides])

rownames(validated_seqs) <- strsplit(rownames(validated_seqs), split = "|", fixed = TRUE) %>% 
  sapply(function(i) i[2]) %>% 
  strsplit(split = "_") %>% 
  sapply(function(i) paste0(i[1L:2], collapse = " "))

# raw_dat[raw_dat[["Name"]] %in% normalized_ngrams[["source"]], 
#         c("Growth.doubling.time..h.", "Growth.rate", 
#           "Min..growth.temp.", "Max..growth.temp.", 
#           "Min..optimal.growth.temp.", "Max..optimal.growth.temp.",  
#           "Min..growth.NaCl", "Max..growth.NaCl", 
#           "Min..optimal.growth.NaCl", "Max..optimal.growth.NaCl", 
#           "Min..growth.pH", "Max..growth.pH", 
#           "Min..optimal.growth.pH", "Max..optimal.growth.pH")]

conditions_dat <- raw_dat[raw_dat[["Name"]] %in% normalized_ngrams[["source"]], 
                          c("Name", 
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




benchmark_ngram_length <- lapply(1L:5, function(ngram_length) {
  
  ngram_matrix <- count_ngrams(validated_seqs, ngram_length, u = c("a", "c", "g", "t"), scale = TRUE)
  
  normalized_ngrams <- ngram_matrix %>% 
    as.matrix %>% 
    data.frame %>% 
    mutate(source = rownames(validated_seqs)) %>% 
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
                        
                        learnerRF <- makeLearner("regr.randomForest")
                        
                        set.seed(1410)
                        benchmark(learnerRF, predict_par, makeResampleDesc("CV", iters = 5L))
                      })
  
  
  lapply(bench_res, function(i) 
    data.frame(i)) %>% 
    do.call(rbind, .) %>% 
    mutate(ngram_length = ngram_length)
})

lapply(1L:5, function(i)
  benchmark_ngram_length[[i]] %>% 
    mutate(ngram_length = i)) %>% 
  do.call(rbind, .) %>% 
  mutate(error = sqrt(mse)) %>% 
  ggplot(aes(x = factor(ngram_length), y = error, color = factor(iter))) +
  geom_point() +
  facet_wrap(~ task.id, scales = "free_y") +
  theme_bw()
 
# group_by(task.id) %>% 
#   summarise(mse = mean(mse)) %>% 
#   mutate(error = sqrt(mse),
#          ngram_length = ngram_length)

