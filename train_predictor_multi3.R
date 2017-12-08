library(dplyr)
library(biogram)
library(seqinr)
library(pbapply)
library(mlr)

configureMlr(show.info = FALSE)

source("./functions/validate_seqs.R")

raw_dat <- read.csv("./data/dump_17-07-23.csv") 

mcra_seqs <- list(read.fasta("./raw_data/McrA_1.txt")) %>% 
  unlist(recursive = FALSE) %>% 
  lapply(function(i) {
    res <- i[i != "-" & i != " "]
    mostattributes(res) <- attributes(i)
    res
  }) %>% 
  validate_seqs

rna_seqs <- list(read.fasta("./raw_data/RNA_1.txt")) %>% 
  unlist(recursive = FALSE) %>% 
  lapply(function(i) {
    res <- i[i != "-" & i != " "]
    mostattributes(res) <- attributes(i)
    res
  }) %>% 
  validate_rna



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
  select(Name, growth_doubl, growth_rate, mean_ogt, mean_ogn, mean_ogp) %>% 
  na.omit %>% 
  filter(Name != "Methanoculleus sediminis")

both_mcra_rna <- intersect(unique(rownames(rna_seqs)), unique(rownames(mcra_seqs)))
both_mcra_conditions <- intersect(as.character(conditions_dat[["Name"]]), unique(rownames(mcra_seqs)))
both_rna_conditions <- intersect(as.character(conditions_dat[["Name"]]), unique(rownames(rna_seqs)))
all_three <- intersect(as.character(conditions_dat[["Name"]]), both_mcra_rna)

chosen_mcra_seqs <- mcra_seqs[rownames(mcra_seqs) %in% all_three, ]
chosen_rna_seqs <- rna_seqs[rownames(rna_seqs) %in% all_three, ]

seq_dat <- list(mcra = lapply(1L:6, function(ith_len) {
  res <- count_ngrams(chosen_mcra_seqs, ith_len, u = c("a", "c", "g", "t"), scale = TRUE)
  colnames(res) <- paste0("mcra_", colnames(res))
  
  lapply(unique(rownames(chosen_mcra_seqs)), function(ith_spec) {
    res[ith_spec == rownames(chosen_mcra_seqs), ] %>% 
      col_means %>% 
      as.list %>% 
      data.frame %>% 
      mutate(species = ith_spec)
  }) %>% 
    do.call(rbind, .)
}), 
rna = lapply(1L:6, function(ith_len) {
  res <- count_ngrams(chosen_rna_seqs, ith_len, u = c("a", "c", "g", "t"), scale = TRUE)
  colnames(res) <- paste0("rna_", colnames(res))
  
  lapply(unique(rownames(chosen_rna_seqs)), function(ith_spec) {
    res[ith_spec == rownames(chosen_rna_seqs), ] %>% 
      col_means %>% 
      as.list %>% 
      data.frame %>% 
      mutate(species = ith_spec)
  }) %>% 
    do.call(rbind, .)
}))

save(seq_dat, file = "./data/seq_dat.RData")

training_data <- expand.grid(type1 = c("mcra", "rna"), len1 = 1L:6, 
                             type2 = c("mcra", "rna"), len2 = 1L:6) %>% 
  split(1L:nrow(.)) 




ngram_dat_list <- lapply(training_data, function(i) {
  if(i[["type1"]] == i[["type2"]] & i[["len1"]] == i[["len2"]]) {
    seq_dat[[i[["type1"]]]][[i[["len1"]]]]
  } else {
    cbind(select(seq_dat[[i[["type1"]]]][[i[["len1"]]]], -species),
          seq_dat[[i[["type2"]]]][[i[["len2"]]]])
  }
})

library("parallelMap")
parallelStartSocket(4)

benchmark_res <- pblapply(c("growth_doubl", "growth_rate", "mean_ogt", 
                          "mean_ogn", "mean_ogp"), function(ith_condition)
                            lapply(c(0.1, 0.25, 0.5), function(feature_frac)
                              lapply(1L:length(ngram_dat_list), function(ith_ngram_dat_list) {
                                try({
                                  ngram_dat <- ngram_dat_list[[ith_ngram_dat_list]]
                                  
                                  dat <- conditions_dat[, c("Name", ith_condition)] %>% 
                                    inner_join(ngram_dat, by = c("Name" = "species")) %>% 
                                    select(-Name)
                                  
                                  predict_ngrams <- makeRegrTask(id = ith_condition, 
                                                                 data = dat, 
                                                                 target = ith_condition)
                                  
                                  filtered_ngrams <- filterFeatures(predict_ngrams, method = "linear.correlation", perc = feature_frac)
                                  n_features <- filtered_ngrams[["task.desc"]][["n.feat"]][["numerics"]]
                                  mtry_possibilities <- round(c(n_features/4, n_features/3, n_features/2), 0)
                                  mtry_possibilities[mtry_possibilities == 0] <- 1
                                  mtry_possibilities <- unique(mtry_possibilities)
                                  
                                  learnerRF <- makeLearner("regr.ranger")
                                  learner_pars <- makeParamSet(
                                    makeDiscreteParam("num.trees", values = c(250, 500, 750, 1000)),
                                    makeDiscreteParam("min.node.size", values = c(3, 5, 7)),
                                    makeDiscreteParam("mtry", values = mtry_possibilities)
                                  )
                                  
                                  set.seed(1410)
                                  
                                  inner <- makeResampleDesc("CV", iters = 5)
                                  outer <- makeResampleDesc("CV", iters = 3)
                                  learnerRF_tuned <- makeTuneWrapper(learnerRF, 
                                                                     resampling = inner, 
                                                                     par.set = learner_pars, 
                                                                     control = makeTuneControlGrid())
                                  
                                  trainedRF <- train(learnerRF_tuned, filtered_ngrams)
                                  nested_cv <- resample(learnerRF_tuned, filtered_ngrams, outer, extract = getTuneResult)
                                  
                                  nested_res <- getNestedTuneResultsOptPathDf(nested_cv) 
                                  
                                  group_by(nested_res, mtry, num.trees, min.node.size) %>% 
                                    summarise(mean_error = mean(mse.test.mean),
                                              sd_error = sd(mse.test.mean)) %>% 
                                    mutate(task.id = ith_condition) %>% 
                                    arrange(mean_error) %>% 
                                    mutate(dat_list = ith_ngram_dat_list,
                                           perc = feature_frac)
                                }, silent = TRUE)
                              })
                            )
)

parallelStop()

save(benchmark_res, file = "./results/ngram_benchmark_full.RData")
