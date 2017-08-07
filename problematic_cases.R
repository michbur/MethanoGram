library(dplyr)
library(biogram)
library(seqinr)
library(pbapply)
library(mlr)

source("./functions/validate_seqs.R")

opt_pars <- read.csv("./results/ngram_benchmark_multi.csv") %>% 
  mutate(mean_error = sqrt(mean_error),
         sd_error = sqrt(sd_error)) %>% 
  left_join(read.csv("./data/full_names.csv")) %>% 
  filter(rna_type == "RNA2", mcra_type == "McrA3", feature_prop == 0.25) %>% 
  group_by(nice) %>% 
  filter(mean_error == min(mean_error))

raw_dat <- read.csv("./data/dump_17-07-23.csv") 

mcra_seqs <- list(read.fasta("./raw_data/McrA_1.txt"),
                  read.fasta("./raw_data/McrA_2.txt"),
                  read.fasta("./raw_data/McrA_3.txt")) %>% 
  lapply(function(j)
    lapply(j, function(i) {
      res <- i[i != "-"]
      mostattributes(res) <- attributes(i)
      res
    })
  ) %>% 
  lapply(validate_seqs)

names(mcra_seqs) <- paste0("McrA", 1L:3)


rna_seqs <- list(read.fasta("./raw_data/RNA_1.txt"),
                 read.fasta("./raw_data/RNA_2.txt")) %>% 
  lapply(function(j)
    lapply(j, function(i) {
      res <- i[i != "-" & i != " "]
      mostattributes(res) <- attributes(i)
      res
    })
  ) %>% 
  lapply(validate_rna)

names(rna_seqs) <- paste0("RNA", 1L:2)



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


train_seqs <- list(rna_seqs = rna_seqs[["RNA2"]], mcra_seqs = mcra_seqs[["McrA3"]])


jackknife_res <- lapply(1L:nrow(opt_pars), function(ith_row) {
  mlr_pars <- opt_pars[ith_row, ]

  ith_seqs_data <- train_seqs[[as.character(mlr_pars[["seq_type"]])]]
  
  ith_seqs_data <- ith_seqs_data[rownames(ith_seqs_data) %in% all_three, ]
  
  ngram_matrix <- count_ngrams(ith_seqs_data, mlr_pars[["ngram_length"]], 
                               u = c("a", "c", "g", "t"), scale = TRUE)
  
  normalized_ngrams <- ngram_matrix %>% 
    as.matrix %>% 
    data.frame %>% 
    mutate(source = rownames(ith_seqs_data)) %>% 
    group_by(source) %>% 
    summarise_all(mean) 
  
  ith_condition <- mlr_pars[["task.id"]]
  dat <- conditions_dat[, c("Name", ith_condition)] %>% 
    inner_join(normalized_ngrams, by = c("Name" = "source")) %>% 
    select(-Name)
  
  dat_names <- conditions_dat[, c("Name", ith_condition)] %>% 
    inner_join(normalized_ngrams, by = c("Name" = "source")) %>% 
    select(Name) %>% 
    unlist()
  
  predict_ngrams <- makeRegrTask(id = ith_condition, 
                                 data = dat, 
                                 target = ith_condition)
  
  filtered_ngrams <- filterFeatures(predict_ngrams, method = "linear.correlation", perc = 0.25)
  
  
  learner_pars <- makeParamSet(
    makeDiscreteParam("num.trees", values = mlr_pars[["num.trees"]]),
    makeDiscreteParam("min.node.size", values = mlr_pars[["min.node.size"]]),
    makeDiscreteParam("mtry", values = mlr_pars[["mtry"]])
  )
  
  
  
  learnerRF <- makeLearner("regr.ranger", par.vals = list(num.trees = mlr_pars[["num.trees"]],
                                                          min.node.size = mlr_pars[["min.node.size"]],
                                                          mtry = mlr_pars[["mtry"]]))
  cv_res <- crossval(learnerRF, filtered_ngrams, nrow(dat))
  
  lapply(getRRPredictionList(cv_res)[["test"]], function(ith_pred)
    data.frame(ith_pred)) %>% 
    do.call(rbind, .) %>% 
    mutate(error = sqrt((response - truth)^2),
           source = dat_names[id],
           task.id = ith_condition) %>% 
    arrange(desc(error)) %>%
    select(task.id, source, truth, error) 
})  


jackknife_res
