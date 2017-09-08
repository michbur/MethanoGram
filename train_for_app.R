library(dplyr)
library(mlr)
library(biogram)
library(seqinr)

source("./functions/validate_seqs.R")

# response -------------------------------
raw_dat <- read.csv("./data/dump_17-07-23.csv") 

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

# training data ------------------------------

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

# tuned classifiers ---------------------------------
tuned_par_all <- read.csv("./results/ngram_benchmark_multi.csv") %>% 
  mutate(mean_error = sqrt(mean_error),
         sd_error = sqrt(sd_error)) %>% 
  left_join(read.csv("./data/full_names.csv")) %>% 
  droplevels() %>% 
  mutate(mtry_nice = round(mtry/((4^ngram_length)*feature_prop*ifelse(seq_type == "both", 2, 1))*12, 0)/12,
         seq_type = factor(seq_type, labels = c("Both", "16 rRNA", "mcrA"))) %>% 
  filter(rna_type == "RNA2", mcra_type == "McrA3", feature_prop == 0.25) 

tuned_par_rna <- filter(tuned_par_all, seq_type == "16 rRNA") %>% 
  group_by(nice) %>% 
  filter(mean_error == min(mean_error)) %>% 
  ungroup

tuned_par_mcra <- filter(tuned_par_all, seq_type == "mcrA") %>% 
  group_by(nice) %>% 
  filter(mean_error == min(mean_error)) %>% 
  ungroup


training_data <- list(rna = list(pars = tuned_par_rna,
                                 seqs = rna_seqs[[2]]),
                      mcra = list(pars = tuned_par_rna,
                                  seqs = mcra_seqs[[3]]))


both_mcra_rna <- intersect(unique(rownames(training_data[["rna"]][["seqs"]])), 
                           unique(rownames(training_data[["mcra"]][["seqs"]])))
both_mcra_conditions <- intersect(as.character(conditions_dat[["Name"]]), 
                                  unique(rownames(training_data[["mcra"]][["seqs"]])))
both_rna_conditions <- intersect(as.character(conditions_dat[["Name"]]), 
                                 unique(rownames(training_data[["rna"]][["seqs"]])))
all_three <- intersect(as.character(conditions_dat[["Name"]]), both_mcra_rna)


pred_list <- lapply(training_data, function(single_data) 
  lapply(single_data[["pars"]][["task.id"]], function(ith_condition) {
    
    opt_pars <- filter(single_data[["pars"]], task.id == ith_condition)
    ith_seqs_data <- single_data[["seqs"]][rownames(single_data[["seqs"]]) %in% all_three, ]
    ngram_matrix <- count_ngrams(ith_seqs_data, opt_pars[["ngram_length"]], 
                                 u = c("a", "c", "g", "t"), scale = TRUE)
    
    normalized_ngrams <- ngram_matrix %>% 
      as.matrix %>% 
      data.frame %>% 
      mutate(source = rownames(ith_seqs_data)) %>% 
      group_by(source) %>% 
      summarise_all(mean) 
    
    dat <- conditions_dat[, c("Name", ith_condition)] %>% 
      inner_join(normalized_ngrams, by = c("Name" = "source")) %>% 
      select(-Name)

    predict_ngrams <- makeRegrTask(id = ith_condition, 
                                   data = dat, 
                                   target = ith_condition)
    
    filtered_ngrams <- filterFeatures(predict_ngrams, method = "linear.correlation", 
                                      perc = opt_pars[["feature_prop"]])
    
    # filtered_ngrams <- filterFeatures(predict_ngrams, method = "linear.correlation", 
    #                                   perc = 1)
    
    learnerRF <- makeLearner("regr.ranger", par.vals = as.list(opt_pars[c("mtry", "num.trees", "min.node.size")]))
    train(learnerRF, filtered_ngrams)
  })
)

save(pred_list, file = "./app/pred_list.RData")


error_df <- rbind(select(tuned_par_rna, Property = task.id, Input.seq = seq_type, Mean.err = mean_error, Sd.err = sd_error),
      select(tuned_par_mcra, Property = task.id, Input.seq = seq_type, Mean.err = mean_error, Sd.err = sd_error)) %>% 
  droplevels %>% 
  mutate(Input.seq = )
  
error_df[["Input.seq"]]
