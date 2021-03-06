# this version was using wrong data types

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
# data for Slawek
# setdiff(both_mcra_rna, as.character(conditions_dat[["Name"]])) %>% 
#   sort() %>% 
#   cat(sep = "\n")

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

configureMlr(show.info = FALSE)

benchmark_ngram_length <- pblapply(2L:7, function(ngram_length) {
  lapply(c(0.25, 0.5, 1), function(feature_prop) {
    lapply(c("both", "mcra_seqs", "rna_seqs"), function(ith_seqs) {
      
      if(ith_seqs == "both") {
        
        ith_seqs_data <- seq_data[["mcra_seqs"]] 
        ith_seqs_data <- ith_seqs_data[rownames(ith_seqs_data) %in% all_three, ]
        
        ngram_matrix <- count_ngrams(ith_seqs_data, ngram_length, u = c("a", "c", "g", "t"), scale = TRUE)
        
        mcra_ngrams <- ngram_matrix %>% 
          as.matrix %>% 
          data.frame %>% 
          mutate(source = rownames(ith_seqs_data)) %>% 
          group_by(source) %>% 
          summarise_all(mean) 
        
        colnames(mcra_ngrams)[-1] <- paste0("mcra_", colnames(mcra_ngrams)[-1])
        
        ith_seqs_data <- seq_data[["rna_seqs"]] 
        ith_seqs_data <- ith_seqs_data[rownames(ith_seqs_data) %in% all_three, ]
        
        ngram_matrix <- count_ngrams(ith_seqs_data, ngram_length, u = c("a", "c", "g", "t"), scale = TRUE)
        
        rna_ngrams <- ngram_matrix %>% 
          as.matrix %>% 
          data.frame %>% 
          mutate(source = rownames(ith_seqs_data)) %>% 
          group_by(source) %>% 
          summarise_all(mean) 
        
        colnames(rna_ngrams)[-1] <- paste0("rna_", colnames(rna_ngrams)[-1])
        
        normalized_ngrams <- cbind(rna_ngrams, mcra_ngrams[, -1])
        
      } else {
        ith_seqs_data <- seq_data[[ith_seqs]] 
        ith_seqs_data <- ith_seqs_data[rownames(ith_seqs_data) %in% all_three, ]
        
        ngram_matrix <- count_ngrams(ith_seqs_data, ngram_length, u = c("a", "c", "g", "t"), scale = TRUE)
        
        normalized_ngrams <- ngram_matrix %>% 
          as.matrix %>% 
          data.frame %>% 
          mutate(source = rownames(ith_seqs_data)) %>% 
          group_by(source) %>% 
          summarise_all(mean) 
      }
      
      bench_res <- lapply(c("growth_doubl", "growth_rate", "mean_ogt", 
                            "mean_ogn", "mean_ogp"),
                          function(ith_condition) {
                            dat <- conditions_dat[, c("Name", ith_condition)] %>% 
                              inner_join(normalized_ngrams, by = c("Name" = "source")) %>% 
                              select(-Name)

                            predict_ngrams <- makeRegrTask(id = ith_condition, 
                                                           data = dat, 
                                                           target = ith_condition)
                            
                            filtered_ngrams <- filterFeatures(predict_ngrams, method = "linear.correlation", perc = feature_prop)
                            n_features <- filtered_ngrams[["task.desc"]][["n.feat"]][["numerics"]]
                            mtry_possibilities <- unique(round(c(n_features/4, n_features/3, n_features/2), 0))
                            
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
                            
                            dat <- getNestedTuneResultsOptPathDf(nested_cv) 
                            group_by(dat, mtry, num.trees, min.node.size) %>% 
                              summarise(mean = mean(mse.test.mean))
                          })
      
      lapply(bench_res, function(i) 
        data.frame(i)) %>% 
        do.call(rbind, .) %>% 
        mutate(ngram_length = ngram_length,
               seq_type = ith_seqs,
               feature_prop = feature_prop)  
    }) %>% 
      do.call(rbind, .)
  }) %>%
    do.call(rbind, .)
}) %>%
  do.call(rbind, .)

write.csv(benchmark_ngram_length, file = "./results/ngram_benchmark_multi.csv", row.names = FALSE)


# group_by(task.id) %>% 
#   summarise(mse = mean(mse)) %>% 
#   mutate(error = sqrt(mse),
#          ngram_length = ngram_length)

