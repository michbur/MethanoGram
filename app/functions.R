extract_ngrams <- function(x, len = 4) {
  lapply(x, tolower) %>% 
    list2matrix %>%
    count_ngrams(len, u = c("a", "c", "g", "t"), scale = TRUE) %>% 
    as.matrix %>% 
    data.frame
}

pred_vals <- function(models, ngrams, seq_names, seq_type) {

  res <- lapply(models, function(single_model) {
    model_ngrams <- extract_ngrams(ngrams, nchar(decode_ngrams(single_model[["features"]])[1]))

    predict(single_model, newdata = model_ngrams[, single_model[["features"]]])
  }) %>% 
    data.frame  
  
  colnames(res) <- c("Growth rate", 
                     "Optimal pH", 
                     "Growth doubling time", 
                     "Optimal temperature",
                     "Optimal NaCl")
  res <- res[, c(1, 3, 4, 2, 5)]
  #res <- res[c(1, 2, 6, 3, 14L:13, 8, 5, 12:11, 7, 4, 10L:9), ]
  rownames(res) <- NULL

  data.frame(Name = seq_names, 
             Input.seq = ifelse(seq_type == "rna", "16S rRNA", "mcrA"),
             res)
}
