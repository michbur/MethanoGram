extract_ngrams <- function(x, len = 4) {
  lapply(x, tolower) %>% 
    list2matrix %>%
    count_ngrams(len, u = c("a", "c", "g", "t"), scale = TRUE) %>% 
    as.matrix %>% 
    data.frame
}

pred_vals <- function(models, ngrams, seq_names, seq_type) {

  res <- lapply(models, function(single_model) {
    raw_names <- single_model[["features"]]
    ngram_lens <- unique(nchar(decode_ngrams(sub("rna_", "", raw_names))))
    
    model_ngrams <- do.call(cbind, lapply(ngram_lens, function(ith_len) {
      extract_ngrams(ngrams, ith_len)
    }))
    
    names(model_ngrams) <- paste0("rna_", names(model_ngrams))
    #model_ngrams <- extract_ngrams(ngrams, nchar(decode_ngrams(sub("rna_", "", raw_names))[1]))


    
    predict(single_model, newdata = model_ngrams[, single_model[["features"]]])
  }) %>% 
    data.frame  

  
  colnames(res) <- c("Growth doubling time",
                     "Optimal temperature",
                     "Growth rate", 
                     "Optimal NaCl",
                     "Optimal pH")
  res <- res[, c(1, 3, 2, 5, 4)]
  res[["Growth doubling time"]] <- exp(res[["Growth doubling time"]])
  res[["Optimal temperature"]] <- exp(res[["Optimal temperature"]])
  
  #res <- res[c(1, 2, 6, 3, 14L:13, 8, 5, 12:11, 7, 4, 10L:9), ]
  rownames(res) <- NULL

  data.frame(Name = seq_names, 
             Input.seq = ifelse(seq_type == "rna", "16S rRNA", "mcrA"),
             res)
}
