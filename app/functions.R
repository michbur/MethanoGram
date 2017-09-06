cont_feats <- c("GrowthDoublingTime", "GrowthRate", "maximalGrowthTemperature", 
                "maximumGrowthNACL", "maximumGrowthPH", "minimalGrowthTemperature", "minimumGrowthNACL", 
                "minimumGrowthPH", "optimalGrowthNACLMaximal", "optimalGrowthNACLMinimal", 
                "optimalGrowthPHMaximal", "optimalGrowthPHMinimal", "optimalGrowthTemperatureMaximal", 
                "optimalGrowthTemperatureMinimal")

extract_ngrams <- function(x, len = 4) {
  lapply(x, tolower) %>% 
    list2matrix %>%
    count_ngrams(len, u = c("a", "c", "g", "t"), scale = TRUE) %>% 
    as.matrix %>% 
    data.frame
}

pred_vals <- function(models, ngrams, seq_names) {

  res <- lapply(models, function(single_model) {
    predict(single_model, newdata = extract_ngrams(ngrams, nchar(decode_ngrams(single_model[["features"]])[1])))
  }) %>% 
    data.frame %>% 
    t %>% 
    data.frame(Property = sapply(models, function(i) i[["task.desc"]][["id"]]), .) 
    
  #res <- res[c(1, 2, 6, 3, 14L:13, 8, 5, 12:11, 7, 4, 10L:9), ]
  rownames(res) <- NULL
  colnames(res)[-1] <- seq_names
  
  res
}
