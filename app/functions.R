cont_feats <- c("GrowthDoublingTime", "GrowthRate", "maximalGrowthTemperature", 
                "maximumGrowthNACL", "maximumGrowthPH", "minimalGrowthTemperature", "minimumGrowthNACL", 
                "minimumGrowthPH", "optimalGrowthNACLMaximal", "optimalGrowthNACLMinimal", 
                "optimalGrowthPHMaximal", "optimalGrowthPHMinimal", "optimalGrowthTemperatureMaximal", 
                "optimalGrowthTemperatureMinimal")

extract_ngrams <- function(x) {
  lapply(x, tolower) %>% 
    list2matrix %>%
    count_ngrams(4, u = c("a", "c", "g", "t"), scale = TRUE) %>% 
    as.matrix %>% 
    data.frame
}

pred_vals <- function(models, ngrams, seq_names) {
  res <- lapply(models, function(single_model)
    predict(single_model[["model"]], newdata = ngrams)
  ) %>% 
    data.frame %>% 
    t %>% 
    data.frame(Property = cont_feats, .) 
    
  res <- res[c(1, 2, 6, 3, 14L:13, 8, 5, 12:11, 7, 4, 10L:9), ]
  rownames(res) <- NULL
  colnames(res)[-1] <- seq_names
  
  res
}
