validate_seqs <- function(seqs) {
  validated_seq_names <- lapply(seqs, function(i) 
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
  only_standard_nucleotides <- sapply(seqs, function(i)
    !any(c("r", "n", "b", "s", "m", "d", "w", "y", "k", "v") %in% i))
  
  validated_seqs <- list2matrix(seqs[validated_seq_names[["is.reference"]] & only_standard_nucleotides])
  
  rownames(validated_seqs) <- strsplit(rownames(validated_seqs), split = "|", fixed = TRUE) %>% 
    sapply(function(i) i[2]) %>% 
    strsplit(split = "_") %>% 
    sapply(function(i) paste0(i[1L:2], collapse = " "))
  
  validated_seqs
}
