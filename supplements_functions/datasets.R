library(dplyr)
library(seqinr)
library(gridExtra)
library(biogram)

source("./functions/validate_seqs.R")

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


both_mcra_rna <- intersect(unique(rownames(rna_seqs)), unique(rownames(mcra_seqs)))
both_mcra_conditions <- intersect(as.character(conditions_dat[["Name"]]), unique(rownames(mcra_seqs)))
both_rna_conditions <- intersect(as.character(conditions_dat[["Name"]]), unique(rownames(rna_seqs)))
all_three <- intersect(as.character(conditions_dat[["Name"]]), both_mcra_rna)

venn_dat <- lapply(names(rna_seqs), function(rna_seq_name)
  lapply(names(mcra_seqs), function(mcra_seq_name)
    lapply(c("mcra_seqs", "rna_seqs"), function(ith_seqs) {
      rna_seq <- rna_seqs[[rna_seq_name]]
      mcra_seq <- mcra_seqs[[mcra_seq_name]]
      
      both_mcra_rna <- intersect(unique(rownames(rna_seq)), unique(rownames(mcra_seq)))
      both_mcra_conditions <- intersect(as.character(conditions_dat[["Name"]]), unique(rownames(mcra_seq)))
      both_rna_conditions <- intersect(as.character(conditions_dat[["Name"]]), unique(rownames(rna_seq)))
      all_three <- intersect(as.character(conditions_dat[["Name"]]), both_mcra_rna)
      
      data.frame(rna = rna_seq_name,
                 mcra = mcra_seq_name,
                 rna_seq = length(unique(rownames(rna_seq))),
                 mcra_seq = length(unique(rownames(mcra_seq))),
                 conditions = length(as.character(conditions_dat[["Name"]])),
                 both_mcra_rna = length(both_mcra_rna),
                 both_mcra_conditions = length(both_mcra_conditions),
                 both_rna_conditions = length(both_rna_conditions),
                 all_three = length(all_three),
                 species = paste0(all_three, collapse = ", ")
                 )
    }) %>% 
      do.call(rbind, .)
  ) %>% 
    do.call(rbind, .)
) %>% 
  do.call(rbind, .) %>% 
  filter(!duplicated(.))

library(VennDiagram)


venn_plots <- lapply(1L:nrow(venn_dat), function(ith_row_id) {
  ith_row <- venn_dat[ith_row_id, ]

  draw.triple.venn(ith_row[["rna_seq"]],
                   ith_row[["mcra_seq"]],
                   ith_row[["conditions"]],
                   ith_row[["both_mcra_rna"]],
                   ith_row[["both_mcra_conditions"]],
                   ith_row[["both_rna_conditions"]],
                   ith_row[["all_three"]],
                   # c(as.character(ith_row[["rna"]]),
                   #   as.character(ith_row[["mcra"]]),
                   #   "Known conditions"),
                   c("16 rRNA",
                     "mcrA",
                     "Known culturing conditions"),
                   fill = c("#f1a340", "#998ec3", "#f7f7f7"),
                   ind = FALSE)
})

# for(i in 1L:length(venn_plots)) {
#   cat("\n\n##", i, "\n\n", sep = " ")
#   grid.draw(venn_plots[[i]])
#   grid.newpage()
#   cat("\n\n\\pagebreak\n\n")
# }


#select(venn_dat, rna, mcra, rna_seq, mcra_seq, all_three)
#filter(select, rna == "RNA1", mcra == "McrA1")