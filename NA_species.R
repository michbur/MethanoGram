library(dplyr)
library(reshape2)

as_numeric <- function(x) as.numeric(as.character(x))

dat <- read.csv(file = "./data/export.csv") %>% 
  select(Name, 
         Growth.doubling.time..h., 
         Growth.rate, 
         Min..optimal.growth.temp., 
         Max..optimal.growth.temp., 
         Min..optimal.growth.NaCl, 
         Max..optimal.growth.NaCl, 
         Min..optimal.growth.pH, 
         Max..optimal.growth.pH) %>% 
  mutate_at(c("Growth.doubling.time..h.", "Growth.rate", "Min..optimal.growth.temp.", 
              "Max..optimal.growth.temp.", "Min..optimal.growth.NaCl", "Max..optimal.growth.NaCl", 
              "Min..optimal.growth.pH", "Max..optimal.growth.pH"), as_numeric)


melt(dat) %>% 
  group_by(Name) %>% 
  summarise(count = sum(is.na(value))) %>% 
  filter(count > 0) %>% 
  select(Name) %>% 
  unlist %>% 
  as.character %>% 
  cat(sep = "\n")
