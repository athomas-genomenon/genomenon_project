#working directory
setwd("~/Desktop/Genomenon")

#csv import
evidence <- read.csv("evidence-7293-2026-01-08.csv", header = TRUE, sep = ",")

#install packages
install.packages("ggplot2")
install.packages("dplyr")
library(ggplot2)
library(dplyr)
library(stringr)

#removing columns we don't want
cleaned_evidence = subset(evidence, select = c("gene_symbol", "mutation", "classification", "all_cats"))

#counting variables in columns
cleaned_evidence %>% count(classification)

#trying to figure out how to have it count all the variables in all_cats that have a specific classification within their string
#for loop that runs through every variable in all_cats and looks for a specific string? unsure

#for (x in cleaned_evidence$all_cats) {
  #str_detect(x, "pp2")
#}

#I think this works to get the counts of the variants that have certain classifications applied?
pm2_count <- cleaned_evidence$all_cats[str_detect(cleaned_evidence$all_cats, "pm2")]
pp2_count <- cleaned_evidence$all_cats[str_detect(cleaned_evidence$all_cats, "pp2")]
length(pm2_count)
length(pp2_count)

#should unweighted cases be removed from this list?
  