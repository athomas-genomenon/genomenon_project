#working directory
setwd("~/Desktop/Genomenon")

#csv import
evidence <- read.csv("evidence-7293-2026-01-08.csv", header = TRUE, sep = ",")

#install packages
install.packages("ggplot2")
install.packages("dplyr")
library(ggplot2)
library(dplyr)
