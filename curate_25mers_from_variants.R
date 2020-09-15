###################################
#   Choose kmers based on founder states
###################################

options(stringsAsFactors = FALSE)
setwd("C:/Users/Kathie/Dropbox (ValdarLab)/MaternalDiet/XCI_paper")


library(tidyverse)
source("masterKmers_source.R")
source("plot.hpd.R")
source("plotting_ratio_regression.R")

### newer w/TSL1 annotation staring file
exonsNew <- read.csv("data/kmer_data/SNPsInExonsTSL1.csv", 
                     colClasses=c("FounderSDP"="character"))
founders <- c("AJ","C57BL6JN","X129S1","NOD","NZOdna","CASTdna","PWKdna","WSBdna")
founderCols <- match(founders,colnames(exonsNew))

deets <- seq(from=grep("ExonID", colnames(exonsNew)), to=grep("FounderSDP", colnames(exonsNew)))
exons <- exonsNew[,c(deets,founderCols, grep("CC", colnames(exonsNew)))]

norm_file <- read.table("data/kmer_data/normalization_factors_founders.csv", sep=",", header=T)
order <- unlist(lapply(unlist(lapply(strsplit(norm_file$founder,"/"), function(y) y[1])), function(x) grep(x, founders)))
norm_file <- norm_file[order(order),]
norm_file$founder = c(founders)



exons_use <- exons %>% 
  filter(Chromosome == "X") %>%
  distinct()
masterSnps <- curateKmers_founders(exons=exons_use, norm_file = norm_file)

## two versions of curateKmers function in masterKmers_source.R
## one is based on founder data, the other is based on CC data
## manuscript included kmers using curateKmers_founders


