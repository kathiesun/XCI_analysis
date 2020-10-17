options(stringsAsFactors = FALSE)
library(tidyverse)

setwd("XCI_analysis")
source("ase_summary_source.R")

C_diet_4 = read.csv("../data/C_matrix_4diets.csv", header=F)
C_diet_2 = read.csv("../data/C_matrix_2diets.csv", header=F)
C_diet_3 = read.csv("../data/C_matrix_3diets.csv", header=F)

indiv_pups = list.files("../data/demographic_info/SP1_haplotypes", pattern="haploBlocks", full.names = T)

phased_CC_haplotype = lapply(indiv_pups, readRDS)
tmp = do.call("rbind", lapply(indiv_pups, function(x) unlist(strsplit(x, "_"))))
names(phased_CC_haplotype) = paste0("Pup.ID_", tmp[,ncol(tmp)-1])

seg_regions = readRDS("../data/demographic_info/SP1_segregating_regions_by_genotyping_perRIX.rds")

## Pup demographic info
allInfo <- read.csv("../data/demographic_info/summary_cegs-mnt_data.csv")
lab = allInfo %>% rename("RIX"="RRIX","CC.1"="CC1","CC.2"="CC2") %>%
  select(one_of("Pup.ID","RIX", "CC.1","CC.2")) %>%
  filter(RIX %in% c(1:4,6:10), !is.na(Pup.ID)) 
pupInfo <- allInfo %>% rename("RIX"="RRIX") %>%
  select(one_of("Pup.ID","RIX", "dir","Diet")) %>%
  mutate(Diet = gsub(" ","", Diet)) %>%
  filter(RIX %in% c(1:4,6:10), !is.na(Pup.ID)) 

data_kmers_lst = ratios_lst = list()

c="X"

## snp info
masterSnps = readRDS("../data/kmer_data/all_kmers_chrX.rds")[[1]]
masterSnps <- do.call("cbind", masterSnps)
masterSnps$seq.consensus <- paste0(masterSnps$seq.end5, masterSnps$seq.end3)

## count data 
ref <- unique(read.csv("../data/kmer_data/chrX_ref_counts_sp1.csv", header = T))
alt <- unique(read.csv("../data/kmer_data/chrX_alt_counts_sp1.csv", header = T))

data_kmers_lst[[paste(c)]] = process_and_plot(chr=c, 
                                              snp_info=masterSnps, 
                                              sample_info=pupInfo, 
                                              RIX_info=lab, 
                                              ref_counts=ref, alt_counts=alt, 
                                              phased_CC_haplotype=phased_CC_haplotype, 
                                              use_gene=F,
                                              problemPups = c(1404, 1716, 1371, 569, 1911, 1951, 1015),
                                              seg_regions=seg_regions)

data_kmers = data_kmers_lst[[paste(c)]]


ratios_lst[[paste(c)]] = run_jags_regress(data_kmers=data_kmers, 
                                          niter=50000, n.thin=5,  
                                          seg_regions=seg_regions,
                                          save_dir=NULL,
                                          STZ=T, use_gene=F,
                                          no_theta=F, alpha=NULL)
