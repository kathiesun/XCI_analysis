options(stringsAsFactors = FALSE)
library(tidyverse)

setwd("XCI_analysis")
source("ase_summary_source.R")
source("jags_kmers_source.R")


#########################################

xce_map = read.csv("../data/demographic_info/xce_info_CC.csv")
cegs_demo = read.csv("../data/demographic_info/summary_cegs-mnt_data.csv")
cegs_demo = cegs_demo[grep("CC", cegs_demo$RRIX),]
cegs_demo$sex = "F"
cegs_demo = cegs_demo %>% rename("sample"="Pup.ID","cross"="CCs", "treatment"="Diet",
                                 "mat_CC"="CC1", "pat_CC"="CC2","ccc"="CC_lab") 

count_files <- list.files("../data/kmer_data/SP2_counts", pattern = ".csv", full.names=T)
rem_string = c("_2020-01-10","_ref","_alt",".csv","XCounts_")
cegs_pups = gsub(paste(rem_string, collapse="|"), "", count_files)
cegs_pups = unlist(lapply(cegs_pups, function(x) strsplit(x, "/")[[1]][length(strsplit(x,"/")[[1]])]))
cegs_strains = do.call("rbind",strsplit(cegs_pups,"_"))
cegs_strains = data.frame(cegs_strains) %>% distinct()
colnames(cegs_strains) = c("iidsht","cross_rna")
cegs_strains$cross_rna = gsub("[A-Z]", "", cegs_strains$cross_rna)

#####  founder blocks from mrcas  #####
hap_pat <- readRDS("../data/demographic_info/SP2_haplotype_data.rds")

#####  X chr snps
snps = readRDS("../data/kmer_data/all_kmers_chrX.rds")[[1]]

ref_counts_f <- count_files[grep("ref", count_files)]
alt_counts_f <- count_files[grep("alt", count_files)]

ref_counts <- lapply(cegs_demo$sample, function(x){
  tmp = read.csv(ref_counts_f[grep(x, ref_counts_f)], header = F)
  colnames(tmp) = c("id", "refseq","dir","count")
  tmp %>% group_by(id, refseq) %>%
    summarize(refcount = sum(count)) -> tmp
  
  tmp$sample = unlist(strsplit(tmp$id,"_"))[c(T,F)]
  tmp %>% left_join(cegs_demo, by="sample") %>%
    left_join(snps$seq, by="refseq")
})
names(ref_counts) = cegs_demo$sample
alt_counts <- lapply(cegs_demo$sample, function(x){
  tmp = read.csv(alt_counts_f[grep(x, alt_counts_f)], header = F)
  colnames(tmp) = c("id", "altseq","dir","count")
  tmp %>% group_by(id, altseq) %>%
    summarize(altcount = sum(count)) -> tmp
  
  tmp$sample = unlist(strsplit(tmp$id,"_"))[c(T,F)]
  tmp %>% left_join(cegs_demo, by="sample") %>%
    left_join(snps$seq, by="altseq")
})
names(alt_counts) = cegs_demo$sample

all_counts = lapply(names(ref_counts), function(x){
  left_join(ref_counts[[x]], alt_counts[[x]], by="ProbeSeq") %>%
    ungroup() %>%
    dplyr::select(-contains(".x")) %>%
    dplyr::rename("CC1" = "mat_CC.y", "CC2" = "pat_CC.y")%>% 
    mutate(ratio = refcount / (refcount + altcount), sum = refcount + altcount) -> tmp
  colnames(tmp) = gsub(".y", "", colnames(tmp))
  tmp = tmp[order(tmp$Position),]
  
  tmp[,c("sample","refcount","altcount","ratio","sum", "treatment","CC1","CC2",
         "Gene","Chromosome","rsId","Position","TranscriptID")]
})

names(all_counts) = names(ref_counts)

lo=3
model_counts = lapply(all_counts, function(x){
  CCs = unlist(strsplit(c(unique(x$CC1), unique(x$CC2)), "/"))
  if(all(CCs %in% names(hap_pat))){
    hap_pat[[CCs[1]]] %>% 
      ungroup %>%
      filter(chr == "chrX") %>% filter(hap != "Z") %>%
      mutate(hap = unlist(lapply(lapply(strsplit(hap, ""), unique), 
                                 function(y) paste(y, collapse="")))) %>%
      mutate(chr = gsub("chr", "", chr)) %>%
      select(-grp) -> mat
    hap_pat[[CCs[2]]] %>% 
      ungroup %>%
      filter(chr == "chrX") %>% filter(hap != "Z") %>%
      mutate(hap = unlist(lapply(lapply(strsplit(hap, ""), unique), 
                                 function(y) paste(y, collapse="")))) %>%
      mutate(chr = gsub("chr", "", chr)) %>%
      select(-grp) -> pat
    if(nrow(pat) > 0 & nrow(mat) > 0){
      mf = lapply(x$Position, function(p) mat$hap[intersect(which(mat$start < p), which(mat$end > p))])
      x$CC1_hap = unlist(lapply(mf, function(y) ifelse(length(y)>0, y, "-")))
      pf = lapply(x$Position, function(p) pat$hap[intersect(which(pat$start < p), which(pat$end > p))])
      x$CC2_hap = unlist(lapply(pf, function(y) ifelse(length(y)>0, y, "-")))
      x %>% filter(CC1_hap != "-", CC2_hap != "-") -> x
      xce_sec = x %>% filter(Position > 99e6, Position < 104e6)
      x %>% filter(CC1_hap != CC2_hap) %>%
        mutate(CC1_xce_hap = paste(unlist(strsplit(unique(xce_sec$CC1_hap),"")), collapse="/"),
               CC2_xce_hap = paste(unlist(strsplit(unique(xce_sec$CC2_hap),"")), collapse="/")) %>%
        ungroup() -> x
      
      dict = data.frame(let = LETTERS[1:8], num = 1:8)
      pos_mat_sdp_1 = apply(x, 1, function(y)
        ifelse(any(snps$sdp[match(as.character(unlist(y["rsId"])), snps$seq$rsId),
                            dict$num[match(unlist(strsplit(unlist(y["CC1_hap"]),"")), dict$let)]] != "0"), 1, 0))
      pos_mat_sdp_2 = apply(x, 1, function(y) 
        ifelse(any(snps$founder_alt[match(as.character(unlist(y["rsId"])),  snps$seq$rsId), 
                                    dict$num[match(unlist(strsplit(unlist(y["CC1_hap"]),"")), dict$let)]] > lo), 1, 0))
      pos_pat_sdp_1 = apply(x, 1, function(y)
        ifelse(any(snps$sdp[match(as.character(unlist(y["rsId"])), snps$seq$rsId),
                            dict$num[match(unlist(strsplit(unlist(y["CC2_hap"]),"")), dict$let)]] != "0"), 1, 0))
      pos_pat_sdp_2 = apply(x, 1, function(y) 
        ifelse(any(snps$founder_alt[match(as.character(unlist(y["rsId"])),  snps$seq$rsId), 
                                    dict$num[match(unlist(strsplit(unlist(y["CC2_hap"]),"")), dict$let)]] > lo), 1, 0))
      flag = union(which(pos_pat_sdp_1 != pos_pat_sdp_2), which(pos_mat_sdp_1 != pos_mat_sdp_2))
      if(length(flag) > 0) x$refcount[flag] = x$altcount[flag] = NA
      
      x$mat_sdp = sapply(1:length(pos_mat_sdp_1), function(n) 
        ifelse(any(pos_mat_sdp_1[n] == 1, pos_mat_sdp_2[n] == 1), 1, 0)) 
      
      x$pat_sdp = sapply(1:length(pos_pat_sdp_1), function(n) 
        ifelse(any(pos_pat_sdp_1[n] == 1, pos_pat_sdp_2[n] == 1), 1, 0))
      
      x %>% filter(mat_sdp != pat_sdp, !is.nan(ratio)) %>%
        mutate(Xist = ifelse(Gene == "Xist", T, F)) -> x
      
      if(nrow(x) > 0){
        for(i in 1:nrow(x)){
          x$fin1[i] = as.numeric(gsub(" ", "", ifelse(x$mat_sdp[i] == 0, ifelse(x$refcount[i] >= lo, x$refcount[i], NA), 
                                                      ifelse(x$altcount[i] >= lo, x$altcount[i], NA)))) 
          x$fin2[i] = as.numeric(gsub(" ", "", ifelse(x$pat_sdp[i] == 0, ifelse(x$refcount[i] >= lo, x$refcount[i], NA), 
                                                      ifelse(x$altcount[i] >= lo, x$altcount[i], NA)))) 
        }
        if(length(which(x$Xist)) > 0) {
          tmp = x$fin2[which(x$Xist)]
          x$fin2[which(x$Xist)] = x$fin1[which(x$Xist)]
          x$fin1[which(x$Xist)] = tmp
        }
        x %>% plyr::rename(c("Gene"="seq.gene", "sample"="Pup.ID", "treatment"="diet")) %>%
          filter(!is.na(fin1), !is.na(fin2)) %>%
          mutate(ratio = fin1/sum) -> x
      }
      
      x            
    }
  }
})


model_counts$C048 = NULL
## removed due to no non-NA ratios -- perhaps only 1 copy of X

xce_strength = c("NOD","a","f","e","NZO","b","c","d")
count_data = list()
rix_summary = list()

for(c in unique(cegs_demo$ccc)){    
  fem = cegs_demo %>% filter(ccc == c) 
  model_rix = lapply(fem$sample, function(x) model_counts[[x]])
  model_rix = do.call("rbind", model_rix)
  
  if(!is.null(model_rix)){
    model_rix$trt = ifelse(model_rix$diet == "placebo", -0.5, 
                           ifelse(model_rix$diet == "drug", 0.5, NA))
    model_rix = model_rix %>% filter(!is.na(fin1), !is.na(fin2))
    CC2_xce_haps = unique(model_rix$CC2_xce_hap)
    CC1_xce_haps = unique(model_rix$CC1_xce_hap)
    xce_founds = c(unique(model_rix$CC1_xce_hap),unique(model_rix$CC2_xce_hap))
    
    if(nchar(CC2_xce_haps) > 2) CC2_xce_haps = unlist(strsplit(unique(model_rix$CC2_xce_hap),"/"))
    if(nchar(CC1_xce_haps) > 2) CC1_xce_haps = unlist(strsplit(unique(model_rix$CC1_xce_hap),"/"))
    CC1_xce_al = factor(unique(xce_map$xce[which(xce_map$founder %in% CC1_xce_haps)]), levels=xce_strength)
    CC2_xce_al = factor(unique(xce_map$xce[which(xce_map$founder %in% CC2_xce_haps)]), levels=xce_strength)
    
    which_numer = ifelse(any(as.numeric(CC1_xce_al) < as.numeric(CC2_xce_al)), "left", "right")
    
    
    model_rix$CC1_xce_al = paste(CC1_xce_al, collapse="/")
    model_rix$CC2_xce_al = paste(CC2_xce_al, collapse="/")
    model_rix$D = apply(model_rix, 1, function(x) ifelse(any(x["CC1_xce_al"]=="n", x["CC2_xce_al"] == "n"), 1, 0))
    
    model_rix = model_rix %>% filter(!is.nan(ratio)) %>%
      arrange(Pup.ID, Position) %>% 
      mutate(RRIX = c, PUP_GENE = paste(Pup.ID, seq.gene, sep="_"),
             Xist = ifelse(seq.gene == "Xist", T, F))
    
    if(which_numer == "right"){
      model_rix$which_numer = "right"
      tmp = model_rix$fin1
      model_rix$fin1 = model_rix$fin2
      model_rix$fin2 = tmp
      model_rix$CC_lab = paste(gsub("[A-Z]|/","",toupper(model_rix$CC2)),gsub("[A-Z]|/","",toupper(model_rix$CC1)), sep="/")
      
      model_rix = model_rix %>% mutate(ratio = fin1 / sum)
    } else {
      model_rix$which_numer = "left"
      model_rix$CC_lab = paste(gsub("[A-Z]|/","",toupper(model_rix$CC1)),gsub("[A-Z]|/","",toupper(model_rix$CC2)), sep="/")
    }
    count_data[[c]] = model_rix 
    count_data[[c]] %>% group_by(Pup.ID, which_numer, RRIX) %>%
      filter(sum > 0) %>%
      summarise(sum=sum(sum), fin1=sum(fin1), fin2=sum(fin2), 
                quan = paste(rat_range = round(quantile(ratio, c(0.25, 0.75), na.rm=T),2), collapse="-"), 
                rat_mean = mean(ratio, na.rm=T)) %>%
      mutate(rat_raw = fin1/sum) -> rix_summary[[c]]
  }
}

count_data_mat = do.call("rbind", count_data)
count_data_mat$RRIX = factor(count_data_mat$RRIX, levels=unique(count_data_mat$RRIX))
count_data_mat = count_data_mat %>% 
  mutate(dir = 1,
         RIX = RRIX, CCs = RRIX, 
         logSum = log(sum)) %>%
  plyr::rename(c("Chromosome" = "seq.Chromosome", "seq.gene" = "seq.Gene",
                 "diet" = "Diet", "Position" = "seq.Position"))

reg = list()

for (i in 1:length(unique(cegs_demo$ccc))){    
  c = unique(cegs_demo$ccc)[i]
  terms = "RRIX"
  quant = NULL
  model_rix = count_data_mat %>% filter(CCs == c) 
  if(nrow(model_rix) > 0) {
    trts = unique(model_rix$trt)  
    if(length(which(is.na(trts))>0)) trts = trts[-which(is.na(trts))]
    if (length(trts) > 1) quant = "TRT"
    niter = 100000 
    reg[[c]] <- jags.genes.run(data=model_rix, mu_g0=0.5, niter=niter, 
                               n.thin = 5, quant=quant,  #quant,  
                               terms=terms, STZ=F, use_gene = F, 
                               add_impute = "\nfor(p in 1:nP){\nD_temp[p] ~ dbern(0.5)\nD[p] <- -0.5 + D_temp[p]\n}\n")
  }
}

rix_summary = do.call("rbind",rix_summary)

