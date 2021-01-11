options(stringsAsFactors = FALSE)
options(scipen=999)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(scales)

setwd("XCI_analysis")
source("ase_summary_source.R")

####################################
##        General datasets        ##
####################################

## updated Xce info
xce_map = read.csv("FileS7_xce_info_CCstrains.csv")

## SP1/2 population information 
xce_pups = read.csv("FileS8_demographic_data_SP1-2.csv")

## X chr genes
all_genes = read.csv("FileS9_X_chr_markers_map.csv")

## load X snps
snps = read.csv("FileS10_all_kmers_chrX.csv")

## SP1 cc haplotypes
indiv_pups = list.files("FileS11_SP1_haplotypes", pattern="haploBlocks", full.names = T)
mnt_haplotypes = lapply(indiv_pups, readRDS)
tmp = do.call("rbind", lapply(indiv_pups, function(x) unlist(strsplit(x, "_"))))
names(mnt_haplotypes) = paste0("Pup.ID_", tmp[,ncol(tmp)-1])

## SP2 cc haplotypes
hap_pat <- readRDS("FileS12_SP2_haplotype_data.rds")

## SP1 regression results
mnt_regSum = readRDS("FileS13_chrX_summary_sp1.rds")
mnt_df    = readRDS("FileS14_chrX_data_sp1.rds")

## SP2 regression results
cegs_regSum = readRDS("FileS15_chrX_summary_sp2.rds")
cegs_df     = readRDS("FileS16_chrX_data_sp2.rds")


## combine SP data
cegs_tmp_df <- rapply(cegs_df, as.character, classes="factor", how="replace")
mnt_tmp_df  <- rapply(mnt_df, as.character, classes="factor", how="replace")
colnames(cegs_tmp_df) = toupper(colnames(cegs_tmp_df))
colnames(mnt_tmp_df) = toupper(colnames(mnt_tmp_df))
keepCols = colnames(mnt_tmp_df)[match(colnames(cegs_tmp_df), colnames(mnt_tmp_df))]
keepCols = keepCols[which(!is.na(keepCols))]

df_comb      = rbind(mnt_tmp_df %>% select(keepCols), cegs_tmp_df %>% select(keepCols))
summary_comb = c(mnt_regSum, cegs_regSum)


## various lookup tables
xce_df = data.frame(founder = LETTERS[1:8],
                    xce_al = c("a","b","a","n","z","c","e","b"))
dict = data.frame(let = LETTERS[1:13], num = 1:13,
                  found = c("AJ","B6","129S1","NOD","NZO","CAST","PWK","WSB", "ALS","BALB","DBA","C3H","FVB"),
                  xce = c("a","b","a","n","z","c","e","b","b","a","b","a","f"))

colorpal = data.frame(violet = "#3b0a53",blue = "#4b77ba",skyblue = "#74affd", green = "#72b57a", slategray = "#2F4F4F",
                      magenta = "#a61f82", pink = "#c987be",yellow = "#f4e53c", gray = "#98999c", jet = "#2B2D30")
##################################
## sample-wide proportion plots ##
##################################

pup_skews = do.call("rbind", lapply(mnt_regSum, function(x) x$summary$mu_p))
pup_skews$Pup.ID = unlist(strsplit(rownames(pup_skews),"[.]"))[c(F,T)]

cegs_skews = do.call("rbind", lapply(cegs_regSum, function(x) x$summary$mu_p))
cegs_skews$Pup.ID = unlist(strsplit(rownames(cegs_skews),"[.]"))[c(F,T)]

rix_mnt_skews = do.call("rbind", lapply(mnt_regSum, function(x) x$summary$mu_r))
rix_mnt_skews$RRIX = unlist(strsplit(rownames(rix_mnt_skews),"_"))[c(F,T)]
rix_cegs_skews = do.call("rbind", lapply(cegs_regSum, function(x) x$summary$mu_r))
rix_cegs_skews$RRIX = rownames(rix_cegs_skews)

rix_skews = rbind(rix_mnt_skews, rix_cegs_skews)
pup_skews = rbind(pup_skews, cegs_skews)

xce_plot_df = left_join(xce_pups, pup_skews, by="Pup.ID")

## summary data for rix-wide boxes
xce_pups %>% 
  ungroup() %>%
  group_by(CC_lab) %>%
  mutate(rix_total = sum(pup_total),n=n()) %>%
  select(-c("Pup.ID","dir","Diet","pup_total")) %>%
  distinct() -> xce_rix

## prepare rix-wide data
rix_plot_df = left_join(xce_rix, rix_skews, by="RRIX")

rix_plot_df$grp = apply(rix_plot_df,1,function(x) 
  ifelse(x["which_numer"] == "right",paste(x["CC2_xce_al"], "vs", x["CC1_xce_al"]),
         paste(x["CC1_xce_al"], "vs", x["CC2_xce_al"])))

grp_order = c("NOD vs NOD","NOD vs NZO","NOD vs a","NOD vs b",   
              "NOD vs e","NZO vs NZO","NZO vs b","NZO vs c",
              "a vs a","a vs e","a vs b","a vs c",  
              "e vs b","e/b vs b","b vs b")
under_50 = c(F,F,F,F,F,F,F,   T,F,T,T,T,T,T,F)
rix_plot_df$grp = factor(rix_plot_df$grp, levels = grp_order, ordered = T)
rix_plot_df = rix_plot_df %>% arrange(grp, Mean)
rix_plot_df$CC_lab = factor(rix_plot_df$CC_lab, levels=unique(rix_plot_df$CC_lab), ordered=T)   
rix_plot_df$under_50 = apply(rix_plot_df, 1, function(x){
  tmp = under_50[which(grp_order == x[["grp"]])]
  ifelse(tmp, ifelse(x[["upper"]] < 0.5, T, F), 
         ifelse(x[["upper"]] < 0.5, F, T))
})

rix_plot_comb = rix_plot_df %>% 
  mutate(rix_sum = 1,
         hap_grp = ifelse(which_numer == "left", paste0(CC1_xce_hap,"/",CC2_xce_hap),
                          paste0(CC2_xce_hap,"/",CC1_xce_hap))) %>%
  #alpha of rix-level data
  select("CC_lab","rix_total","hap_grp","Mode","Median","Mean","lower","upper",
         "rix_sum","grp","under_50") %>%
  rename("pup_total" = "rix_total") %>%
  arrange(CC_lab)
rix_plot_comb$under_50[which(rix_plot_comb$CC_lab %in% c("CC028/CC025"))] = F
## determine facet groups

## combine rix-wide and sample-wide data to plot in same figure
xce_plot_df[,which(colnames(xce_plot_df) == "Mode"):which(colnames(xce_plot_df) == "upper.1")] = 
  apply(xce_plot_df[,which(colnames(xce_plot_df) == "Mode"):which(colnames(xce_plot_df) == "upper.1")], 2, as.numeric)
xce_plot_df$CC_lab = factor(xce_plot_df$CC_lab, levels=levels(rix_plot_comb$CC_lab))
xce_plot_df$hap_grp = rix_plot_comb$hap_grp[match(xce_plot_df$CC_lab, rix_plot_comb$CC_lab)]
xce_plot_df$grp = rix_plot_comb$grp[match(xce_plot_df$CC_lab, rix_plot_comb$CC_lab)]
xce_plot_df$under_50 = rix_plot_comb$under_50[match(xce_plot_df$CC_lab, rix_plot_comb$CC_lab)]

xce_plot_df = xce_plot_df %>% arrange(CC_lab, Mean)
xce_plot_df$rix_sum = 0

## set plot order, a and b panels
xce_plot_df$plot = 1
xce_plot_df$plot[grep("NOD|NZO",xce_plot_df$grp)] = 0
xce_plot_df$x_val = 1:nrow(xce_plot_df)


rix_plot_comb$plot = xce_plot_df$plot[match(rix_plot_comb$CC_lab, xce_plot_df$CC_lab)]        
rix_plot_comb$x_val = unlist(lapply(as.character(levels(rix_plot_comb$CC_lab)), function(x) 
  mean((xce_plot_df %>% filter(CC_lab == x))$x_val)
))

combine_plot_df = data.frame(xce_plot_df) %>% 
  dplyr::select("CC_lab","pup_total","hap_grp","Mode","Median","Mean",
         "lower","upper","rix_sum","grp","plot","x_val","under_50") %>%
  rbind(data.frame(rix_plot_comb)) %>%
  mutate(shape = ifelse(rix_sum == 0, F, T))

combine_plot_df = combine_plot_df %>% arrange(CC_lab)
combine_plot_df$xlab = paste(combine_plot_df$CC_lab, combine_plot_df$hap_grp,
                             sep="::")

combine_plot_df$CC_lab = gsub("CC","",combine_plot_df$CC_lab)

plot_data = combine_plot_df %>% filter(plot == 0)

cols = c("TRUE" = colorpal$slategray,"FALSE" = colorpal$magenta)

pfin = ggplot(data=plot_data, aes(y=Mean, x=x_val, col=under_50, 
                                  alpha=rix_sum, shape=shape)) +  
  scale_color_manual(values = cols, name = "As expected?") + 
  geom_point(aes(size=as.numeric(pup_total))) + 
  geom_hline(yintercept=0.5, linetype="dotted") + 
  geom_segment(aes(y=lower, yend=upper, x=x_val, xend=x_val)) + 
  scale_x_continuous(name="CC-RIX",
                   breaks=plot_data$x_val[which(plot_data$rix_sum == 1)],
                   labels=as.character(plot_data$xlab[which(plot_data$rix_sum == 1)])) +   #breaks=seq(1,1000)
  scale_y_continuous(limits = c(0, 0.8), breaks=seq(0.1,0.8, by=0.1), labels = seq(0.1,0.8, by=0.1)) + 
  scale_shape_manual(name="Estimate", labels=c("Sample","RIX"), values=c(15, 19)) + 
  scale_alpha_continuous(name=NULL, NULL, range=c(0.3,1)) + 
  theme_classic() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 6.5)) +
  facet_grid(grp~., scales="free", space="free") + 
  scale_size_continuous(name=c("Total counts"),labels = comma) + 
  coord_flip()


pfin


##################################
##    heterozygous ideograms    ##
##################################


## ordering for CC-RIX in plot
xce_guide = xce_guide %>% arrange(grp)
nod_rix = xce_guide$CCs[which(xce_guide$grp %in% c("NOD vs stronger","NOD vs NOD"))]
nzo_rix = xce_guide$CCs[which(xce_guide$grp %in% c("NZO vs c","NZO vs b","NZO vs NZO"))]

## allele probabilities for mnt and cegs
alleleprobs_minimuga <- readRDS("FileS17_interp_alleleprobs_sp1.rds")
alleleprobs_cegs = readRDS("FileS18_interp_alleleprobs_sp2.rds")

## all markers with probabilities
X_markers = dimnames(alleleprobs_minimuga)[[3]]

## make probability objects into dataframes and combine
alleleprobs_mnt = list()
for(i in 1:dim(alleleprobs_minimuga)[1]){
  tmp = t(alleleprobs_minimuga[i,,])
  tmp = data.frame(marker = rownames(tmp), tmp) %>% left_join(all_genes, by="marker")
  alleleprobs_mnt[[dimnames(alleleprobs_minimuga)[[1]][i]]] = tmp
}

alleleprobs_cegs_df = do.call("rbind", alleleprobs_cegs)
alleleprobs_cegs_df$Pup.ID = unlist(strsplit(rownames(alleleprobs_cegs_df),"[.]"))[c(T,F)]
alleleprobs_cegs_df$CC.1 = xce_pups$CC1[match(alleleprobs_cegs_df$Pup.ID, xce_pups$Pup.ID)]
alleleprobs_cegs_df$CC.2 = xce_pups$CC2[match(alleleprobs_cegs_df$Pup.ID, xce_pups$Pup.ID)]
alleleprobs_cegs_df$RRIX = apply(alleleprobs_cegs_df, 1, function(x) paste0(x["CC.1"],"/", x["CC.2"]))

alleleprobs_mnt_df = do.call("rbind", alleleprobs_mnt)
alleleprobs_mnt_df$Pup.ID = unlist(strsplit(rownames(alleleprobs_mnt_df),"[.]"))[c(T,F)]
alleleprobs_mnt_df = alleleprobs_mnt_df %>% 
  left_join(xce_pups[,-which(colnames(xce_pups) == "D")], by="Pup.ID") %>%
  rename("CC.1"="CC1","CC.2"="CC2")


alleleprobs_comb = rbind(alleleprobs_cegs_df[, intersect(colnames(alleleprobs_cegs_df), colnames(alleleprobs_mnt_df))], 
                         alleleprobs_mnt_df[, intersect(colnames(alleleprobs_cegs_df), colnames(alleleprobs_mnt_df))])

## sort through dataframe and remove missing data
rem = which(is.nan(alleleprobs_comb$A))
alleleprobs_comb[-rem, ] %>% arrange(RRIX, Pup.ID, pos) %>%
  distinct() -> alleleprobs_comb
alleleprobs_comb %>%
  filter(!Pup.ID %in% problemPups) %>%
  arrange(RRIX, pos, Pup.ID) %>% distinct() -> alleleprobs_comb
alleleprobs_comb$RIX = paste0(alleleprobs_comb$CC.1,"/",alleleprobs_comb$CC.2)

## round out probs
alleleprobs_comb[,LETTERS[1:8]] = round(alleleprobs_comb[,LETTERS[1:8]], 3)
alleleprobs_comb$RIX = factor(alleleprobs_comb$RIX, levels=xce_guide$CCs, ordered = T)
alleleprobs_comb = alleleprobs_comb %>% arrange(RIX, pos)
keep_round = alleleprobs_comb
keep_round[,LETTERS[1:8]] = round(keep_round[,LETTERS[1:8]], 1)
keep_round = keep_round %>% select(-"Pup.ID") %>% 
  arrange(RIX, pos) %>% distinct() %>%
  mutate(RIX_pos = paste0(pos, "_", RIX))

## only keep positions with agreement among all samples in a CC-RIX
keep_pos = keep_round$RIX_pos
alleleprobs_comb_nod = alleleprobs_comb %>% select(-"Pup.ID") %>% 
  arrange(RIX, pos) %>% distinct() %>%
  mutate(RIX_pos = paste0(pos, "_", RIX)) %>%
  filter(RIX_pos %in% keep_pos) %>%
  group_by(RRIX, RIX, RIX_pos, CC.1, CC.2, marker, pos) %>%
  dplyr::summarize(A = mean(A),B = mean(B),C = mean(C),D = mean(D),
                   E = mean(E),F = mean(F),G = mean(G),H = mean(H))

## identify homozygous regions 
homog = apply(alleleprobs_comb_nod, 1, function(x) ifelse(any(x[LETTERS[1:8]] > 0.8), T, F))
plot_hets=data.frame(alleleprobs_comb_nod)
plot_hets$homog = homog

## determine top two founders
het_founds = data.frame(do.call("rbind",apply(plot_hets, 1, function(x){
  sort(names(sort(x[LETTERS[1:8]][which(x[LETTERS[1:8]] > 0.05)],decreasing=T)[1:2]))
})))

colnames(het_founds) = c("founder1", "founder2")
plot_hets = cbind(plot_hets, het_founds)
plot_hets$founders = apply(plot_hets, 1, function(x) 
  paste(sort(c(paste(x["founder1"]), paste(x["founder2"]))), collapse = ""))
plot_hets %>% 
  arrange(RIX, pos) %>%
  group_by(founders) -> plot_hets

## wide-to-long format
plot_hets_points = plot_hets %>% select("RIX","pos","founder1","founder2","homog") %>%  
  gather(which_f, founder, founder1:founder2, factor_key = T) %>% 
  rename("Founder"="founder") %>%
  arrange(RIX, pos, which_f) 

## order in plot
plot_hets_points$NOD = 3
plot_hets_points$NOD[which(plot_hets_points$RIX %in% nod_rix)] = 1
plot_hets_points$NOD[which(plot_hets_points$RIX %in% setdiff(nzo_rix, nod_rix))] = 2
plot_hets_points %>% arrange(NOD) -> plot_hets_points
plot_hets_points$RIX = factor(plot_hets_points$RIX, levels=unique(plot_hets_points$RIX), ordered=T)

## convert x-axis to numeric such that each CC-RIX gets two points
plot_hets_points$Pup = as.numeric(factor(plot_hets_points$RIX)) +  
  as.numeric(paste0("0.",gsub("founder","", as.character(plot_hets_points$which_f))))*2
plot_hets_points$Founder = factor(plot_hets_points$Founder, levels = LETTERS[1:8])
npup = length(unique(plot_hets_points$RIX))


## lower transparency of homozygous regions
plot_hets_points$alpha = 1
plot_hets_points$alpha[which(plot_hets_points$homog)] = 0.0001

## separate nod, nzo groups
nzo_x = range(plot_hets_points$Pup[which(plot_hets_points$NOD == 2)])
area_bounds = c(nzo_x[1]-0.2, nzo_x[2]+0.2)


## determine consistent region for Xce
low = (plot_hets_points %>% filter(RIX == "CC035/CC062", Founder == "D", pos < 120, homog) %>% arrange(-pos))$pos[1]
hi1 = (plot_hets_points %>% filter(RIX == "CC041/CC051", Founder == "D", pos < 120) %>% arrange(-pos))$pos[1]
hi2 = (plot_hets_points %>% filter(RIX == "CC041/CC012", founders == "DD", pos < 120) %>% arrange(-pos))$pos[1]
ribbon = data.frame(x=c(-Inf, Inf),
                    lo=low,  hi=min(hi1, hi2))

## classic CC colors
cols = c("yellow", "gray", "pink", "blue", "deepskyblue", "forestgreen", "red", "purple")

p = ggplot() + 
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=area_bounds[1], xmax=area_bounds[2]),
            fill="gray", alpha=0.8) + 
  geom_ribbon(data=ribbon, aes(ymin=lo, ymax=hi, x=x), fill="#a3752e", alpha=0.2) + 
  geom_point(data=plot_hets_points, aes(x=Pup, y=pos, col=Founder, alpha=alpha), size=0.5,shape=15) + 
  scale_colour_manual(values=cols) + 
  scale_alpha_continuous(NULL, NULL) + 
  scale_x_continuous(name ="CC-RIX", breaks=c(1:npup)+0.3, #96
                     labels = as.character(unique(plot_hets_points$RIX)),#Pup.ID
                     limits=c(1,npup+0.5)) +
  scale_y_continuous(name ="Position", breaks=seq(0,170,5)) +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45,hjust=1, size=8),
        axis.text.y=element_text(size=8))

p

##################################
##          cnv plots           ##
##################################
kmer_bounds = read.table("FileS19_positions_duplications.txt",header=T)

## meta data for samples providing DNA-seq kmer counts
meta = read.csv("FileS20_dna_45kmer_meta_full.csv")
meta_use = meta %>% filter(use)

## all unfiltered kmer counts for CC, inbreds, and sister strains
kmersAll = read.csv("FileS21_dna_45kmer_counts_allCCs.csv")

## normalize counts by mode or mean
counts = kmersAll[,grep("_",colnames(kmersAll))]
colnames(counts) = gsub("^X", "", colnames(counts))
counts = counts[,which(colnames(counts) %in% meta_use$Sample)]
normed = sapply(1:ncol(counts), function(i)
  counts[,i] / meta_use$mode[which(meta_use$Sample == colnames(counts)[i])], simplify=F)
normed = do.call("cbind", normed)
colnames(normed) = colnames(counts)

## determine number of copies in B6 strains
B6_CC = intersect(which(meta_use$xce_hap == "B"), grep("CC", meta_use$CC))
tmp = setdiff(which(meta_use$xce_hap == "B"), B6_CC)
B6_sis = intersect(which(is.na(meta_use$sex)), which(meta_use$xce_hap == "B"))
B6_inbred = setdiff(tmp, B6_sis)
B6_use = sort(c(B6_CC, B6_inbred))
B6_inbred_ave = apply(normed[,B6_inbred], 1, mean)
B6_CC_ave = apply(normed[,B6_CC], 1, mean)
delta = B6_inbred_ave - B6_CC_ave
B6_all_ave = apply(normed[,B6_use], 1, mean)  #which(meta_use$xce_hap == "B")
kmersAll$B6_copies = round(B6_all_ave)

## remove overrepresented kmers
remove_Pos = sort(unique(c(which(abs(delta) > 2),
                           which(kmersAll$Sequence %in% names(which(sort(table(kmersAll$Sequence)) > 1))) )))
if(length(remove_Pos) > 0){
  counts = counts[-remove_Pos, which(colnames(counts) %in% meta_use$Sample)]
  normed = normed[-remove_Pos,]
  kmer_deets = kmer_deets[-remove_Pos,]
  kmersAll = kmersAll[-remove_Pos, 
                      c(which(colnames(kmersAll) %in% c("Position","Sequence","CN","B6_copies")), 
                        which(colnames(kmersAll) %in% meta_use$Sample))]
}

## remove sister strains
rem = grep("sis", meta_use$grp)
meta_use = meta_use[-rem,]
normed = normed[,-rem]


## calculate averages and deltas for each founder
aves = lapply(sort(unique(meta_use$xce_hap)), function(i){    # use grp to distinguish CC and inbred
  if(i == "B"){
    return(data.frame("B_inbred" = apply(normed[,which(meta_use$grp == "B_inbred")],1,mean),
                      "B_CC" = apply(normed[,which(meta_use$grp == "B_CC")],1,mean)))
  } else if(length(which(meta_use$xce_hap == i)) > 1){
    cc = paste(i)
    return(data.frame(cc = apply(normed[,which(meta_use$xce_hap == i)],1,mean)))
  } else {
    return(data.frame(cc = normed[,which(meta_use$xce_hap == i)]))
  }
})
aves = data.frame(do.call("cbind", aves))
colnames(aves) = paste0(c("A", "B_inbred","B_CC", sort(unique(meta_use$xce_hap))[3:8]), "_ave")
aves = aves[,c("B_inbred_ave",colnames(aves)[-which(colnames(aves) %in% c("B_inbred_ave"))])]

ave_delt = data.frame(lapply(colnames(aves)[-1], function(i)             
  apply(aves,1,function(y)  y[i]-y["B_inbred_ave"])))
ave_delt = do.call("cbind", ave_delt)
colnames(ave_delt) = paste0(colnames(aves)[-1],"_delt")

## duplicate info in kmer_deets, in case
## determine boundaries of SD's
normed45_delt = cbind(kmer_deets, ave_delt)
in_sd = apply(kmer_bounds, 1, function(x) 
  which(normed45_delt$Position > as.numeric(x["start"]) & normed45_delt$Position < as.numeric(x["end"])))
normed45_delt$in_sd = 0
for(i in 1:length(in_sd)){
  normed45_delt$in_sd[in_sd[[i]]] = i
}

## make data wide from long
normed45_delt = gather(normed45_delt, key="founder", value="diff", 
                       colnames(ave_delt)[1]:colnames(ave_delt)[ncol(ave_delt)]) %>% 
  mutate(type = "delta")
Ddup = data.frame(start=102.802501, end=102.839301) ## inclusive

xmin=102.7
xmax=102.9
  
make_cnv_plots = function(founder_label = NULL, kmer_len = 45){
  Ddup = data.frame(start=102.802501, end=102.839301) ## inclusive
  f = founder_label
  same_range = above_range = below_range = NULL
  centers = c(-1,0,1)
  if(length(grep("F|G", f)) > 0) centers = c(-2,-1,0,1)
  if(length(grep("E", f)) > 0) centers = c(-1,0,1,2,3)
  
  k=length(centers)
  plot = normed45_delt %>% filter(founder == f, 
                                  Position > xmin, Position < xmax)#, #overlap=="none")#,
  plot$in_dup = F
  in_dup = which(plot$Position > Ddup$start & plot$Position < Ddup$end)
  plot$in_dup[in_dup] = T
  
  km_stats = kmeans(data.frame(plot$diff), centers = centers)
  centers = as.vector(km_stats$centers)
  above_0 = which.max(centers)
  below_0 = which.min(centers)
  same = setdiff(1:k, c(above_0, below_0))
  other=NA
  if(length(same) > 1) 
    same = which.min(abs((0 - centers)))
  other = setdiff(setdiff(1:k, c(above_0, below_0)), same)
  
  centers = centers[c(same, below_0, other, above_0)]
  
  ## remove centers from duplicated kmers
  plot$plot_diff = plot$diff
  if(length(grep("^[A|C|D]", f)) > 0){
    plot$plot_diff[grep("_with_5",plot$overlap)] = plot$diff[grep("_with_5",plot$overlap)] - centers[length(centers)]
    plot$plot_diff[grep("6_with_4",plot$overlap)] = plot$diff[grep("6_with_4",plot$overlap)] - centers[length(centers)]
  } else if(f == "H_ave_delt"){
    plot$plot_diff[grep("_with_1",plot$overlap)] = plot$diff[grep("_with_1",plot$overlap)] - centers[length(centers)]
  } else if(f == "E_ave_delt"){
    plot$plot_diff[grep("_with_5",plot$overlap)] = plot$diff[grep("_with_5",plot$overlap)] - centers[length(centers)]
    plot$plot_diff[grep("6_with_4",plot$overlap)] = plot$diff[grep("6_with_4",plot$overlap)] - centers[length(centers)-1]
  }
 
  ## colors for k-mean clusters
  ## other comments below for coloring in clusters too
  colors = brewer.pal(2*(k+1), "Paired")  #palette()
  colors = c(c("gray60","gray40"), colors)  
  dot_colors = colors[(1:(k+1))*2-1]  
  lin_colors = colors[(1:(k+1))*2]  
  
  plot$cpn = factor(km_stats$cluster, levels=c(same, below_0, other, above_0, (length(centers)+1)), ordered=T)
  same_range = range(plot$diff[which(plot$cpn == same)])
  below_range = range(plot$diff[which(plot$cpn == below_0)])
  above_range = range(plot$diff[which(plot$cpn == above_0)])
  if(f %in% c("B_CC_ave_delt", "H_ave_delt")){
    same_range = above_range = below_range = c(0,0)
  }
  plot$col = factor(plot$CN, levels=c("0","1","2","3","4+"), ordered=T)   #B6 copies
  plot$col[is.na(plot$col)] = "4+"
  nsnps = length(which(plot$col == "1" & plot$diff < below_range[2]))
  nsnps_inDup = length(intersect(which(plot$col == "1" & plot$diff < below_range[2]),
                                 which(plot$in_dup==T)))
  perc_snp = nsnps/nrow(plot)
  perc_snp_inDup = nsnps_inDup / length(which(plot$in_dup==T))
  plot$alpha = 0.1
  plot$alpha[which(plot$sd != 0)] = plot$alpha[which(plot$sd != 0)]+0.99
  ymax = round(max(centers)+3)
  ymin = round(min(centers)-1.5)
  kmer_bounds$y = rep(c(ymax-0.25, ymax-0.35),ceiling(nrow(kmer_bounds)/2))[1:nrow(kmer_bounds)]

  Ddup$y = ymax-0.5
  if(f == "E_ave_delt") Ddup$end = kmer_bounds$end[6]
  if(f == "H_ave_delt"){
    Ddup$start = 102.7302 
    Ddup$end = kmer_bounds$end[1]
  }
  
  title = ifelse(length(grep("sis",f))>0, unique(meta_use$CC[which(meta_use$grp == gsub("_ave_delt","",f))]),
                 unlist(strsplit(f,"_"))[c(T,F)])
  title_long = data.frame(let = LETTERS[1:8], CC = c("A/J","C57BL/6J","129S1/SvImJ","NOD","NZO","CAST/EiJ","PWK/PhJ","WSB/EiJ"))
  title = title_long$CC[which(title_long$let == title)]
  sub_nsnps = ifelse(length(grep("F|G", f)) > 0, "", paste0("\n# SNPs: ", nsnps,", ", round(perc_snp,3), "%"))
  sub_nsnps_inDup = ifelse(length(grep("F|G", f)) > 0, "", paste0("\n# SNPs in dup: ", nsnps_inDup,", ", round(perc_snp_inDup,3),"%"))
  break_pos = seq(xmin,xmax, by = 0.01)
  break_lab = break_pos
  break_lab[grep("[.]", paste(break_pos/0.05))] = ""
  p = ggplot(data=plot) +
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=same_range[1], ymax=same_range[2]),fill="gray",col=NA,alpha=0.8) + 
    geom_point(aes(x=Position, y=plot_diff, col=col, alpha=alpha)) +
    geom_segment(data=kmer_bounds, 
                 aes(x=start, xend=end, y=y, yend=y),
                 lineend = "round", linejoin = "bevel",
                 arrow = arrow(length = unit(0.3, "cm"))) + 
    geom_segment(data=Ddup, col="mediumvioletred", size=2,
                 aes(x=start, xend=end, y=y, yend=y)) + 
    scale_y_continuous(limits = c(ymin, ymax), 
                       breaks=seq(ymin, ymax, 1)) + 
    scale_x_continuous(breaks = break_pos, labels = break_lab) + 
    theme_classic() + 
    scale_color_viridis(discrete = TRUE, name ="# copies in ref") + 
    scale_alpha_continuous(NULL,NULL,NULL) +
    geom_vline(xintercept = c(Ddup$start, Ddup$end),col="mediumvioletred", linetype = 2) +
    geom_vline(alpha=0.3, col="mediumvioletred", linetype = "dotted",
               xintercept = c(kmer_bounds$start, kmer_bounds$end[-1])) +
    theme(plot.title = element_text(size=10), 
          plot.subtitle=element_text(size=8, face="italic", color="gray20")) + #,
    labs(y = paste("Inbred C57BL/6J counts -", title)) #,            

  return(list(plot = p, data = plot))
}
  
plts = list()
pdf("all_CC_CNV_plots.pdf",width=7, height=4)
for(f in unique(normed45_delt$founder)){
  plts[[f]] = make_cnv_plots(f)
  print(plts[[f]]$plot)
}
dev.off()

##############################################
##        sample-specific proportions       ##
##############################################

#df_comb
#summary_comb

full_df = list()
for(i in 1:length(unique(df_comb$CC_LAB))){
  cc <- unique(df_comb$CC_LAB)[i]
  #rix <- as.numeric(unique(gsub("a|b", "", unique(df$RIX))))[i]
  chr <- unique(df_comb$SEQ.CHROMOSOME)
  rixdat <- df_comb %>% filter(CC_LAB == cc) %>%
    mutate(SEQ.POSITION = as.numeric(SEQ.POSITION)) %>%
    arrange(DIR, PUP.ID, SEQ.POSITION) 
  rixdat$PUP.ID <- factor(rixdat$PUP.ID, levels=unique(rixdat$PUP.ID))
  use_sum_ind <- unique(match(rixdat$RRIX ,gsub("RIX_", "", toupper(names(summary_comb)))))
  use_sum <- summary_comb[[use_sum_ind]]$summary$mu_p
  use_sum$Pup.ID = rownames(use_sum)
  
  rixdat$kRat_g <- as.numeric(summary_comb[[use_sum_ind]]$summary$mu_g$Mean[match(rixdat$PUP_GENE,
                                                                           rownames(summary_comb[[use_sum_ind]]$summary$mu_g))])
  rixdat$kLB_g <- as.numeric(summary_comb[[use_sum_ind]]$summary$mu_g$lower[match(paste(rixdat$PUP.ID, rixdat$SEQ.GENE, sep="_"),
                                                                           rownames(summary_comb[[use_sum_ind]]$summary$mu_g))])
  rixdat$kUB_g <- as.numeric(summary_comb[[use_sum_ind]]$summary$mu_g$upper[match(paste(rixdat$PUP.ID, rixdat$SEQ.GENE, sep="_"),
                                                                           rownames(summary_comb[[use_sum_ind]]$summary$mu_g))])
  
  colnames(use_sum) = tolower(colnames(use_sum))
  rixdat$pup_lower <- as.numeric(use_sum$lower[match(rixdat$PUP.ID, use_sum$pup.id)])
  rixdat$pup_upper <- as.numeric(use_sum$upper[match(rixdat$PUP.ID, use_sum$pup.id)])
  rixdat$pup_mean <- as.numeric(use_sum$mean[match(rixdat$PUP.ID, use_sum$pup.id)])
  full_df[[cc]] = rixdat
}

full_df = do.call("rbind",full_df)
full_df$XIST = ifelse(full_df$XIST, "Yes","No")

mvec = c(2422,2138,1656)
### 2422 in CC041/CC051, D vs H
### 2138 in CC006/CC026, D vs D
### 1656 in CC023/CC047, D vs E

rixdat = full_df %>% filter(PUP.ID %in% mvec)

for(c in unique(full_df$CC_LAB)){
  rixdat = full_df %>% filter(CC_LAB == c)
  
  p = ggplot(data=rixdat, aes(x=SEQ.POSITION, alpha=LOGSUM, color=XIST, shape=factor(DIR))) + 
    geom_point(aes(y = kRat_g)) + #
    geom_segment(aes(x=SEQ.POSITION, xend=SEQ.POSITION, y = kLB_g, yend=kUB_g), size=0.5) + 
    ylab("Ratio") + xlab("Position") +
    ggtitle(c) + 
    theme_bw() + 
    scale_shape_discrete(name=NULL, NULL) + 
    scale_alpha_continuous(name=NULL, NULL) + 
    scale_x_continuous(labels=function(n){paste(n/1000000, "M")}) +
    geom_hline(yintercept=0.5, col=colorpal$blue, linetype= "dotted", size=1) + 
    scale_colour_manual(values=c(colorpal$slategray, colorpal$magenta)) +
    geom_ribbon(aes(x=SEQ.POSITION, ymin=pup_lower, ymax=pup_upper), inherit.aes = T, 
                alpha=0.2, col="gray") + 
    geom_hline(aes(yintercept = pup_mean), col="#292d57") + 
    facet_wrap(. ~ PUP.ID, scales="free_x") + ylim(c(0,1))
  
  print(p)
}

##############################################
##         XCI proportions from lit         ##
##############################################

library(reshape2)
cor = data.frame(a = c(0.5,0.33,0.4,0.24), 
           e = c(0.73,0.5,0.56,0.34),
           b = c(0.61,0.57,0.5,0.35),
           c = c(0.79,0.79,0.67, 0.5))
rownames(cor) = colnames(cor)
cor_plot = melt(cor)
cor_plot$Paternal = factor(rep_len(rownames(cor), nrow(cor_plot)), levels=colnames(cor))
cor_plot = cor_plot %>% rename("Maternal" = "variable")  
cor_plot$color = "white"
cor_plot$color[which(cor_plot$value > 0.6)] = "black"
ggplot(data = cor_plot, aes(Paternal, Maternal, fill = value))+
  geom_tile() + 
  scale_fill_viridis(name="Maternal Xa\nproportion") + 
  geom_text(aes(label = value, color = color), size = 4) +
  scale_color_manual(values = c("black","white"), NULL, NULL) + 
  theme_bw()+ 
  coord_fixed()



###############################
###         alpha0          ###
###############################
plot_gams = read.csv("../tables_figures/TableS2_alpha0.csv")
plot_gams = rbind(plot_gams, plot_gams[which(plot_gams$rix == "grand_sum"),])

for(i in 1:nrow(plot_gams)){
  if(i == 1){
    lwduse = ifelse(plot_gams$rix[i] == "sum", 3,
                    ifelse(plot_gams$rix[i] == "grand_sum", 2, 0.5))
    coluse = ifelse(plot_gams$sp[i] == 1, colorpal$magenta,
                    ifelse(plot_gams$rix[i] == "grand_sum", "white", colorpal$green))
    ltyuse = ifelse(plot_gams$rix[i] == "grand_sum", 0, 1)
    plot(x = seq(1,1000), y = dgamma(seq(1,1000),
                                     shape = plot_gams$shape[i], rate=plot_gams$rate[i]), 
         type="l", col=coluse, log="x", lwd=lwduse, lty=ltyuse, 
         xlab = expression(paste("Number of brain precursor cells in the E5.5 epiblast (" *alpha[0]*")")), 
         ylab = "Posterior density")   #, main = expression(paste("Combined"~alpha[0]~"estimates")))
  } else {
    lwduse = ifelse(plot_gams$rix[i] == "sum", 3,
                    ifelse(plot_gams$rix[i] == "grand_sum", 2, 0.5))
    coluse = ifelse(plot_gams$sp[i] == 1, colorpal$magenta,
                    ifelse(plot_gams$rix[i] == "grand_sum", colorpal$skyblue, colorpal$green))
    ltyuse = ifelse(plot_gams$rix[i] == "grand_sum", 2, 1)
    lines(x = seq(1,1000), y = dgamma(seq(1,1000),
                                      shape = plot_gams$shape[i], rate=plot_gams$rate[i]), 
          lwd=lwduse, type="l", lty=ltyuse, col=coluse)
  }
  
}

legend("topright", inset=0.02, legend=c("SP1", "SP2","Combined"),
       col=as.character(colorpal[c(6,4,3)]),lty=c(1,1,2), lwd=2)
