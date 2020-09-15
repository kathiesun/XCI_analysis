######## Paternity test
#setwd("/nas/depts/006/valdar-lab/users/sunk")
#setwd("C:/DB Mount/Dropbox\ (ValdarLab)/outputs/")
#### matnut_outputs

library(tidyverse)
#library(MASS)


xce <- as.data.frame(rbind(rix1 = c("CC011", "CC001", "e", "a"),
                    rix2 = c("CC041", "CC051", "b", "b"),
                    rix3 = c("CC004", "CC017", "b", "e"),
                    rix4 = c("CC023", "CC047", "b", "b"),
                    rix6 = c("CC026", "CC006", "b", "b"),
                    rix7 = c("CC003", "CC014", "e", "a"),
                    rix8 = c("CC035", "CC062", "a", "a"),
                    rix9 = c("CC032", "CC042", "c", "a"),
                    rix10 = c("CC005", "CC040", "b", "b")))
colnames(xce) = c("num_cc", "dem_cc", "num_let", "dem_let")


xce_df = data.frame(founder = LETTERS[1:8],
                    xce_al = c("a","b","a","n","z","c","e","b"))

curateKmers_founders <- function(exons, founders=c("AJ", "C57BL6J", "129", "NOD","NZO","CAST","PWK","WSB"), 
                                 norm_file = NULL, up_bound = 150, lo_bound = 3){
  chrAll <- unique(exons$Chromosome)  
  fin <- list()
  for(j in 1:length(chrAll)) {
    chr <- paste(chrAll[j])
    ### read in data
    exonsUse <- exons[which(exons$Chromosome %in% chr),] %>% 
      mutate(ProbeSeq = as.character(ProbeSeq)) 
    
    ## remove duplicates across genes
    # keep none
    exonsUse %>%
      dplyr::select(one_of("rsId","Position","ProbeSeq")) %>% 
      group_by(ProbeSeq) %>%
      distinct() %>%
      mutate(n = n()) %>%
      filter(n>1) %>%
      arrange(desc(n), ProbeSeq) %>%
      dplyr::select(ProbeSeq) -> remove
    
    ## remove duplicates within gene
    # keep one
    exonsUse <- filter(exonsUse, !ProbeSeq %in% t(remove)) %>% 
      group_by(ProbeSeq) %>%
      filter(row_number() == 1)
    
    ### Get rid of variants within 25 bp to exon start/stop
    avoid <- cbind(exonsUse$start + 13, exonsUse$end - 13)	
    tooClose <- sapply(1:length(exonsUse$Position), function(x) 
      intersect(which(exonsUse$Position[x] < avoid[x,1]), which(exonsUse$Position[x] > avoid[x,2])), simplify=F) 
    del <- which(unlist(lapply(tooClose, length)) > 0)
    if(length(del) > 0) exonsUse <- exonsUse[-del,]
    
    ### Separate out ref and alt alleles in founder counts
    founderMiceNum <- unlist(lapply(founders, function(x) grep(x, colnames(exonsUse))[1]))
    founderMice <- colnames(exonsUse)[founderMiceNum]
    founderCounts <- exonsUse[,founderMiceNum]
    founderCounts_ref <- apply(founderCounts,2,function(x) as.numeric(paste(unlist(strsplit(x, "/"))[c(T,F)])))
    founderCounts_alt <- apply(founderCounts,2,function(x) as.numeric(paste(unlist(strsplit(x, "/"))[c(F,T)])))
    
    if(!is.null(norm_file)){
      founderCounts_ref <- sapply(1:ncol(founderCounts_ref), function(x){
        n = which(colnames(founderCounts_ref)[x] == norm_file$founder)
        founderCounts_ref[,x] * norm_file$mult_fact[n]
      })
      founderCounts_alt <- sapply(1:ncol(founderCounts_alt), function(x){
        n = which(colnames(founderCounts_alt)[x] == norm_file$founder)
        founderCounts_alt[,x] * norm_file$mult_fact[n]
      })
      colnames(founderCounts_ref) = colnames(founderCounts_alt) = founderMice
    }
    
    snp <- exonsUse[, c("Gene","Chromosome","rsId","Position","ProbeSeq","TranscriptID")]
    
    seq <- as.data.frame(do.call("rbind",strsplit(as.character(exonsUse$ProbeSeq), "\\[|\\]|/")))
    colnames(seq) <- c("end5","ref","alt","end3")
    seq$refseq <- paste0(seq$end5, seq$ref, seq$end3)
    seq$altseq <- paste0(seq$end5, seq$alt, seq$end3)
    
    comb <- data.frame(snp, seq)
    
    sdp <- do.call("rbind", strsplit(exonsUse$FounderSDP,""))
    isAlt <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 1), which(founderCounts_ref[x,] <= 3 )))
    isRef <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 0), which(founderCounts_ref[x,] > 3)))
    isAlt2 <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 1), which(founderCounts_alt[x,] > 3)))
    isRef2 <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 0), which(founderCounts_alt[x,] <= 3)))
    altTrue <- sapply(1:nrow(sdp), function(x) ifelse(identical(isAlt[[x]], isAlt2[[x]]), 1, 0), simplify = T)
    refTrue <- sapply(1:nrow(sdp), function(x) ifelse(identical(isRef[[x]], isRef2[[x]]), 1, 0), simplify = T)
    bothTrue <- intersect(which(altTrue == 1), which(refTrue == 1))
    
    ### trim out variants that disagree with founder sdp
    var_trim <- list()
    var_trim[["seq"]] <- comb[bothTrue,] 
    var_trim[["sdp"]] <- sdp[bothTrue,]
    var_trim[["founder_ref"]] <- founderCounts_ref[bothTrue,]
    var_trim[["founder_alt"]] <- founderCounts_alt[bothTrue,]
    
    hetsOrZeros <- sapply(1:nrow(var_trim[["founder_ref"]]), function(x){
      sapply(1:ncol(var_trim[["founder_ref"]]), function(y){
        ifelse(var_trim[["founder_ref"]][x,y] > up_bound | 
                 var_trim[["founder_alt"]][x,y] > up_bound, "M", 
               ifelse(var_trim[["founder_ref"]][x,y] > lo_bound &
                        var_trim[["founder_alt"]][x,y] > lo_bound, "H", 
                      ifelse(var_trim[["founder_ref"]][x,y] <= lo_bound & 
                               var_trim[["founder_alt"]][x,y] <= lo_bound, "Z", 
                             ifelse(var_trim[["founder_ref"]][x,y] > var_trim[["founder_alt"]][x,y], "R", "A"))))
      })
    })
    
    hetsOrZeros <- t(hetsOrZeros)	
    posCounts <- var_trim[["founder_ref"]] + var_trim[["founder_alt"]]
    
    keep <- which(apply(posCounts, 1, function(x) ifelse(sum(!is.na(x)) == length(x), T, F)))
    posCounts <- posCounts[keep,]
    kmer_mu <- apply(posCounts, 1, function(x) mean(x[which(!is.na(x))]))	
    kmer_med <- apply(posCounts, 1, function(x) median(x[which(!is.na(x))]))
    kmer_length <- apply(posCounts, 1, function(x) length(x[which(!is.na(x))]))	
    kmer_var <- apply(posCounts, 1, function(x) var(x[which(!is.na(x))]))
    kmer_range <- apply(posCounts, 1, function(x) range(x[which(!is.na(x))]))
    kmer_quant <- apply(posCounts, 1, function(x) quantile(x[which(!is.na(x))]))
    lo_bounds_p <- unlist(lapply(kmer_mu, function(x) qpois(0.0025, x)))
    hi_bounds_p <- unlist(lapply(kmer_mu, function(x) qpois(0.9975, x)))
    lo_bounds_n <- unlist(sapply(1:length(kmer_mu), function(x) qnorm(0.0005, kmer_mu[x], sqrt(kmer_var[x]))))
    hi_bounds_n <- unlist(sapply(1:length(kmer_mu), function(x) qnorm(0.9995, kmer_mu[x], sqrt(kmer_var[x]))))
    
    hi_bounds = hi_bounds_p
    lo_bounds = rep(lo_bound*2, length(hi_bounds))    
    
    kmer_stats <- data.frame(mu = kmer_mu, med = kmer_med, var = kmer_var, 
                             sd = apply(posCounts,1, function(x) sd(x, na.rm=T)),
                             min = apply(posCounts,1, function(x) min(x[which(!is.na(x))])),
                             max = apply(posCounts,1, function(x) max(x[which(!is.na(x))])), 
                             l_bound_p = lo_bounds_p, h_bound_p = hi_bounds_p,
                             l_bound_n = lo_bounds_n, h_bound_n = hi_bounds_n)
    
    keep <- seq(1,nrow(posCounts))
    
    remove <- sapply(keep, function(x)  
      union(which(posCounts[x,] > hi_bounds[x]), which(posCounts[x,] < lo_bounds[x])), simplify = T)
    
    kmer_stats$num_removed <- 0
    kmer_stats$num_removed[keep] <- unlist(lapply(remove, length))
    ###
    remove_rows = which(kmer_stats$num_removed > 0)
    
    var_trim$hetsOrZeros <- hetsOrZeros
    var_trim$posCounts <- posCounts
    var_trim = lapply(var_trim, function(x) x[-remove_rows,])
    fin[[j]] <- var_trim
  }
  return(fin)
}


sort_counts <- function(ref, alt){
  colnames(ref) = colnames(alt) <- c("file", "k.mer", "dir", "counts")
  
  ref$Pup <- unlist(strsplit(as.character(ref$file), "/"))[c(F,F,F,F,T)]
  rem <- grep("old", ref$Pup)
  if(length(rem) > 0) ref <- ref[-rem,]
  ref$Pup.ID <- as.numeric(paste(gsub("[^0-9]", "", ref$Pup)))
  
  alt$Pup <- unlist(strsplit(as.character(alt$file), "/"))[c(F,F,F,F,T)]
  rem <- grep("old", alt$Pup)
  if(length(rem) > 0) alt <- alt[-rem,]
  alt$Pup.ID <- as.numeric(paste(gsub("[^0-9]", "", alt$Pup)))
  
  tmp <- do.call("rbind", strsplit(as.character(alt$k.mer), ""))
  tmp <- tmp[,-((ncol(tmp)/2)+1)]
  consensus <- apply(tmp, 1, function(x) paste(x, collapse=""))
  alt$consensus = consensus
  alt$consensus_dir = paste(alt$consensus, alt$dir, sep="_")
  
  tmp <- do.call("rbind", strsplit(as.character(ref$k.mer), ""))
  tmp <- tmp[,-((ncol(tmp)/2)+1)]
  consensus <- apply(tmp, 1, function(x) paste(x, collapse=""))
  ref$consensus = consensus
  ref$consensus_dir = paste(ref$consensus, ref$dir, sep="_")
  
  ref$allele = 0
  alt$allele = 1
  
  if(identical(ref$consensus_dir, alt$consensus_dir) && identical(ref$file, alt$file)){
    ref <- ref %>% dplyr::select(-one_of("file","consensus_dir","allele")) %>%
      rename("counts_ref"="counts","k.mer_ref"="k.mer")
    alt <- alt %>% dplyr::select(one_of("counts","k.mer")) %>%
      rename("counts_alt"="counts","k.mer_alt"="k.mer")
  }
  
  
  all <- cbind(ref, alt) %>% distinct()
  
  all %>% mutate(counts = counts_ref+counts_alt) %>%
    group_by(consensus, dir) %>% 
    summarise(sum = sum(counts)) %>% 
    arrange(desc(sum)) %>% 
    dplyr::slice(1) %>%
    mutate(filt = paste(consensus, dir, sep="_")) -> keepDir
  
  all %>% 
    mutate(filt = paste(consensus, dir, sep="_")) %>%
    filter(filt %in% keepDir$filt) %>%
    dplyr::select(-one_of("filt")) %>%
    mutate(counts = counts_ref+counts_alt) -> counts
  return(counts)
}



curateKmers <- function(exons, myCC=NULL, founders=c("AJ", "C57BL6J", "129", "NOD","NZO","CAST","PWK","WSB"), 
                        up_bound = 150, lo_bound = 3, 
                        founderHaps=NULL, check_CC_sdp=T, impute=F) {
	chrAll <- unique(exons$Chromosome)  
	fin <- list()
	for(j in 1:length(chrAll)) {
		chr <- as.factor(chrAll[j])
		### read in data
		exonsUse <- exons[which(exons$Chromosome %in% chr),] %>% 
		  mutate(ProbeSeq = as.character(ProbeSeq)) 
		founderHapsUse <- lapply(founderHaps, function(x) na.omit(x[which(gsub("chr","",x$chr) == chr),]))

		## remove duplicates across genes
		# keep none
		exonsUse %>%
		  dplyr::select(one_of("rsId","Position","ProbeSeq")) %>% 
		  group_by(ProbeSeq) %>%
		  distinct() %>%
		  mutate(n = n()) %>%
		  filter(n>1) %>%
		  arrange(desc(n), ProbeSeq) %>%
		  dplyr::select(ProbeSeq) -> remove

		## remove duplicates within gene
		# keep one
		exonsUse <- filter(exonsUse, !ProbeSeq %in% t(remove)) %>% 
		  group_by(ProbeSeq) %>%
		  filter(row_number() == 1)
		
		avoid <- cbind(exonsUse$start + 13, exonsUse$end - 13)	
	  tooClose <- sapply(1:length(exonsUse$Position), function(x) 
			intersect(which(exonsUse$Position[x] < avoid[x,1]), which(exonsUse$Position[x] > avoid[x,2])), simplify=F) 
		del <- which(unlist(lapply(tooClose, length)) > 0)
		if(length(del) > 0) exonsUse <- exonsUse[-del,]
	
	  ### Separate out ref and alt alleles, CC and founder counts
		ccMiceNum <- grep("CC", colnames(exonsUse)) 
		ccMice <- unique(unlist(strsplit(colnames(exonsUse)[ccMiceNum], "[.]"))[c(T,F)])
		if(is.null(myCC)) myCC=ccMice

		keepccMiceNum <- ccMiceNum#[which(ccMice %in% myCC)]
		impCC <- myCC[which(!myCC %in% ccMice)]
		ccMiceCounts <- exonsUse[,keepccMiceNum]
		ccMiceCounts_ref <- apply(ccMiceCounts,2,function(x) as.numeric(paste(unlist(strsplit(x, "/"))[c(T,F)])))
		ccMiceCounts_alt <- apply(ccMiceCounts,2,function(x) as.numeric(paste(unlist(strsplit(x, "/"))[c(F,T)])))

		founderMiceNum <- unlist(lapply(founders, function(x) grep(x, colnames(exonsUse))[1]))
		founderMice <- colnames(exonsUse)[founderMiceNum]
		founderCounts <- exonsUse[,founderMiceNum]
		founderCounts_ref <- apply(founderCounts,2,function(x) as.numeric(paste(unlist(strsplit(x, "/"))[c(T,F)])))
		founderCounts_alt <- apply(founderCounts,2,function(x) as.numeric(paste(unlist(strsplit(x, "/"))[c(F,T)])))
		
		
		

		snp <- exonsUse[, c("Gene","Chromosome","rsId","Position","ProbeSeq","TranscriptID")]
		
		seq <- as.data.frame(do.call("rbind",strsplit(as.character(exonsUse$ProbeSeq), "\\[|\\]|/")))
		colnames(seq) <- c("end5","ref","alt","end3")
		seq$refseq <- paste0(seq$end5, seq$ref, seq$end3)
		seq$altseq <- paste0(seq$end5, seq$alt, seq$end3)
		
		comb <- data.frame(snp, seq)

		sdp <- do.call("rbind", strsplit(exonsUse$FounderSDP,""))
		isAlt <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 1), which(founderCounts_ref[x,] <= 3 )))
		isRef <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 0), which(founderCounts_ref[x,] > 3)))
		isAlt2 <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 1), which(founderCounts_alt[x,] > 3)))
		isRef2 <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 0), which(founderCounts_alt[x,] <= 3)))
		altTrue <- sapply(1:nrow(sdp), function(x) ifelse(identical(isAlt[[x]], isAlt2[[x]]), 1, 0), simplify = T)
		refTrue <- sapply(1:nrow(sdp), function(x) ifelse(identical(isRef[[x]], isRef2[[x]]), 1, 0), simplify = T)
		bothTrue <- intersect(which(altTrue == 1), which(refTrue == 1))

		### trim out variants that disagree with founder sdp
		var_trim <- list()
		var_trim[["seq"]] <- comb[bothTrue,] 
		var_trim[["sdp"]] <- sdp[bothTrue,]
		var_trim[["founder_ref"]] <- founderCounts_ref[bothTrue,]
		var_trim[["founder_alt"]] <- founderCounts_alt[bothTrue,]
		var_trim[["CC_ref"]] <- ccMiceCounts_ref[bothTrue,]
		var_trim[["CC_alt"]] <- ccMiceCounts_alt[bothTrue,]
		
	  # impute for CC014 - no sequence data	
		if((!is.null(founderHaps)) && length(impCC) > 0 && impute){
		  for(c in 1:length(impCC)){
		    tmp_r <- rep(NA, nrow(var_trim[["seq"]]))
		    tmp_a <- rep(NA, nrow(var_trim[["seq"]]))
		    ind <- grep(impCC[c], names(founderHapsUse))
		    hap <- founderHapsUse[[ind]] %>% dplyr::filter(hap != 0) %>% dplyr::mutate(start = as.numeric(paste(start)), end = as.numeric(paste(end)))
		    hap$hap <- as.factor(hap$hap)
		    
		    for(k in 1:nrow(hap)) {
		      assign <- intersect(which(var_trim[["seq"]]$Position > hap$start[k]),
		                          which(var_trim[["seq"]]$Position < hap$end[k]))
		      hapUse <- unique(unlist(strsplit(as.character(hap$hap[k]), "")))
		      if(length(hapUse) == 1 && hapUse != "Z"){
		        tmp_r[assign] <- var_trim[["founder_ref"]][assign, as.numeric(as.factor(hapUse))]	
		        tmp_a[assign] <- var_trim[["founder_alt"]][assign, as.numeric(as.factor(hapUse))]
		      }
		    }
		    var_trim[["CC_ref"]] <- cbind(var_trim[["CC_ref"]], tmp_r)
		    var_trim[["CC_alt"]] <- cbind(var_trim[["CC_alt"]], tmp_a)
		    colnames(var_trim[["CC_ref"]])[ncol(var_trim[["CC_ref"]])] <- paste0(impCC[c],".impute")
		    colnames(var_trim[["CC_alt"]])[ncol(var_trim[["CC_alt"]])] <- paste0(impCC[c],".impute")
		  }
		}
		
    # flag hets (CC is alt and ref) or zeros (CC is neither) from set
		hetsOrZeros <- sapply(1:nrow(var_trim[["CC_ref"]]), function(x){
		  sapply(1:ncol(var_trim[["CC_ref"]]), function(y){
			ifelse(var_trim[["CC_ref"]][x,y] > up_bound | 
			       var_trim[["CC_alt"]][x,y] > up_bound, "M", 
		   	ifelse(var_trim[["CC_ref"]][x,y] > lo_bound &
			         var_trim[["CC_alt"]][x,y] > lo_bound, "H", 
		        ifelse(var_trim[["CC_ref"]][x,y] <= lo_bound & 
			             var_trim[["CC_alt"]][x,y] <= lo_bound, "Z", 
			        ifelse(var_trim[["CC_ref"]][x,y] > var_trim[["CC_alt"]][x,y], "R", "A"))))
		    })
		})
	
		# 0 is flag, 1 is keep
		# Z = zero, H = het, M = multiple
		# B = out of bounds, S = segregating
		# eB = extreme bounds, mB = my bounds
		# r = ref, a = alt
    #	 keep positive counts for each CC, i.e., only the allele that should be present
    #  unless are hets
		hetsOrZeros <- t(hetsOrZeros)	
		posCounts <- sapply(1:nrow(var_trim[["CC_ref"]]), function(x){
		  sapply(1:ncol(var_trim[["CC_ref"]]), function(y){
			ifelse(hetsOrZeros[x,y] == "H", 
			       var_trim[["CC_ref"]][x,y] + var_trim[["CC_alt"]][x,y],
			  ifelse(hetsOrZeros[x,y] %in% c("M","Z"), NA, 
			    ifelse(var_trim[["CC_ref"]][x,y] > var_trim[["CC_alt"]][x,y],
			           var_trim[["CC_ref"]][x,y], var_trim[["CC_alt"]][x,y])))
			})
	  })

		posCounts <- t(posCounts)
	
  	# check distribution of counts in CC ################
		
		keep <- which(apply(posCounts, 2, function(x) ifelse(sum(is.na(x)) == length(x), F, T)))
		posCounts <- posCounts[,keep]
		hetsOrZeros <- hetsOrZeros[,keep]

		kmer_mu <- apply(posCounts, 2, function(x) mean(x[which(!is.na(x))]))	
		kmer_med <- apply(posCounts, 2, function(x) median(x[which(!is.na(x))]))
		kmer_length <- apply(posCounts, 2, function(x) length(x[which(!is.na(x))]))	
	  kmer_var <- apply(posCounts, 2, function(x) var(x[which(!is.na(x))]))
	  kmer_range <- apply(posCounts, 2, function(x) range(x[which(!is.na(x))]))
	  kmer_quant <- apply(posCounts, 2, function(x) quantile(x[which(!is.na(x))]))
		lo_bounds_p <- unlist(lapply(kmer_mu, function(x) qpois(0.005, x)))
		hi_bounds_p <- unlist(lapply(kmer_mu, function(x) qpois(0.995, x)))
		lo_bounds_n <- unlist(sapply(1:length(kmer_mu), function(x) qnorm(0.16, kmer_mu[x], sqrt(kmer_var[x]))))
		hi_bounds_n <- unlist(sapply(1:length(kmer_mu), function(x) qnorm(0.84, kmer_mu[x], sqrt(kmer_var[x]))))
		
		hi_bounds = hi_bounds_p    
		lo_bounds = rep(lo_bound*2, length(hi_bounds))    
		 
		kmer_stats <- data.frame(mu = kmer_mu, med = kmer_med, var = kmer_var, 
		                         sd = apply(posCounts, 2, function(x) sd(x, na.rm=T)),
		                         min = apply(posCounts, 2, function(x) min(x[which(!is.na(x))])),
		                         max = apply(posCounts, 2, function(x) max(x[which(!is.na(x))])), 
		                         l_bound_p = lo_bounds_p, h_bound_p = hi_bounds_p,
		                         l_bound_n = lo_bounds_n, h_bound_n = hi_bounds_n)
	
		keep <- seq(1,ncol(posCounts))
		
		remove <- sapply(keep, function(x)  
		  union(which(posCounts[,x] > hi_bounds[x]), which(posCounts[,x] < lo_bounds[x])))
		remove_extreme <- sort(unique(unlist(remove)))
		
		kmer_stats$num_removed <- 0
		kmer_stats$num_removed[keep] <- unlist(lapply(remove, length))
		###
		hetsOrZeros[remove_extreme, ] <- paste("eB",hetsOrZeros[remove_extreme,],sep=",")
		
		for(e in 1:length(keep)){
		  hetsOrZeros[,keep[e]][remove[[e]]] <- "mB"
		}

		
		##############
	  ### make sure counts match with expected founders #################3
    var_trim[["founder_map"]] <- matrix(NA, nrow=nrow(hetsOrZeros), ncol=ncol(hetsOrZeros))
		
		if(check_CC_sdp){
		  for(d in 1:ncol(hetsOrZeros)){
		    CC=unlist(strsplit(colnames(hetsOrZeros)[d],"[.]"))[c(T,F)]
		    ind <- grep(CC, names(founderHapsUse))
		    hap <- founderHapsUse[[ind]] %>% dplyr::filter(hap != 0) %>%
		      dplyr::mutate(start = as.numeric(paste(start)), end = as.numeric(paste(end))) %>% arrange(start)
		    hap$hap <- as.character(hap$hap)
		    var_trim[["founder_map"]][,d] <- sapply(1:nrow(hetsOrZeros), function(y){
		      lo=which(var_trim[["seq"]]$Position[y] > hap$start)
		      hi=which(var_trim[["seq"]]$Position[y] < hap$end)
		      groupInd=ifelse(length(lo) == 0 | length(hi) == 0, NA, 
		                      intersect(which(var_trim[["seq"]][y,]$Position > hap$start), 
		                                which(var_trim[["seq"]][y,]$Position < hap$end)))
		      ifelse(is.na(groupInd), NA, hap$hap[groupInd])
		    }, simplify=T)
		    
		    hetsOrZeros[,d] <- sapply(1:nrow(hetsOrZeros), function(y){
		      found_lets <- unlist(strsplit(var_trim[["founder_map"]][y,d],""))
		      found_nums <- as.numeric(factor(found_lets, levels=c("A","B","C","D","E","F","G","H")))
		      orig <- hetsOrZeros[y,d]
		      if (any(is.na(found_nums)) | length(unique(found_nums)) > 1) orig <- paste("fS_r",orig,sep=",")
		                                                                      # founder segregating in region fS_r
		      orig <- ifelse(length(grep("A|R", hetsOrZeros[y,d])) == 0, orig, 
            ifelse(all(var_trim$founder_ref[y,found_nums] > lo_bound), ## if both founders have ref
                ifelse((var_trim[["CC_ref"]][y,colnames(hetsOrZeros)[d]] > lo_bound), orig, # and CC is ref
                        paste("fD",orig,sep=",")), # if CC is not ref, fD for disagree w/ founder
		            ifelse(all(var_trim$founder_ref[y,found_nums] <= lo_bound), ## or if both founders are alt
		                ifelse(var_trim[["CC_ref"]][y,colnames(hetsOrZeros)[d]] <= lo_bound,orig, # and CC is alt
		                       paste("fD", hetsOrZeros[y,d],sep=",")), paste("fS_l",orig,sep=","))))
		      orig
			  }, simplify=T)
		  }

		
		}
		
		var_trim$CC_alt <- var_trim$CC_alt[,match(colnames(posCounts), colnames(var_trim$CC_alt))]
		var_trim$CC_ref <- var_trim$CC_ref[,match(colnames(posCounts), colnames(var_trim$CC_ref))]
		
	  var_trim$hetsOrZeros <- hetsOrZeros
	  var_trim$posCounts <- posCounts
	
	  fin[[j]] <- var_trim
	  }
    return(fin)
}



##############

### get patterns for kmers to keep and which to disregard in each CC; i.e. var_trim$hetsOrZeros
plot_counts <- function(snps, counts, pupInfo, CC_labels, lo_bound = 3,
                        keepEB=T, keepFS=T, keepFD=T){
  ref_counts = alt_counts = 
    matrix(NA, nrow=length(unique(counts$k.mer_ref)), ncol=length(unique(counts$Pup)))
  pups <- unique(counts$Pup)
  temp_ref <- sapply(1:ncol(ref_counts), function(x){
    counts %>% filter(Pup.ID == pups[x]) -> temp
    temp[match(as.character(snps$seq.refseq), as.character(temp$k.mer_ref)),]$counts_ref
  }, simplify=T)
  
  temp_alt <- sapply(1:ncol(alt_counts), function(x){
    counts %>% filter(Pup.ID == pups[x]) -> temp
    temp[match(snps$seq.refseq, temp$k.mer_ref),]$counts_alt
  }, simplify=T)
  
  consensus <- snps$seq.consensus
  
  colnames(temp_alt) = colnames(temp_ref) = pups
  rownames(temp_alt) = rownames(temp_ref) = consensus
  
  if(length(grep("hetsOrZeros|keep",colnames(snps))) > 0){
    CC_keep <- data.frame(snps[,grep("keep|hetsOrZeros",colnames(snps))])
    rownames(CC_keep) <- snps$seq.consensus
    seq = snps[,grep("seq", colnames(snps))]
    xc <- list(ref = temp_ref, alt = temp_alt, 
               seq = seq,
               CC_keep = CC_keep)
  } else{
    xc <- list(ref = temp_ref, alt = temp_alt, 
               seq = snps[,grep("seq", colnames(snps))],
               CC_keep = snps[,grep("sdp", colnames(snps))])
  }
  
  
  remove <- intersect(which(apply(xc$ref, 1, sum) == 0), which(apply(xc$alt, 1, sum) == 0))
  xc <- lapply(xc, function(x) x[-remove,])
  
  pupn <- as.numeric(unlist(strsplit(as.character(colnames(xc$ref)),"_"))[c(F,F,T)])
  
  tmpList <- apply(CC_labels[,grep("CC",colnames(CC_labels))], 2, function(x){
    unlist(lapply(x, function(y) grep(y, colnames(xc$CC_keep))[1]))
  })
  
  CC_labels[,grep("CC", colnames(CC_labels))] <- apply(CC_labels[,grep("CC", colnames(CC_labels))], 2, as.character)
  CC_labels$CC1_col <- tmpList[,1]
  CC_labels$CC2_col <- tmpList[,2]
  CC_labels <- CC_labels[match(pupn, CC_labels$Pup.ID),]
  switch <- which(!as.character(CC_labels$CC.1) %in% as.character(xce$num_cc))
  switch = c()
  CC_labels$CC_n <- CC_labels$CC.1
  CC_labels$CC_d <- CC_labels$CC.2
  CC_labels$n_col <- CC_labels$CC1_col
  CC_labels$d_col <- CC_labels$CC2_col
  CC_labels$CC_n[switch] <- CC_labels$CC.2[switch]
  CC_labels$CC_d[switch] <- CC_labels$CC.1[switch]
  CC_labels$n_col[switch] <- CC_labels$CC2_col[switch]
  CC_labels$d_col[switch] <- CC_labels$CC1_col[switch]
  
  remove <- grep("[0-9]", colnames(CC_labels))
  CC_labels <- CC_labels[,-remove]
  xc$pups_use <- CC_labels#[-remove,]
  
  tossOut <- c("mB", "Z", "H","M","NA")
  if(!keepEB){
    tossOut <- c(tossOut, "eB")
  }
  if(!keepFS){
    tossOut <- c(tossOut, "fS")
  }
  if(!keepFD){
    tossOut <- c(tossOut, "fD")
  }
  
  ratios <- list()
  metrics <- list()
  for(i in 1:ncol(xc$ref)){
    cols <- grep("col", colnames(xc$pups_use))
    alts <- list()
    refs <- list()
    ok <- list()
    metrics[[paste(xc$pups_use$Pup.ID[i])]] <- setNames(data.frame(matrix(ncol = 4, nrow = 3)), c("n_problems", "n_refs", "n_alts", "n_ok"))
    rownames(metrics[[paste(xc$pups_use$Pup.ID[i])]]) <- unlist(c(xc$pups_use[i, c("CC_n", "CC_d")], "both"))
    
    for(j in c(1:2)){
      useCol <- xc$pups_use[i,cols[j]]
      if(!is.na(useCol)){
        if(length(grep("keep|posCounts",colnames(snps))) > 0){
          keep <- as.character(paste(xc$CC_keep[,useCol]))
        } else {
          keep <- ifelse(xc$CC_keep[,useCol] == 1,"A", 
                         ifelse(xc$CC_keep[,useCol] == 0, "R", xc$CC_keep[,useCol])) 
        }
        problems <- grep(paste(tossOut, collapse="|"), keep)
        problems <- c(problems, which(is.na(keep)))
        
        if(length(problems)>0) keep[problems] <- ""
        
        refs[[j]] <- grep("R", keep)
        alts[[j]] <- grep("A", keep)
        ok[[j]] <- which(keep != "")
        
        metrics[[paste(xc$pups_use$Pup.ID[i])]][j,] <- c(length(problems), length(refs[[j]]), length(alts[[j]]), length(ok[[j]]))
        
      } else {
        refs[[j]] <- ""
        alts[[j]] <- ""
        ok[[j]] <- ""
      }
    } 
    usable = NULL
    if(!any(unlist(lapply(refs, is.na)))){
      ok_both <- intersect(ok[[1]], ok[[2]])
      usable <- sort(intersect(ok_both, union(refs[[1]][which(refs[[1]] %in% alts[[2]])], alts[[1]][which(alts[[1]] %in% refs[[2]])])))
      metrics[[paste(xc$pups_use$Pup.ID[i])]][3,1:2] <- c(length(ok_both), length(usable))
    }
    
    if(length(usable) > 0){
      sums <- xc$alt[usable,i] + xc$ref[usable,i]
      usable <- usable[which(sums > lo_bound*2)]
      col1_is_ref <- grep("R", xc$CC_keep[usable, xc$pups_use[i,cols[1]]])
      col2_is_ref <- grep("R", xc$CC_keep[usable, xc$pups_use[i,cols[2]]])
      col1_is_alt <- grep("A", xc$CC_keep[usable, xc$pups_use[i,cols[1]]])
      col2_is_alt <- grep("A", xc$CC_keep[usable, xc$pups_use[i,cols[2]]])
      # check if col1_is_ref = col2_is_alt, col2_is_ref = col1_is_alt
      CC_numer = CC_denom = rep(0, length(usable))
      if(length(col1_is_ref) > 0){
        CC_numer[col1_is_ref] <- xc$ref[usable[col1_is_ref], i]
        CC_denom[col2_is_alt] <- xc$alt[usable[col2_is_alt], i]
      }
      if(length(col1_is_alt) > 0){
        CC_numer[col1_is_alt] <- xc$alt[usable[col1_is_alt], i]
        CC_denom[col2_is_ref] <- xc$ref[usable[col2_is_ref], i]
      }
      CC_numer2 <- sapply(1:length(usable), function(y){
        ifelse(length(grep("A", keep[usable[y]])) == 1, xc$ref[usable[y],i], xc$alt[usable[y],i])})
      CC_denom2 <- sapply(1:length(usable), function(y){
        ifelse(length(grep("A", keep[usable[y]])) == 1, xc$alt[usable[y],i], xc$ref[usable[y],i])})
      # check if CC_numer = CC_numer, CC_denom = CC_denom
      ratio <- matrix(NA, nrow=length(keep), ncol=2) 
      ratio[usable,] <- c(CC_numer, CC_denom)                   
      colnames(ratio) = c("CC_1", "CC_2")
      
      ratios[[i]] <- ratio
      metrics[[paste(xc$pups_use$Pup.ID[i])]][3,3:4] <- c(length(usable), length(ratio))
    } else {
      ratios[[i]] <- NA
    }
  }
  
  plotList <- list()
  plotList <- sapply(1:length(ratios), function(x){ 
    if(!is.null(ratios[[x]])){
      if(!all(is.na(as.matrix(ratios[[x]])))){
        s1 <- cbind(xc$seq, as.matrix(ratios[[x]]), 
                    sum = (xc$ref[,x] + xc$alt[,x]))
        s1 <- s1[intersect(which(!is.na(s1$CC_1)), which(s1$sum > lo_bound)),]
        if(nrow(s1) > 0){
          s1$Pup <- rep(colnames(xc$ref)[x])
          s1$CCs <- rep(paste(xc$pups_use$CC_n[x], xc$pups_use$CC_d[x], sep="/"))
          
          s1$kRat <- s1$CC_1 / s1$sum
          s1$kLB = qbeta(0.05,s1$CC_1 + 1, s1$CC_2 + 1)
          s1$kUB = qbeta(0.95,s1$CC_1 + 1, s1$CC_2 + 1)
          s1$kLB[which(is.nan(s1$kRat))] = s1$kUB[which(is.nan(s1$kRat))] <- 0
          s1
        } else{
          NULL
        }
      } 
    }
  }, simplify=F)
  
  
  plotGenes <- lapply(plotList, function(x){
    if(!is.null(x)){
      x %>% 
        group_by(seq.Gene) %>% 
        dplyr::select(-one_of("seq.TranscriptID","seq.end5","seq.end3","seq.consensus")) %>%
        mutate(minPos = min(seq.Position), maxPos = max(seq.Position)) %>% 
        group_by(Pup, CCs, seq.Chromosome, seq.Gene, minPos, maxPos) %>% 
        summarize(sumTot = sum(sum), sum1 = sum(CC_1), sum2 = sum(CC_2)) %>%
        mutate(gRat = sum1/sumTot, 
               gLB = qbeta(0.05, sum1 + 1, sum2 + 1),
               gUB = qbeta(0.95, sum1 + 1, sum2 + 1))
    } else {
      NULL
    }
  })
  
  plotGenes <- unique(do.call("rbind", plotGenes))
  plotGenes$logSumTot <- log(plotGenes$sumTot)
  plotGenes$Pup.ID <- as.numeric(unlist(strsplit(plotGenes$Pup,"_"))[c(F,F,T)])
  
  ### add pup demographic info
  pupInfo$Diet <- gsub(" $","", pupInfo$Diet, perl=T)
  
  
  plotGenes <- left_join(plotGenes, pupInfo, by="Pup.ID") %>%
    mutate(Xist = if_else(seq.Gene == "Xist", T, F))
  
  plotMat <- unique(do.call("rbind", plotList))
  plotMat$logSum <- log(plotMat$sum)
  plotMat$Pup.ID <- as.numeric(unlist(strsplit(plotMat$Pup,"_"))[c(F,F,T)])
  
  plotMat <- left_join(plotMat, pupInfo, by="Pup.ID") %>%
    mutate(Xist = if_else(seq.Gene == "Xist", T, F))
  
  
  return(list(plotMat=plotMat, plotGenes=plotGenes, metrics=metrics))
}

#####################

