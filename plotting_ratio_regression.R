source("plot.hpd.R")
source("jags_kmers_source.R")

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

CC_pairs <- c("CC001/011", "CC041/051", "CC017/004", "CC023/047", "NA/NA", "CC026/006", "CC014/003", "CC035/062", "CC032/042", "CC005/040")
CC_pairs <- paste(xce$num_cc, xce$dem_cc, sep = "/")

plot_simple_ratios <- function(plotMat, save_path="./plots.pdf", plotOn=T){
  p <- list()
  for(i in 1:length(unique(plotMat$CC_lab))){
    cc <- unique(plotMat$CC_lab)[i]
    chr <- unique(plotMat$seq.Chromosome)
    title <- paste0("Chr ",chr,"-wide proportion of parental alleles in ",cc)
    if(length(grep("seq.Position", colnames(plotMat))) == 0) {
      plotMat %>% mutate(seq.Position = mean(minPos, maxPos)) -> plotMat
    }
    rixdat <- plotMat %>% filter(CC_lab == cc) %>%
      mutate(seq.Position = as.numeric(seq.Position)) %>%
      arrange(dir, Pup.ID, seq.Position) 
    #rixdat$ratio <- as.numeric(rixdat$ratio)
    rixdat$Pup <- factor(rixdat$Pup, levels=unique(rixdat$Pup))
    
    p[[i]] = ggplot(data=rixdat, aes(x=seq.Position, y=kRat, alpha=logSum, 
                                     color=Xist, shape=dir)) + 
      geom_point() + #
      geom_segment(aes(x=seq.Position, xend=seq.Position, y = kLB, yend=kUB), size=0.5) + 
      ylab("Ratio") + xlab("Position") +
      theme_bw() + 
      ggtitle(title) + 
      geom_hline(yintercept=0.5, col="gray") + 
      theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
      facet_wrap( ~ Pup, scales = "fixed") + ylim(c(0,1)) + 
      scale_colour_manual(values=c("royalblue", "sandybrown"))
  }
  
  if(plotOn){
    pdf(file.path(save_path), width=8.5, height=11)
    print(p)
    dev.off()
  }
  
  return(p)
}

plot_ratios <- function(plotMat, summary=NULL, 
                        gene=T, byPup=T, ylim=NULL, xist_opp=F){
  p <- list()
  if(!is.null(summary) & length(grep("RRIX", toupper(names(summary)))) > 0){
    plotMat <- plotMat %>% filter(RRIX %in% gsub("[A-Z]|_|-", "", toupper(names(summary))))
  } 
  plotGene <- plotMat %>% group_by(Pup.ID, seq.Gene) %>% 
    mutate(minPos = min(seq.Position), maxPos = max(seq.Position), 
           RIX = RRIX) %>% 
    group_by(seq.Chromosome, Pup.ID, CC_lab, seq.Gene, minPos, maxPos, dir, RIX, Xist, Diet) %>% 
    dplyr::summarize(sumTot = sum(sum), sum1 = sum(fin1), sum2 = sum(fin2)) %>%
    dplyr::mutate(kRat = sum1/sumTot, 
           kLB = qbeta(0.05, sum1 + 1, sum2 + 1),
           kUB = qbeta(0.95, sum1 + 1, sum2 + 1),
           logSum = log(sumTot),
           seq.Position = mean(minPos, maxPos))
  if(is.null(summary)){
    plotGene %>% group_by(Pup.ID) %>% 
      summarize(sumTotp = sum(sumTot), sum1p = sum(sum1), sum2p = sum(sum2)) %>%
      mutate(mean = sum1p/sumTotp, 
             lower = qbeta(0.05, sum1p + 1, sum2p + 1),
             upper = qbeta(0.95, sum1p + 1, sum2p + 1)) -> use_sum_all
    use_sum_all = plotGene %>% group_by(Pup.ID) %>%
      mutate(mean = weighted.mean(kRat, logSum), lower = weighted.mean(kLB, logSum), upper = weighted.mean(kUB, logSum))
  }
  if(gene){
    plotMat <- plotGene
  }
  for(i in 1:length(unique(plotMat$CC_lab))){
    cc <- unique(plotMat$CC_lab)[i]
    if(length(grep("a|b", unique(plotMat$RIX))) > 0){
      rix <- as.numeric(unique(gsub("a|b", "", unique(plotMat$RIX))))[i]
    } else {
      rix = unique(plotMat$RIX)[i]
    }
    chr <- unique(plotMat$seq.Chromosome)

    if(length(grep("seq.Position", colnames(plotMat))) == 0) {
      plotMat %>% mutate(seq.Position = mean(minPos, maxPos)) -> plotMat
    }
    
    rixdat <- plotMat %>% filter(CC_lab == cc) %>%
      mutate(seq.Position = as.numeric(seq.Position)) %>%
      arrange(dir, Pup.ID, seq.Position) 
    rixdat$Pup.ID <- factor(rixdat$Pup.ID, levels=unique(rixdat$Pup.ID))
    if(is.null(summary)){
      use_sum <- use_sum_all %>% filter(Pup.ID %in% unique(rixdat$Pup.ID))
    } else {
      if(length(grep("RIX", toupper(names(summary)))) > 0){
        use_sum_ind <- match(rix ,gsub("[A-Z]|_|-", "", toupper(names(summary))))
      } else {
        use_sum_ind <- match(rix, toupper(names(summary)))
        if(is.na(use_sum_ind)){
          sep_rix = unlist(strsplit(as.character(rix), "_|/"))
          sum_names = do.call("rbind",lapply(names(summary), function(x) unlist(strsplit(as.character(x),"_|/"))))
          use_sum_ind = which(apply(sum_names, 1, function(x) ifelse(all(sep_rix %in% x), T, F)))
        }
      }
      use_sum <- summary[[use_sum_ind]]$summary$mu_p
      use_sum$Pup.ID = rownames(use_sum)
      
      if (!gene){
        rixdat$kRat <- as.numeric(summary[[use_sum_ind]]$summary$mu_g$Mean[match(paste(rixdat$Pup.ID, rixdat$seq.Gene, sep="_"),
                                                                                 rownames(summary[[use_sum_ind]]$summary$mu_g))])
        rixdat$kLB <- as.numeric(summary[[use_sum_ind]]$summary$mu_g$lower[match(paste(rixdat$Pup.ID, rixdat$seq.Gene, sep="_"),
                                                                                 rownames(summary[[use_sum_ind]]$summary$mu_g))])
        rixdat$kUB <- as.numeric(summary[[use_sum_ind]]$summary$mu_g$upper[match(paste(rixdat$Pup.ID, rixdat$seq.Gene, sep="_"),
                                                                                 rownames(summary[[use_sum_ind]]$summary$mu_g))])
      }
    }
    
    colnames(use_sum) = tolower(colnames(use_sum))
    rixdat$pup_lower <- as.numeric(use_sum$lower[match(rixdat$Pup.ID, use_sum$pup.id)])
    rixdat$pup_upper <- as.numeric(use_sum$upper[match(rixdat$Pup.ID, use_sum$pup.id)])
    rixdat$pup_mean <- as.numeric(use_sum$mean[match(rixdat$Pup.ID, use_sum$pup.id)])
    
    if(xist_opp){
      rixdat$kRat = apply(rixdat, 1, function(x) as.numeric(ifelse(x["Xist"], 1- as.numeric(x["kRat"]), x["kRat"])))
      rixdat$kLB = apply(rixdat, 1, function(x) as.numeric(ifelse(x["Xist"], 1- as.numeric(x["kLB"]),  x["kLB"])))
      rixdat$kUB = apply(rixdat, 1, function(x) as.numeric(ifelse(x["Xist"], 1- as.numeric(x["kUB"]), x["kUB"])))
    }
    
    if(byPup){
      mvec <- unique(rixdat$Pup.ID)
    } else {
      mvec <- list(mvec = unique(rixdat$Pup.ID))
      
    } 
    p[[cc]] <- list()
    for(m in mvec){
      pup <- ifelse(byPup, paste0("Pup ", m, " "), "")
      title <- paste0("Chr ",chr,"-wide proportion of parental alleles in ",pup, cc)
      
      rixdat_use <- rixdat %>% filter(Pup.ID == m)
      p[[cc]][[m]] = ggplot(data=rixdat_use, aes(x=seq.Position, alpha=logSum, 
                                                color=Xist, shape=factor(dir))) + 
        geom_point(aes(y = kRat)) + #
        geom_segment(aes(x=seq.Position, xend=seq.Position, y = kLB, yend=kUB), size=0.5) + 
        ylab("Ratio") + xlab("Position") +
        theme_bw() + 
        ggtitle(title) + 
        scale_x_continuous(labels=function(n){paste(n/1000000, "M")}) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
        geom_hline(yintercept=0.5, col="azure", linetype= "dotted") + 
        scale_colour_manual(values=c("royalblue", "sandybrown")) +
        geom_ribbon(aes(x=seq.Position, ymin=pup_lower, ymax=pup_upper), inherit.aes = T, 
                    alpha=0.2, col="gray") + 
        geom_hline(aes(yintercept = pup_mean), col="navy")
      if(!is.null(ylim)){
        p[[cc]][[m]] <- p[[cc]][[m]] + ylim(ylim) 
      }
    }
  }
    
  return(p)
}

plot_gene_ratios <- function(plotGenes){
  p <- list()
  dat <- list()
  for(i in 1:length(unique(plotGenes$CC_lab))){
    cc <- unique(plotGenes$CC_lab)[i]
    title <- paste(cc)
    
    rixdat <- plotGenes %>% filter(CC_lab == cc) %>%
      group_by() %>%
      mutate(minPos = as.numeric(minPos), maxPos = as.numeric(maxPos),
             meanPos = (minPos+maxPos)/2) %>%
      arrange(dir, Pup, minPos) 
    rixdat$Pup <- factor(rixdat$Pup, levels=unique(rixdat$Pup))
    p[[i]] = ggplot(data=rixdat, aes(x=meanPos)) + #, size=logSumTot
      geom_point(aes(y=gRat, alpha=logSumTot, color=Xist)) + 
      geom_segment(aes(xend=meanPos, y = gLB, yend=gUB, alpha=logSumTot, color=Xist), size=1) + 
      ylab("Ratio") + xlab("Position") +
      theme_bw() + 
      ggtitle(title) + 
      theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
      facet_wrap( ~ Pup, scales = "fixed") + ylim(c(0,1)) + 
      scale_colour_manual(values=c("royalblue", "sandybrown"))
    
    if("X.Mean" %in% colnames(plotGenes)){
      rixdat[which(is.na(rixdat$X.mean)), c("X.mean", "X.Low", "X.Up")] <- 0
      
      p[[i]] = p[[i]] + geom_ribbon(aes(ymin=X.Low, ymax=X.Up))
    }
    dat[[i]] <- rixdat
  }

  return(list(p=p, dat=dat))
}

plot_reg_results <- function(df, kmer_counts,ylim=NULL){
  plotMu <- list()
  plotDiet <- list()
  plotDietPORIX <- list()
  
  if(length(grep("RIX", toupper(names(df)))) > 0){
    rixes <- as.numeric(paste(factor(gsub("[A-Z]|_|-","",toupper(names(df))), levels=levels(kmer_counts$RRIX))))
  } else {
    rixes <- toupper(names(df))
  }
  
  chr <- unique(kmer_counts$seq.Chromosome)
  CC_lab <- unique(kmer_counts$CC_lab)[sort.list(as.numeric(unique(kmer_counts$RRIX)))]

  for(i in 1:length(df)){
    if(length(grep("RIX", toupper(names(df)))) > 0){
      CC <- CC_lab[sort.list(as.numeric(gsub("[A-Z]|_|-", "", toupper(names(df)))))[i]]
    } else {
      CC <- CC_lab[i]
    }
    
    plotMu[[i]] <- plot.mine.hpd(coda.object=df[[i]]$summary$mu_p,
                                 wanted=df[[i]]$codes$Level[which(df[[i]]$codes$Variable == "PUP.ID")],
                                 prob.wide=0.95,
                                 yint=0.5,
                                 xlab="HPD interval",
                                 title = paste0(CC, " pup-wide skew estimates on chr ",chr),
                                 grouped = df[[i]]$indices$DIR, 
                                 names=NULL,
                                 ylim=ylim)
    if(!is.null(df[[i]]$summary$b_diet)){
      plotDiet[[i]] <- plot.mine.hpd(coda.object=df[[i]]$summary$b_diet,
                                     wanted=df[[i]]$codes$Level[which(df[[i]]$codes$Variable == "DIET")],
                                     prob.wide=0.95,
                                     xlab="HPD interval",
                                     title = paste0(CC, " diet-effects on chr ",chr),
                                     names=NULL)
    }

  
    if(!is.null(df[[i]]$summary$b_dietPORIX)){
      plotDietPORIX[[i]] <- plot.mine.hpd(coda.object=df[[i]]$summary$b_dietPORIX,
                                          wanted=df[[i]]$codes$Level[which(df[[i]]$codes$Variable == "DIETRIX")],
                                          prob.wide=0.95,
                                          xlab="HPD interval",
                                          title = paste0(CC, " diet-by-PO-effects on chr ",chr),
                                          names=NULL)
    }
    
  }
  
  allGmu <- c()
  allBetaPO <- c()
  plotGmu = plotBetaPO = NULL
  if(length(df) > 1){
    for(i in 1:length(df)){
      if("b0" %in% names(df[[i]]$summary)){
        arr <- unlist(exp(df[[i]]$summary$b0)/(1+exp(df[[i]]$summary$b0)))
      } else if("mu_r" %in% names(df[[i]]$summary)){
        arr <- df[[i]]$summary$mu_r
        
      }
      allGmu <- rbind(allGmu, arr)
      arr <- unlist(df[[i]]$summary$b_PORIX)
      hasPO = ifelse(is.null(arr), F, T)
      if(hasPO){
        allBetaPO <- rbind(allBetaPO, arr)
      }
    }
    rownames(allGmu) = gsub("rix_","", names(df))
    if(hasPO) rownames(allBetaPO) = gsub("rix_","", names(df))
    getCodes = c("RRIX", "Diet", "Pup.ID")
    if(hasPO) getCodes = c(getCodes, "DietRIX", "dir")
    bigEncoded <- getEncoding(kmer_counts, getCodes)  
    plotGmu <- plot.mine.hpd(coda.object=allGmu,        
                             wanted=bigEncoded$Level[which(bigEncoded$Variable == "RRIX")],   
                             prob.wide=0.95,
                             yint = 0.5,
                             xlab="HPD interval",
                             title = paste0("Overall estimates of skew per RIX on chr ",chr),
                             names=CC_lab[1:length(rixes)],
                             ylim=ylim)
    if(hasPO){
      plotBetaPO <- plot.mine.hpd(coda.object=allBetaPO,
                                  wanted=bigEncoded$Level[which(bigEncoded$Variable == "RRIX")],
                                  prob.wide=0.95,
                                  yint = 0,
                                  title = paste0("Overall estimates of PO-effects per RIX on chr ",chr),
                                  names=CC_lab[1:length(rixes)])
    }
  }
  
  return(list(mu=plotMu, b_diet=plotDiet, b_dietPORIX=plotDietPORIX,
         mu_r=plotGmu, b_PO=plotBetaPO))
}

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}
