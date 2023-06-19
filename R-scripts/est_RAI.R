est_RAI <- function(n, outputtype) {
  # function to plot estimated rank abundance distribution:
  
  # n is the abundance vector (number of individuals per species)
  # outputtype: 'fit' provides the goodness-of-fit;
  #             'plot' provides the SAD figure;
  #             'estimates' provides the METE estimates of the abundances
  #             'differences' provides the differences per rank between the 
  #                           actual abundance and the METE estimate
  
  df2 <- data.frame(spp = 1:length(n), 
                    abund = sort(n, decreasing = T),
                    rank = 1:length(n)) 
  
  esf1 <- meteESF(spp=df2$spp, abund=df2$abund)
  sad1 <- sad(esf1)
  sad1 # print function returns useful summary
  
  n2 <- sort(sad1$r(length(n)), decreasing = T)
  for (i in 2:1000) {
    n2 <- rbind(n2, sort(sad1$r(length(n)), decreasing = T))
  }
  n2 <- round(colMeans(n2))
  
  df2$estRank <- 1:length(n)
  df2$estAbund <- n2
  
  df2$dif <- df2$abund[order(df2$rank)] - df2$estAbund[order(df2$estRank)]
  
  if (outputtype == 'plot'){
    p2 <- ggplot(df2, aes(x=rank, y=abund )) +
      geom_line(aes(x=estRank, y=estAbund), color='mediumseagreen', linewidth=1.1) +
      geom_point() + 
      scale_y_continuous(trans='log10') + 
      xlab('Rank') + 
      ylab('Abundance')
    print(p2)
  }
  
  if (outputtype == 'fit'){
    return(sum(df2$dif^2))
  }
  
  if (outputtype == 'estimates'){
    return(df2$estAbund)
  }
  
  if (outputtype == 'differences'){
    return(df2$dif)
  }
}
