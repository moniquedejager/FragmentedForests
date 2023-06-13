# I am curious where the specific Rank abundance distribution that follows 
# the maximum entropy estimation comes from, and why we think that this 
# corresponds to a healthy ecosystem. 

# if we start with an equal amount of individuals per species and assume 
# neutrality, what will happen if we randomly draw new offspring from the 
# growing pool of individuals? 

# The effect of the first random selection will be much larger if there are 
# fewer individuals per species to start with... 

# Should we also include speciation?
# And death of individuals?
# what happens if we would include space? 

# function to plot estimated rank abundance distribution:
plot_est_RAI <- function(n) {
  df2 <- data.frame(spp = 1:length(n), 
                    abund = n,
                    rank = length(n) + 1 - rank(n + runif(length(n), 0, 0.1))) 
  
  esf1 <- meteESF(spp=df2$spp, abund=df2$abund)
  sad1 <- sad(esf1)
  sad1 # print function returns useful summary
  
  n2 <- sort(sad1$r(length(n)), decreasing = T)
  for (i in 2:1000) {
    n2 <- rbind(n2, sort(sad1$r(length(n)), decreasing = T))
  }
  n2 <- round(colMeans(n2))
  #lines(n2, col='blue')
  
  df2$estRank <- 1:length(n)
  df2$estAbund <- n2
  df2$relAbund <- df2$abund / sum(df2$abund)
  df2$relEstAbund <- df2$estAbund / sum(df2$estAbund)
  
  p2 <- ggplot(df2, aes(x=rank, y=relAbund )) +
    geom_line(aes(x=estRank, y=relEstAbund), color='mediumseagreen', linewidth=1.1) +
    geom_point() + 
    scale_y_continuous(trans='log10') + 
    xlab('Rank') + 
    ylab('Relative abundance')
  
  df2$dif <- df2$relAbund[order(df2$rank)] - df2$relEstAbund[order(df2$estRank)]
  df2$direction <- df2$dif > 0
  p3 <- ggplot(df2, aes(x=rank[order(rank)], 
                        y=abs(dif),
                        color=direction)) + 
    geom_point() + 
    #scale_y_continuous(trans='log10') + 
    xlab('Rank') + 
    ylab(expression(paste(Delta, 'Relative aAbundance', sep=' '))) + 
    scale_y_continuous(trans='log10') + 
    theme(legend.position = 'none')
  
  windows(height= 5, width=10)
  print(ggarrange(p2, p3, nrow=1, ncol=2))
  return(mean(abs(1/df2$dif)))
}

# Model without space and speciation:
i <- 1:200 # 200 species, with 1 individual per species
for (j in 1:10000) {
  i <- c(i, sample(i, 1))
}
n <- tapply(i, as.factor(i), length)
plot_est_RAI(n)

# Model without space, with speciation:
spec_rate <- 0.1 # speciation rate
i <- 1:100 # 200 species, with 1 individual per species
for (j in 1:10000) {
  if (runif(1, 0, 1) <= spec_rate){
    i <- c(i, max(i) + 1)
  } else {
    i <- c(i, sample(i, 1))
  }
}
n <- tapply(i, as.factor(i), length)
plot_est_RAI(n)

# the higher the speciation rate, the better the fit! 
# the number of species to start with, also changes the fit...
# it seems we need to examine the combination of these two parameters... 
# as well as the number of iterations... 

