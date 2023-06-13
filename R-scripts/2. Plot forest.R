# 2. Plot the forest
library(ggplot2)
library(meteR) # rominger and merow, 2016
library(ggpubr)

# function to plot METEs estimate of the SAD ####
plotSAD <- function(df){
  p1 <- ggplot(df, aes(x=x, y=y, color=species)) + 
    geom_point() + 
    xlab('X-coordinate (m)') + 
    ylab('Y-coordinate (m)') + 
    theme_minimal() +
    theme(legend.position="none")
  
  n <- tapply(df$species, df$species, length)
  df2 <- data.frame(spp = unique(df$species), 
                    abund = n,
                    rank = length(n) + 1 - rank(n + runif(length(n), 0, 0.1))) 
  
  esf1 <- meteESF(spp=df2$spp, abund=df2$abund)
  sad1 <- sad(esf1)
  sad1 # print function returns useful summary
  
  #plot(sad1, ptype='rad',log='y') 
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
  
  p3 <- ggplot(df2, aes(x=rank[order(rank)], 
                        y=(relAbund[order(rank)] - relEstAbund[order(estRank)]))) + 
    geom_point() + 
    #scale_y_continuous(trans='log10') + 
    xlab('Rank') + 
    ylab(expression(paste(Delta, 'Relative aAbundance', sep=' ')))
  
  windows(height= 5, width=15)
  print(ggarrange(p1, p2, p3, nrow=1, ncol=3))
}

for (iForest in 1:10) {
  # dataframes are too large for the meteR package to handle, 
  # therefore, we split the area into 1 km^2 areas:
  
  df <- read.table(paste('./results/created forest/Forest', iForest, '.txt'))
  df <- data.frame(x=df$V1, y=df$V2, species=df$V4)
  
  ix <- 0:10 * 1000
  for (i in 1:10) {
    df2 <- df[(df$y > ix[i])&(df$y <= ix[i+1]),]
    plotSAD(df2)
  }
}

# RAI of all data together:
df <- read.table(paste('./results/created forest/Forest', 1, '.txt'))
for (i in 2:10){
  df2 <- read.table(paste('./results/created forest/Forest', i, '.txt'))
  df  <- rbind(df, df2)
}
df    <- data.frame(x=df$V1, y=df$V2, species=df$V4)
n     <- tapply(df$x, df$species, length)
df2   <- data.frame(spp = sort(unique(df$species)), 
                  abund = n,
                  rank = length(n) + 1 - rank(n + runif(length(n), 0, 0.1))) 

ggplot(df2, aes(x=rank, y=abund )) +
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  xlab('Rank') + 
  ylab('Relative abundance')

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

ggplot(df2, aes(x=rank, y=abund )) +
  geom_line(aes(x=estRank, y=estAbund), color='mediumseagreen', linewidth=1.1) +
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  xlab('Rank') + 
  ylab('Abundance')

