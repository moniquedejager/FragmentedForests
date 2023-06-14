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

library(ggplot2)
library(meteR) # rominger and merow, 2016
library(ggpubr)

# function to plot estimated rank abundance distribution:
plot_est_RAI <- function(n) {
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

  p2 <- ggplot(df2, aes(x=rank, y=abund )) +
    geom_line(aes(x=estRank, y=estAbund), color='mediumseagreen', linewidth=1.1) +
    geom_point() + 
    scale_y_continuous(trans='log10') + 
    xlab('Rank') + 
    ylab('Relative abundance')
  p2
  
  df2$dif <- df2$abund[order(df2$rank)] - df2$estAbund[order(df2$estRank)]
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
  
  #print(ggarrange(p2, p3, nrow=1, ncol=2))
  print(p2)
  return(mean(df2$dif^2))
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
i <- 1:100 # 100 species, with 1 individual per species
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

# the speciation rate (which is actually the probability that seeds enter 
# the community from outside) should depend on the number of seeds produced
# within the community

# keep track of changes in a community, from the point where there are 
# 1000 individuals at the site. From there on, a random individual gets
# replaced by a random seedling. 

# to do: spec_rate omschrijven naar een dispersal parameter
# goed opschrijven wat je allemaal gedaan hebt, inclusief alle assumpties
# Verschillende beginwaardes proberen van dispersal parameter, 
# in combinatie met verschillende veranderingen in deze parameter 
# (bijvoorbeeld toenemende fragmentatie, of verschillende fragmentatieniveaus)

spec_rate <- 50
i         <- 1
for (j in 2:1000){
  if (runif(1, 0, 1) <= spec_rate/length(i)){
    i <- c(i, max(i) + 1)
  } else {
    i <- c(i, sample(i, 1))
  }
}

n <- tapply(i, as.factor(i), length)
n <- sort(n, decreasing = T)
n <- c(n, rep(0, 1000 - length(n)))
df <- data.frame(generation = 0, abundance = n, rank = 1:length(n))

for (i_gen in 1:100){
  for (j in 1:1000) {
    ix <- sample(1:1000, 1)
    if (runif(1, 0, 1) <= spec_rate/length(i)){
      i[ix] <- max(i) + 1
    } else {
      i[ix] <-sample(i, 1)
    }
  }
  n <- tapply(i, as.factor(i), length)
  n <- sort(n, decreasing = T)
  n <- c(n, rep(0, 1000 - length(n)))
  df2 <- data.frame(generation = i_gen, abundance = n, rank = 1:length(n))
  df <- rbind(df, df2)
}

ggplot(df, aes(x=rank, y=abundance, color=generation)) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') + 
  scale_color_continuous(type='viridis') + 
  xlab('Rank') + 
  ylab('Abundance')
  
df3 <- data.frame(generation = 0, 
                  fit = plot_est_RAI(df$abundance[(df$generation == 0)
                                                  &(df$abundance > 0)]))
for (j in 1:100){
  df4 <- data.frame(generation = j, 
                    fit = plot_est_RAI(df$abundance[(df$generation == j)
                                                    &(df$abundance > 0)]))
  df3 <- rbind(df3, df4)
}

ggplot(df3, aes(x=generation, y=fit)) + geom_point()

# what happens if we run it another 100 generations, but without 
# input from outside (no connections with other areas)?
df$spec_rate <- spec_rate
spec_rate <- 10

for (i_gen in 1:100){
  for (j in 1:1000) {
    ix <- sample(1:1000, 1)
    if (runif(1, 0, 1) <= spec_rate/length(i)){
      i[ix] <- max(i) + 1
    } else {
      i[ix] <-sample(i, 1)
    }
  }
  n <- tapply(i, as.factor(i), length)
  n <- sort(n, decreasing = T)
  n <- c(n, rep(0, 1000 - length(n)))
  df2 <- data.frame(generation = i_gen, 
                    abundance = n, 
                    rank = 1:length(n),
                    spec_rate = spec_rate)
  df <- rbind(df, df2)
}

df3$type <- 'withMigration'
for (j in 1:100){
  df4 <- data.frame(generation = j, 
                    fit = plot_est_RAI(df$abundance[(df$generation == j)
                                                    &(df$abundance > 0)&
                                                      (df$spec_rate == spec_rate)]),
                    type = 'withoutMigration')
  df3 <- rbind(df3, df4)
}

ggplot(df3, aes(x=generation, y=fit)) + geom_point() + 
  facet_wrap(vars(type))

