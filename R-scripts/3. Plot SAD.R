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

# Parameters:
tau    <- 1000 # number of replacements per generation
P_m1   <- 0.05 # Initial probability of dispersal from metacommunity 
P_m2   <- 0.01 # Pm after fragmentation
S_meta <- 1000 # number of species in the metacommunity
X      <- 1000 # number of cells in the local community that 
               #can each hold one individual
nGens1 <- 100  # number of generations with the initial Pm
nGens2 <- 300  # number of generations with the second Pm

library(ggplot2)
library(meteR) # rominger and merow, 2016
library(ggpubr)
source('./R-scripts/est_RAI.R')
source('./R-scripts/create_SAD.R')

df <- create_SAD(tau, P_m1, P_m2, S_meta, X, nGens1, nGens2)

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

# to do: 
# goed opschrijven wat je allemaal gedaan hebt, inclusief alle assumpties
# Verschillende beginwaardes proberen van dispersal parameter, 
# in combinatie met verschillende veranderingen in deze parameter 
# (bijvoorbeeld toenemende fragmentatie, of verschillende fragmentatieniveaus)


df3 <- data.frame(generation = 0, 
                  fit = est_RAI(df$abundance[(df$generation == 0)], 'fit'),
                  n_species = length(df$abundance[(df$generation == 0)]),
                  type = P_m1)
for (j in 1:(nGens1 + nGens2)){
  df4 <- data.frame(generation = j, 
                    fit = est_RAI(df$abundance[(df$generation == j)], 'fit'),
                    n_species = length(df$abundance[(df$generation == j)]),
                    type = df$Pm[(df$generation == j)])
  df3 <- rbind(df3, df4)
}

windows(height=4, width=4)
ggplot(df3, aes(x=n_species, y=log(1/fit), color=type)) +
  geom_point() + 
  xlab('Number of species') + 
  ylab('Fit') + 
  theme(legend.position = 'none')

df4a <- data.frame(generation=df3$generation,
                  y = log(1/df3$fit),
                  type = df3$type, 
                  type2 = 'GOF')
df4b <- data.frame(generation=df3$generation,
                   y = df3$n_species,
                   type = df3$type, 
                   type2 = 'Number of Species')
df4 <- rbind(df4a, df4b)

windows(height=5, width=8)
ggplot(df4, aes(x=generation, y=y, color=type)) + 
  geom_line(linewidth=1.2) + 
  xlab('Generation') + 
  ylab('') + 
  facet_grid(rows=vars(type2), scales='free_y', 
             switch = "y") + 
  theme(legend.position = 'none') + 
  theme(strip.placement = "outside")

# plot SAD's of generations 100 and 200 to show the difference:
# with predictions of METE:

df5 <- data.frame(rank = c(1:length(df$generation[df$generation == 100]),
                           1:length(df$generation[df$generation == 200])),
                  abundance = c(sort(df$abundance[df$generation == 100],
                                     decreasing=T),
                                sort(df$abundance[df$generation == 200], 
                                     decreasing = T)),
                  estAbundance = c(est_RAI(df$abundance[(df$generation==100)], 
                                           'estimates'),
                                   est_RAI(df$abundance[(df$generation==200)], 
                                             'estimates')),
                  generation = paste('Generation', 
                                     df$generation[df$generation %in% c(100, 200)]),
                  difference = c(est_RAI(df$abundance[(df$generation==100)], 
                                         'differences'),
                                 est_RAI(df$abundance[(df$generation==200)], 
                                         'differences')))

windows(height=4, width=8)  
ggplot(df5, 
       aes(x=rank, y=abundance)) + 
  geom_point(size=2) + 
  geom_line(linewidth=1.2) + 
  geom_line(aes(x=rank, y=estAbundance), 
            linetype='dashed',
            color='red',
            size=1.2)+
  geom_point(aes(x=rank, y=estAbundance), 
             color='red',
             shape=1, size=2) + 
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(legend.position = 'none') + 
  xlab('Rank') + 
  ylab('Abundance') + 
  facet_wrap(vars(as.factor(generation)))

windows(height=6, width=6) 
p1 <- ggplot(df, 
       aes(x=generation, y=abundance, color=as.factor(species))) + 
  geom_line(linewidth=1.2) + 
  theme(legend.position = 'none') + 
  xlab('Generation') + 
  ylab('Abundance') + 
  xlim(c(50, 100)) + 
  ylim(c(0, 100))
p2 <- ggplot(df, 
       aes(x=generation, y=abundance, color=as.factor(species))) + 
  geom_line(linewidth=1.2) + 
  theme(legend.position = 'none') + 
  xlab('Generation') + 
  ylab('Abundance') + 
  xlim(c(200, 250))
ggarrange(p1+ rremove("xlab"), p2, nrow=2, ncol=1)

# is there a pattern in the deviation in abundance per rank, per value of Pm?
#ggplot(df5, aes(x=rank, y=difference)) + 
#  geom_point(size=2) + 
#  geom_line(linewidth=1.2) + 
#  #scale_y_continuous(trans='log10') +
#  #scale_x_continuous(trans='log10') +
#  theme(legend.position = 'none') + 
#  xlab('Rank') + 
#  ylab('Difference') + 
#  facet_wrap(vars(as.factor(generation)))



