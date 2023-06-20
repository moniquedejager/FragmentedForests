# create a large database with the differences in abundances per rank between
# simulated data and the METE fit to this data, per generation and dispersal
# probability P_m2.

# The first part of the simulation will be kept constant; for the second part, 
# we will vary P_m2. For generations 150 - 200, we will record the difference
# between observed and expected abundance per rank. 

# if we use 10 different values of P_m2, we can run the simulations in parallel. 

SAD_diff_per_Pm <- function(P_m2){
  # Parameters:
  P_m1   <- 0.05 # Initial probability of dispersal from metacommunity 
  #P_m2   <- 0.01 # Pm after fragmentation
  S_meta <- 1000 # number of species in the metacommunity
  X      <- 100 # number of cells in the local community 
  tau    <- X    # number of replacements per generation
  nGens1 <- 100  # number of generations with the initial Pm
  nGens2 <- 100  # number of generations with the second Pm
  
  library(ggplot2)
  library(meteR) # rominger and merow, 2016
  library(ggpubr)
  source('./R-scripts/est_RAI.R')
  source('./R-scripts/create_SAD.R')
  
  df <- create_SAD(tau, P_m1, P_m2, S_meta, X, nGens1, nGens2)
  
  df2 <- data.frame(rank = 1:length(df$generation[df$generation == 150]),
                    generation = df$generation[df$generation == 150],
                    difference = est_RAI(df$abundance[(df$generation==150)], 
                                         'differences'),
                    abundance = sort(df$abundance[df$generation == 150],
                                       decreasing=T))
  
  for (gen in 151:200){
    df2a <- data.frame(rank = 1:length(df$generation[df$generation == gen]),
                       generation = df$generation[df$generation == gen],
                       difference = est_RAI(df$abundance[(df$generation==gen)], 
                                            'differences'),
                       abundance = sort(df$abundance[df$generation == gen],
                                        decreasing=T))
    df2 <- rbind(df2, df2a)
  }
  
  df2$Pm <- P_m2
  
  filename <- paste('./results/SAD differences per Pm/SAD_dif_Pm1=', P_m1,
                    '_Pm2=', P_m2, 
                    '_X=', X, '.txt', sep='')
  write.table(df2, 
              filename,
              append=F, col.names = T, row.names = F)
}

library(parallel)
nCores  <- 11
cl      <- makeCluster(nCores)
P_m2    <- 0:10/1000
results <- parSapply(cl, P_m2, SAD_diff_per_Pm)
stopCluster(cl)
