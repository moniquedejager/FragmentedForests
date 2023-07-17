# create a large database with the abundances per rank 
# of the first 20 ranks (no use of mete here) of generation 200 
# per combination of dispersal probability P_m2, local community size X, and 
# metacommunity size S_meta. This database will be used to explain
# P_m2 from the abundances using a deep neural network approach.  

create_SAD_df <- function(ix){
  
  source('./R-scripts/create_SAD.R')
  source('./R-scripts/est_RAI.R')
  library(meteR) # rominger and merow, 2016

  # create a full orthagonal sampling of the parameters P_m2, X, and S_meta:
  params <- expand.grid(P_m2 = 0:20/1000, 
                        X = 1:20*100, 
                        S_meta = 1:20*100)
  
  # filenames to write the output to: 
  filenames <- c(paste('./results/SADs for deep neural network/SAD_first20_df', 
                    ix,'.txt', sep=''),
                 paste('./results/SADs for deep neural network/SAD_binned_df', 
                       ix,'.txt', sep=''),
                 paste('./results/SADs for deep neural network/SAD_diff_df', 
                       ix,'.txt', sep=''))
  
  for (i in 1:length(params$P_m2)) {
    
    # Parameters:
    P_m1   <- 0.05             # Initial probability of dispersal from metacommunity 
    P_m2   <- params$P_m2[i]   # Pm after fragmentation
    S_meta <- params$S_meta[i] # number of species in the metacommunity
    X      <- params$X[i]      # number of cells in the local community 
    tau    <- X                # number of replacements per generation
    nGens1 <- 100              # number of generations with the initial Pm
    nGens2 <- 100              # number of generations with the second Pm

    df  <- create_SAD(tau, P_m1, P_m2, S_meta, X, nGens1, nGens2)
    df2 <- data.frame(rank = 1:length(df$generation[df$generation == 200]),
                    abundance = sort(df$abundance[df$generation == 200],
                                     decreasing=T))
    df3 <- data.frame(rank = 1:20,
                      abundance = c(df2$abundance[df2$rank <= 20],
                                    rep(0, 
                                        20 - length(df2$abundance
                                                    [df2$rank <= 20]))))
    
    # if there are less than 20 species, we do not calculate METE's SAD, and
    # we don't make 20 bins. 
    record <- 1
    if (max(df2$rank) >= 20){
      df2$difference <- est_RAI(df$abundance[(df$generation==200)], 
                                'differences')
      
      group <- ceiling(df2$rank / max(df2$rank) * 20)
      df3$binned_abundance <- tapply(df2$abundance, group, mean)
      df3$binned_difference <- tapply(df2$difference, group, mean)
      record <- 3
    } 
    
    for (irec in 1:record){
      m <- as.matrix(df3)
      df4 <- data.frame(r1 = m[1,irec+1],
                        r2 = m[2,irec+1],
                        r3 = m[3,irec+1],
                        r4 = m[4,irec+1],
                        r5 = m[5,irec+1],
                        r6 = m[6,irec+1],
                        r7 = m[7,irec+1],
                        r8 = m[8,irec+1],
                        r9 = m[9,irec+1],
                        r10 = m[10,irec+1],
                        r11 = m[11,irec+1],
                        r12 = m[12,irec+1],
                        r13 = m[13,irec+1],
                        r14 = m[14,irec+1],
                        r15 = m[15,irec+1],
                        r16 = m[16,irec+1],
                        r17 = m[17,irec+1],
                        r18 = m[18,irec+1],
                        r19 = m[19,irec+1],
                        r20 = m[20,irec+1],
                        P_m2 = P_m2,
                        X = X,
                        S_meta = S_meta)
      
      filename <- filenames[irec]
      if (file.exists(filename)) {
        write.table(df4, 
                    filename,
                    append=T, col.names = F, row.names = F)
      } else {
        write.table(df4, 
                    filename,
                    append=T, col.names = T, row.names = F)
      }
    }
  }
}

library(parallel)
nCores  <- 10
cl      <- makeCluster(nCores)
results <- parSapply(cl, 1:10, create_SAD_df)
stopCluster(cl)
