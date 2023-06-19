create_SAD <- function(tau, P_m1, P_m2, S_meta, X, nGens1, nGens2)
{
  # Function to produce species rank abundance distributions, starting with
  # a single individual of species 1 in the local community and S_meta species
  # in the metacommunity. First, all X cells are filled, either by offspring 
  # from the metacommunity or from the local community. At full capacity, 
  # the probability of dispersal from the metacommunity equals P_m1. This 
  # value is kept constant for the next nGens1 generations (with tau random
  # replacements per generation). After nGens1 generations, the dispersal 
  # probability changes to P_m2, which is used for nGens2 generations. 
  # Every generation, the rank abundance distribution is recorded. 
  # The output of this function is a dataframe that contains the following 
  # columns: generation, abundance, species, and dispersal probability (Pm). 
  
  # Parameters:
  # tau: number of replacements per generation
  # P_m1: Initial probability of dispersal from metacommunity 
  # P_m2: Pm after fragmentation
  # S_meta: number of species in the metacommunity
  # X: number of cells in the local community that can each hold one individual
  # nGens1: number of generations with the initial Pm
  # nGens2: number of generations with the second Pm
  
  i  <- 1
  for (j in 2:X){
    if (runif(1, 0, 1) <= P_m1^(length(i)/X)){
      i <- c(i, sample(1:S_meta, 1))
    } else {
      i <- c(i, sample(i, 1))
    }
  }
  
  n  <- tapply(i, as.factor(i), length)
  s  <- sort(unique(i))
  df <- data.frame(generation = 0, 
                   abundance = n, 
                   species = s,
                   Pm = P_m1)
  
  for (i_gen in 1:nGens1){
    for (j in 1:tau) {
      ix <- sample(1:X, 1)
      if (runif(1, 0, 1) <= P_m1){
        i[ix] <- sample(1:S_meta, 1)
      } else {
        i[ix] <-sample(i, 1)
      }
    }
    n <- tapply(i, as.factor(i), length)
    s <- sort(unique(i))
    df2 <- data.frame(generation = i_gen, 
                      abundance = n, 
                      species = s,
                      Pm = P_m1)
    df <- rbind(df, df2)
  }
  
  for (i_gen in (nGens1+1):(nGens1+nGens2)){
    for (j in 1:tau) {
      ix <- sample(1:1000, 1)
      if (runif(1, 0, 1) <= P_m2){
        i[ix] <- sample(1:S_meta, 1)
      } else {
        i[ix] <-sample(i, 1)
      }
    }
    n <- tapply(i, as.factor(i), length)
    s <- sort(unique(i))
    df2 <- data.frame(generation = i_gen, 
                      abundance = n, 
                      species = s,
                      Pm = P_m2)
    df <- rbind(df, df2)
  }
  return(df)
}
