# R-script to create a spatially explicit representation of a forest, 
# based on rank abundance data provided at specific locations. 
# The modelled area size can be adjusted with the parameters area_length
# and area_width (both in m). Each 1 x 1 m patch can hold up to one tree. 

# We start with an empty space containing area_length x area_width patches.
# At specified locations (specified in the input file 'plot_locations.txt'), 
# sources of tree seeds can be found. These sources contain a number of 
# individuals of each species, that is specified in the input file 
# 'RAI_per_plot.txt'. 
# Each patch can become occupied by a tree with probability (p_occ). 
# The probability that a seed comes from a specific source population
# depends on the relative distance of the source to the patch (i.e. the 
# probability is highest that a seed comes from the nearest source population).
# Every tree in a source population has an equal change of colonizing the
# empty patch. 

# The created forest is saved to 10 files ('forest xx.txt'), with x, y, 
# and species data. 

# to speed up, we need to fill in the forest using parallel runs:
library(parallel)
nCores <- 10

fillForest <- function(ix) {
  # 1. Set the parameters (area size and Probability that a patch becomes 
  # occupied): 
  area_length <- 10000
  area_width <- 10000
  
  # 2. Load in the locations and numbers of individuals per species of 
  # the source populations:
  data_per_plot <- read.table('./data/RAI_per_plot.txt', header=T)
  location_per_plot <-read.table('./data/plot_locations.txt', header=T)
  
  # what is the density at the plots? 
  # the plots are 50 x 100 m each
  p_occ  <- length(data_per_plot$site) / (length(location_per_plot$site) * 50 * 100) 
  
  x1     <- 0:10 * (area_length / 10)
  df     <- expand.grid(x = (x1[ix]+1):x1[ix+1], 
                    y = 1:area_width)
  df     <- data.frame(x = df$x, y = df$y)
  df$occ <- runif(length(df$x), 0, 1)
  df     <- df[df$occ <= p_occ,]
  
  for (i in 1:length(df$x)) {
    # calculate distances to source populations and select a source population:
    dists  <- sqrt((df$x[i] - location_per_plot$x_location)^2 + 
                    (df$y[i] - location_per_plot$y_location)^2)
    prob   <- 0.99^dists
    source <- sample(location_per_plot$site, 1, prob = prob)
    
    # from the source population, select a tree:
    source_ppn <- data_per_plot[data_per_plot$site == source, ]
    
    DATA <- df[i,]
    DATA$species <- source_ppn[sample(1:length(source_ppn$site), 1),]$taxon
    
    write.table(DATA, paste('Forest', ix, '.txt'), 
                append=T, 
                row.names = F, 
                col.names = F)
  }
}


# run in parallel:
cl   = makeCluster(nCores)
results = parSapply(cl, 1:nCores, fillForest)
stopCluster(cl)





