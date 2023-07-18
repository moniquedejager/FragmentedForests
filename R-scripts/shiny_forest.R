# Shiny forest:

# te doen: bij het simuleren van veel subcommunities, moet er parallel gerekend 
# worden!!

library(shiny)
# Launches an app, with the app's source code included
#runExample("05_sliders")

# Lists more prepackaged examples
#runExample()

# Define UI for my Shiny Forest app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Shiny forest"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Control what happens:
      "Controls",
      actionButton("gogobutt","Go"),
      actionButton("stopbutt","Stop"),
      
      # br() element to introduce extra vertical spacing ----
      br(),
      
      # Input: Slider for the number of subpopulations to generate ----
      sliderInput("n",
                  "Number of subpopulations:",
                  value = 25,
                  min = 1,
                  max = 1000),
      
      br(),
      
      # Input: Slider for the dispersal probability ----
      sliderInput("Pm",
                  "Dispersal probability:",
                  value = 0.1,
                  min = 0,
                  max = 1),
      
      br(),
      
      # Input: Slider for the number of individuals per subcommunity ----
      sliderInput("n_individuals",
                  "Number of individuals per subcommunity:",
                  value = 500,
                  min = 100,
                  max = 1000),
      
      br(),
      
      # Input: Slider for the number of species in the metacommunity ----
      sliderInput("S_meta",
                  "Number of species in the metacommunity:",
                  value = 5000,
                  min = 100,
                  max = 10000)
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      "Map of the shiny forest",
      plotOutput("map"),
      
      br(),
      
      "Species' rank abundance distributions",
      plotOutput("SAD"),
      
      br(),
      
      "Species-area relation",
      plotOutput("SpeciesArea")
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output, session) {
  library(ggplot2)
  library(matrixStats)
  
  rv <- reactiveValues(run = F)
  
  observe({
    rv$n     <- input$n 
    rv$Pm    <- input$Pm
    rv$n_ind <- input$n_individuals
    rv$S_meta<- input$S_meta
    rv$t     <- 0
    
    isolate({
      # create a lattice that approximates n cells:
      rv$nx <- ceiling(sqrt(rv$n))
      rv$ny <- floor(sqrt(rv$n))
      
      # calculate the actual number of cells:
      rv$n <- rv$nx * rv$ny
      
      # position the cells in the lattice:
      rv$x <- rep(1:rv$nx, rv$ny)
      rv$y <- sort(rep(1:rv$ny, rv$nx))
      
      # create the meta community around the forest:
      rv$x <- c(rv$x, rep(0, rv$ny), rep(rv$nx+1, rv$ny), 1:rv$nx, 1:rv$nx)
      rv$y <- c(rv$y, 1:rv$ny, 1:rv$ny, rep(0, rv$nx), rep(rv$ny+1, rv$nx))
      
      # each cell can hold a maximum of n_ind individuals.
      # start with 0 individuals in all cells:
      rv$species <- matrix(0, rv$n_ind, rv$n)
      rv$nspecies <- rep(0, rv$n)
      
      # the first individual always comes from the metacommunity:
      # (the number represents the species)
      for (j in 1:rv$n){
        rv$species[1, j] <- sample(1:rv$S_meta, 1)
        rv$nspecies[j] <- length(unique(rv$species[,j])) - 1
      }
      
      # 3 subcommunities to follow in SAD-plot:
      # first one = bottom-left
      # center one
      # last one = top right
      rv$plot_comm <- c(1, round(rv$n/2), rv$n)
      rv$plot_name <- c('Bottom-left', 'Center', 'Top right')
    })
  })
  
  autoInvalidate <- reactiveTimer(intervalMs=200, session)
  
  observe({
    autoInvalidate()
    isolate({
      if (rv$run) {
        rv$species_new <- rv$species
        
        rv$t <- rv$t + 1
        
        # per cell, select one of the subcommunities (or the metacommunity) 
        # closest to this distance
        if (sum(rv$species[,1] > 0) == rv$n_ind){
          
          for (j in 1:rv$n) {
            
            # using the cumulative probability distribution, we calculate the 
            # distance from which the new individual should approximately come,
            # per old individual that will be replaced:
            r    <- runif(rv$n_ind)
            d <- log(r)/log(rv$Pm)
            
            # use matrices to randomly choose one of the cells closest to d: 
            dist      <- sqrt((rv$x[j] - rv$x)^2 + (rv$y[j] - rv$y)^2)
            dist_mat  <- matrix(dist, ncol=rv$n_ind, nrow=length(dist))
            d_mat     <- t(matrix(d, nrow=rv$n_ind, ncol=length(dist)))
            dist_mat2 <- abs(dist_mat - d_mat) + 
              matrix(runif(rv$n_ind*length(dist), 0, 0.001), 
                     ncol=rv$n_ind, nrow=length(dist))
            min_mat <- t(matrix(colMins(dist_mat2), 
                                nrow=rv$n_ind, ncol=length(dist)))
            
            chosen <- matrix(1:length(dist), 
                             ncol=rv$n_ind, 
                             nrow=length(dist))[dist_mat2 - min_mat == 0]
            
            # select individuals from the metacommunity for those individuals
            # with chosen > rv$n:
            rv$species_new[chosen > rv$n, j] <- sample(1:rv$S_meta, 
                                                       sum(chosen > rv$n), 
                                                       replace = T)
            
            # for the others, select a random individual from the chosen 
            # subcommunity:
            chosen_communities <- rv$species[,chosen[chosen <= rv$n]]
            nr <- nrow(chosen_communities)
            nc <- ncol(chosen_communities)
            rand_individuals   <- matrix(runif(nr*nc),nr,nc) * 
              (chosen_communities > 0)
            max_rand <- t(matrix(colMaxs(rand_individuals), nc, nr))
            specs <- chosen_communities[rand_individuals == max_rand]
            rv$species_new[chosen <= rv$n, j] <- specs
            
            spec <- unique(rv$species_new[,j])
            rv$nspecies[j] <- length(spec[spec > 0]) 
          }
          rv$species <- rv$species_new
        } else {
          # per cell, select one of the subcommunities (or the metacommunity) 
          # closest to this distance
          for (i in 2:rv$n_ind) {
            rv$species_new <- rv$species
            r    <- runif(rv$n)
            d <- log(r)/log(rv$Pm)
            
            for (j in 1:rv$n) {
              dist <- sqrt((rv$x[j] - rv$x)^2 + (rv$y[j] - rv$y)^2)
              closest <- (1:length(rv$x))[abs(d[j] - dist) == min(abs(d[j] - dist))]
              chosen <- sample(closest, 1)
              
              # does the chosen subcommunity exist, or is it part of the 
              # meta community?
              if (chosen <= rv$n) {
                # it is a subcommunity, randomly select an individual to reproduce
                specs <- rv$species[,chosen]
                specs <- rep(specs[specs > 0], 2)
                rv$species_new[i,j] <- sample(specs, 1)
              } else {
                # choose an individual from the metacommunity:
                rv$species_new[i, j] <- sample(1:rv$S_meta, 1)
              }
              rv$nspecies[j] <- length(unique(rv$species[,j])) - 1
            }
            rv$species <- rv$species_new
          }
        }
      }
    })
  })
  
  observeEvent(input$gogobutt, { isolate({ rv$run=T  }) })
  observeEvent(input$stopbutt, { isolate({ rv$run=F  }) })
  
  # Generate a map of the data ----
  output$map <- renderPlot({
    # Show the number of different species per cell:
    df <- data.frame(x = rv$x[1:rv$n], y = rv$y[1:rv$n], nspecies = rv$nspecies)
    ggplot(df, aes(x=x, y=y, fill=nspecies)) + 
      geom_raster() + 
      scale_fill_gradient(low = "white", high = "darkgreen", 
                          name='Number of species') + 
      theme_void() + 
      theme(legend.position = 'right') + 
      labs(caption = paste('gen = ', rv$t, sep=''), 
           x='', y='')
  })
  
  output$SAD <- renderPlot({
    # first, create a rank abundance distribution of all subcommunities 
    # combined:
    spec <- as.vector(rv$species)
    spec      <- spec[spec > 0]
    abundance <- sort(tapply(spec, as.factor(spec), length), decreasing = T)
    df <- data.frame(Rank = 1:length(abundance), 
                     Abundance = abundance, 
                     Subcommunity = 'All', 
                     Type = 'All')
    
    # rank abundance plot of four random subcommunities:
    for (i in 1:3){
      spec      <- rv$species[,rv$plot_comm[i]]
      spec      <- spec[spec > 0]
      abundance <- sort(tapply(spec, as.factor(spec), length), decreasing = T)
      df2 <- data.frame(Rank = 1:length(abundance), 
                        Abundance = abundance, 
                        Subcommunity = rv$plot_name[i], 
                        Type = 'Subcommunities')
      df <- rbind(df, df2)
    }
    if (length(df2$Rank) > 1) {
      ggplot(df, aes(x=Rank, y=Abundance, color=Subcommunity)) + 
        geom_line() + 
        facet_wrap(vars(Type), ncol=2, scales='free') + 
        scale_y_continuous(trans='log10') + 
        scale_color_discrete(type='viridis', name='Subcommunity')
    }
  })
  
  output$SpeciesArea <- renderPlot({
    # create a species-area plot:
    # number of species per area size
    df <- data.frame(area_size = 1, nspecies = rv$nspecies) 
    
    for (i in 2:10){
      x2 <- round(rv$x[1:rv$n]/i)
      y2 <- round(rv$y[1:rv$n]/i)
      
      areas <- paste(x2, y2, sep='-')
      area_size <- tapply(areas, areas, length)
      
      nspecies <- vector(length=0)
      for (j in sort(unique(areas))){
        spec <- rv$species[,areas == j]
        nspecies <- c(nspecies, length(unique(spec[spec > 0])))
      }
      df2 <- data.frame(area_size = area_size,
                        nspecies = nspecies)
      df <- rbind(df, df2)
    }
    
    ggplot(df, aes(x=area_size, y=nspecies)) + 
      geom_point() + 
      scale_x_continuous(trans='log10') + 
      scale_y_continuous(trans='log10') + 
      xlab('Area size (number of cells)') + 
      ylab('Number of species')
  })
}

# Create Shiny app ----
shinyApp(ui, server)
