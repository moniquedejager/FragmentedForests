# Shiny forest:

# te doen: SAD plots toevoegen...

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
                  value = 100,
                  min = 1,
                  max = 1000),
      
      # br() element to introduce extra vertical spacing ----
      br(),
      
      # Input: Slider for the dispersal probability ----
      sliderInput("Pm",
                  "Dispersal probability:",
                  value = 0.5,
                  min = 0,
                  max = 1),
      
      # br() element to introduce extra vertical spacing ----
      br(),
      
      # Input: Slider for the number of individuals per subcommunity ----
      sliderInput("n_individuals",
                  "Number of individuals per subcommunity:",
                  value = 500,
                  min = 100,
                  max = 1000)
    ),
    
    # Main panel for displaying outputs ----
    mainPanel("Map of the shiny forest",
              plotOutput("map")
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output, session) {
  library(ggplot2)
  
  rv <- reactiveValues(run = F)
  
  observe({
    rv$n     <- input$n 
    rv$Pm    <- input$Pm
    rv$n_ind <- input$n_individuals
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
      
      rv$S_meta <- nrow(rv$species)
      
      # the first individual always comes from the metacommunity:
      # (the number represents the species)
      for (j in 1:rv$n){
        rv$species[1, j] <- sample(1:rv$S_meta, 1)
        rv$nspecies[j] <- length(unique(rv$species[,j])) - 1
      }
    })
  })
  
  autoInvalidate <- reactiveTimer(intervalMs=200, session)
  
  observe({
    autoInvalidate()
    isolate({
      if (rv$run) {
        rv$t <- rv$t + 1
        r    <- runif(rv$n)
        
        # using the cumulative probability distribution, we calculate the 
        # distance from which the new individual should approximately come:
        d <- log(r)/log(rv$Pm)
        
        # per cell, select one of the subcommunities (or the metacommunity) 
        # closest to this distance
        if (sum(rv$species[,1] > 0) == rv$n_ind){
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
              rv$species[sample(1:rv$n_ind, 1),j] <- sample(specs, 1)
            } else {
              # choose an individual from the metacommunity:
              rv$species[sample(1:rv$n_ind, 1), j] <- sample(1:rv$S_meta, 1)
            }
            rv$nspecies[j] <- length(unique(rv$species[,j])) - 1
          }
        } else {
          # per cell, select one of the subcommunities (or the metacommunity) 
          # closest to this distance 
          i <- sum(rv$species[,1] > 0) + 1
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
              rv$species[i,j] <- sample(specs, 1)
            } else {
              # choose an individual from the metacommunity:
              rv$species[i, j] <- sample(1:rv$S_meta, 1)
            }
            rv$nspecies[j] <- length(unique(rv$species[,j])) - 1
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
      labs(caption = paste('t = ', rv$t, sep=''), 
           x='', y='')
  })
}

# Create Shiny app ----
shinyApp(ui, server)
