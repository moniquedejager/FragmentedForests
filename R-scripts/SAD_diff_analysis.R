# Analyse patterns in the differences in abudance per rank between simulated
# data and METE estimations

library(ggplot2)

list = list.files(path = './results/SAD differences per Pm' ,
                 full.names=TRUE,recursive=TRUE)

df <- read.table(list[1], header =T)
lx <- unlist(
  strsplit(
    unlist(
      strsplit(
        unlist(
          strsplit(list[1], '_')), '=')), '.txt'))
df$Pm1 <- as.numeric(lx[4])
df$X   <- as.numeric(lx[8])

for (i in 2:length(list)){
  dfa <- read.table(list[i], header =T)
  lx <- unlist(
    strsplit(
      unlist(
        strsplit(
          unlist(
            strsplit(list[i], '_')), '=')), '.txt'))
  dfa$Pm1 <- as.numeric(lx[4])
  dfa$X   <- as.numeric(lx[8])
  
  df <- rbind(df, dfa)
}

ggplot(df[df$generation == 200,], 
       aes(x=rank, 
           y=difference/abundance, 
           color=as.factor(X))) +
  #geom_point() + 
  geom_line() + 
  scale_x_continuous(trans='log10') + 
  facet_wrap(vars(Pm)) 
  xlim(c(0, 20))  
  theme(legend.position = 'none')

windows(height=10, width=15)
ggplot(df[df$Pm == 0,], 
       aes(x=rank, 
           y=abs(difference), 
           color=as.factor(generation))) +
  geom_point() + 
  geom_line() + 
  facet_wrap(vars(generation)) +
  xlim(c(0, 20)) + 
  theme(legend.position = 'none')

windows(height=10, width=15)
ggplot(df[df$Pm == 0.03,], 
       aes(x=rank, 
           y=abs(difference), 
           color=as.factor(generation))) +
  geom_point() + 
  geom_line() + 
  facet_wrap(vars(generation)) +
  xlim(c(0, 20)) + 
  theme(legend.position = 'none')


df$diff_transformed <- log(abs(df$difference)+1) - 1
df$diff_transformed[df$difference < 0] <- 
  df$diff_transformed[df$difference < 0] * -1

ggplot(df[df$generation == 175,], 
       aes(x=rank, 
           y=Pm, 
           fill=diff_transformed)) +
  geom_tile() + 
  scale_fill_continuous(type = 'viridis') + 
  xlim(c(0, 20))

new = T
for (iPm in unique(df$Pm)){
  for (iGen in unique(df$generation)){
    df2 <- df[(df$Pm == iPm)&(df$generation == iGen),]
    df3a <- df2[2:length(df2$rank),]
    df3b <- df2[1:(length(df2$rank)-1),]
    df3b$ddif <- df3a$difference - df3b$difference
    if (new){ 
      df3 <- df3b 
      new <- F
    }else{
      df3 <- rbind(df3, df3b) 
    }
  }
}

ggplot(df3[df3$generation %in% 175:185,], 
       aes(x=rank, 
           y=abs(ddif), 
           color=as.factor(generation))) +
  geom_point() + 
  geom_line() + 
  facet_wrap(vars(Pm)) +
  xlim(c(0, 20)) + 
    theme(legend.position = 'none')


ggplot(df3[df3$generation > 0,], 
       aes(x=rank, 
           y=log(abs(ddif)+1), 
           color=as.factor(generation))) +
  geom_point() + 
  geom_line() + 
  facet_wrap(vars(Pm)) +
  xlim(c(0, 100))

df3$y <- abs(df3$ddif)
mod <- glm(y~generation + rank + Pm + 
             generation:rank + rank:Pm,
           data=df3[(df3$rank <= 200),], family = 'poisson')
summary(mod)
 
mod0 <- glm(y~1, data=df3[(df3$rank <= 200),],family='poisson')
# McFadden’s R-Squared = 1 – (log likelihoodmodel / log likelihoodnull)
r2 <- 1 - (logLik(mod) / logLik(mod0))


