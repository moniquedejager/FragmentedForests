# Analyse patterns in the differences in abudance per rank between simulated
# data and METE estimations

library(ggplot2)
library(lme4)

df <- read.table(paste('./results/SAD differences per Pm/SAD_dif_', 0, '.txt'), 
                 header =T)
for (i in 1:10/100){
  df <- rbind(df, 
              read.table(paste('./results/SAD differences per Pm/SAD_dif_', i, '.txt'), 
                    header =T))
}

ggplot(df[df$generation == 175,], 
       aes(x=rank, 
           y=abs(difference), 
           color=as.factor(generation))) +
  geom_point() + 
  geom_line() + 
  facet_wrap(vars(Pm)) +
  xlim(c(0, 20))

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
