# deep learning neural network to predict dispersal capacity from SAD

# see 'keras-numerixY-YT_example.R', and
# see youtube: 
# https://www.youtube.com/watch?v=SrQw_fWo4lw&list=PL34t5iLfZddtC6LqEfalIBhQGSZX77bOn&index=5

# Install devtools
#install.packages('devtools')

# Install tensorflow:
#devtools::install_github('rstudio/tensorflow')

# install reticulate:
#install.packages('reticulate')
#install_keras()

# Libraries
library(reticulate)
use_condaenv("r-reticulate")
library(keras)
library(mlbench) 
library(dplyr)
library(magrittr)
library(neuralnet)
library(tensorflow)

# Data
ix <- 1
filename <- paste('./results/SADs for deep neural network/SAD_diff_df', 
                  ix,'.txt', sep='')
data <- read.table(filename, header=T)
for (ix in (2:10)){
  filename <- paste('./results/SADs for deep neural network/SAD_diff_df', 
                    ix,'.txt', sep='')
  data <- rbind(data, read.table(filename, header=T))
}
str(data)

# because P_m2 is restricted between 0 and 1, 
# we need to logit-transform the data...
# P_m2 = 0 cannot be logit-transformed, so we remove the data corresponding 
# with this value:
data <- data[data$P_m2 > 0,]
data$Pm <- log(data$P_m2 / (1 - data$P_m2))

# examine the data: what combinations of P_m2, X, and S_meta were included?
tab <- table(data$X, data$S_meta, data$P_m2)
tab[,,1]
# X > 1000 is always in... 
# -> for small populations (fewer than 1,100 individuals), not all parameter 
# combinations are represented! 

# remove all data from X <= 1000:
#data <- data[data$X > 1000,]

# for now, select only the data with X = 1500:
#data <- data[(data$S_meta == 1500),]

# Matrix
data <- as.matrix(data)
dimnames(data) <- NULL

# can we see differences in the patterns ourselves?
par(mfrow=c(3,3))
plot(as.vector(data[1, 1:20]), main=paste(data[1,21], data[1, 22], data[1, 23], sep='; '))
plot(as.vector(data[3, 1:20]), main=paste(data[3,21], data[3, 22], data[3, 23], sep='; '))
plot(as.vector(data[5, 1:20]), main=paste(data[5,21], data[5, 22], data[5, 23], sep='; '))
plot(as.vector(data[7, 1:20]), main=paste(data[7,21], data[7, 22], data[7, 23], sep='; '))
plot(as.vector(data[9, 1:20]), main=paste(data[9,21], data[9, 22], data[9, 23], sep='; '))
plot(as.vector(data[11, 1:20]), main=paste(data[11,21], data[11, 22], data[11, 23], sep='; '))
plot(as.vector(data[13, 1:20]), main=paste(data[13,21], data[13, 22], data[13, 23], sep='; '))
plot(as.vector(data[15, 1:20]), main=paste(data[15,21], data[15, 22], data[15, 23], sep='; '))
plot(as.vector(data[20, 1:20]), main=paste(data[20,21], data[20, 22], data[20, 23], sep='; '))

# Normalize
data2 <- data
m <- colMeans(data2)
s <- apply(data2, 2, sd)
data2 <- scale(data2, center = m, scale = s)

df <- data.frame(SSQ = rowSums((data2[,1:20])^2),
                 Pm = data2[,24],
                 X  = data2[,22],
                 Smeta = data2[,23])

library(ggplot2)

ggplot(df, aes(x=Pm, y=X, fill=SSQ)) + 
  geom_raster() + 
  facet_wrap(vars(Smeta)) + 
  scale_fill_continuous(type='viridis', trans='log10')

ggplot(df, aes(x=(Pm), y=SSQ, color=as.factor(X))) +
  geom_line() + 
  facet_wrap(vars(Smeta)) + 
  scale_y_continuous(trans='log10') + 
  #scale_x_continuous(trans='log10') + 
  scale_color_discrete(type='viridis')

df$y <- log(df$SSQ)
mod1 <- lm(y~Pm + I(Pm^2) + X + I(X^2) + Smeta + I(Smeta^2) + 
             Pm:X + X:Smeta + Pm:Smeta + Pm:X:Smeta + 
             I(X^2):I(Smeta^2) + I(Pm^2):I(X^2):I(Smeta^2) + 
             I(Pm^3) + I(X^3) + I(Smeta^3) + 
             I(Pm^3):I(Smeta^3), data=df)
summary(mod1)
mod2 <- lm(y~Pm + I(Pm^2) + X + I(X^2) + Smeta + I(Smeta^2) + 
             Pm:X + X:Smeta + Pm:Smeta + Pm:X:Smeta + 
             I(X^2):I(Smeta^2) + I(Pm^2):I(X^2):I(Smeta^2) + 
             I(Pm^3) + I(X^3) + I(Smeta^3) + 
             I(Pm^3):I(Smeta^3) + I(X^3):I(Smeta^3), data=df)
summary(mod2)
anova(mod1, mod2)


# Partition
set.seed(1234)
ind <- sample(2, nrow(data), replace = T, prob = c(.7, .3))
training <- abs(data[ind==1,c(1:20, 22, 23)])
test <- abs(data[ind==2, c(1:20, 22, 23)])
trainingtarget <- data[ind==1, c(24)]
testtarget <- data[ind==2, c(24)]
#trainingtarget <- df$SSQ[ind==1]
#testtarget <- df$SSQ[ind==2]

# Normalize
m <- colMeans(training)
s <- apply(training, 2, sd)
training <- scale(training, center = m, scale = s)
test <- scale(test, center = m, scale = s)

m <- mean(trainingtarget)
s <- sd(trainingtarget)
trainingtarget <- scale(trainingtarget, center = m, scale = s)
testtarget <- scale(testtarget, center = m, scale = s)

# Create Model
model <- keras_model_sequential()
model %>% 
  layer_dense(units = 5, activation = 'relu', input_shape = c(22)) %>%
  layer_dense(units = 1)

# Compile
model %>% compile(loss = 'mse',
                  optimizer = 'rmsprop',
                  metrics = 'mae')

# Fit Model
mymodel <- model %>%
  fit(training,
      trainingtarget,
      epochs = 200,
      batch_size = 32,
      validation_split = 0.2)

# Evaluate
model %>% evaluate(test, testtarget)
pred <- model %>% predict(test)

plot(testtarget, pred)
lines(c(min(pred, testtarget), max(pred, testtarget)), c(min(pred, testtarget), max(pred, testtarget)), lwd=2, lty=2)

# more extensive neural network, with extra hidden layer:
model <- keras_model_sequential()
model %>% 
  layer_dense(units = 100, activation = 'relu', input_shape = c(22)) %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 50, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 50, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 20, activation = 'relu') %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1)
summary(model)


# neural network more suited to calculate the sum of squares:
model <- keras_model_sequential()
model %>% 
  layer_dense(units = 200, activation = 'relu', input_shape = c(20)) %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 50, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 10, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 1)
summary(model)

# Compile
model %>% compile(loss = 'mse',
                  optimizer = optimizer_rmsprop(learning_rate = 0.002),
                  #optimizer = 'rmsprop',
                  metrics = 'mae')

# Fit Model
mymodel <- model %>%
  fit(training,
      trainingtarget,
      epochs = 100,
      batch_size = 32,
      validation_split = 0.2)

# Evaluate
model %>% evaluate(test, testtarget)
pred <- model %>% predict(test)

plot(testtarget, pred)
lines(c(min(pred, testtarget), max(pred, testtarget)), 
      c(min(pred, testtarget), max(pred, testtarget)), lwd=2, lty=2)

testtarget2 <- 1/(1 + exp(-1*testtarget))
pred2 <- 1/(1 + exp(-1*pred))
plot(testtarget2, pred2)
lines(c(min(pred2, testtarget2), max(pred2, testtarget2)), 
      c(min(pred2, testtarget2), max(pred2, testtarget2)), lwd=2, lty=2)



par(mfrow=c(1,2))
plot(testtarget[,1], pred[,1])
lines(c(0, max(pred[,1], testtarget[,1])), c(0, max(pred[,1], testtarget[,1])), lwd=2, lty=2)
plot(testtarget[,2], pred[,2])
lines(c(0, max(pred[,2], testtarget[,2])), c(0, max(pred[,2], testtarget[,2])), lwd=2, lty=2)
