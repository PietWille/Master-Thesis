## Install and load required packages
install.packages("VGAM")
install.packages("mvtnorm")

library(VGAM)
library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(scales)
library(ggthemes)

set.seed(35)


### Daten simulieren
## Simulierte Markov Chain

N = 2
T = 500
states = rep(NA, times = T)
delta = rep(1 / N, times = N)
Gamma = matrix(c(0.90, 0.10, 0.10, 0.90), ncol = N)
states[1] = sample(1:N, size = 1, prob = delta)
for(t in 2:T) {
  states[t] = sample(1:N, size = 1, prob = Gamma[states[t - 1], ])
}


## Simulierte Datenpunkte

x = runif(n = T, min = -2, max = 3.4)
y = matrix(NA, nrow = T, ncol = 2)
mu = cbind(c(7.3, 2), c(-5.4, 2))
for(t in 1:T) {
  y[t, 1] = rnorm(n = 1, mean = mu[states[t], 1] * x[t], sd = 1)
  y[t, 2] = rnorm(n = 1, mean = mu[states[t], 2] * x[t], sd = 1)
}

### Funktion, welche das Markov-switching Modell mit der VGLM Funktion (also MS-VGLM) mit einem Regressor schätzt
ms_vgam = function(x, y, fitted_0, N = 2, max_iter = 100, conv_tol = 1e-03) { 
  delta = rep(0.5, times = 2)
  gamma = matrix(c(0.95, 0.05, 0.05, 0.95), ncol = 2)
  fitted = fitted_0
  mod = list()
  term = FALSE
  old = 0
  sigma = matrix(0, ncol = ncol(y), nrow = ncol(y))
  for (i in 1:ncol(y)) {  
    diag(sigma)[i] = var(y[,i])
  }
  
  while(term == FALSE) {
    for(i in 1:max_iter) {
      delta_next = delta
      gamma_next = gamma
      fitted_next = fitted
      allprobs = matrix(NA, T, N)
      lalpha = lbeta = matrix(NA, nrow = N, ncol = T)
      
      for(j in 1:N) {
        for(t in 1:T) {
          allprobs[t, j] = dmvnorm(y[t,], mean = fitted_next[[j]][t,],
                                   sigma = sigma)                               
        }
      }
      
      allprobs = ifelse(!is.na(allprobs), allprobs, 1)
      foo = delta * allprobs[1,]
      sumfoo = sum(foo)
      lscale = log(sumfoo)
      foo = foo / sumfoo                                                        
      lalpha[, 1] = log(foo) + lscale
      
      for(j in 2:nrow(y)) {                                                     
        foo = foo %*% gamma * allprobs[j,]                                          
        sumfoo = sum(foo) 
        lscale = lscale + log(sumfoo)                                          
        foo = foo / sumfoo                                                      
        lalpha[, j] = log(foo) + lscale                                         
      }
      
      foo = rep(1 / N, times = N)
      lbeta[, nrow(y)] = rep(0, times = N)
      lscale = log(N)
      foo = foo / sum(foo)
      
      for(j in (nrow(y) - 1):1) {                                               
        foo = gamma %*% (allprobs[j + 1,] * foo)
        lbeta[, j] = log(foo) + lscale
        sumfoo = sum(foo)
        foo = foo / sumfoo
        lscale = lscale + log(sumfoo)                                          
      }
      
      lallprobs = log(allprobs)
      llh = max(lalpha[, nrow(y)]) + log(sum(exp(lalpha[, nrow(y)] - max(lalpha[, nrow(y)])))) 
      weights = matrix(NA, N, nrow(y))
      
      for(j in 1:nrow(y)) {
        weights[, j] = exp(lalpha[, j] + lbeta[, j] - llh)                      
      }
      
      for(j in 1:N) {
        for(k in 1:N) {
          gamma_next[j, k] = gamma[j, k] * sum(exp(lalpha[j, 1:(nrow(y) - 1)] + lallprobs[2:nrow(y), k] + lbeta[k, 2:nrow(y)] - llh))
        }                                                                       
      }
      
      gamma_next = gamma_next / apply(gamma_next, 1, sum)
      delta_next = exp(lalpha[, 1] + lbeta[, 1] - llh)
      delta_next = delta_next / sum(delta_next)
      conv_crit = sum(abs(gamma - gamma_next)) + sum(abs(delta - delta_next))
      ind = weights
      
      for(j in 1:N){                                                            
        mod[[j]] = vglm(y ~ x, weights = ind[j,], family = uninormal)
        fitted_next[[j]][, 1] = as.vector(fitted(mod[[j]])[, 1])
        fitted_next[[j]][, 2] = as.vector(fitted(mod[[j]])[, 2])
      }
      
      conv_crit = abs(llh - old)
      cat("Iteration = ", i, ", criterion = ", round(conv_crit, 3), "\r", sep = "") 
      
      if(conv_crit < conv_tol | i == max_iter) {
        if(i == max_iter) {
          print(paste("No convergence within", max_iter, "iterations"))
        }
      else{
          print(paste("Convergence after", i, "iterations"))
        }
        term = TRUE
        break
      }  
      delta = delta_next                                                        
      gamma = gamma_next                                                        
      fitted = fitted_next                                                     
      old = llh                                                              
    }
  }
  return(list(x = x, y = y, delta = delta_next, gamma = gamma_next, mod = mod, llh = llh, state_probs = weights))
}

## Model fitting
fitted_0 = list() # Startwerte für den Initialstep setzen
fitted_0[[1]] = cbind(mu[1, 1] * x, mu[1, 2] * x)
fitted_0[[2]] = cbind(mu[2, 1] * x, mu[2, 2] * x)

## MSVGLM schätzen

mod = ms_vgam(x = x, y = y, fitted_0 = fitted_0)


### Geschätzte Daten plotten

## Farben für die States definieren und lokales State Decoding durch die stateprobs durchführen

Farbe1 = "green"
Farbe2 = "red"

# lokales State Decoding

Farbe1 = "green"
Farbe2 = "red"

WZustand = apply(mod$state_probs, MARGIN = 2, FUN = which.max)
n = length(WZustand)
colour = rep(NA, n)
for (i in 1:n){                                                              
  if (WZustand[i]==1)
  {
    colour[i] = Farbe1
  }
  else 
  {
    colour[i] = Farbe2
  }
}


## Plots mit geschätzten Werten

GGdata = data.frame(y,x)

Ploty1 = ggplot(data = GGdata, aes(x,y[,1])) + geom_point(col ="grey", na.rm = TRUE)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 1]), na.rm = TRUE, colour = Farbe1, size = 1)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 1]), na.rm = TRUE, colour = Farbe2, size = 1)+
  xlab("X") + ylab("Responsevariable y_1") + ggtitle("Geschätze Modelle Responsevariable 1")+
  theme_bw()

Ploty2 = ggplot(data = GGdata, aes(x,y[,2])) + geom_point(col = "grey", na.rm = TRUE)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 2]), na.rm = TRUE, colour = Farbe1, size = 1)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 2]), na.rm = TRUE, colour = Farbe2, size = 1)+
  xlab("X") + ylab("Responsevariable y_2") + ggtitle("Geschätze Modelle Responsevariable 2")+
  theme_bw()

grid.arrange(Ploty1, Ploty2, ncol = 2)


## Plots mit geschätzten Werten und zuordnung der States durch das lokale State Decoding

GGdata = data.frame(y,x)

FPloty1 = ggplot(data = GGdata, aes(x,y[,1])) + geom_point(col = colour, na.rm = TRUE, alpha = 0.3)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 1]), na.rm = TRUE, colour = "black", size = 1.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 1]), na.rm = TRUE, colour = "black", size = 1.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 1]), na.rm = TRUE, colour = Farbe1, size = 1)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 1]), na.rm = TRUE, colour = Farbe2, size = 1)+
  xlab("X") + ylab("Responsevariable y_1") + ggtitle("Geschätze Modelle Responsevariable 1")+
  theme_bw()

FPloty2 = ggplot(data = GGdata, aes(x,y[,2])) + geom_point(col = colour, na.rm = TRUE, alpha = 0.3)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 2]), na.rm = TRUE, colour = "black", size = 1.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 2]), na.rm = TRUE, colour = "black", size = 1.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 2]), na.rm = TRUE, colour = Farbe1, size = 1)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 2]), na.rm = TRUE, colour = Farbe2, size = 1)+
  xlab("X") + ylab("Responsevariable y_2") + ggtitle("Geschätze Modelle Responsevariable 2")+
  theme_bw()

grid.arrange(FPloty1, FPloty2, ncol = 2)

grid.arrange(Ploty1, Ploty2, FPloty1, FPloty2, ncol = 2)

