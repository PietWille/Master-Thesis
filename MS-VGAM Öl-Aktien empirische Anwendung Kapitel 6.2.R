## Install and load required packages
install.packages("VGAM")
install.packages("ggthemes", dependencies = TRUE)
install.packages("scales")
install.packages("mvtnorm")
install.packages("ggplot2")
install.packages("fHMM")
install.packages("gridExtra")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("functional")
install.packages("plotly")

library(VGAM)
library(mvtnorm)
library(fHMM)
library(tidyverse)
library(functional)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales)
library(ggthemes)
library(plotly)

set.seed(35)

## Daten ziehen 

Boundary ="1990-01-01"
Boundaryup ="2024-08-09"

# Shell
Shell = download_data(symbol ="SHEL", from = Boundary, to = Boundaryup)
Shell = select(Shell, Date, Close)
Shell = Shell %>% rename(CloseShell = Close)


# BP
BP = download_data(symbol ="BP", from = Boundary, to = Boundaryup)
BP = select(BP, Date, Close)
BP = BP %>% rename(CloseBP = Close)


# Total Energie
TE = download_data(symbol ="TTE", from = Boundary, to = Boundaryup)
TE = select(TE, Date, Close)
TE = TE %>% rename(CloseTE = Close)


# Exxon Mobil
EXX = download_data(symbol ="XOM", from = Boundary, to = Boundaryup)
EXX = select(EXX, Date, Close)
EXX = EXX %>% rename(CloseEXX = Close)

# Öl Preise
WTI = download_data(symbol = "CL=F", from = Boundary, to = Boundaryup)
WTI = select(WTI, Date, Close)
WTI = WTI %>% rename(CloseWTI = Close)

Brent = download_data(symbol = "BZ=F", from = Boundary, to = Boundaryup)
Brent = select(Brent, Date, Close)
Brent = Brent %>% rename(CloseBrent = Close)


# data nur mit WTI, da die Zeitreihen dann länger sind.

data = merge(Shell, BP, by = "Date", all = TRUE)
data = merge(data, TTE, by = "Date", all = TRUE)
data = merge(data, EXX, by = "Date", all = TRUE)
data = merge(data, WTI, by = "Date", all = TRUE)


## Datenaufbereitung

data = na.omit(data)
data = subset(x = data, subset = CloseWTI>0)
data$Date = as.Date(data$Date)


### Variablen für MS-VGAM


## Regressionsvariblen

x = as.matrix(data$CloseWTI)
y = as.matrix(data[,2:5])

xdim = ncol(x)
ycol = ncol(y)
T = nrow(y)


## Properties Markov Chain

N = 2
delta = rep(1 / N, times = N)
Gamma = matrix(c(0.95, 0.05, 0.05, 0.95), ncol = N)


## weitere Eigenschaften

mu = matrix(NA, nrow = N, ncol = ycol)
for(i in 1:ycol) {
  mu[,i] = cbind(c(mean(y[,i]*0.85),mean(y[,i]*1.15)))}

Farbe1 = "green"
Farbe2 = "red"


### Funktion, welche das Markov-switching Modell mit der VGAM Funktion (also MS-VGAM) mit einem Regressor schätzt
ms_vgam = function(x, y, fitted_0, delta_0, gamma_0, N = 2, max_iter = 100, conv_tol = 1e-03) { 
  delta = delta_0
  gamma = gamma_0
  fitted = fitted_0
  mod = list()
  term = FALSE
  old = 0
  sigma = matrix(0, ncol = ycol, nrow = ycol)
  for (i in 1:ycol) {  
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
      
      for(j in 2:T) {                                                     
        foo = foo %*% gamma * allprobs[j,]                                          
        sumfoo = sum(foo) 
        lscale = lscale + log(sumfoo)                                           
        foo = foo / sumfoo                                                     
        lalpha[, j] = log(foo) + lscale                                         
      }
      
      foo = rep(1 / N, times = N)
      lbeta[, T] = rep(0, times = N)
      lscale = log(N)
      foo = foo / sum(foo)
      
      for(j in (T - 1):1) {                                              
        foo = gamma %*% (allprobs[j + 1,] * foo)
        lbeta[, j] = log(foo) + lscale
        sumfoo = sum(foo)
        foo = foo / sumfoo
        lscale = lscale + log(sumfoo)                                           
      }
      
      lallprobs = log(allprobs)
      llh = max(lalpha[, T]) + log(sum(exp(lalpha[, T] - max(lalpha[, T])))) 
      weights = matrix(NA, N, T)
      
      for(j in 1:T) {
        weights[, j] = exp(lalpha[, j] + lbeta[, j] - llh)                      
      }
      
      for(j in 1:N) {
        for(k in 1:N) {
          gamma_next[j, k] = gamma[j, k] * sum(exp(lalpha[j, 1:(T - 1)] + lallprobs[2:T, k] + lbeta[k, 2:T] - llh))
        }                                                                       
      }
      
      gamma_next = gamma_next / apply(gamma_next, 1, sum)
      delta_next = exp(lalpha[, 1] + lbeta[, 1] - llh)
      delta_next = delta_next / sum(delta_next)
      conv_crit = sum(abs(gamma - gamma_next)) + sum(abs(delta - delta_next))
      ind = weights
      
      for(j in 1:N){                                                            
        mod[[j]] = vglm(y ~ x, weights = ind[j,], family = uninormal)        
        for (m in 1:xdim) {
          fitted_next[[j]][, m] = as.vector(fitted(mod[[j]])[, m])          
        }
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


### Model fitting
fitted_0 = list()
fitted_0[[1]] = cbind(rep(mu[1,1], T), rep(mu[1,2], T),rep(mu[1,3], T), rep(mu[1,4], T))
fitted_0[[2]] = cbind(rep(mu[2,1], T), rep(mu[2,2], T),rep(mu[2,3], T), rep(mu[2,4], T))
mod = ms_vgam(x = x, y = as.matrix(y), fitted_0 = fitted_0, delta_0 = delta, gamma_0 = Gamma)


### Zeitreihen und States plotten

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

## Plots für eine Regressorvariable x

GGdata = data.frame(x,y)

Shellplot = ggplot(data = GGdata, aes(x,y[,1])) + geom_point(col ="grey", na.rm = TRUE)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 1]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 1]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("Aktienschlusskurs") + ggtitle("Geschätze Modelle der Aktie Shell der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

BPplot = ggplot(data = GGdata, aes(x,y[,2])) + geom_point(col ="grey", na.rm = TRUE)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 2]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 2]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("Aktienschlusskurs") + ggtitle("Geschätze Modelle der Aktie BP der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

TEplot = ggplot(data = GGdata, aes(x,y[,3])) + geom_point(col ="grey", na.rm = TRUE)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 3]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 3]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("Aktienschlusskurs") + ggtitle("Geschätze Modelle der Aktie TE der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

EXXplot = ggplot(data = GGdata, aes(x,y[,4])) + geom_point(col ="grey", na.rm = TRUE)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 4]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 4]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("Aktienschlusskurs") + ggtitle("Geschätze Modelle der Aktie EXX der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

grid.arrange(Shellplot, BPplot, TEplot, EXXplot, ncol = 2)

## Plots mit zuordnung der States

CShellplot = ggplot(data = GGdata, aes(x,y[,1])) + geom_point(col =colour, na.rm = TRUE, alpha = 0.05)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 1]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 1]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 1]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 1]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("Aktienschlusskurs") + ggtitle("Geschätze Modelle der Aktie Shell der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

CBPplot = ggplot(data = GGdata, aes(x,y[,2])) + geom_point(col =colour, na.rm = TRUE, alpha = 0.05)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 2]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 2]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 2]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 2]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("Aktienschlusskurs") + ggtitle("Geschätze Modelle der Aktie BP der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

CTEplot = ggplot(data = GGdata, aes(x,y[,3])) + geom_point(col =colour, na.rm = TRUE, alpha = 0.05)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 3]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 3]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 3]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 3]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("Aktienschlusskurs") + ggtitle("Geschätze Modelle der Aktie TE der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

CEXXplot = ggplot(data = GGdata, aes(x,y[,4])) + geom_point(col =colour, na.rm = TRUE, alpha = 0.05)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 4]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 4]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 4]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 4]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("Aktienschlusskurs") + ggtitle("Geschätze Modelle der Aktie EXX der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

grid.arrange(CShellplot, CBPplot, CTEplot, CEXXplot, ncol = 2)




### Pseudoresiduen

Zustandswahrscheinlichkeiten = mod$state_probs
x_vals = normal_density = standardisierte_pseudo_res = pseudo_res = PMF1 = PMF2 = matrix(NA, nrow = T, ncol = ncol(y))
means1 = coef(mod$mod[[1]], matrix = TRUE)[1,c(1,3,5,7)]
sds1 = coef(mod$mod[[1]], matrix = TRUE)[1,c(2,4,6,8)]
means2 = coef(mod$mod[[2]], matrix = TRUE)[1,c(1,3,5,7)]
sds2 = coef(mod$mod[[2]], matrix = TRUE)[1,c(2,4,6,8)] 
mean_pseudo_res = rep(NA, times = ncol(y))
sd_pseudo_res = rep(NA, times = ncol(y))
par(mfrow = c(sqrt(ncol(y)), sqrt(ncol(y))))
data_new = data[2:5]
histname = c("Shell", "BP", "TE", "EXX")
for (i in 1:ncol(y)) {
  
  PMF1[,i] = pnorm(data_new[,i], 
                mean = means1[i], 
                sd = exp(sds1[i])^2)
  
  PMF2[,i] = pnorm(data_new[,i], 
                    mean = means2[i], 
                    sd = exp(sds2[i])^2)
  
  # Berechnung der Pseudoresiduen
  pseudo_res[,i] = qnorm(Zustandswahrscheinlichkeiten[1,] * PMF1[,i] + Zustandswahrscheinlichkeiten[2,] * PMF2[,i]) 
  
  ## Standartisieren der Pseudoresiduen
  
  # Berechne Mittelwert und Standardabweichung der Pseudoresiduen
  mean_pseudo_res[i] = mean(pseudo_res[,i])
  sd_pseudo_res[i] = sd(pseudo_res[,i])
  
  # Standardisiere die Pseudoresiduen
  standardisierte_pseudo_res[,i] = (pseudo_res[,i] - mean_pseudo_res[i]) / sd_pseudo_res[i]
  
  ## Histo mit Standardnormalverteilung
  
  # Berechne die Normalverteilungskurve für den Bereich der standardisierten Pseudoresiduen
  x_vals[,i] = seq(min(standardisierte_pseudo_res[,i]), max(standardisierte_pseudo_res[,i]), length = T)
  normal_density[,i] = dnorm(x_vals[,i])
  
  # Histogramm der standardisierten Pseudoresiduen
  hist(standardisierte_pseudo_res[,i], freq = FALSE, breaks = 30, 
       main = paste("Histogram der Pseudoresiduen von", histname[i]), 
       xlab = "Standardisierte Pseudoresiduen")
  
  # Zeichne die Normalverteilungskurve in das Histogramm
  lines(x_vals[,i], normal_density[,i], col = "red", lwd = 2)
  
}

# QQ-Plot der Pseudoresiduen

par(mfrow = c(sqrt(ncol(y)), sqrt(ncol(y))))

for (i in 1:ncol(y)) {
  
  qqnorm(pseudo_res[,i], main = paste("Normal QQ-Plot von", histname[i]))
  qqline(pseudo_res[,i], col = "red", lwd = 2)

}







