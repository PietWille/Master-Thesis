### Installiere und lade die Packete
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

set.seed(35)

### Datenimport via fHMM

# Zeitgrenzen damit immer der selbe Datensatz importiert wird 

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

# Öl Preis WTI

WTI = download_data(symbol = "CL=F", from = Boundary, to = Boundaryup)
WTI = select(WTI, Date, Close)
WTI = WTI %>% rename(CloseWTI = Close)


## Datan zusammenfügen

data = merge(Shell, BP, by = "Date", all = TRUE)
data = merge(data, TE, by = "Date", all = TRUE)
data = merge(data, WTI, by = "Date", all = TRUE)
data = merge(data, EXX, by = "Date", all = TRUE)


## Datenaufbereitung

data = na.omit(data)
data = subset(x = data, subset = CloseWTI>0)
data$Date = as.Date(data$Date)

## Modlelparameter definieren

x = as.matrix(data$CloseWTI)
y = cbind(data$CloseShell, data$CloseBP, data$CloseTE, data$CloseEXX)
xdim = ncol(x)
ycol = ncol(y)
T = nrow(y)

## Modellparameter der Markov Chain

N = 2
delta = rep(1 / N, times = N)
Gamma = matrix(c(0.95, 0.05, 0.05, 0.95), ncol = N)

## Weitere Parameter

mu = matrix(NA, nrow = N, ncol = ycol)
for(i in 1:ycol) {
  mu[,i] = cbind(c(mean(y[,i]*0.85),mean(y[,i]*1.15)))}





### Veranschaulichung der Daten um ein Verständis über die Zeitreihen und evtl. vorliegenden States zu erlangen

# BP und WTI

plotBP = ggplot(data = data, mapping = aes(Date, CloseBP))+
  geom_line(na.rm = TRUE)+
  ggtitle("Aktienkurs von BP")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  (scale_x_date(breaks = date_breaks("2 years"),labels = date_format("%Y")))+
  theme_bw()

SP_BP = ggplot(data = data, mapping = aes(CloseWTI, CloseBP))+
  geom_point(na.rm = TRUE, col = "darkgrey")+
  ggtitle("Scatterplot Öl-Preis und Aktienkurs von BP")+
  xlab("Tagesschlusskurs WTI") + ylab("Tagesschlusskurs")+
  xlim(0,150)+
  theme_bw()


# Shell und WTI

plotShell = ggplot(data = data, mapping = aes(Date, CloseShell))+
  geom_line(na.rm = TRUE)+
  ggtitle("Aktienkurs von Shell")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  (scale_x_date(breaks = date_breaks("2 years"),labels = date_format("%Y")))+
  theme_bw()

SP_Shell = ggplot(data = data, mapping = aes(CloseWTI, CloseShell))+
  geom_point(na.rm = TRUE, col = "darkgrey")+
  ggtitle("Scatterplot Öl-Preis und Aktienkurs von Shell")+
  xlab("Tagesschlusskurs WTI") + ylab("Tagesschlusskurs")+
  xlim(0,150)+
  theme_bw()


# TE und WTI

plotTE = ggplot(data = data, mapping = aes(Date, CloseTE))+
  geom_line(na.rm = TRUE)+
  ggtitle("Aktienkurs von TE")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  (scale_x_date(breaks = date_breaks("2 years"),labels = date_format("%Y")))+
  theme_bw()

SP_TE = ggplot(data = data, mapping = aes(CloseWTI, CloseTE))+
  geom_point(na.rm = TRUE, col = "darkgrey")+
  ggtitle("Scatterplot Öl-Preis und Aktienkurs von TE")+
  xlab("Tagesschlusskurs WTI") + ylab("Tagesschlusskurs")+
  xlim(0,150)+
  theme_bw()


# EXX und WTI

plotEXX = ggplot(data = data, mapping = aes(Date, CloseEXX))+
  geom_line(na.rm = TRUE)+
  ggtitle("Aktienkurs von EXX")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  (scale_x_date(breaks = date_breaks("2 years"),labels = date_format("%Y")))+
  theme_bw()

SP_EXX = ggplot(data = data, mapping = aes(CloseWTI, CloseEXX))+
  geom_point(na.rm = TRUE, col = "darkgrey")+
  ggtitle("Scatterplot Öl-Preis und Aktienkurs von EXX")+
  xlab("Tagesschlusskurs WTI") + ylab("Tagesschlusskurs")+
  xlim(0,150)+
  theme_bw()


# Scatterplot und Zeitreihen kombiniert

grid.arrange(plotShell, SP_Shell, plotBP, SP_BP, plotTE, SP_TE, plotEXX, SP_EXX,ncol = 2)


### Modellselektion

colnames(y) = c("Shell", "BP", "TE", "EXX")

LM = vglm(y ~ x, family = uninormal)
LM = rbind(AIC(LM), BIC(LM))

vglm1 = vglm(y ~ x + I(x^2), family = uninormal)
VGLM1 = rbind(AIC(vglm1), BIC(vglm1))

vglm2 = vglm(y ~ x + I(x^2) + I(x^3), family = uninormal)
VGLM2 = rbind(AIC(vglm2), BIC(vglm2))

vgam1 = vgam(y ~ x + s(x), family = uninormal)
VGAM1 = rbind(AIC(vgam1), BIC(vgam1))

vgam2 = vgam(y ~ s(x), family = uninormal)

round(coef(vgam2, matrix = TRUE)[,1:4],4)
round(coef(vgam2, matrix = TRUE)[,5:8],4)
VGAM2 = rbind(AIC(vgam2), BIC(vgam2))

Modellselection = cbind(LM, VGLM1, VGLM2, VGAM1, VGAM2)
row.names(Modellselection) = c("AIC","BIC")
colnames(Modellselection) = c("LM", "VGLM1", "VGLM2", "VGAM1", "VGAM2")
Modellselection


### Funktion, welche das the Markov-switching Modell mit der VGAM Funktion (also MS-VGAM) mit einem Regressor schätzt
ms_vgam = function(x, y, fitted_0, delta_0, gamma_0, N = 2, max_iter = 100, conv_tol = 1e-03) { 
  delta = delta_0
  gamma = gamma_0
  fitted = fitted_0
  mod = list()
  term = FALSE
  old = 0
  sigma = matrix(0, ncol = ycol, nrow = ycol) # Festlegen der Kovarianzmatrix der PMF
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
      
      for(j in 1:N) { # Berechnung alle Wahrscheinlichkeiten aller y_i und Beobachtungen durch die PMF
        for(t in 1:T) {
          allprobs[t, j] = dmvnorm(y[t,], mean = fitted_next[[j]][t,],
                                   sigma = sigma)                               
        }
      }
      
      allprobs = ifelse(!is.na(allprobs), allprobs, 1) 
      foo = delta * allprobs[1,]  # Berechnung der initial forw. Probs
      sumfoo = sum(foo)
      lscale = log(sumfoo) # Log Variante verwenden
      foo = foo / sumfoo  # Scale Paramater \Phi                                                     
      lalpha[, 1] = log(foo) + lscale
      
      for(j in 2:T) { # Iterative Berechnung der forw. Probs                                                    
        foo = foo %*% gamma * allprobs[j,]                                          
        sumfoo = sum(foo) 
        lscale = lscale + log(sumfoo)                                           
        foo = foo / sumfoo                                                     
        lalpha[, j] = log(foo) + lscale                                         
      }
      
      foo = rep(1 / N, times = N)
      lbeta[, T] = rep(0, times = N)
      lscale = log(N)
      foo = foo / sum(foo) # Scale Paramater \Phi 
      
      for(j in (T - 1):1) { # Iterative Berechnung der backw. Probs                                            
        foo = gamma %*% (allprobs[j + 1,] * foo)
        lbeta[, j] = log(foo) + lscale # Log Variante verwenden
        sumfoo = sum(foo)
        foo = foo / sumfoo
        lscale = lscale + log(sumfoo)                                           
      }
      
      lallprobs = log(allprobs)
      llh = max(lalpha[, T]) + log(sum(exp(lalpha[, T] - max(lalpha[, T])))) # Berechnung der LLH
      weights = matrix(NA, N, T)
      
      for(j in 1:T) { # u's mit exp weil oben der log angewendet wurde  
        weights[, j] = exp(lalpha[, j] + lbeta[, j] - llh)                   
      }
      
      for(j in 1:N) {
        for(k in 1:N) { # v's mit exp weil oben der log angewendet wurde  
          gamma_next[j, k] = gamma[j, k] * sum(exp(lalpha[j, 1:(T - 1)] + lallprobs[2:T, k] + lbeta[k, 2:T] - llh))
        }                                                                       
      }
      
      gamma_next = gamma_next / apply(gamma_next, 1, sum) # Berechnung der neuen Parameter
      delta_next = exp(lalpha[, 1] + lbeta[, 1] - llh) # Berechnung der neuen Parameter
      delta_next = delta_next / sum(delta_next)
      conv_crit = sum(abs(gamma - gamma_next)) + sum(abs(delta - delta_next))
      ind = weights # Setzen der Zustandswahrscheinlichkeiten als Gewichtung
      
      for(j in 1:N){ # Schätzen der Parameter durch das gewählt VGAM. Dieses muss je nach gewähltem Modell angepasst werden                                                          
        mod[[j]] = vgam(y ~ s(x), weights = ind[j,], family = uninormal)        
        for (m in 1:xdim) {
          fitted_next[[j]][, m] = as.vector(fitted(mod[[j]])[, m])          
        }
      }
      
      conv_crit = abs(llh - old) # Berechnung des Differenzschritt für das Konvergenzkriteriums
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
      delta = delta_next  # Neue Paramter für den neuen E-Step setzten                                                     
      gamma = gamma_next  # Neue Paramter für den neuen E-Step setzten                                                      
      fitted = fitted_next # Neue Paramter für den neuen E-Step setzten durch das VGAM setzen (Summand 3)                                                 
      old = llh   
    }
  }
  return(list(x = x, y = y, delta = delta_next, gamma = gamma_next, mod = mod, llh = llh, state_probs = weights))
}


### Model fitting
fitted_0 = list() # Startwerte für den Initialstep setzen
fitted_0[[1]] = cbind(rep(mu[1,1], T), rep(mu[1,2], T),rep(mu[1,3], T), rep(mu[1,4], T))
fitted_0[[2]] = cbind(rep(mu[2,1], T), rep(mu[2,2], T),rep(mu[2,3], T), rep(mu[2,4], T))

## MSVGAM schätzen

mod = ms_vgam(x = x, y = as.matrix(y), fitted_0 = fitted_0, delta_0 = delta, gamma_0 = Gamma)


## Koeffizientenmatrix der geschätzten Modelle

round(coef(mod$mod[[1]], matrix = TRUE), 3)
round(coef(mod$mod[[2]], matrix = TRUE), 3)



### Geschätzte Daten plotten


## Farben für die States definieren und lokales State Decoding durch die stateprobs durchführen

Farbe1 = "green"
Farbe2 = "red"

# lokales State Decoding

WZustand = apply(mod$state_probs, MARGIN = 2, FUN = which.max)
n = length(WZustand)
colour = rep(NA, n)
for (i in 1:n){                                                                 # Farbe für die weigths definieren
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

## Plots mit geschätzten Werten und zuordnung der States durch das lokale State Decoding

CShellplot = ggplot(data = GGdata, aes(x,y[,1])) + geom_point(col =colour, na.rm = TRUE, alpha = 0.05)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 1]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 1]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 1]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 1]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("fitted Aktienkurs") + ggtitle("Geschätze Modelle der Aktie Shell der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

CBPplot = ggplot(data = GGdata, aes(x,y[,2])) + geom_point(col =colour, na.rm = TRUE, alpha = 0.05)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 2]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 2]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 2]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 2]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("fitted Aktienkurs") + ggtitle("Geschätze Modelle der Aktie BP der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

CTEplot = ggplot(data = GGdata, aes(x,y[,3])) + geom_point(col =colour, na.rm = TRUE, alpha = 0.05)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 3]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 3]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 3]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 3]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("fitted Aktienkurs") + ggtitle("Geschätze Modelle der Aktie TE der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

CEXXplot = ggplot(data = GGdata, aes(x,y[,4])) + geom_point(col =colour, na.rm = TRUE, alpha = 0.05)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 4]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 4]), na.rm = TRUE, colour = "black", size = 2.5)+
  geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 4]), na.rm = TRUE, colour = Farbe1, size = 2)+
  geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 4]), na.rm = TRUE, colour = Farbe2, size = 2)+
  xlab("Öl Preis") + ylab("fitted Aktienkurs") + ggtitle("Geschätze Modelle der Aktie EXX der 2 Zustände")+
  expand_limits(x=c(-5,200), y=c(20,100))+
  theme_bw()

grid.arrange(CShellplot, CBPplot, CTEplot, CEXXplot, ncol = 2)

### Modellchecking via Pseudoresiduan

Zustandswahrscheinlichkeiten = mod$state_probs
x_vals = normal_density = standardized_pseudo_res = pseudo_res = PMF1 = PMF2 = matrix(NA, nrow = T, ncol = ncol(y))
means1 = coef(mod$mod[[1]], matrix = TRUE)[1,c(1,3,5,7)]
sds1 = coef(mod$mod[[1]], matrix = TRUE)[1,c(2,4,6,8)]
means2 = coef(mod$mod[[2]], matrix = TRUE)[1,c(1,3,5,7)]
sds2 = coef(mod$mod[[2]], matrix = TRUE)[1,c(2,4,6,8)] 
mean_pseudo_res = rep(NA, times = ncol(y))
sd_pseudo_res = rep(NA, times = ncol(y))
data_new = data[,c(2:4,6)]
histname = c("Shell", "BP", "TE", "EXX")

par(mfrow = c(sqrt(ncol(y)), sqrt(ncol(y))))
for (i in 1:ncol(y)) {
  
  ## Verteilungsfunktion der Beobachtungen berechnen
  
  # State 1
  
  PMF1[,i] = pnorm(data_new[,i], 
                   mean = means1[i], 
                   sd = exp(sds1[i])^2)
  
  # State 1
  
  PMF2[,i] = pnorm(data_new[,i], 
                   mean = means2[i], 
                   sd = exp(sds2[i])^2)
  
  
  # Berechnung der Pseudoresiduan
  
  pseudo_res[,i] = qnorm(Zustandswahrscheinlichkeiten[1,] * PMF1[,i] + Zustandswahrscheinlichkeiten[2,] * PMF2[,i]) 
  
  ## Standartisieren der Pseudoresiduan
  
  
  # Berechne Mittelwert und Standardabweichung der Pseudoresiduan
  
  mean_pseudo_res[i] = mean(pseudo_res[,i])
  sd_pseudo_res[i] = sd(pseudo_res[,i])
  
  
  # Standardisieren
  
  standardized_pseudo_res[,i] = (pseudo_res[,i] - mean_pseudo_res[i]) / sd_pseudo_res[i]
  
  ## Histo mit Standardnormalverteilung
  
  # Berechne die Normalverteilungskurve für den Bereich der standardisierten Pseudoresiduan
  
  x_vals[,i] = seq(min(standardized_pseudo_res[,i]), max(standardized_pseudo_res[,i]), length = T)
  normal_density[,i] = dnorm(x_vals[,i])
  
  
  # Histogramm der standardisierten Pseudoresiduan
  
  hist(standardized_pseudo_res[,i], freq = FALSE, breaks = 30, 
       main = paste("Histogram der Pseudoresiduan von", histname[i]), 
       xlab = "Standardisierte Pseudoresiduan")
  
  
  # Zeichne die Normalverteilungskurve in das Histogramm
  
  lines(x_vals[,i], normal_density[,i], col = "red", lwd = 2)
  
}

## QQ-Plot der Pseudoresiduan

par(mfrow = c(sqrt(ncol(y)), sqrt(ncol(y))))

for (i in 1:ncol(y)) {
  
  qqnorm(pseudo_res[,i], main = paste("Normal QQ-Plot von", histname[i]))
  qqline(pseudo_res[,i], col = "red", lwd = 2)
  
}


### Zeitreihen der Beobachtungen

## Zeitreihe über den ganzen Zeitraum und durch lokales State Decoding eingefärbt

MSBP = ggplot(data = data, mapping = aes(Date, CloseBP))+
  geom_line(na.rm = TRUE, col = colour)+
  ggtitle("Aktienkurs von BP")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  (scale_x_date(breaks = date_breaks("3 years"),labels = date_format("%Y")))+
  theme_bw()

MSShell = ggplot(data = data, mapping = aes(Date, CloseShell))+
  geom_line(na.rm = TRUE, col = colour)+
  ggtitle("Aktienkurs von Shell")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  (scale_x_date(breaks = date_breaks("3 years"),labels = date_format("%Y")))+
  theme_bw()

MSTE = ggplot(data = data, mapping = aes(Date, CloseTE))+
  geom_line(na.rm = TRUE, col = colour)+
  ggtitle("Aktienkurs von TE")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  (scale_x_date(breaks = date_breaks("3 years"),labels = date_format("%Y")))+
  theme_bw()

MSEXX = ggplot(data = data, mapping = aes(Date, CloseEXX))+
  geom_line(na.rm = TRUE, col = colour)+
  ggtitle("Aktienkurs von EXX")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  (scale_x_date(breaks = date_breaks("3 years"),labels = date_format("%Y")))+
  theme_bw()


grid.arrange(MSShell, MSBP,  MSTE, MSEXX, ncol = 2)


## Verkürzte Zeitreihe um mit dem multiplen MSVGAM zuvergleichen

datav = subset(x = data, subset = Date >= "2007-07-30")
tail = length(datav$Date)
colourv = tail(colour, tail)

VMSBP = ggplot(data = datav %>% filter(Date >= as.Date("2007-07-30")), mapping = aes(Date, CloseBP))+
  geom_line(na.rm = TRUE, col = colourv)+
  ggtitle("Aktienkurs von BP")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  scale_x_date(breaks = date_breaks("3 years"), labels = date_format("%Y"))+
  theme_bw()

VMSShell = ggplot(data = datav %>% filter(Date >= as.Date("2007-07-30")), mapping = aes(Date, CloseShell))+
  geom_line(na.rm = TRUE, col = colourv)+
  ggtitle("Aktienkurs von Shell")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  scale_x_date(breaks = date_breaks("3 years"), labels = date_format("%Y"))+
  theme_bw()

VMSTE = ggplot(data = datav %>% filter(Date >= as.Date("2007-07-30")), mapping = aes(Date, CloseTE))+
  geom_line(na.rm = TRUE, col = colourv)+
  ggtitle("Aktienkurs von TE")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  scale_x_date(breaks = date_breaks("3 years"), labels = date_format("%Y"))+
  theme_bw()

VMSEXX = ggplot(data = datav %>% filter(Date >= as.Date("2007-07-30")), mapping = aes(Date, CloseEXX))+
  geom_line(na.rm = TRUE, col = colourv)+
  ggtitle("Aktienkurs von EXX")+
  xlab("Jahr") + ylab("Tagesschlusskurs")+
  scale_x_date(breaks = date_breaks("3 years"), labels = date_format("%Y"))+
  theme_bw()

grid.arrange(VMSBP, VMSShell, VMSTE, VMSEXX, ncol = 2)
