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

### Funktion, welche das the Markov-switching Modell mit der VGAM Funktion (also MS-VGAM) mit multiplen Regressoren schätzt
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
                                   sigma = sigma)                               # Wir könnten hier zu Verallgemeinerung noch eine Variable einfügen, die man in die Funktion einfügen muss, um die Density zuberechnen
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
        mod[[j]] = vgam(y ~ s(x1) + s(x2), weights = ind[j,], family = uninormal)        # können wir diesen Schritt irgendwie automatisieren?
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

# Öl Preis Brent

Brent = download_data(symbol = "BZ=F", from = Boundary, to = Boundaryup)
Brent = select(Brent, Date, Close)
Brent = Brent %>% rename(CloseBrent = Close)


## Daten zusammenfügen

data = merge(Shell, BP, by = "Date", all = TRUE)
data = merge(data, TE, by = "Date", all = TRUE)
data = merge(data, EXX, by = "Date", all = TRUE)
data = merge(data, WTI, by = "Date", all = TRUE)
data = merge(data, Brent, by = "Date", all = TRUE)

## Datenaufbereitung

data = na.omit(data)
data = subset(x = data, subset = CloseWTI>0)
data$Date = as.Date(data$Date)

min(data$Date)
max(data$Date)



## Modlelparameter definieren

x = as.matrix(data[,6:7])
x1 = x[,1]
x2 = x[,2]
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
  mu[,i] = cbind(c(mean(y[,i]*0.85),mean(y[,i]*1.15)))  
}


### Model fitting
fitted_0 = list() # Startwerte für den Initialstep setzen
fitted_0[[1]] = cbind(rep(mu[1,1], T), rep(mu[1,2], T),rep(mu[1,3], T), rep(mu[1,4], T))
fitted_0[[2]] = cbind(rep(mu[2,1], T), rep(mu[2,2], T),rep(mu[2,3], T), rep(mu[2,4], T))

## multiples MSVGAM schätzen

mod = ms_vgam(x = x, y = as.matrix(y), fitted_0 = fitted_0, delta_0 = delta, gamma_0 = Gamma)

## Koeffizientenmatrix der geschätzten Modelle

round(coef(mod$mod[[1]], matrix = TRUE), 3)
round(coef(mod$mod[[2]], matrix = TRUE), 3)


### 3D Plots erstellen mit gefitteten Werten

Farbe1 = "green"
Farbe2 = "red"

Shell3D <- plot_ly(x = x[, 1], y = x[, 2], z = fitted(mod$mod[[1]])[,1], type = "scatter3d", mode = 'markers', marker = list(color = Farbe1), name = 'State 1') %>%
  add_trace(x = x[, 1], y = x[, 2], z = fitted(mod$mod[[2]])[,1], type = "scatter3d", mode = 'markers', marker = list(color = Farbe2), name = 'State 2') %>%
  layout(scene = list(
    xaxis = list(title = 'WTI'),
    yaxis = list(title = 'Brent'),
    zaxis = list(title = 'Aktienkurs Shell')
  )) %>%
  layout(title = 'Aktienkurs von Shell in Abhängigkeit beider Ölsorten (fitted Values)')

BP3D <- plot_ly(x = x[, 1], y = x[, 2], z = fitted(mod$mod[[1]])[,2], type = "scatter3d", mode = 'markers', marker = list(color = Farbe1), name = 'State 1') %>%
  add_trace(x = x[, 1], y = x[, 2], z = fitted(mod$mod[[2]])[,2], type = "scatter3d", mode = 'markers', marker = list(color = Farbe2), name = 'State 2') %>%
  layout(scene = list(
    xaxis = list(title = 'WTI'),
    yaxis = list(title = 'Brent'),
    zaxis = list(title = 'Aktienkurs BP')
  )) %>%
  layout(title = 'Aktienkurs von BP in Abhängigkeit beider Ölsorten (fitted Values)')

TE3D <- plot_ly(x = x[, 1], y = x[, 2], z = fitted(mod$mod[[1]])[,3], type = "scatter3d", mode = 'markers', marker = list(color = Farbe1), name = 'State 1') %>%
  add_trace(x = x[, 1], y = x[, 2], z = fitted(mod$mod[[2]])[,3], type = "scatter3d", mode = 'markers', marker = list(color = Farbe2), name = 'State 2') %>%
  layout(scene = list(
    xaxis = list(title = 'WTI'),
    yaxis = list(title = 'Brent'),
    zaxis = list(title = 'Aktienkurs TE')
  )) %>%
  layout(title = 'Aktienkurs von TE in Abhängigkeit beider Ölsorten (fitted Values)')

EXX3D <- plot_ly(x = x[, 1], y = x[, 2], z = fitted(mod$mod[[1]])[,4], type = "scatter3d", mode = 'markers', marker = list(color = Farbe1), name = 'State 1') %>%
  add_trace(x = x[, 1], y = x[, 2], z = fitted(mod$mod[[2]])[,4], type = "scatter3d", mode = 'markers', marker = list(color = Farbe2), name = 'State 2') %>%
  layout(scene = list(
    xaxis = list(title = 'WTI'),
    yaxis = list(title = 'Brent'),
    zaxis = list(title = 'Aktienkurs EXX')
  )) %>%
  layout(title = 'Aktienkurs von EXX in Abhängigkeit beider Ölsorten (fitted Values)')

Shell3D
BP3D
TE3D
EXX3D

### 3D Plots erstellen mit kompletter Ebene, welche die Vorhersagen angibt

## Datenaufbereitung

x1_range = seq(min(x1), max(x1), length.out = 100)
x2_range = seq(min(x2), max(x2), length.out = 100)
grid = expand.grid(x1 = x1_range, x2 = x2_range)

# Berechne die Vorhersagen für die Gitter

pred_grid_shell1 = predict(mod$mod[[1]], newdata = grid, type = "link")[,1:2]
pred_grid_shell2 = predict(mod$mod[[2]], newdata = grid, type = "link")[,1:2]

pred_grid_BP1 = predict(mod$mod[[1]], newdata = grid, type = "link")[,3:4]
pred_grid_BP2 = predict(mod$mod[[2]], newdata = grid, type = "link")[,3:4]

pred_grid_TE1 = predict(mod$mod[[1]], newdata = grid, type = "link")[,5:6]
pred_grid_TE2 = predict(mod$mod[[2]], newdata = grid, type = "link")[,5:6]

pred_grid_EXX1 = predict(mod$mod[[1]], newdata = grid, type = "link")[,7:8]
pred_grid_EXX2 = predict(mod$mod[[2]], newdata = grid, type = "link")[,7:8]

# Setze alle negativen Werte auf 0. da hier keine negativen Aktienkurse inbetracht gezogen werden

pred_grid_shell1 = pmax(pred_grid_shell1, 0)
pred_grid_shell2 = pmax(pred_grid_shell2, 0)

pred_grid_BP1 = pmax(pred_grid_BP1, 0)
pred_grid_BP2 = pmax(pred_grid_BP2, 0)

pred_grid_TE1 = pmax(pred_grid_TE1, 0)
pred_grid_TE2 = pmax(pred_grid_TE2, 0)

pred_grid_EXX1 = pmax(pred_grid_EXX1, 0)
pred_grid_EXX2 = pmax(pred_grid_EXX2, 0)


## 3D Plots erstellen

Shell = plot_ly() %>%
  add_surface(
    x = matrix(grid$x1, nrow = length(x1_range)),
    y = matrix(grid$x2, ncol = length(x2_range)),
    z = matrix(pred_grid_shell1, nrow = length(x1_range)),
    colorscale = list(c(0, 1), c('green', 'yellow')),
    showscale = TRUE,
    name = 'State 1'
  ) %>%
  add_surface(
    x = matrix(grid$x1, nrow = length(x1_range)),
    y = matrix(grid$x2, ncol = length(x2_range)),
    z = matrix(pred_grid_shell2, nrow = length(x1_range)),
    colorscale = list(c(0, 1), c('blue', 'red')),
    showscale = TRUE,
    name = 'State 2'
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = 'Ölpreis WTI'),
      yaxis = list(title = 'Ölpreis Brent'),
      zaxis = list(
        title = 'Aktienschlusskurs Shell',
        range = c(0, max(c(pred_grid_shell1, pred_grid_shell2)))
      )
    ),
    title = 'Aktienkurs von Shell in Abhänigkeit beider Ölsorten'
  )

BP = plot_ly() %>%
  add_surface(
    x = matrix(grid$x1, nrow = length(x1_range)),
    y = matrix(grid$x2, ncol = length(x2_range)),
    z = matrix(pred_grid_BP1, nrow = length(x1_range)),
    colorscale = list(c(0, 1), c('green', 'yellow')),
    showscale = TRUE,
    name = 'State 1'
  ) %>%
  add_surface(
    x = matrix(grid$x1, nrow = length(x1_range)),
    y = matrix(grid$x2, ncol = length(x2_range)),
    z = matrix(pred_grid_BP2, nrow = length(x1_range)),
    colorscale = list(c(0, 1), c('blue', 'red')),
    showscale = TRUE,
    name = 'State 2'
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = 'Ölpreis WTI'),
      yaxis = list(title = 'Ölpreis Brent'),
      zaxis = list(
        title = 'Aktienschlusskurs BP',
        range = c(0, max(c(pred_grid_BP1, pred_grid_BP2)))
      )
    ),
    title = 'Aktienkurs von BP in Abhänigkeit beider Ölsorten'
  )

TE = plot_ly() %>%
  add_surface(
    x = matrix(grid$x1, nrow = length(x1_range)),
    y = matrix(grid$x2, ncol = length(x2_range)),
    z = matrix(pred_grid_TE1, nrow = length(x1_range)),
    colorscale = list(c(0, 1), c('green', 'yellow')),
    showscale = TRUE,
    name = 'State 1'
  ) %>%
  add_surface(
    x = matrix(grid$x1, nrow = length(x1_range)),
    y = matrix(grid$x2, ncol = length(x2_range)),
    z = matrix(pred_grid_TE2, nrow = length(x1_range)),
    colorscale = list(c(0, 1), c('blue', 'red')),
    showscale = TRUE,
    name = 'State 2'
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = 'Ölpreis WTI'),
      yaxis = list(title = 'Ölpreis Brent'),
      zaxis = list(
        title = 'Aktienschlusskurs TE',
        range = c(0, max(c(pred_grid_TE1, pred_grid_TE2)))
      )
    ),
    title = 'Aktienkurs von TE in Abhänigkeit beider Ölsorten'
  )

EXX = plot_ly() %>%
  add_surface(
    x = matrix(grid$x1, nrow = length(x1_range)),
    y = matrix(grid$x2, ncol = length(x2_range)),
    z = matrix(pred_grid_EXX1, nrow = length(x1_range)),
    colorscale = list(c(0, 1), c('green', 'yellow')),
    showscale = TRUE,
    name = 'State 1'
  ) %>%
  add_surface(
    x = matrix(grid$x1, nrow = length(x1_range)),
    y = matrix(grid$x2, ncol = length(x2_range)),
    z = matrix(pred_grid_EXX2, nrow = length(x1_range)),
    colorscale = list(c(0, 1), c('blue', 'red')),
    showscale = TRUE,
    name = 'State 2'
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = 'Ölpreis WTI'),
      yaxis = list(title = 'Ölpreis Brent'),
      zaxis = list(
        title = 'Aktienschlusskurs EXX',
        range = c(0, max(c(pred_grid_EXX1, pred_grid_EXX2)))
      )
    ),
    title = 'Aktienkurs von EXX in Abhänigkeit beider Ölsorten'
  )

Shell
BP
TE
EXX


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
  
  PMF1[,i] = pnorm(data_new[,i], 
                   mean = means1[i], 
                   sd = exp(sds1[i])^2)
  
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

# QQ-Plot der Pseudoresiduan

par(mfrow = c(sqrt(ncol(y)), sqrt(ncol(y))))

for (i in 1:ncol(y)) {
  
  qqnorm(pseudo_res[,i], main = paste("Normal QQ-Plot von", histname[i]))
  qqline(pseudo_res[,i], col = "red", lwd = 2)
  
}




### Zeitreihen und States plotten

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

grid.arrange(MSShell, MSBP, MSTE, MSEXX, ncol = 2)


