### Installiere und lade die Packete
install.packages("VGAM")
install.packages("mvtnorm")
install.packages("dplyr")
install.packages("pryr")
install.packages("xtable")



library(VGAM)
library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(scales)
library(ggthemes)
library(dplyr)
library(pryr)
library(xtable)

set.seed(5000)

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

x = rep(NA, times = T)
y = matrix(NA, nrow = T, ncol = 2)
x = rnorm(n = T, mean = 5, sd = 1) 
mu = cbind(c(7.3, 2), c(-5.4, 4))
for(t in 1:T) {
  y[t, 1] = sin(x[t]) - 2*cos(x[t])* rnorm(n = 1, mean = mu[states[t], 1]*x[t], sd = 4)
  y[t, 2] = -2*sin(x[t]) + 2*cos(x[t])* rnorm(n = 1, mean = mu[states[t], 2]*x[t], sd = 4)
}


### Funktion, welche das Markov-switching Modell mit der VGAM Funktion (also MS-VGAM) mit einem Regressor schätzt
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
        mod[[j]] = vgam(y ~ s(x), weights = ind[j,], family = uninormal)
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
  return(list(x = x, y = y, delta = delta_next, gamma = gamma_next, mod = mod, llh = llh, state_probs = weights, allprobs = allprobs, fitted = fitted))
}


### Simulation von M MSVGAMs mit erster Plot Initialdaten von Oberhalb der Funktion


M = 4
yM = matrix(NA, nrow = T, ncol = 2)
xM = matrix(NA, nrow = T, ncol = M)
Zustandsvektor = rep(NA, times = T)
maxZustandsmatrix = matrix(NA, nrow = T, ncol = M)
Zustandsmatrix_next = matrix(0, nrow = 2, ncol = T)

#start_time = Sys.time() #Falls ein Benchmark durchgeführt wird
#start_memory = mem_used() #Falls ein Benchmark durchgeführt wird

for (i in 1:M) {
  set.seed(i)
  
    if (i == 1) {
      
      ## Initial Mod und Definitionen für die Plots
      
      fitted_0 = list()
      fitted_0[[1]] = cbind(mu[1, 1] * x, mu[1, 2] * x)
      fitted_0[[2]] = cbind(mu[1, 2] * x, mu[2, 2] * x)
      mod = ms_vgam(x = x, y = y, fitted_0 = fitted_0)
      
      GGdata = data.frame(y,x)
      
      ### Plot der y und ihrer States des Initialsteps mit Scatterplot, um die Weiteren Simulationen später einzuzeichnen
      
      # Erstellen einer Sequenz von x-Werten für die Berechnung der Funktionen und SE
      x_seq = seq(min(GGdata$x), max(GGdata$x), length.out = T)
      
      # Berechnung der y-Werte und Dummy SE für die erste Funktion
      y_fun1 = sin(x_seq) - 2 * cos(x_seq) * (mu[1, 1] * x_seq)
      se_fun1 = 0.1 * abs(y_fun1) # Dummy SE-Werte, passen Sie diese entsprechend Ihrer Daten an
      
      # Berechnung der y-Werte und Dummy SE für die zweite Funktion
      y_fun2 = 2 * sin(x_seq) - 2 * cos(x_seq) * (mu[2, 1] * x_seq)
      se_fun2 = 0.1 * abs(y_fun2) # Dummy SE-Werte, passen Sie diese entsprechend Ihrer Daten an
      
      # Berechnung der y-Werte und Dummy SE für die dritte Funktion
      y_fun3 = -sin(x_seq) + 2 * cos(x_seq) * (mu[1, 2] * x_seq)
      se_fun3 = 0.1 * abs(y_fun3) # Dummy SE-Werte, passen Sie diese entsprechend Ihrer Daten an
      
      # Berechnung der y-Werte und Dummy SE für die vierte Funktion
      y_fun4 = -2 * sin(x_seq) + 2 * cos(x_seq) * (mu[2, 2] * x_seq)
      se_fun4 = 0.1 * abs(y_fun4) # Dummy SE-Werte, passen Sie diese entsprechend Ihrer Daten an
      
      
      # Plot für Responsevariable y_1 mit gestrichelten SE-Linien
      GGdata = data.frame(y,x,x_seq,y_fun1,se_fun1,y_fun2,se_fun2,y_fun3,se_fun3,y_fun4,se_fun4)
      
      Zustandsvektor = apply(mod$state_probs, MARGIN = 2, FUN = which.max)

      if (Zustandsvektor[T] == 1){
        Ploty1 = ggplot(data = GGdata, aes(x, y = y[,1])) +
          geom_point(aes(color = state), col ="grey", na.rm = TRUE, alpha = 0.5) +
          geom_line(aes(x = x, y = fitted(mod$mod[[1]])[, 1], color = "State 1"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
          geom_line(aes(x = x, y = fitted(mod$mod[[2]])[, 1], color = "State 2"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
          xlab("X") + ylab("Responsevariable y_1") + ggtitle("Geschätzte Modelle Responsevariable 1") +
          scale_color_manual(values = c("State 2" = "red", "State 1" = "green", "Wahre Verteilung" = "black"), name = "Legende") +
          theme_bw() + theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0))
        
        Ploty2 = ggplot(data = GGdata, aes(x, y[,2])) + 
          geom_point(col = "grey", na.rm = TRUE, alpha = 0.5) +
          geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 2], colour = "State 1"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
          geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 2], colour = "State 2"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
          xlab("X") + ylab("Responsevariable y_2") + ggtitle("Geschätzte Modelle Responsevariable 2") +
          scale_color_manual(values = c("State 2" = "red", "State 1" = "green", "Wahre Verteilung" = "black"), name = "Legende") +
          theme_bw() + theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0))  
        
        
      }
      
      else {
        
        Ploty1 = ggplot(data = GGdata, aes(x, y = y[,1])) + 
          geom_point(aes(color = state), col ="grey", na.rm = TRUE, alpha = 0.5) +
          geom_line(aes(x = x, y = fitted(mod$mod[[1]])[, 1], color = "State 2"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
          geom_line(aes(x = x, y = fitted(mod$mod[[2]])[, 1], color = "State 1"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
          xlab("X") + ylab("Responsevariable y_1") + ggtitle("Geschätzte Modelle Responsevariable 1") +
          scale_color_manual(values = c("State 1" = "red", "State 2" = "green", "Wahre Verteilung" = "black"), name = "Legende") +
          theme_bw() + theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0))
        
        
        Ploty2 = ggplot(data = GGdata, aes(x, y[,2])) + 
          geom_point(col = "grey", na.rm = TRUE, alpha = 0.5) +
          geom_line(mapping = aes(x, fitted(mod$mod[[1]])[, 2], colour = "State 2"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
          geom_line(mapping = aes(x, fitted(mod$mod[[2]])[, 2], colour = "State 1"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
          xlab("X") + ylab("Responsevariable y_2") + ggtitle("Geschätzte Modelle Responsevariable 2") +
          scale_color_manual(values = c("State 1" = "red", "State 2" = "green", "Wahre Verteilung" = "black"), name = "Legende") +
          theme_bw() + theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0)) 
        
      }
    }
  
  # Simulations neuer Datenpunkte für jede Schleife
  
  xM[,i] = rnorm(n = T, mean = 5, sd = 1)
  for (t in 1:T) {
    yM[t, 1] = sin(xM[t,i]) - 2 * cos(xM[t,i]) * rnorm(n = 1, mean = mu[states[t], 1] * xM[t,i], sd = 4)
    yM[t, 2] = -2 * sin(xM[t,i]) + 2 * cos(xM[t,i]) * rnorm(n = 1, mean = mu[states[t], 2] * xM[t,i], sd = 4)
  }
  
  # Startwerte für den Initialstep für jede Schleife setzen
  fitted_0 = list()
  fitted_0[[1]] = cbind(mu[1, 1] * xM[,i], mu[1, 2] * xM[,i])
  fitted_0[[2]] = cbind(mu[1, 2] * xM[,i], mu[2, 2] * xM[,i])
  
  # MSVGAM schätzen für jede Schleife 
  modM = ms_vgam(x = xM[,i], y = yM, fitted_0 = fitted_0)
  
  # dataframes für die einzelnen Plots
  data_a = data.frame(
    x = xM[,i],
    a = fitted(modM$mod[[1]])[, 1])
    
  data_b = data.frame(
    x = xM[,i],
    b = fitted(modM$mod[[2]])[, 1])
  
  data_c = data.frame(
    x = xM[,i],
    c = fitted(modM$mod[[1]])[, 2])
  
  data_d = data.frame(
    x = xM[,i],
    d = fitted(modM$mod[[2]])[, 2])
  
  # lokales State decoding
  maxZustandsmatrix[,i] = apply(modM$state_probs, MARGIN = 2, FUN = which.max)


  # Plots innerhalb des Loops plotten und so in die initial Plots einzeichnen
  
  if (maxZustandsmatrix[T,i] == 1){
    
    Ploty1 = Ploty1 + 
      geom_line(data = data_a, aes(x = x, y = a, color = "State 1"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
      geom_line(data = data_b, aes(x = x, y = b, color = "State 2"), na.rm = TRUE, size = 0.7, alpha = 0.25)
    
    Ploty2 = Ploty2 +
      geom_line(data = data_c, aes(x, c, colour = "State 1"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
      geom_line(data = data_d, aes(x, d, colour = "State 2"), na.rm = TRUE, size = 0.7, alpha = 0.25)
    
    
  }
  
  else {
    Ploty1 = Ploty1 + 
      geom_line(data = data_a, aes(x = x, y = a, color = "State 2"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
      geom_line(data = data_b, aes(x = x, y = b, color = "State 1"), na.rm = TRUE, size = 0.7, alpha = 0.25)
    
    Ploty2 = Ploty2 +
      geom_line(data = data_c, aes(x, c, colour = "State 2"), na.rm = TRUE, size = 0.7, alpha = 0.25) +
      geom_line(data = data_d, aes(x, d, colour = "State 1"), na.rm = TRUE, size = 0.7, alpha = 0.25)
      
    
  }
  
  if (i == M){
    Ploty1 = Ploty1 + 
      geom_line(data = data.frame(x = x_seq, y = y_fun1 + se_fun1), aes(x = x, y = y), linetype = "dashed", color = "black") +
      geom_line(data = data.frame(x = x_seq, y = y_fun1 - se_fun1), aes(x = x, y = y), linetype = "dashed", color = "black") +
      stat_function(fun = function(x) sin(x) - 2 * cos(x) * (mu[1, 1] * x), aes(color = "Wahre Verteilung"), size = 1, linetype = "solid") +
      geom_line(data = data.frame(x = x_seq, y = y_fun2 + se_fun2), aes(x = x, y = y), linetype = "dashed", color = "black") +
      geom_line(data = data.frame(x = x_seq, y = y_fun2 - se_fun2), aes(x = x, y = y), linetype = "dashed", color = "black") +
      stat_function(fun = function(x) 2 * sin(x) - 2 * cos(x) * (mu[2, 1] * x), color = "black", size = 1, linetype = "solid")
      
    Ploty2 = Ploty2 + 
      geom_line(data = data.frame(x = x_seq, y = y_fun3 + se_fun3), aes(x = x, y = y), linetype = "dashed", color = "black") +
      geom_line(data = data.frame(x = x_seq, y = y_fun3 - se_fun3), aes(x = x, y = y), linetype = "dashed", color = "black") +
      stat_function(fun = function(x) -sin(x) + 2 * cos(x) * (mu[1, 2] * x), aes(color = "Wahre Verteilung"), size = 1, linetype = "solid") +
      geom_line(data = data.frame(x = x_seq, y = y_fun4 + se_fun4), aes(x = x, y = y), linetype = "dashed", color = "black") +
      geom_line(data = data.frame(x = x_seq, y = y_fun4 - se_fun4), aes(x = x, y = y), linetype = "dashed", color = "black") +
      stat_function(fun = function(x) -2 * sin(x) + 2 * cos(x) * (mu[2, 2] * x), color = "black", size = 1, linetype = "solid")
      
  }
  
  # lokales State Decoding nutzen, um Modellgenauigkeit (Anzahl der richtig geschätzten States) festzustellen
    
  Zustandsmatrix = Zustandsmatrix_next + modM$state_probs
  Zustandsmatrix_next = Zustandsmatrix
  
  if (i == M) {
    Zustandsmatrix = Zustandsmatrix / (M+1)
    Zustandsmatrix = apply(Zustandsmatrix, MARGIN = 2, FUN = which.max)
    anzahl_nullen = sum((states - Zustandsmatrix) == 0)
    grid.arrange(Ploty1, Ploty2, ncol = 2)
    print(paste("Dies ist die Anzahl der richtig geschätzten States. Anzahl:", anzahl_nullen))
  }
  if (i != M){print(paste("iterations", i, "of", M))}

}

#end_time = Sys.time()
#end_memory = mem_used()

## Ergebnisse des Benchmarks

#print(paste("Gesamte Rechenzeit:", as.numeric(difftime(end_time, start_time, units = "secs")), "Sekunden"))
#print(paste("Gesamter zusätzlicher Speicherverbrauch:", end_memory - start_memory, "bytes"))







### Simulation von M MSVGAMs in K Schleifen um Genaugkeit über M Simulationen K mal zuteste

#start_time = Sys.time() #Falls ein Benchmark durchgeführt wird
#start_memory = mem_used() #Falls ein Benchmark durchgeführt wird

K= 2
M = 3

yM = matrix(NA, nrow = T, ncol = 2)
xM = matrix(NA, nrow = T, ncol = M)
Zustandsmatrix_next = matrix(0, nrow = 2, ncol = T)
Genauigkeit = matrix(NA, nrow = 1, ncol = K)


for (k in 1:K) {

for (i in 1:M) {
  set.seed(k * i + M)
  
  # Modelldaten generieren
  
  xM[,i] = rnorm(n = T, mean = 5, sd = 1)
  for (t in 1:T) {
    yM[t, 1] = sin(xM[t,i]) - 2 * cos(xM[t,i]) * rnorm(n = 1, mean = mu[states[t], 1] * xM[t,i], sd = 4)
    yM[t, 2] = -2 * sin(xM[t,i]) + 2 * cos(xM[t,i]) * rnorm(n = 1, mean = mu[states[t], 2] * xM[t,i], sd = 4)
  }
  
  # Modellschätzen (MSVGAM)
  
  fitted_0 = list()
  fitted_0[[1]] = cbind(mu[1, 1] * xM[,i], mu[1, 2] * xM[,i])
  fitted_0[[2]] = cbind(mu[1, 2] * xM[,i], mu[2, 2] * xM[,i])
  
  modM = ms_vgam(x = xM[,i], y = yM, fitted_0 = fitted_0)

  # Durchschnittsbildung der Stateprobs über unsere M Simulationen
  
  Zustandsmatrix = Zustandsmatrix_next + modM$state_probs
  Zustandsmatrix_next = Zustandsmatrix
  
  # max Schritt des lokalen States Decoding über den Durchschnitt der Stateprobs
  
  if (i == M) {
    Zustandsmatrix = Zustandsmatrix / (M+1)
    Zustandsmatrix = apply(Zustandsmatrix, MARGIN = 2, FUN = which.max)
    anzahl_nullen = sum((states - Zustandsmatrix) == 0)
  }
  print(paste("Iteration", i, "von", M))
}
  
  # Genauigkeitsmatrix über die K loops
  
  Genauigkeit[,k] = anzahl_nullen
  
  yM = matrix(NA, nrow = T, ncol = 2)
  xM = matrix(NA, nrow = T, ncol = M)
  Zustandsmatrix_next = matrix(0, nrow = 2, ncol = T)
  
  if (k == K) {

      print(paste("Dies ist die Anzahl der richtig geschätzten States des jeweiligen Durchlaufs"))
      print(Genauigkeit)

  }
  if (k != K){print(paste("Schleife", k, "von", K))}

}


#end_time = Sys.time()
#end_memory = mem_used()



