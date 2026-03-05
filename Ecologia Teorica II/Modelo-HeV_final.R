#Barreras térmicas: Modelado mecanístico de la supervivencia y transmisión del 
#Hendra virus (HeV) basado en la termotolerancia

#PARTE 1: CURVAS DE SUPERVIVIENCIA VIRUS----------------------------------------
#Instalar paquetería y librerías 
install.packages("doParallel")
install.packages("deSolve")
install.packages("terra")
install.packages("foreach")
library(doParallel)
library(deSolve)
library(terra)
library(ggplot2)
library(lattice)
library(foreach) 
library(tidyr)

#Bases de datos
datos <- read.csv("HeV-survival.csv")

d.4 <- subset(datos, Temp == 4)
d.22 <- subset(datos, Temp == 22)
d.56 <- subset(datos, Temp == 56)

#Gráficos 

ggplot(d.4) + geom_point(aes(x = Time.h, y = ln.S)) +
  geom_smooth(aes(x = Time.h, y = ln.S), method = "lm")

ggplot(d.22) + geom_point(aes(x = Time.h, y = ln.S)) +
  geom_smooth(aes(x = Time.h, y = ln.S), method = "loess")

ggplot(d.56) + geom_point(aes(x = Time.h, y = ln.S)) +
  geom_smooth(aes(x = Time.h, y = ln.S), method = "loess")

# Weibull
m4 <- nls(ln.S ~ -( rho * Time.h )^ (kappa), 
          data = d.4, 
          start = list(rho = 1, kappa = 0.9), #definir parámetros 
          lower = c(0.00001, 0.1), 
          upper = c(1, 1),
          algorithm =  "port")

m22 <- nls(ln.S ~ -( rho * Time.h )^ (kappa), 
           data = d.22, 
           start = list(rho = 1, kappa = 0.9),
           lower = c(0.0001, 0.1),
           upper = c(1, 1),
           algorithm =  "port")

m56 <- nls(ln.S ~ -( rho * Time.h )^ (kappa), 
           data = d.56, 
           start = list(rho = 1, kappa = 0.9),
           lower = c(0.1, 0.1),
           upper = c(8000, 1.5),
           algorithm =  "port")

summary (m4)
summary(m22)
summary(m56)

#Gráficar cambios en los coeficientes respecto a la temperatura 
par.estim <- data.frame(rbind(coef(m4), coef(m22), coef(m56)))
par.estim$Temp <-c(4, 22, 56)

ggplot(par.estim) + geom_point(aes(x = Temp, y = log(rho))) +
  geom_smooth(aes(x = Temp, y = log(rho)), method = "lm")

ggplot(par.estim) + geom_point(aes(x = Temp, y = log(kappa))) +
  geom_smooth(aes(x = Temp, y = log(kappa)), method = "lm")


m.rho <- lm(log(rho) ~ Temp, data = par.estim)
m.kappa <- lm(log(kappa) ~ Temp + Temp, data = par.estim)

rho.coef <- coef(m.rho) 
kappa.coef <- coef(m.kappa)

coef.pk <- data.frame(par = c("ap", "Bp", "ak", "Bk"), 
                      valor = c(rho.coef, kappa.coef)) #condiciones iniciales 

#Ajustar el modelo 
pars <- as.list(coef.pk$valor)
names(pars) <- coef.pk$par
mod <- nls(ln.S ~ -(exp(ap  + Bp * Temp) * Time.h)^exp(ak + Bk * Temp), 
           data = datos,
           start = pars,
           lower = c(-12, 0, -1, -0.05),
           upper = c(-3, 0.5, 0.1, 0.1),
           algorithm = "port")
summary(mod)

#Renombrar coeficientes 
coef.mod <- data.frame(par = names(pars),
                       valor = coef(mod))

write.csv(coef.mod, "Coeficientes.csv")  

#Nueva base de datos 
datos.nuevos <- expand.grid(Tiempo = seq(0, 12, len = 50),
                            Temp = seq(4, 56, len = 50))

knitr::kable(head(datos.nuevos))

#Función que hará los calculos: 
weib <- function(Tiempo = NA, Temp = NA, pars = NA){
  ap <- pars$ap
  Bp <- pars$Bp
  ak <- pars$ak
  Bk <- pars$Bk
  
  p <- exp(ap + Bp * Temp)
  k <- exp(ak + Bk * Temp)
  
  S <- exp(- (p * Tiempo) ^ k)
  return(S)
}

#Gráfico: Efecto de la temperatura

pars.1 <- as.list(coef(mod))

Sup <- weib(Tiempo = datos.nuevos$Tiempo,
            Temp = datos.nuevos$Temp,
            pars = pars.1)

datos.nuevos$Sup <- Sup #comlumna de supervivencia 

#Gráficar 
wireframe(Sup ~ Tiempo + Temp, 
          data = datos.nuevos,
          drape = T,
          screen = list(z = -135, x = -70, y = 3))

#Mundo Geográfico
#Simulaciones espaciales con el modelo de supervivencia

#1) ESCENARIO SIMPLE 
bio1 <- rast("Bio1.tif")
plot(bio1) 

bio1.df <- as.data.frame(bio1, xy = T) 

pars <- as.list(coef.pk$valor)
names(pars) <- coef.pk$par

Sup.bio1 <- weib(Tiempo = 12,
                 Temp = bio1.df$Bio1,
                 pars = pars)

#Añadir las predicciones de supervivencia a bio1.df y transformar a raster
bio1.df$Sup <- Sup.bio1 
Sup.r <- rast(bio1.df[, -3]) 
plot(Sup.r)

#2) ESCENARIO REALISTA 
pars.df <- read.csv("Coeficientes.csv")

weib.ode <- function(t, y, params){
  with(params, {
    S <- y
    
    Temp <- Tmax - Tdif * cos(pi * t/24)^2
    
    rho <- exp(ap + Bp * Temp)
    
    kappa <- exp(ak + Bk * Temp)
    
    dS <- - rho * kappa * (-log(S))^(1 - 1/kappa) * S
    
    list(c(dS))
  })
}

Tmax <- rast("Tmax-01.tif")
Tmin <- rast("Tmin-01.tif")

Tdif <- Tmax - Tmin

Tmax.df <- as.data.frame(Tmax, xy = T) 
Tdif.df <- as.data.frame(Tdif, xy = F)

Temps.df <- data.frame(Tmax.df, Tdif = Tdif.df$`Tmax-01`)
names(Temps.df) <- c("x", "y", "Tmax", "Tdif")

y <- 0.9999 
t <- seq(0, 12, length(100))

#simulaciones 
registerDoParallel(cores = 2)

pars <- as.list(pars.df$valor)
names(pars) <- pars.df$par


sims <- foreach(i = 1:nrow(Temps.df), .combine = c,
                .packages = "deSolve") %dopar% {
                  params <- pars
                  params$Tmax <- Temps.df$Tmax[i]
                  params$Tdif <- Temps.df$Tdif[i]
                  out <- lsoda(y = y, times = t,
                               func = weib.ode,
                               parms = params)
                  return(out[nrow(out), 2])
                }

#Convertir a raster          
Sup.fluc <- data.frame(Temps.df[, c("x", "y")], Sup = sims)
Sup.fluc.r <- rast(Sup.fluc) 
plot(Sup.fluc.r)

#PARTE 2: MODELO DE LA DISTRIBUCIÓN DEL VIRUS-----------------------------------
#Subir archivos 
archivos <- list.files(".", "Var",
                       full.names = T,
                       recursive = F)
r <- rast(archivos)
plot(r[[1]])
class(r)

#Funciones de formateo
source("imFromStack.R")
r.im <- imFromStack(r)
class(r.im) 

#Imagen binaria 
source("winFromRaster.R")
w <- winFromRaster(r)
class(w) 
plot(w)

#Simulación 
set.seed(984573)
puntos <- data.frame(crds(r)[sample(1:840, 200),])
puntos$x <- puntos$x + rnorm(200, 0, 0.05)
puntos$y <- puntos$y + rnorm(200, 0, 0.05) 
plot(r$`Var-1`)
points(puntos, col = "black", pch = 20)

#Formateo :Transformar los puntos a un PP
puntos.ppp <- ppp(x = puntos$x,
                  y = puntos$y,
                  window = w,
                  check = F)
class(puntos.ppp)
plot(puntos.ppp)

#Análisis exploratorio
K <- envelope(puntos.ppp, fun = Kest, nsim = 39)
plot(K) 

#Análisis de respuesta o covariables 
Q <- pixelquad(X = puntos.ppp, W = as.owin(w)) 

source("plotQuantIntens.R")
plotQuantIntens(imList = r.im, 
                noCuts = 5, 
                Quad = Q, 
                p.pp = puntos.ppp, 
                dir = "",  
                name = "Respuestas") 

#Ajuste del modelo
names(r.im) <- c("Var.1", "Var.2", "Var.3")
m1 <- ppm(Q = puntos.ppp,
          trend = ~ Var.1,
          covariates = r.im)

summary(m1)
par(mar = c(1,1,1,1))
plot()

m2 <- ppm(Q = puntos.ppp,
          trend = ~ Var.2,
          covariates = r.im)

summary(m2)

#FIN :)