#PRÁCTICA I: CURVAS DE SUPERVIVIENCIA VIRUS 
install.packages("doParallel")
install.packages("deSolve")
install.packages("terra")

library(doParallel)
library(deSolve)
library(terra)

#Bases de datos
datos <- read.csv("HeV-survival.csv")

d.4 <- subset(datos, Temp == 4)
d.22 <- subset(datos, Temp == 22)
d.56 <- subset(datos, Temp == 56)

#Gráficos 
library(ggplot2)

ggplot(d.4) + geom_point(aes(x = Time.h, y = ln.S)) +
  geom_smooth(aes(x = Time.h, y = ln.S), method = "lm")

ggplot(d.22) + geom_point(aes(x = Time.h, y = ln.S)) +
  geom_smooth(aes(x = Time.h, y = ln.S), method = "loess")

ggplot(d.56) + geom_point(aes(x = Time.h, y = ln.S)) +
  geom_smooth(aes(x = Time.h, y = ln.S), method = "loess")

# Weibull 
# S( t ) = exp( − ( ρ t)k)
# log ( t ) = − ( ρ t)k
#Ajustar modelos respecto a sus mínimos y máximos 

m4 <- nls(ln.S ~ -( rho * Time.h )^ (kappa), #tipo 2 k=1
          data = d.4, 
          start = list(rho = 1, kappa = 0.9), #definimos cuánto valen los parametros 
          lower = c(0.00001, 0.1), #para rho y k no busques valores menores a...
          upper = c(1, 1), #para rho y k no busques valores mayores a... 
          algorithm =  "port")

m22 <- nls(ln.S ~ -( rho * Time.h )^ (kappa), 
           data = d.22, 
           start = list(rho = 1, kappa = 0.9),#tipo 1 k<1
           lower = c(0.0001, 0.1),
           upper = c(1, 1),
           algorithm =  "port")

m56 <- nls(ln.S ~ -( rho * Time.h )^ (kappa), 
           data = d.56, 
           start = list(rho = 1, kappa = 0.9),#tipo 1 k<1
           lower = c(0.1, 0.1),
           upper = c(8000, 1.5),
           algorithm =  "port")

summary (m4)
summary(m22)
summary(m56)

#Gráficar cómo cambian los coeficientes respecto a la temperatura 
#La escala logaritmica nos ayuda a hacer valores positivos 
#Kappa y rho siempre positivos

par.estim <- data.frame(rbind(coef(m4), coef(m22), coef(m56)))

par.estim$Temp <-c(4, 22, 56)

#rho
ggplot(par.estim) + geom_point(aes(x = Temp, y = log(rho))) +
  geom_smooth(aes(x = Temp, y = log(rho)), method = "lm")

#kappa
#conforme aumenta la tempertarura se vuelve más tipo 1 
#inactivación esta más en los primeros momentois de activación 
#menor temp, mayor activación 
ggplot(par.estim) + geom_point(aes(x = Temp, y = log(kappa))) +
  geom_smooth(aes(x = Temp, y = log(kappa)), method = "lm")

#
m.rho <- lm(log(rho) ~ Temp, data = par.estim)
m.kappa <- lm(log(kappa) ~ Temp + Temp, data = par.estim)

rho.coef <- coef(m.rho) #si la temperatura aumenta 1°, Rho aumenta 0.28 unidades 
kappa.coef <- coef(m.kappa)

coef.pk <- data.frame(par = c("ap", "Bp", "ak", "Bk"), 
                      valor = c(rho.coef, kappa.coef)) #condiciones iniciales 

#Ajustando el modelo 
pars <- as.list(coef.pk$valor)
names(pars) <- coef.pk$par
mod <- nls(ln.S ~ -(exp(ap  + Bp * Temp) * Time.h)^exp(ak + Bk * Temp), 
           data = datos,
           start = pars,
           lower = c(-12, 0, -1, -0.05),
           upper = c(-3, 0.5, 0.1, 0.1),
           algorithm = "port")
summary(mod)

#siempre es más facil encontrar rangos en valores logaritmicos 
#si tengo valores muy grandes o pequeños es mejor transformarlos a valores log
#¿Qué pasa si mis valores son negativos????


coef.mod <- data.frame(par = names(pars),
                       valor = coef(mod))

write.csv(coef.mod, "Coeficientes.csv") #guardar coeficientes 

#Nueva base de datos 
#expand.grid: divide el intervalo en x(len) elementos 
datos.nuevos <- expand.grid(Tiempo = seq(0, 12, len = 50),
                            Temp = seq(4, 56, len = 50))

knitr::kable(head(datos.nuevos))

#Función que hará los calculos 
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

#gráfico que nos permita ver dicho efecto de la temperatura:

pars.1 <- as.list(coef(mod)) #parametros creados 

Sup <- weib(Tiempo = datos.nuevos$Tiempo,
            Temp = datos.nuevos$Temp,
            pars = pars.1)

datos.nuevos$Sup <- Sup #agregamos la comlumna de supervivencia 

#Gráficar :)))) 
library(lattice)

#Lado     izq    #derecha
wireframe(Sup ~ Tiempo + Temp, 
          data = datos.nuevos,
          drape = T,
          screen = list(z = -135, x = -70, y = 3))

#Modelo exponencial:mecanistico o correlativo puede ser 

#Mundo Geográfico--------------------------------------
#Simulaciones espaciales con el modelo de supervivencia


#ESCENARIO SIMPLE <<<<<<<
bio1 <- rast("Bio1.tif")
plot(bio1) #recortamos la distribución a la del murcielago (vector)

bio1.df <- as.data.frame(bio1, xy = T) #tempertaura 

pars <- as.list(coef.pk$valor)
names(pars) <- coef.pk$par

Sup.bio1 <- weib(Tiempo = 12,
                 Temp = bio1.df$Bio1,
                 pars = pars)

#Añadir las predicciones de supervivencia a bio1.df y transformar a raster:
bio1.df$Sup <- Sup.bio1 #primeras dos colunmnas son las coordenaDAS (X,y)
Sup.r <- rast(bio1.df[, -3]) #capa raster 
plot(Sup.r)

#ESCENARIO REALISTA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#Murcis tiene buna memoria: recuerdan todos los arboles que visitan :0
#cos (0) = 1
#cos (180) = -1
# cos va de -1 a 1 -> solo valores + -> || o al cuadrado (representar fluctuación de T°) 

# Tmin * Cos (pi*t)
# Tmin + TdifCs2(pi*t)

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

Tmax.df <- as.data.frame(Tmax, xy = T) #(capa, xy = T es reusltado incluye coordenadas de cada pixel)
Tdif.df <- as.data.frame(Tdif, xy = F)

Temps.df <- data.frame(Tmax.df, Tdif = Tdif.df$`Tmax-01`)
names(Temps.df) <- c("x", "y", "Tmax", "Tdif")

#condiciones 
#1/0 -> infinito 

y <- 0.9999 #iniciales 
t <- seq(0, 12, length(100))

#simulaciones 
install.packages("foreach")
library(foreach) #usa bucles for más efectivamente 
library(tidyr)
library(doParallel) #extensión de foreach 
library(terra)

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

#Convertir a capa raster          coordenadas   supervivencia
Sup.fluc <- data.frame(Temps.df[, c("x", "y")], Sup = sims)
Sup.fluc.r <- rast(Sup.fluc) #rast para transformar a raster 
plot(Sup.fluc.r)

#Limitacion: se le asigna el mismo peso a todos los pixeles 
#

#FIN :) 

