#Práctica 2
install.packages("spatstat")
library(spatstat)
library(terra)
library(foreach)

#Subir archivos 
archivos <- list.files(".", "Var",
                       full.names = T,
                       recursive = F)
r <- rast(archivos)
plot(r[[1]])
class(r)

#funciones de formateo
source("imFromStack.R")

#uso de la función
r.im <- imFromStack(r)
class(r.im) #lista 

#imagen binaria 
source("winFromRaster.R")
w <- winFromRaster(r)
class(w) #tipo de objeto 
plot(w)

#simulación 
set.seed(984573)
puntos <- data.frame(crds(r)[sample(1:840, 200),])
puntos$x <- puntos$x + rnorm(200, 0, 0.05)
puntos$y <- puntos$y + rnorm(200, 0, 0.05) #200 puntos de presencia 
plot(r$`Var-1`)
points(puntos, col = "black", pch = 20)

#formateo :Transformar los puntos a un PP
#Que puntos se usan y cuales no: coordenadas, ventana de trabajo 
#y checar si los puntos estan dentro
puntos.ppp <- ppp(x = puntos$x,
                  y = puntos$y,
                  window = w,
                  check = F)
class(puntos.ppp)
plot(puntos.ppp)

#Análisis exploratorio---------------------------------------------------------
#keystimate y nsim:numero de simulaciones 
K <- envelope(puntos.ppp, fun = Kest, nsim = 39)
plot(K) #es aleatorio 

#Análisis de respuesta o covariables 
#cuadrantes 
Q <- pixelquad(X = puntos.ppp, W = as.owin(w))
#patron de cuadrantura: cuenta en numero de puntos por pixel 

source("plotQuantIntens.R")
plotQuantIntens(imList = r.im, #lista de imagenes 
                 noCuts = 5, #numero de cuartiles cortes
                 Quad = Q, #los cuadrantes
                 p.pp = puntos.ppp, #patron de puntos
                 dir = "", #donde lo guarda 
                 name = "Respuestas") #nombre

#detectar como el patron de puntos como responde a las covariables 

#Ajuste del modelo
names(r.im) <- c("Var.1", "Var.2", "Var.3")
m1 <- ppm(Q = puntos.ppp,
          trend = ~ Var.1,
          covariates = r.im)

summary(m1)
par(mar = c(1,1,1,1))
plot()

#reflexión en fomato de abstract sobre modelo correlativo y mecanisticos 

m2 <- ppm(Q = puntos.ppp,
          trend = ~ Var.2,
          covariates = r.im)

summary(m2)
