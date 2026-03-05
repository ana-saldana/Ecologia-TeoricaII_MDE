### Instalación de  spatstat
install.packages("spatstat")

### Paquetes adicionales

library(spatstat)
library(rgdal)
library(raster)
library(foreach)

#Abrir archivos

archivos <- list.files("Datos-ejemplos/", "tif", 
                       full.names = T, 
                       recursive = F)
r <- stack(archivos)
plot(r[[1]])
class(r) #Verificar tipo de objeto

### Funciones de formateo a spatstat

source("Funciones-spatstat/imFromStack.R")

#Uso de la función
r.im <- imFromStack(r)
class(r.im)

### Ventana de trabajo

source("Funciones-spatstat/winFromRaster.R")
w <- winFromRaster(r)
class(w)

### Puntos de ocurrencia

#Simulación
set.seed(984573)
puntos <- data.frame(coordinates(r)[sample(1:840, 200),])
puntos$x <- puntos$x + rnorm(200, 0, 0.05)
puntos$y <- puntos$y + rnorm(200, 0, 0.05)

#Formateo
puntos.ppp <- ppp(x = puntos$x,
                  y = puntos$y,
                  window = w,
                  check = F)
class(puntos.ppp)

plot(puntos.ppp)

### Análisis exploratorio

## Detección/medición de autocorrelación

K <- envelope(puntos.ppp, fun = Kest, nsim = 39)
plot(K)

### Análisis de respuesta a covariables

#Cuadrantes
Q <- pixelquad(X = puntos.ppp, W = as.owin(w))

source("Funciones-spatstat/plotQuantIntens.R") #La función

plotQuantIntens(imList = r.im, # Lista de imágenes
                noCuts = 5, # No. de cortes
                Quad = Q, # Los cuadrantes
                p.pp = puntos.ppp, # El proceso de puntos
                dir = "", # Directorio donde se guardará
                name = "Respuestas") # Nombre del archivo
# Ajuste del modelo

m1 <- ppm(Q = puntos.ppp,
          trend = ~ Var.1,
          covariates = r.im)

# Resultados
summary(m1)

### Comparación de modelos

m2 <- ppm(Q = puntos.ppp,
          trend = ~ Var.2,
          covariates = r.im)

m3 <- ppm(Q = puntos.ppp,
          trend = ~ Var.3,
          covariates = r.im)

AIC(m1); AIC(m2); AIC(m3)

### El modelo "mejor"
plot(m3, trend = T, se = F)
