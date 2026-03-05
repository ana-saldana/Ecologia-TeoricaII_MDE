datos <- read.csv("Nicho/Datos-HeV/HeV-survival.csv")

d.4 <- subset(datos, Temp == 4)
d.22 <- subset(datos, Temp == 22)
d.56 <- subset(datos, Temp == 56)

library(ggplot2)

ggplot(d.4) + geom_point(aes(x = Time.h, y = ln.S)) +
  geom_smooth(aes(x = Time.h, y = ln.S), method = "lm")

ggplot(d.22) + geom_point(aes(x = Time.h, y = ln.S)) +
  geom_smooth(aes(x = Time.h, y = ln.S), method = "loess")

ggplot(d.56) + geom_point(aes(x = Time.h, y = ln.S)) +
  geom_smooth(aes(x = Time.h, y = ln.S), method = "loess")

# 

# Regresiones

m4 <- nls(ln.S ~ -(rho * Time.h)^(kappa), 
          data = d.4, 
          start = list(rho = 1, kappa = 0.9),
          lower = c(0.00001, 0.1),
          upper = c(1, 1),
          algorithm =  "port")
summary(m4)

m22 <- nls(ln.S ~ -( rho * Time.h )^ (kappa), 
          data = d.22, 
          start = list(rho = 1, kappa = 0.9),
          lower = c(0.0001, 0.1),
          upper = c(1, 1),
          algorithm =  "port")
summary(m22)

m56 <- nls(ln.S ~ -( rho * Time )^ (kappa), 
          data = d.56, 
          start = list(rho = 1, kappa = 0.9),
          lower = c(0.1, 0.1),
          upper = c(105, 1.5),
          algorithm =  "port")
summary(m56)

# Identificando efecto de temperatura sobre rho y kappa

par.estim <- data.frame(rbind(coef(m4), coef(m22), coef(m56)))

par.estim$Temp <-c(4, 22, 56)

ggplot(par.estim) + geom_point(aes(x = Temp, y = log(rho))) +
  geom_smooth(aes(x = Temp, y = log(rho)), method = "lm")
ggplot(par.estim) + geom_point(aes(x = Temp, y = (kappa))) +
  geom_smooth(aes(x = Temp, y = (kappa)), method = "lm")

# Modelos para los parámetros

m.rho <- lm(log(rho) ~ Temp, data = par.estim)
m.kappa <- lm(log(kappa) ~ Temp, data = par.estim)

par(mfrow = c(2,2))
plot(m.rho)
plot(m.kappa)

rho.coef <- coef(m.rho)
kappa.coef <- coef(m.kappa)

coef.pk <- data.frame(par = c("ap", "Bp", "ak", "Bk"), 
                      valor = c(rho.coef, kappa.coef))

write.csv(coef.pk, "Nicho/Coeficientes-pk-lm.csv", row.names = F)

# Modelo completo
pars <- as.list(coef.pk$valor)
names(pars) <- coef.pk$par
mod <- nls(ln.S ~ -(exp(ap  + Bp * Temp) * Time.h)^exp(ak + Bk * Temp), 
           data = datos,
           start = pars,
           lower = c(-12, 0, -1, -0.05),
           upper = c(-3, 0.5, 0.1, 0.1),
           algorithm = "port")
summary(mod)

coef.mod <- data.frame(par = names(pars),
                       valor = coef(mod))

write.csv(coef.mod, "Nicho/Coeficientes-pk.csv")
