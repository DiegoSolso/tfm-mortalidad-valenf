# Weibull frente a RSF (y Cox) con las variables del modelo final y las mismas particiones

library(survival)
library(randomForestSRC)

ruta <- "datos/capacidad_funcional_1VALENF_LARGO_con_constantes.csv"

variables <- c("perfil", "frec_cardiaca", "ta_sistolica", "ta_diastolica", "sat_o2",
               "tipo_ingreso", "tipoproceso")

# RSF con la mejor combinación de la búsqueda en rejilla del estudio preliminar
arboles <- 500
nodesize <- 60
mtry     <- 3


# Datos (como en modelo_perfiles.R)

datos <- read.csv(ruta, sep = ",", dec = ",", stringsAsFactors = FALSE,
                  na.strings = c("", "NA"))
datos <- datos[datos$orden_valenf == 1 & datos$ctes_4vars_en_24h == 1, ]

valoracion <- as.POSIXct(datos$fecha_valoracion)
alta       <- as.POSIXct(datos$fecha_alta)
exitus     <- as.POSIXct(datos$fecha_exitus)
murio      <- as.logical(datos$hay_exitus_episodio)
dia <- as.numeric(as.POSIXct(format(exitus, "%Y-%m-%d")))
ini <- pmax(as.numeric(valoracion), dia)
fin <- ifelse(murio, (ini + dia + 86400) / 2, as.numeric(alta))
datos$tiempo <- (fin - as.numeric(valoracion)) / 86400
datos$evento <- as.integer(murio)
datos <- datos[datos$evento == 1 | datos$tiempo > 0, ]

categorica <- function(x, niveles) factor(x, levels = niveles, labels = paste0("=", niveles))
datos$p1 <- categorica(datos$respuesta_pregunta1, c(15, 10, 5, 0))
datos$p2 <- categorica(datos$respuesta_pregunta2, 4:1)
datos$p3 <- categorica(datos$respuesta_pregunta3, 4:1)
datos$p4 <- categorica(datos$respuesta_pregunta4, 4:1)
datos$p5 <- categorica(datos$respuesta_pregunta5, 2:1)
datos$p6 <- categorica(datos$respuesta_pregunta6, 0:1)
datos$p7 <- categorica(datos$respuesta_pregunta7, 0:1)
datos$sexo         <- categorica(datos$sexo, c("H", "M"))
datos$perfil       <- categorica(datos$primer_perfil, 1:5)
datos$tipo_ingreso <- categorica(datos$tipo_ingreso, c("P", "U"))
datos$tipoproceso  <- categorica(datos$tipoproceso, c("M", "Q"))

todas <- c(paste0("p", 1:7), "edad", "sexo", "valor_charlson", "valor_cf", "valor_rLPP", "valor_rC",
           "perfil", "frec_cardiaca", "ta_sistolica", "ta_diastolica", "sat_o2", "tipo_ingreso", "tipoproceso")
datos <- datos[complete.cases(datos[, todas]), ]
datos <- datos[, c("tiempo", "evento", variables)]

cat("N =", nrow(datos), " exitus =", sum(datos$evento), "\n")
cat("Variables:", variables, "\n")


# 1. Ajuste de los tres modelos en cada partición

formula_modelo <- as.formula(paste("Surv(tiempo, evento) ~", paste(variables, collapse = " + ")))

# riesgo acumulado de cada paciente en su tiempo t (H: pacientes x tiempos)
acumulado_en <- function(tiempos, H, t) {
  j <- findInterval(t, tiempos)
  ifelse(j == 0, 0, H[cbind(seq_along(t), pmax(j, 1))])
}

set.seed(1)
particion <- sample(rep(1:5, length.out = nrow(datos)))

modelos <- c("Weibull", "Cox", "RSF")
riesgo <- H7 <- H14 <- matrix(NA, nrow(datos), 3, dimnames = list(NULL, modelos))
indice_c <- matrix(NA, 5, 3, dimnames = list(paste("particion", 1:5), modelos))

for (k in 1:5) {
  cat("Particion", k, "de 5...\n")
  ent <- datos[particion != k, ]
  pru <- datos[particion == k, ]
  i   <- which(particion == k)
  t7  <- pmin(pru$tiempo, 7)
  t14 <- pmin(pru$tiempo, 14)

  # Weibull
  w  <- survreg(formula_modelo, data = ent, dist = "weibull", control = survreg.control(maxiter = 100))
  mu <- predict(w, newdata = pru, type = "lp")
  riesgo[i, "Weibull"] <- -mu
  H7[i, "Weibull"]  <- (t7 / exp(mu))^(1 / w$scale)
  H14[i, "Weibull"] <- (t14 / exp(mu))^(1 / w$scale)

  # Cox
  cx <- coxph(formula_modelo, data = ent)
  sf <- survfit(cx, newdata = pru, censor = FALSE)
  Hc <- t(-log(sf$surv))
  riesgo[i, "Cox"] <- predict(cx, newdata = pru, type = "lp")
  H7[i, "Cox"]  <- acumulado_en(sf$time, Hc, t7)
  H14[i, "Cox"] <- acumulado_en(sf$time, Hc, t14)

  # RSF
  rf <- rfsrc(formula_modelo, data = ent, ntree = arboles, nodesize = nodesize, mtry = mtry,
              splitrule = "logrank", samptype = "swor", seed = -k)
  pr <- predict(rf, newdata = pru)
  riesgo[i, "RSF"] <- pr$predicted
  H7[i, "RSF"]  <- acumulado_en(pr$time.interest, pr$chf, t7)
  H14[i, "RSF"] <- acumulado_en(pr$time.interest, pr$chf, t14)

  for (m in modelos) {
    r <- riesgo[i, m]
    indice_c[k, m] <- concordance(Surv(pru$tiempo, pru$evento) ~ r, reverse = TRUE)$concordance
  }
}


# 2. Índice C

cat("\n\n INDICE C CON DATOS NUEVOS (por particion) \n")
print(round(indice_c, 3))
cat("\nMedia:      ", sprintf("%-9s %.3f", modelos, colMeans(indice_c)), "\n")
cat("Desviacion: ", sprintf("%-9s %.3f", modelos, apply(indice_c, 2, sd)), "\n")
for (m in c("Weibull", "Cox"))
  cat(m, "supera al RSF en", sum(indice_c[, m] > indice_c[, "RSF"]), "de 5 particiones",
      " (diferencia media", sprintf("%+.3f", mean(indice_c[, m] - indice_c[, "RSF"])), ")\n")


# 3. Muertes observadas frente a esperadas a los 7 y 14 días

grupos <- list("Todos" = rep(TRUE, nrow(datos)))
for (p in levels(datos$perfil)) grupos[[paste("Perfil", sub("=", "", p))]] <- datos$perfil == p

tabla_oe <- function(d, H) do.call(rbind, lapply(names(grupos), function(g) {
  i <- grupos[[g]]
  O <- sum(datos$evento[i] == 1 & datos$tiempo[i] <= d)
  do.call(rbind, lapply(modelos, function(m) {
    E  <- sum(H[i, m])
    ic <- poisson.test(O, E)$conf.int
    data.frame(grupo = g, modelo = m, reales = O, esperadas = round(E, 1),
               O_E = round(O / E, 2), IC_inf = round(ic[1], 2), IC_sup = round(ic[2], 2))
  }))
}))
oe <- list("7" = tabla_oe(7, H7), "14" = tabla_oe(14, H14))

cat("\n\n REALIDAD FRENTE A MODELO \n")
for (d in names(oe)) {
  cat("\n Muertes en los primeros", d, "dias\n")
  print(oe[[d]], row.names = FALSE)
}


# Gráficos

colores <- c(Weibull = "#2a78d6", Cox = "grey45", RSF = "#eb6834"); tinta <- "grey20"
carpeta <- "graficos_comparacion"
dir.create(carpeta, showWarnings = FALSE)

# índice C de cada modelo en cada partición
png(file.path(carpeta, "1_indice_c.png"), width = 900, height = 650, res = 120)
par(mar = c(4, 4.5, 3.5, 1))
rango <- range(indice_c) + c(-0.01, 0.01)
plot(NA, xlim = c(0.6, 3.4), ylim = rango, xaxt = "n", xlab = "", ylab = "indice C con datos nuevos",
     main = "Capacidad de ordenar pacientes en cada particion")
axis(1, 1:3, modelos)
for (k in 1:5) lines(1:3, indice_c[k, ], col = "grey80")
for (j in 1:3) {
  points(rep(j, 5), indice_c[, j], pch = 19, col = colores[j])
  segments(j - 0.18, mean(indice_c[, j]), j + 0.18, mean(indice_c[, j]), lwd = 3, col = colores[j])
  text(j + 0.2, mean(indice_c[, j]), sprintf("%.3f", mean(indice_c[, j])), pos = 4, cex = 0.85, col = tinta)
}
mtext("punto = una particion; raya = media", side = 3, line = 0.3, cex = 0.8, col = "grey40")
invisible(dev.off())

# O/E de cada modelo por perfil
png(file.path(carpeta, "2_O_E.png"), width = 1300, height = 650, res = 120)
par(mfrow = c(1, 2), mar = c(4.5, 7, 4, 1))
filas <- length(grupos):1
for (d in names(oe)) {
  r <- oe[[d]]
  tope <- max(r$IC_sup) * 1.08
  plot(NA, xlim = c(0, tope), ylim = c(0.5, length(grupos) + 0.5), yaxt = "n", ylab = "",
       xlab = "O/E = muertes reales / esperadas", main = paste("Primeros", d, "dias"))
  axis(2, filas, names(grupos), las = 1)
  abline(h = length(grupos) - 0.5, col = "grey85")
  abline(v = 1, lty = 2, col = "red")
  for (j in 1:3) {
    y <- filas + c(0.22, 0, -0.22)[j]
    rj <- r[r$modelo == modelos[j], ]
    segments(rj$IC_inf, y, rj$IC_sup, y, lwd = 2, col = colores[j])
    points(rj$O_E, y, pch = 19, col = colores[j])
  }
  if (d == "7") legend("bottomright", modelos, col = colores, pch = 19, lwd = 2, bty = "n", cex = 0.85)
}
invisible(dev.off())

cat("\nGraficos en la carpeta '", carpeta, "'\n", sep = "")
