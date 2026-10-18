# Elección del modelo: riesgos proporcionales, criterio de la recta y otras distribuciones

library(survival)

ruta    <- "datos/capacidad_funcional_1VALENF_LARGO_con_constantes.csv"
carpeta <- "graficos_eleccion"
dir.create(carpeta, showWarnings = FALSE)


# Datos (como en modelo_perfiles.R)

datos <- read.csv(ruta, sep = ",", dec = ",", stringsAsFactors = FALSE, na.strings = c("", "NA"))
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
datos$perfil       <- categorica(datos$primer_perfil, c(5, 1, 2, 3, 4))   # referencia: perfil 5
datos$tipo_ingreso <- categorica(datos$tipo_ingreso, c("P", "U"))
datos$tipoproceso  <- categorica(datos$tipoproceso, c("M", "Q"))
candidatas <- c(paste0("p", 1:7), "edad", "sexo", "valor_charlson", "valor_cf", "valor_rLPP",
                "valor_rC", "perfil", "frec_cardiaca", "ta_sistolica", "ta_diastolica", "sat_o2",
                "tipo_ingreso", "tipoproceso")
datos <- datos[complete.cases(datos[, candidatas]), ]

variables <- c("perfil", "frec_cardiaca", "ta_sistolica", "ta_diastolica", "sat_o2",
               "tipo_ingreso", "tipoproceso")
derecha   <- paste(variables, collapse = " + ")
cat("N =", nrow(datos), " exitus =", sum(datos$evento), "\n")

ajustar <- function(v = variables, d = datos, dist = "weibull") {
  survreg(as.formula(paste("Surv(tiempo, evento) ~", paste(v, collapse = " + "))), data = d,
          dist = dist, control = survreg.control(maxiter = 100))
}
indiceC <- function(lp, d = datos) {
  c <- concordance(Surv(d$tiempo, d$evento) ~ lp)$concordance
  if (c < 0.5) 1 - c else c
}
set.seed(1)
particion <- sample(rep(1:5, length.out = nrow(datos)))   # mismas particiones que modelo_perfiles.R


# 1. Riesgos proporcionales

cox_perfil <- coxph(Surv(tiempo, evento) ~ perfil, data = datos)
cox_final  <- coxph(as.formula(paste("Surv(tiempo, evento) ~", derecha)), data = datos)
cat("\n RIESGOS PROPORCIONALES (cox.zph) \n")
cat("Solo el perfil:\n");      print(cox.zph(cox_perfil)$table)
cat("Modelo completo:\n");     print(cox.zph(cox_final)$table)


# 2. Criterio de la recta: log(-log S) frente a log t

km     <- survfit(Surv(tiempo, evento) ~ 1, data = datos)
puntos <- data.frame(t = km$time, S = km$surv)[km$n.event > 0 & km$surv > 0, ]
puntos$x <- log(puntos$t); puntos$y <- log(-log(puntos$S))
recta  <- lm(y ~ x, data = puntos)
curva  <- lm(y ~ x + I(x^2), data = puntos)
cat("\n CRITERIO DE LA RECTA \n")
cat("Puntos:", nrow(puntos), " R2 =", round(summary(recta)$r.squared, 4),
    " pendiente (alfa) =", round(coef(recta)[2], 3),
    " termino cuadratico =", round(coef(curva)[3], 4),
    " p =", signif(summary(curva)$coefficients[3, 4], 2), "\n")


# 3. Otras distribuciones con las variables del modelo final

c_por_particion <- function(dist) sapply(1:5, function(k) {
  m <- ajustar(d = datos[particion != k, ], dist = dist)
  indiceC(predict(m, newdata = datos[particion == k, ], type = "lp"), datos[particion == k, ])
})
c_cox_particion <- sapply(1:5, function(k) {
  m <- coxph(as.formula(paste("Surv(tiempo, evento) ~", derecha)), data = datos[particion != k, ])
  indiceC(predict(m, newdata = datos[particion == k, ], type = "lp"), datos[particion == k, ])
})

cat("\n OTRAS DISTRIBUCIONES \n")
familias <- c(Exponencial = "exponential", Weibull = "weibull", Lognormal = "lognormal",
              Loglogistica = "loglogistic")
comparacion <- do.call(rbind, lapply(names(familias), function(f) {
  m <- ajustar(dist = familias[f])
  data.frame(familia = f, AIC = round(AIC(m), 1), C_vc = round(mean(c_por_particion(familias[f])), 3))
}))
comparacion$dAIC <- comparacion$AIC - min(comparacion$AIC)
print(comparacion[order(comparacion$AIC), ], row.names = FALSE)
cat("Cox: C con validacion cruzada", round(mean(c_cox_particion), 3), "\n")

modelo <- ajustar()
expo   <- ajustar(dist = "exponential")
x <- 2 * (as.numeric(logLik(modelo)) - as.numeric(logLik(expo)))
cat("Weibull frente a exponencial: chi2 =", round(x, 1), " gl = 1  p =",
    signif(pchisq(x, 1, lower.tail = FALSE), 3), "\n")


# Gráfico: recta de Weibull con todos los pacientes y según el perfil

coma    <- function(x, d = 1) formatC(x, format = "f", digits = d, decimal.mark = ",")
colores <- c("#2b83ba", "#5ab4ac", "#e0a526", "#e0592b", "#9e1b32")   # perfil 1 a 5
tinta   <- "grey15"; azul <- "#2b6cb0"

pdf(file.path(carpeta, "recta_weibull.pdf"), width = 7.2, height = 3.8)
par(mfrow = c(1, 2), mar = c(4.3, 4.6, 2.2, 1), mgp = c(2.8, 0.7, 0), las = 1, cex.axis = 0.85)
ejes <- function() {
  marcas <- c(0.5, 1, 2, 5, 10, 20, 50)
  axis(1, at = log(marcas), labels = sapply(marcas, function(m) coma(m, ifelse(m < 1, 1, 0))))
  abline(v = log(marcas), col = "grey92")
}
plot(puntos$x, puntos$y, type = "n", xlab = "Días (escala logarítmica)",
     ylab = expression(log(-log~hat(S)(t))), xaxt = "n", bty = "l", main = "(a) Todos los pacientes",
     cex.main = 0.95, font.main = 1)
ejes(); abline(h = pretty(puntos$y), col = "grey92")
points(puntos$x, puntos$y, pch = 16, cex = 0.45, col = tinta)
abline(recta, col = azul, lwd = 2)
legend("bottomright", c("Kaplan-Meier", sprintf("Recta (pendiente %s)", coma(coef(recta)[2], 2))),
       pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), col = c(tinta, azul), bty = "n", cex = 0.8)

km_p <- survfit(Surv(tiempo, evento) ~ primer_perfil, data = datos)
tramos <- rep(1:5, km_p$strata)
lm_p <- data.frame(t = km_p$time, S = km_p$surv, e = km_p$n.event, g = tramos)
lm_p <- lm_p[lm_p$e > 0 & lm_p$S < 1 & lm_p$S > 0, ]
plot(log(lm_p$t), log(-log(lm_p$S)), type = "n", xlab = "Días (escala logarítmica)",
     ylab = expression(log(-log~hat(S)(t))), xaxt = "n", bty = "l", main = "(b) Según el perfil",
     cex.main = 0.95, font.main = 1)
ejes(); abline(h = pretty(log(-log(lm_p$S))), col = "grey92")
for (g in 1:5) {
  s <- lm_p[lm_p$g == g, ]
  lines(log(s$t), log(-log(s$S)), type = "s", col = colores[g], lwd = 1.8)
}
legend("topleft", paste("Perfil", 1:5), col = colores, lwd = 2.5, seg.len = 2.5, bty = "n", cex = 0.8)
invisible(dev.off())

cat("\nGrafico guardado en", file.path(carpeta, "recta_weibull.pdf"), "\n")
