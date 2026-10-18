# Modelo loglineal de Weibull según el perfil de la primera valoración

library(survival)

ruta <- "datos/capacidad_funcional_1VALENF_LARGO_con_constantes.csv"


# Variables del modelo (TRUE = entra)

usar <- c(
  p1 = FALSE, p2 = FALSE, p3 = FALSE, p4 = FALSE,     # preguntas de VALENF
  p5 = FALSE, p6 = FALSE, p7 = FALSE,
  edad           = FALSE,
  sexo           = FALSE,
  valor_charlson = FALSE,
  valor_cf       = FALSE,     # Barthel
  valor_rLPP     = FALSE,     # Braden
  valor_rC       = FALSE,     # Downton
  perfil         = TRUE,
  frec_cardiaca  = TRUE,
  ta_sistolica   = TRUE,
  ta_diastolica  = TRUE,
  sat_o2         = TRUE,
  tipo_ingreso   = TRUE,
  tipoproceso    = TRUE
)
variables <- names(usar)[usar]


# Datos

datos <- read.csv(ruta, sep = ",", dec = ",", stringsAsFactors = FALSE,
                  na.strings = c("", "NA"))

# primera valoración de cada episodio, con las 4 constantes a 24 h o menos
datos <- datos[datos$orden_valenf == 1 & datos$ctes_4vars_en_24h == 1, ]

valoracion <- as.POSIXct(datos$fecha_valoracion)
alta       <- as.POSIXct(datos$fecha_alta)
exitus     <- as.POSIXct(datos$fecha_exitus)
murio      <- as.logical(datos$hay_exitus_episodio)

# del exitus solo se sabe el día: se toma el punto medio del intervalo posible
dia <- as.numeric(as.POSIXct(format(exitus, "%Y-%m-%d")))
ini <- pmax(as.numeric(valoracion), dia)
fin <- ifelse(murio, (ini + dia + 86400) / 2, as.numeric(alta))

datos$tiempo <- (fin - as.numeric(valoracion)) / 86400      # en días
datos$evento <- as.integer(murio)

# fuera las altas anteriores a la valoración o a la misma hora
excluidos <- datos[datos$evento == 0 & datos$tiempo <= 0, ]
datos     <- datos[datos$evento == 1 | datos$tiempo > 0, ]

# la primera categoría de cada variable es la de referencia
categorica <- function(x, niveles) factor(x, levels = niveles, labels = paste0("=", niveles))
datos$p1 <- categorica(datos$respuesta_pregunta1, c(15, 10, 5, 0))
datos$p2 <- categorica(datos$respuesta_pregunta2, 4:1)
datos$p3 <- categorica(datos$respuesta_pregunta3, 4:1)
datos$p4 <- categorica(datos$respuesta_pregunta4, 4:1)
datos$p5 <- categorica(datos$respuesta_pregunta5, 2:1)
datos$p6 <- categorica(datos$respuesta_pregunta6, 0:1)
datos$p7 <- categorica(datos$respuesta_pregunta7, 0:1)
datos$sexo         <- categorica(datos$sexo, c("H", "M"))
datos$perfil       <- categorica(datos$primer_perfil, c(5, 1, 2, 3, 4))
datos$tipo_ingreso <- categorica(datos$tipo_ingreso, c("P", "U"))
datos$tipoproceso  <- categorica(datos$tipoproceso, c("M", "Q"))

# mismos pacientes sea cual sea la combinación de variables
datos <- datos[complete.cases(datos[, names(usar)]), ]

cat("N =", nrow(datos), " exitus =", sum(datos$evento),
    " censurados =", sum(datos$evento == 0), " excluidos =", nrow(excluidos), "\n")
cat("Variables en el modelo:", variables, "\n")


# Funciones

ajustar <- function(v, d = datos) {
  derecha <- if (length(v) == 0) "1" else paste(v, collapse = " + ")
  survreg(as.formula(paste("Surv(tiempo, evento) ~", derecha)), data = d, dist = "weibull",
          control = survreg.control(maxiter = 100))
}

# razón de verosimilitudes entre dos modelos anidados
prueba_lr <- function(grande, pequeno) {
  x  <- 2 * (as.numeric(logLik(grande)) - as.numeric(logLik(pequeno)))
  gl <- sum(!is.na(coef(grande))) - sum(!is.na(coef(pequeno)))
  c(chi2 = x, gl = gl, p = pchisq(x, gl, lower.tail = FALSE),
    menos_log10_p = -pchisq(x, gl, lower.tail = FALSE, log.p = TRUE) / log(10))
}

indiceC <- function(m, d = datos) {
  lp <- predict(m, newdata = d, type = "lp")
  c  <- concordance(Surv(d$tiempo, d$evento) ~ lp)$concordance
  if (c < 0.5) 1 - c else c
}

# factor sobre el tiempo con su IC del 95 %, y cociente de riesgos (HR = factor^-alfa)
tabla_coeficientes <- function(m) {
  b  <- coef(m)
  ee <- sqrt(diag(vcov(m)))[1:length(b)]
  data.frame(
    beta   = round(b, 4),
    EE     = round(ee, 4),
    z      = round(b / ee, 2),
    p      = signif(2 * pnorm(-abs(b / ee)), 3),
    factor_tiempo = round(exp(b), 3),
    IC_inf = round(exp(b - 1.96 * ee), 3),
    IC_sup = round(exp(b + 1.96 * ee), 3),
    HR     = round(exp(-b / m$scale), 3)
  )
}


# 1. Modelo

modelo <- ajustar(variables)

cat("\n MODELO \n")
print(tabla_coeficientes(modelo))


# 2. Aportación de cada variable (se quita del modelo y se compara)

aporte <- do.call(rbind, lapply(variables, function(v) {
  q <- prueba_lr(modelo, ajustar(setdiff(variables, v)))
  data.frame(variable = v, gl = q["gl"], chi2 = round(q["chi2"], 1), p = signif(q["p"], 3),
             mlp = q["menos_log10_p"])
}))
aporte <- aporte[order(aporte$p), ]

cat("\n\n APORTACIÓN POR VARIABLE \n")
print(aporte[, 1:4], row.names = FALSE)


# 3. Resumen e índice C con validación cruzada (5 particiones)

set.seed(1)
particion   <- sample(rep(1:5, length.out = nrow(datos)))
c_particion <- sapply(1:5, function(k) {
  indiceC(ajustar(variables, datos[particion != k, ]), datos[particion == k, ])
})
c_nuevos    <- mean(c_particion)
c_modelo    <- indiceC(modelo)

cat("\n\n RESUMEN \n")
cat("parametros =", sum(!is.na(coef(modelo))) + 1,
    "  logLik =", round(as.numeric(logLik(modelo)), 2),
    "  AIC =", round(AIC(modelo), 1), "\n")
cat("C =", round(c_modelo, 3), "  C con datos nuevos =", round(c_nuevos, 3),
    "  (por particion:", round(c_particion, 3), ")\n")
cat("intercepto =", round(coef(modelo)[1], 4), "  sigma =", round(modelo$scale, 4),
    "  forma =", round(1 / modelo$scale, 4), "\n")


# 4. Ejemplo: riesgo de morir en los primeros 7 días según el perfil (con las variables del modelo final)

ejemplo <- data.frame(perfil = levels(datos$perfil), frec_cardiaca = 80, ta_sistolica = 130,
                      ta_diastolica = 75, sat_o2 = 95, tipo_ingreso = "=U", tipoproceso = "=M")
if (all(variables %in% names(ejemplo))) {
  mu <- predict(modelo, newdata = ejemplo, type = "lp")
  ejemplo$riesgo_7d <- round(100 * (1 - exp(-(7 / exp(mu))^(1 / modelo$scale))), 1)
  cat("\n\n EJEMPLO (FC 80, TAS 130, TAD 75, SatO2 95, urgente, medico): riesgo a 7 dias (%) \n")
  print(ejemplo[order(ejemplo$perfil), c("perfil", "riesgo_7d")], row.names = FALSE)
}


# 5. Correlación entre las variables del modelo

codigos <- model.matrix(modelo)[, -1, drop = FALSE]
correlacion <- cor(codigos, method = "spearman")

cat("\n\n CORRELACION ENTRE VARIABLES \n")
print(round(correlacion, 2))


# 6. Muertes observadas frente a esperadas, con las predicciones de la validación cruzada

mu_nuevo <- sigma_nuevo <- numeric(nrow(datos))
for (k in 1:5) {
  m <- ajustar(variables, datos[particion != k, ])
  mu_nuevo[particion == k]    <- predict(m, newdata = datos[particion == k, ], type = "lp")
  sigma_nuevo[particion == k] <- m$scale
}
esperadas <- function(dias) (pmin(datos$tiempo, dias) / exp(mu_nuevo))^(1 / sigma_nuevo)
reales    <- function(dias) datos$evento == 1 & datos$tiempo <= dias

grupos <- list("Todos" = rep(TRUE, nrow(datos)))
for (p in 1:5) grupos[[paste("Perfil", p)]] <- datos$primer_perfil == p

tabla_realidad <- function(d) do.call(rbind, lapply(names(grupos), function(g) {
  i  <- grupos[[g]]
  O  <- sum(reales(d)[i])
  E  <- sum(esperadas(d)[i])
  ic <- poisson.test(O, E)$conf.int
  data.frame(grupo = g, n = sum(i), reales = O, esperadas = round(E, 1),
             O_E = round(O / E, 2), IC_inf = round(ic[1], 2), IC_sup = round(ic[2], 2))
}))
realidad <- list("7" = tabla_realidad(7), "14" = tabla_realidad(14))

cat("\n\n REALIDAD FRENTE A MODELO \n")
for (d in names(realidad)) {
  cat("\n Muertes en los primeros", d, "dias\n")
  print(realidad[[d]], row.names = FALSE)
}


# Gráficos

azul <- "#2a78d6"; gris <- "grey75"; tinta <- "grey20"
carpeta <- "graficos_perfiles"
dir.create(carpeta, showWarnings = FALSE)

# p valor de cada variable (en azul las que aportan)
texto_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

dibujar_p <- function(mlp, p, titulo) {
  tope <- max(mlp, 2) * 1.35
  bp <- barplot(rev(mlp), horiz = TRUE, names.arg = rev(aporte$variable), las = 1, cex.names = 0.85,
                xlim = c(0, tope), border = NA, col = rev(ifelse(p < 0.05, azul, gris)),
                xaxt = "n", xlab = "p valor", main = titulo)
  marcas <- if (tope <= 6) c(0.05, 0.01, 0.001, 1e-5) else if (tope <= 25) c(0.05, 0.001, 1e-5, 1e-10, 1e-20)
            else c(0.05, 1e-10, 1e-20, 1e-40, 1e-80, 1e-160)
  marcas <- marcas[-log10(marcas) <= tope]
  axis(1, -log10(marcas), ifelse(marcas >= 0.001, as.character(marcas), paste0("1e-", round(-log10(marcas)))))
  abline(v = -log10(0.05), lty = 2, col = "red")
  text(pmax(rev(mlp), -log10(0.05)), bp, rev(texto_p(p)), pos = 4, cex = 0.7, col = tinta)
}

png(file.path(carpeta, "1_p_valor.png"), width = 900, height = 200 + 38 * nrow(aporte), res = 125)
par(mar = c(4.5, 9, 3, 1))
dibujar_p(aporte$mlp, aporte$p, "Lo que aporta cada variable con las demas dentro")
invisible(dev.off())

# índice C
png(file.path(carpeta, "2_indice_c.png"), width = 900, height = 420, res = 120)
par(mar = c(4.5, 10, 4, 2))
valores <- c("con datos nuevos" = c_nuevos, "C" = c_modelo)
bp <- barplot(valores, horiz = TRUE, las = 1, xlim = c(0.5, 1), xpd = FALSE, border = NA,
              col = c("grey60", azul), xlab = "indice C", main = "Capacidad de ordenar pacientes")
abline(v = c(0.7, 0.8), lty = 3, col = "grey50")
text(valores, bp, sprintf("%.3f", valores), pos = 4, cex = 0.9, col = tinta)
invisible(dev.off())

# correlación
nv <- ncol(correlacion)
png(file.path(carpeta, "3_correlacion.png"), width = 1100, height = 1000, res = 120)
par(mar = c(8, 8, 3, 2))
image(1:nv, 1:nv, correlacion[, nv:1, drop = FALSE], zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
      col = colorRampPalette(c("firebrick", "white", "steelblue"))(41),
      main = "Correlacion entre las variables del modelo")
axis(1, 1:nv, colnames(correlacion), las = 2, cex.axis = 0.8)
axis(2, 1:nv, rev(colnames(correlacion)), las = 1, cex.axis = 0.8)
text(rep(1:nv, nv), rep(nv:1, each = nv), sprintf("%.2f", correlacion), cex = 0.6)
invisible(dev.off())

# muertes observadas y esperadas acumuladas en los primeros 14 días (figura de la memoria)
dias <- seq(0, 14, by = 0.25)
pdf(file.path(carpeta, "4_realidad_vs_modelo.pdf"), width = 7.2, height = 5)
par(mfrow = c(2, 3), mar = c(4, 4.2, 2.2, 0.8), mgp = c(2.5, 0.6, 0), las = 1, cex.axis = 0.85)
for (g in names(grupos)) {
  i   <- grupos[[g]]
  obs <- sapply(dias, function(d) sum(reales(d)[i]))
  esp <- sapply(dias, function(d) sum(esperadas(d)[i]))
  plot(dias, obs, type = "n", ylim = c(0, max(obs, esp) * 1.08), xaxt = "n", bty = "l",
       xlab = "Días desde la valoración", ylab = "Muertes acumuladas",
       main = sprintf("%s (n = %s)", g, formatC(sum(i), big.mark = ".", decimal.mark = ",")),
       cex.main = 0.95, font.main = 1)
  axis(1, at = c(0, 7, 14)); abline(v = c(7, 14), col = "grey90")
  lines(dias, esp, lwd = 2.2, col = "#2b6cb0")
  lines(dias, obs, type = "s", lwd = 2, col = "grey15")
  if (g == "Todos") legend("topleft", c("Observadas", "Esperadas por el modelo"), col = c("grey15", "#2b6cb0"),
                           lwd = 2, bty = "n", cex = 0.85)
}
invisible(dev.off())

# O/E por perfil con su IC del 95 %, a los 7 y a los 14 días
naranja  <- "#eb6834"
colores  <- c("7" = azul, "14" = naranja)
desplaza <- c("7" = 0.15, "14" = -0.15)
filas    <- length(grupos):1
tope     <- max(sapply(realidad, function(r) r$IC_sup)) * 1.12
png(file.path(carpeta, "5_O_E.png"), width = 1000, height = 650, res = 120)
par(mar = c(4.5, 7, 4.5, 2))
plot(NA, xlim = c(0, tope), ylim = c(0.5, length(grupos) + 0.5), yaxt = "n", ylab = "",
     xlab = "O/E = muertes reales / muertes esperadas por el modelo",
     main = "Muertes reales frente a esperadas, por perfil")
axis(2, filas, names(grupos), las = 1)
abline(h = length(grupos) - 0.5, col = "grey85")
abline(v = 1, lty = 2, col = "red")
mtext(c("el modelo esperaba mas  ", "  el modelo esperaba menos"), side = 3, line = 0.3,
      at = 1, adj = c(1, 0), cex = 0.8, col = "grey40")
for (d in names(realidad)) {
  r <- realidad[[d]]
  y <- filas + desplaza[d]
  segments(r$IC_inf, y, r$IC_sup, y, lwd = 2, col = colores[d])
  points(r$O_E, y, pch = 19, col = colores[d])
  text(r$IC_sup, y, sprintf("%.2f", r$O_E), pos = 4, cex = 0.75, col = colores[d])
}
legend("topright", paste(names(realidad), "dias"), col = colores, pch = 19, lwd = 2,
       bty = "n", cex = 0.9)
invisible(dev.off())


cat("\nGraficos en la carpeta '", carpeta, "'\n", sep = "")
