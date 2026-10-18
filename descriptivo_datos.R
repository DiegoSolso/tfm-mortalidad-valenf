# Capítulo de datos: filtros, tabla por perfil, mortalidad, Kaplan-Meier y correlaciones

library(survival)

ruta    <- "datos/capacidad_funcional_1VALENF_LARGO_con_constantes.csv"
carpeta <- "graficos_datos"
dir.create(carpeta, showWarnings = FALSE)

dia_de <- function(x) substr(x, 1, 10)     # solo la fecha
ctes   <- c("frec_cardiaca", "ta_sistolica", "ta_diastolica", "sat_o2")


# 1. Fichero

todo <- read.csv(ruta, sep = ",", dec = ",", stringsAsFactors = FALSE, na.strings = c("", "NA"))
valoraciones_ep <- table(todo$id_anonim_episodio)

cat("\n FICHERO \n")
cat("Valoraciones:", nrow(todo), " episodios:", length(valoraciones_ep),
    " pacientes:", length(unique(todo$id_anonim_paciente)), "\n")
cat("Valoraciones por episodio: de", min(valoraciones_ep), "a", max(valoraciones_ep),
    "| episodios con una sola:", sprintf("%.1f%%", 100 * mean(valoraciones_ep == 1)), "\n")
cat("Edad: de", min(todo$edad), "a", max(todo$edad), "\n")
cat("Ingresos: del", min(dia_de(todo$fecha_ingreso)), "al", max(dia_de(todo$fecha_ingreso)), "\n")


# 2. Filtros

# primera valoración de cada episodio
primeras <- todo[which(todo$orden_valenf == 1), ]
horas    <- as.numeric(difftime(as.POSIXct(primeras$fecha_valoracion), as.POSIXct(primeras$fecha_ingreso),
                                units = "hours"))
con_ingreso <- length(unique(todo$id_anonim_episodio[which(todo$motivo_FINAL == "Ingreso")]))

cat("\n FILTROS \n")
cat("Primera valoracion de cada episodio:", nrow(primeras), "\n")
cat("  horas desde el ingreso: mediana", round(median(horas), 1),
    "| en las primeras 24 h:", sprintf("%.1f%%", 100 * mean(horas <= 24)), "\n")
cat("  episodios con alguna valoracion etiquetada como Ingreso:", con_ingreso, "\n")

murio1 <- as.logical(primeras$hay_exitus_episodio)
cat("  exitus el mismo dia que la primera valoracion:",
    sum(murio1 & dia_de(primeras$fecha_exitus) == dia_de(primeras$fecha_valoracion), na.rm = TRUE), "\n")

ultima <- tapply(todo$fecha_valoracion, todo$id_anonim_episodio, max)
ultima <- dia_de(ultima[as.character(primeras$id_anonim_episodio)])
cat("  exitus con fecha anterior a la ultima valoracion:",
    sum(murio1 & dia_de(primeras$fecha_exitus) < ultima, na.rm = TRUE), "\n")

# las 4 constantes a 24 h o menos de la valoración
sin_ctes <- primeras[which(primeras$ctes_4vars_en_24h != 1 | is.na(primeras$ctes_4vars_en_24h)), ]
datos    <- primeras[which(primeras$ctes_4vars_en_24h == 1), ]
falta    <- !complete.cases(sin_ctes[, ctes])
cat("Con las 4 constantes a 24 h o menos:", nrow(datos), " excluidos:", nrow(sin_ctes),
    "(falta alguna:", sum(falta), "| alguna a mas de 24 h:", sum(!falta), ")\n")
cat("  entre los excluidos: programados", sprintf("%.1f%%", 100 * mean(sin_ctes$tipo_ingreso == "P")),
    "| quirurgicos", sprintf("%.1f%%", 100 * mean(sin_ctes$tipoproceso == "Q")),
    "| exitus", sum(as.logical(sin_ctes$hay_exitus_episodio)), "\n")
distancia <- abs(as.matrix(datos[, paste0(ctes, "_delta_h")]))
cat("  distancia a la valoracion: mediana", round(60 * median(distancia), 0), "minutos",
    "| las 4 a menos de 6 h:", sprintf("%.1f%%", 100 * mean(apply(distancia, 1, max) < 6)), "\n")

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

# altas anteriores a la valoración o a la misma hora
excluidos <- datos[datos$evento == 0 & datos$tiempo <= 0, ]
datos     <- datos[datos$evento == 1 | datos$tiempo > 0, ]
cat("Altas anteriores o iguales a la valoracion:", nrow(excluidos),
    "(anteriores:", sum(excluidos$tiempo < 0), "| a la misma hora:", sum(excluidos$tiempo == 0), ")",
    " exitus con tiempo <= 0:", sum(datos$evento == 1 & datos$tiempo <= 0), "\n")

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

# datos completos en todas las variables candidatas
candidatas <- c(paste0("p", 1:7), "edad", "sexo", "valor_charlson", "valor_cf", "valor_rLPP",
                "valor_rC", "perfil", "frec_cardiaca", "ta_sistolica", "ta_diastolica", "sat_o2",
                "tipo_ingreso", "tipoproceso")
n_antes <- nrow(datos)
datos   <- datos[complete.cases(datos[, candidatas]), ]
episodios_paciente <- table(datos$id_anonim_paciente)

cat("Sin datos completos:", n_antes - nrow(datos), "\n")
cat("MUESTRA FINAL: episodios =", nrow(datos), " pacientes =", length(episodios_paciente),
    " exitus =", sum(datos$evento), " altas =", sum(datos$evento == 0), "\n")
cat("  pacientes con mas de un episodio:", sum(episodios_paciente > 1),
    "| exitus el mismo dia que la valoracion:",
    sum(datos$evento == 1 & dia_de(datos$fecha_exitus) == dia_de(datos$fecha_valoracion)), "\n")
cat("  altas con vida con fecha de exitus (murieron despues del alta):",
    sum(datos$evento == 0 & !is.na(datos$fecha_exitus)), "\n")


# 3. Tabla descriptiva por perfil

media <- function(x) sprintf("%.1f (%.1f)", mean(x), sd(x))
pct   <- function(x) sprintf("%.1f", 100 * mean(x))
med   <- function(x) sprintf("%.1f (%.1f-%.1f)", median(x), quantile(x, 0.25), quantile(x, 0.75))

grupos <- c(list(Todos = rep(TRUE, nrow(datos))),
            setNames(lapply(levels(datos$perfil), function(p) datos$perfil == p), paste("Perfil", 1:5)))
descriptivo <- sapply(grupos, function(i) {
  d <- datos[i, ]
  c("Episodios (%)"          = sprintf("%d (%.1f)", nrow(d), 100 * nrow(d) / nrow(datos)),
    "Edad"                   = media(d$edad),
    "Mujeres (%)"            = pct(d$sexo == "=M"),
    "Charlson"               = media(d$valor_charlson),
    "Barthel"                = media(d$valor_cf),
    "Braden"                 = media(d$valor_rLPP),
    "Downton"                = media(d$valor_rC),
    "FC"                     = media(d$frec_cardiaca),
    "TAS"                    = media(d$ta_sistolica),
    "TAD"                    = media(d$ta_diastolica),
    "SatO2"                  = media(d$sat_o2),
    "Ingreso urgente (%)"    = pct(d$tipo_ingreso == "=U"),
    "Proceso quirurgico (%)" = pct(d$tipoproceso == "=Q"),
    "Seguimiento (dias)"     = med(d$tiempo),
    "Exitus (%)"             = sprintf("%d (%.1f)", sum(d$evento), 100 * mean(d$evento)))
})
cat("\n TABLA DESCRIPTIVA \n")
print(noquote(descriptivo))


# 4. Mortalidad

m <- datos[datos$evento == 1, ]
cat("\n MORTALIDAD \n")
cat("Dias hasta el exitus: mediana", round(median(m$tiempo), 1),
    "| exitus en 7 dias:", sum(m$tiempo <= 7), "| en 14 dias:", sum(m$tiempo <= 14), "\n")
cat("Siguen ingresados a los 14 dias:", sum(datos$tiempo > 14),
    sprintf("(%.1f%%)", 100 * mean(datos$tiempo > 14)), "\n")
cat("Mortalidad (%) por tipo de ingreso:\n"); print(round(100 * tapply(datos$evento, datos$tipo_ingreso, mean), 1))
cat("Mortalidad (%) por tipo de proceso:\n"); print(round(100 * tapply(datos$evento, datos$tipoproceso, mean), 1))
cat("Ingresos programados que son quirurgicos:",
    sprintf("%.1f%%", 100 * mean(datos$tipoproceso[datos$tipo_ingreso == "=P"] == "=Q")), "\n")


# 5. Kaplan-Meier por perfil

km <- survfit(Surv(tiempo, evento) ~ perfil, data = datos)
cat("\n SUPERVIVENCIA ESTIMADA A LOS 7, 14 Y 30 DIAS \n")
s <- summary(km, times = c(7, 14, 30), extend = TRUE)
print(data.frame(perfil = sub("perfil==", "", s$strata), dias = s$time, en_riesgo = s$n.risk,
                 S = round(s$surv, 3)))
lr <- survdiff(Surv(tiempo, evento) ~ perfil, data = datos)
cat("\nLog-rank: chi2 =", round(lr$chisq, 1), " gl =", length(lr$n) - 1,
    " p =", signif(pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE), 3), "\n")

colores <- c("#2b83ba", "#5ab4ac", "#e0a526", "#e0592b", "#9e1b32")   # perfil 1 a 5
coma    <- function(x, d = 1) formatC(x, format = "f", digits = d, decimal.mark = ",")
miles   <- function(x) formatC(x, format = "d", big.mark = ".", decimal.mark = ",")
tiempos <- seq(0, 30, by = 5)
riesgo  <- summary(km, times = tiempos, extend = TRUE)
riesgo  <- split(riesgo$n.risk, riesgo$strata)

pdf(file.path(carpeta, "km_perfil.pdf"), width = 7, height = 5.2)
par(mar = c(10, 6.5, 1, 1.5), mgp = c(3, 0.7, 0), las = 1, cex.axis = 0.9)
plot(NA, xlim = c(0, 30), ylim = c(0.4, 1), xaxs = "i", xaxt = "n", yaxt = "n", bty = "l",
     xlab = "", ylab = "Supervivencia estimada")
abline(h = seq(0.4, 1, 0.1), col = "grey90")
abline(v = tiempos, col = "grey90")
axis(1, at = tiempos)
axis(2, at = seq(0.4, 1, 0.1), labels = coma(seq(0.4, 1, 0.1)))
mtext("Días desde la primera valoración", side = 1, line = 2.3)
lines(km, col = colores, lwd = 2.2, mark.time = FALSE)
legend("bottomleft", sprintf("Perfil %d (n = %s)", 1:5, miles(as.numeric(table(datos$perfil)))),
       col = colores, lwd = 2.2, bty = "n", cex = 0.85, inset = c(0.02, 0.02))
# pacientes en riesgo debajo del gráfico
mtext("Pacientes que siguen ingresados", side = 1, line = 4, at = 15, cex = 0.8, font = 2)
for (k in 1:5) {
  mtext(paste("Perfil", k), side = 1, line = 4 + k, at = -1.8, adj = 1, cex = 0.78, col = colores[k])
  mtext(miles(riesgo[[k]]), side = 1, line = 4 + k, at = tiempos, cex = 0.78)
}
invisible(dev.off())


# 6. Correlación de Spearman (urgente y quirúrgico = 1)

finales <- data.frame(Perfil = as.numeric(datos$perfil), FC = datos$frec_cardiaca,
                      TAS = datos$ta_sistolica, TAD = datos$ta_diastolica, SatO2 = datos$sat_o2,
                      Urgente = as.numeric(datos$tipo_ingreso == "=U"),
                      Quirurgico = as.numeric(datos$tipoproceso == "=Q"))
cat("\n CORRELACION DE SPEARMAN (predictores del modelo final) \n")
print(round(cor(finales, method = "spearman"), 2))

# perfil frente a las variables con las que se definieron los perfiles
componentes <- data.frame(Perfil = as.numeric(datos$perfil), Edad = datos$edad,
                          Mujer = as.numeric(datos$sexo == "=M"), Charlson = datos$valor_charlson,
                          Barthel = datos$valor_cf, Braden = datos$valor_rLPP, Downton = datos$valor_rC)
cat("\n CORRELACION DE SPEARMAN (perfil y variables de los perfiles) \n")
print(round(cor(componentes, method = "spearman"), 2))

cat("\nGrafico guardado en", file.path(carpeta, "km_perfil.pdf"), "\n")
