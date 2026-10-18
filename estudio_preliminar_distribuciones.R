# Estudio preliminar: comparación de distribuciones con el modelo reducido de 10 variables
# (primera valoración de cada episodio y tiempo desde el ingreso)

library(survival)

ruta <- "datos/capacidad_funcional_1VALENF_LARGO_con_constantes.csv"

datos <- read.csv(ruta, sep = ",", dec = ",", stringsAsFactors = FALSE)
datos <- datos[datos$orden_valenf == 1, ]
datos <- data.frame(evento = as.integer(as.logical(datos$hay_exitus_episodio)),
                    tiempo = as.numeric(datos$dias_estancia), edad_ingreso = as.numeric(datos$edad_ingreso),
                    sexo = as.factor(datos$sexo), valor_charlson = as.numeric(datos$valor_charlson),
                    grupo = as.factor(datos$primer_perfil),
                    datos[, paste0("respuesta_pregunta", 1:7)])
datos <- na.omit(datos[datos$tiempo > 0, ])

f10 <- as.formula(paste("Surv(tiempo, evento) ~",
                        paste(c("edad_ingreso", "sexo", "valor_charlson", paste0("respuesta_pregunta", 1:7)),
                              collapse = " + ")))

cat("N =", nrow(datos), " exitus =", sum(datos$evento), "\n")
for (dist in c("lognormal", "loglogistic", "weibull", "exponential")) {
  m <- survreg(f10, data = datos, dist = dist)
  cat(sprintf("  %-12s AIC = %.1f   C (datos del ajuste) = %.4f\n", dist, AIC(m),
              concordance(Surv(tiempo, evento) ~ predict(m, type = "lp"), data = datos)$concordance))
}
cat(sprintf("  %-12s C (datos del ajuste) = %.4f\n", "Cox", summary(coxph(f10, data = datos))$concordance[1]))
