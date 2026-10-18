# Comprobaciones que se citan en la discusión

library(survival)

ruta <- "datos/capacidad_funcional_1VALENF_LARGO_con_constantes.csv"


# Datos y modelo final (como en modelo_perfiles.R)

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

modelo <- ajustar()
c_weib <- sapply(1:5, function(k) {
  m <- ajustar(d = datos[particion != k, ])
  indiceC(predict(m, newdata = datos[particion == k, ], type = "lp"), datos[particion == k, ])
})

grupos <- c(list(Todos = rep(TRUE, nrow(datos))),
            setNames(lapply(1:5, function(p) datos$primer_perfil == p), paste("Perfil", 1:5)))
reales <- function(d) datos$evento == 1 & datos$tiempo <= d


# 1. Muertes en los 30 días siguientes al alta

vivos     <- datos$evento == 0
tras_alta <- as.numeric(as.Date(datos$fecha_exitus) - as.Date(datos$fecha_alta))
muere_30d <- vivos & !is.na(tras_alta) & tras_alta >= 0 & tras_alta <= 30
cat("\nAltas con vida:", sum(vivos), " muertes en los 30 dias siguientes al alta:", sum(muere_30d),
    sprintf("(%.1f %%)", 100 * sum(muere_30d) / sum(vivos)), "\n")
print(do.call(rbind, lapply(1:5, function(p) {
  i <- datos$primer_perfil == p
  data.frame(perfil = p, muertes_hospital = sum(datos$evento[i]), altas_con_vida = sum(vivos & i),
             muertes_30d_tras_alta = sum(muere_30d & i),
             porcentaje = round(100 * sum(muere_30d & i) / sum(vivos & i), 1))
})), row.names = FALSE)


# 2. Episodios excluidos por las constantes

todo <- read.csv(ruta, sep = ",", dec = ",", stringsAsFactors = FALSE, na.strings = c("", "NA"))
todo <- todo[todo$orden_valenf == 1, ]
excluidos <- todo[is.na(todo$ctes_4vars_en_24h) | todo$ctes_4vars_en_24h != 1, ]
cat("\nExcluidos por las constantes:", nrow(excluidos), "de", nrow(todo), "episodios\n")
print(table(perfil = excluidos$primer_perfil, exitus = as.logical(excluidos$hay_exitus_episodio)))
cat("Episodios del perfil 5 antes de excluir:", sum(todo$primer_perfil == 5), "\n")
rm(todo)


# 3. Errores estándar robustos agrupados por paciente

robusto <- survreg(as.formula(paste("Surv(tiempo, evento) ~", derecha)), data = datos, dist = "weibull",
                   cluster = id_anonim_paciente, control = survreg.control(maxiter = 100))
coef_r <- coef(robusto)
ee_r   <- sqrt(diag(vcov(robusto)))[seq_along(coef_r)]
ee_n   <- sqrt(diag(vcov(modelo)))[seq_along(coef_r)]
cat("\nErrores robustos agrupados por paciente (aumento = error robusto / error del modelo)\n")
print(signif(data.frame(factor = exp(coef_r), inf = exp(coef_r - 1.96 * ee_r), sup = exp(coef_r + 1.96 * ee_r),
                        p = 2 * pnorm(-abs(coef_r / ee_r)), aumento = ee_r / ee_n), 4))


# 4. Validación con particiones por paciente y con la lognormal

validar <- function(part, dist = "weibull") {
  lp <- escala <- numeric(nrow(datos)); C <- numeric(5)
  for (k in 1:5) {
    m <- ajustar(d = datos[part != k, ], dist = dist)
    lp[part == k]     <- predict(m, newdata = datos[part == k, ], type = "lp")
    escala[part == k] <- m$scale
    C[k] <- indiceC(lp[part == k], datos[part == k, ])
  }
  H <- function(d) {   # riesgo acumulado hasta el alta, la muerte o el día d
    z <- (log(pmin(datos$tiempo, d)) - lp) / escala
    if (dist == "weibull") exp(z) else -pnorm(z, lower.tail = FALSE, log.p = TRUE)
  }
  cat("Indice C:", sprintf("%.3f", mean(C)), "\n")
  print(do.call(rbind, lapply(c(7, 14), function(d) do.call(rbind, lapply(names(grupos), function(g) {
    i <- grupos[[g]]; O <- sum(reales(d)[i]); E <- sum(H(d)[i])
    data.frame(dias = d, grupo = g, O = O, E = round(E, 1), OE = round(O / E, 2))
  })))), row.names = FALSE)
}
pacientes <- unique(datos$id_anonim_paciente)
set.seed(1)
particion_paciente <- sample(rep(1:5, length.out = length(pacientes)))[match(datos$id_anonim_paciente, pacientes)]
cat("\nWeibull con las particiones hechas por paciente\n");  validar(particion_paciente)
cat("\nLognormal con las particiones de siempre\n");         validar(particion, "lognormal")


# 5. Perfiles 1 y 2: diferencia entre ellos y modelo con los dos juntos

m_ref2 <- ajustar(d = transform(datos, perfil = relevel(perfil, ref = "=2")))
b1 <- coef(m_ref2)["perfil=1"]; e1 <- sqrt(vcov(m_ref2)["perfil=1", "perfil=1"])
cat("\nPerfil 1 frente al 2: factor =", round(exp(b1), 3), " HR =", round(exp(-b1 / m_ref2$scale), 2),
    " p =", signif(2 * pnorm(-abs(b1 / e1)), 2), "\n")
datos$perfil12 <- factor(ifelse(datos$primer_perfil %in% 1:2, "1-2", datos$primer_perfil),
                         levels = c("5", "1-2", "3", "4"))
v12 <- c("perfil12", setdiff(variables, "perfil"))
m12 <- ajustar(v12)
x12 <- 2 * (as.numeric(logLik(modelo)) - as.numeric(logLik(m12)))
c12 <- sapply(1:5, function(k) {
  m <- ajustar(v12, d = datos[particion != k, ])
  indiceC(predict(m, newdata = datos[particion == k, ], type = "lp"), datos[particion == k, ])
})
cat("Perfiles 1 y 2 juntos: chi2 =", round(x12, 2), " gl = 1  p =", signif(pchisq(x12, 1, lower.tail = FALSE), 2),
    "| indice C con validacion cruzada:", sprintf("%.3f", mean(c12)), "(con los cinco perfiles:", sprintf("%.3f", mean(c_weib)), ")\n")
