# Estudio preliminar: Cox frente a RSF para la mortalidad intrahospitalaria,
# con búsqueda en rejilla de los hiperparámetros del RSF

# temporales en el home (más espacio que /tmp)
tmp_dir <- path.expand("~/Rtmp")
if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_dir)

# randomForestSRC en un solo hilo para no saturar la memoria
options(rf.cores = 1, mc.cores = 1)

library(DBI)
library(RPostgres)
library(dplyr)
library(tidyr)
library(lubridate)
library(survival)
library(randomForestSRC)
library(ggplot2)
library(broom)

set.seed(456)

output_dir <- "graficos_v6"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

df_final_cache <- file.path(output_dir, "df_final.rds")

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# 1-3. Datos (se guardan en df_final.rds para no repetir la consulta)
if (file.exists(df_final_cache)) {
  log_msg("Cache detectado: cargando df_final desde disco")
  df_final <- readRDS(df_final_cache)
} else {
  log_msg("[1] Cargando valenf.csv")
  df_enfermeria <- read.csv(
    "datos/valenf.csv",
    sep = ","
  )
  colnames(df_enfermeria) <- tolower(colnames(df_enfermeria))
  df_enfermeria <- df_enfermeria %>%
    filter(motivo_final == "['Ingreso']") %>%
    distinct(id_anonim_episodio, .keep_all = TRUE)

  log_msg("[2] Conectando a la BD")
  # credenciales en variables de entorno (ver README.md)
  con <- dbConnect(RPostgres::Postgres(),
                   dbname   = Sys.getenv("VALENF_BD"),
                   host     = Sys.getenv("VALENF_HOST"),
                   port     = 5432,
                   user     = Sys.getenv("VALENF_USUARIO"),
                   password = Sys.getenv("VALENF_CLAVE"))
  df_hosp <- dbGetQuery(con, "
    SELECT ID_ANONIM_EPISODIO, EDAD_INGRESO, SEXO,
           TIPO_INGRESO, TIPOPROCESO,
           FECHA_INGRESO, FECHA_ALTA,
           FECHA_EXITUS, MOTIVO_ALTA_PK
    FROM FTR_HOSP
  ")
  dbDisconnect(con)
  colnames(df_hosp) <- tolower(colnames(df_hosp))

  log_msg("[3] Construyendo dataset final")
  df_final <- df_hosp %>%
    inner_join(df_enfermeria, by = "id_anonim_episodio") %>%
    mutate(
      fecha_ingreso_dt = ymd_hms(fecha_ingreso),
      fecha_alta_dt    = ymd_hms(fecha_alta),
      evento_mortalidad = ifelse(motivo_alta_pk == 4, 1, 0),
      tiempo_mortalidad = as.numeric(
        difftime(fecha_alta_dt, fecha_ingreso_dt, units = "days")
      ),
      sexo         = as.factor(sexo),
      tipo_ingreso = as.factor(tipo_ingreso),
      tipoproceso  = as.factor(tipoproceso),
      cat_braden   = as.factor(categoria_braden),
      cat_downton  = as.factor(categoria_downton),
      valor_barthel = as.numeric(gsub(",", ".", as.character(valor_barthel))),
      respuesta_pregunta1 = as.numeric(respuesta_pregunta1),
      respuesta_pregunta2 = as.numeric(respuesta_pregunta2),
      respuesta_pregunta3 = as.numeric(respuesta_pregunta3),
      respuesta_pregunta4 = as.numeric(respuesta_pregunta4),
      respuesta_pregunta5 = as.numeric(respuesta_pregunta5),
      respuesta_pregunta6 = as.numeric(respuesta_pregunta6),
      respuesta_pregunta7 = as.numeric(respuesta_pregunta7)
    ) %>%
    filter(tiempo_mortalidad > 0) %>%
    select(tiempo_mortalidad, evento_mortalidad,
           edad_ingreso, sexo, tipo_ingreso, tipoproceso,
           cat_braden, cat_downton, valor_barthel,
           respuesta_pregunta1:respuesta_pregunta7) %>%
    na.omit()

  saveRDS(df_final, df_final_cache)
  rm(df_hosp, df_enfermeria); gc()
}

log_msg(sprintf("    N=%d, eventos=%d (%.2f%%)",
                nrow(df_final),
                sum(df_final$evento_mortalidad),
                100 * mean(df_final$evento_mortalidad)))

# 4. Conjuntos de predictores y particiones
log_msg("[4] Configurando modelos y folds")

basicos   <- c("edad_ingreso", "tipo_ingreso", "tipoproceso",
               "cat_braden", "cat_downton", "valor_barthel")
preguntas <- paste0("respuesta_pregunta", 1:7)

modelos <- list(
  "M1_baseline_conSexo"  = c(basicos, "sexo"),
  "M2_baseline_sinSexo"  = basicos,
  "M3_preguntas_conSexo" = c(basicos, "sexo", preguntas),
  "M4_preguntas_sinSexo" = c(basicos, preguntas)
)

modelos_label <- c(
  "M1_baseline_conSexo"  = "M1: básicas y sexo",
  "M2_baseline_sinSexo"  = "M2: básicas",
  "M3_preguntas_conSexo" = "M3: básicas, sexo y preguntas",
  "M4_preguntas_sinSexo" = "M4: básicas y preguntas"
)

# folds estratificados por el evento
K <- 5
set.seed(123)
folds <- integer(nrow(df_final))
idx_event   <- which(df_final$evento_mortalidad == 1)
idx_noevent <- which(df_final$evento_mortalidad == 0)
folds[idx_event]   <- sample(rep(1:K, length.out = length(idx_event)))
folds[idx_noevent] <- sample(rep(1:K, length.out = length(idx_noevent)))

log_msg(sprintf("    Eventos por fold: %s",
                paste(sapply(1:K, function(k)
                  sum(df_final$evento_mortalidad[folds == k])),
                  collapse = " ")))

# índice C en el fold de prueba
cindex_cox <- function(train, test, formula_i) {
  fit <- coxph(formula_i, data = train, x = TRUE)
  test$lp <- predict(fit, newdata = test, type = "lp")
  res <- concordance(Surv(tiempo_mortalidad, evento_mortalidad) ~ lp,
                     data = test, reverse = TRUE)$concordance
  rm(fit); gc(verbose = FALSE)
  res
}

cindex_rsf <- function(train, test, formula_i,
                       ntree = 500, nodesize = 15, mtry = NULL,
                       splitrule = "logrank") {
  args <- list(formula = formula_i, data = train,
               ntree = ntree, nodesize = nodesize,
               splitrule = splitrule,
               importance = "none", samptype = "swor")
  if (!is.null(mtry)) args$mtry <- mtry
  fit  <- do.call(rfsrc, args)
  pred <- predict(fit, newdata = test)
  test$risk <- pred$predicted
  res <- concordance(Surv(tiempo_mortalidad, evento_mortalidad) ~ risk,
                     data = test, reverse = TRUE)$concordance
  rm(fit, pred); gc(verbose = FALSE)
  res
}

# 5. Validación cruzada con el RSF por defecto
log_msg("[5] CV inicial: Cox vs RSF (default)")

cv_default <- list()
for (nombre in names(modelos)) {
  predictores <- modelos[[nombre]]
  formula_i <- as.formula(
    paste("Surv(tiempo_mortalidad, evento_mortalidad) ~",
          paste(predictores, collapse = " + "))
  )
  c_cox <- numeric(K); c_rsf <- numeric(K)
  cat("    -", nombre, "Fold:")
  for (k in 1:K) {
    cat(" ", k)
    train <- df_final[folds != k, ]
    test  <- df_final[folds == k, ]
    c_cox[k] <- cindex_cox(train, test, formula_i)
    c_rsf[k] <- cindex_rsf(train, test, formula_i, ntree = 500)
  }
  cat("\n")
  cv_default[[nombre]] <- tibble(Modelo = nombre, Fold = 1:K,
                                 Cox = c_cox, RSF = c_rsf)
}
cv_default_df <- bind_rows(cv_default)

tabla_default <- cv_default_df %>%
  group_by(Modelo) %>%
  summarise(
    Cox_mean = mean(Cox), Cox_sd = sd(Cox),
    RSF_mean = mean(RSF), RSF_sd = sd(RSF),
    Diferencia_RSF_Cox = mean(RSF) - mean(Cox),
    p_Wilcoxon = wilcox.test(RSF, Cox, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(Modelo = modelos_label[Modelo]) %>%
  mutate_if(is.numeric, ~round(., 4))

cat("\n===== RSF con hiperparámetros por defecto =====\n")
print(tabla_default)
write.csv(tabla_default,
          file.path(output_dir, "tabla_default.csv"), row.names = FALSE)
write.csv(cv_default_df,
          file.path(output_dir, "tabla_default_porfold.csv"), row.names = FALSE)
gc()


# 6. Búsqueda en rejilla sobre M3
log_msg("[6] Grid search RSF sobre M3")

grid <- expand.grid(
  ntree    = 500,
  nodesize = c(15, 30, 60),
  mtry     = c(3, 5, 7),
  stringsAsFactors = FALSE
) %>% as_tibble()
log_msg(sprintf("    %d combinaciones x %d folds = %d entrenamientos",
                nrow(grid), K, nrow(grid) * K))

formula_M3 <- as.formula(
  paste("Surv(tiempo_mortalidad, evento_mortalidad) ~",
        paste(modelos$M3_preguntas_conSexo, collapse = " + "))
)

grid_results <- list()
for (i in seq_len(nrow(grid))) {
  params <- grid[i, ]
  cat(sprintf("    [%2d/%d] ntree=%d nodesize=%d mtry=%d",
              i, nrow(grid), params$ntree, params$nodesize, params$mtry))
  c_vals <- numeric(K)
  for (k in 1:K) {
    train <- df_final[folds != k, ]
    test  <- df_final[folds == k, ]
    c_vals[k] <- cindex_rsf(train, test, formula_M3,
                            ntree    = params$ntree,
                            nodesize = params$nodesize,
                            mtry     = params$mtry)
  }
  grid_results[[i]] <- tibble(
    ntree = params$ntree, nodesize = params$nodesize, mtry = params$mtry,
    cindex_mean = mean(c_vals), cindex_sd = sd(c_vals)
  )
  cat(sprintf("  C = %.4f ± %.4f\n", mean(c_vals), sd(c_vals)))

  # guardado intermedio (recuperable si peta)
  saveRDS(grid_results, file.path(output_dir, "grid_intermedio.rds"))
  gc()
}

grid_df <- bind_rows(grid_results) %>% arrange(desc(cindex_mean))
cat("\n===== TOP 5 combinaciones =====\n")
print(head(grid_df, 5))
write.csv(grid_df, file.path(output_dir, "grid_search.csv"),
          row.names = FALSE)

best <- grid_df[1, ]
log_msg(sprintf("    Mejor: ntree=%d nodesize=%d mtry=%d, C = %.4f",
                best$ntree, best$nodesize, best$mtry, best$cindex_mean))

# 7. Los cuatro conjuntos con la mejor combinación
log_msg("[7] Re-evaluando 4 modelos con hiperparámetros tuned")

cv_tuned <- list()
for (nombre in names(modelos)) {
  predictores <- modelos[[nombre]]
  formula_i <- as.formula(
    paste("Surv(tiempo_mortalidad, evento_mortalidad) ~",
          paste(predictores, collapse = " + "))
  )
  c_cox <- numeric(K); c_rsf <- numeric(K)
  cat("    -", nombre, "Fold:")
  for (k in 1:K) {
    cat(" ", k)
    train <- df_final[folds != k, ]
    test  <- df_final[folds == k, ]
    c_cox[k] <- cindex_cox(train, test, formula_i)
    c_rsf[k] <- cindex_rsf(train, test, formula_i,
                           ntree    = best$ntree,
                           nodesize = best$nodesize,
                           mtry     = best$mtry)
  }
  cat("\n")
  cv_tuned[[nombre]] <- tibble(Modelo = nombre, Fold = 1:K,
                               Cox = c_cox, RSF_tuned = c_rsf)
}
cv_tuned_df <- bind_rows(cv_tuned)

tabla_tuned <- cv_tuned_df %>%
  group_by(Modelo) %>%
  summarise(
    Cox_mean = mean(Cox), Cox_sd = sd(Cox),
    RSF_tuned_mean = mean(RSF_tuned), RSF_tuned_sd = sd(RSF_tuned),
    Diferencia = mean(RSF_tuned) - mean(Cox),
    p_Wilcoxon = wilcox.test(RSF_tuned, Cox, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(Modelo = modelos_label[Modelo]) %>%
  mutate_if(is.numeric, ~round(., 4))

cat("\n===== RSF tuned vs Cox =====\n")
print(tabla_tuned)
write.csv(tabla_tuned, file.path(output_dir, "tabla_tuned.csv"),
          row.names = FALSE)
write.csv(cv_tuned_df, file.path(output_dir, "tabla_tuned_porfold.csv"),
          row.names = FALSE)
gc()


# 8. Gráficos
log_msg("[8] Generando gráficos")

tema_tfm <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position  = "top",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold")
  )

# índice C: Cox, RSF por defecto y RSF optimizado
df_comp <- bind_rows(
  cv_default_df %>%
    select(Modelo, Cox, RSF) %>%
    pivot_longer(c(Cox, RSF), names_to = "Algoritmo", values_to = "CIndex") %>%
    mutate(Config = ifelse(Algoritmo == "RSF", "RSF default", Algoritmo)),
  cv_tuned_df %>%
    select(Modelo, RSF_tuned) %>%
    pivot_longer(RSF_tuned, names_to = "Algoritmo", values_to = "CIndex") %>%
    mutate(Config = "RSF tuned")
) %>%
  mutate(Config = factor(Config, levels = c("Cox", "RSF default", "RSF tuned")),
         Modelo = modelos_label[Modelo]) %>%
  group_by(Modelo, Config) %>%
  summarise(media = mean(CIndex), sd = sd(CIndex), .groups = "drop")

p_comp <- ggplot(df_comp, aes(x = Modelo, y = media, fill = Config)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd),
                position = position_dodge(width = 0.75),
                width = 0.18, color = "grey30") +
  geom_text(aes(label = sprintf("%.3f", media)),
            position = position_dodge(width = 0.75),
            vjust = -1.6, size = 3.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Cox" = "#1F77B4",
                               "RSF default" = "#FFB347",
                               "RSF tuned"   = "#D95F0E")) +
  coord_cartesian(ylim = c(0.7, max(df_comp$media + df_comp$sd) + 0.04)) +
  labs(
    title    = "C-Index: Cox vs RSF (default vs tuned)",
    subtitle = "Mortalidad intrahospitalaria, media ± SD entre 5 folds",
    x = NULL, y = "C-Index"
  ) +
  tema_tfm +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(output_dir, "01_cindex_default_vs_tuned.png"),
       p_comp, width = 11, height = 6, dpi = 300)

# búsqueda en rejilla
p_grid <- ggplot(grid_df,
                 aes(x = factor(nodesize), y = cindex_mean,
                     color = factor(mtry), group = mtry)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = cindex_mean - cindex_sd,
                    ymax = cindex_mean + cindex_sd),
                width = 0.2, alpha = 0.5) +
  labs(
    title    = "Grid search RSF, M3 (mortalidad intrahospitalaria)",
    subtitle = "C-Index medio ± SD entre 5 folds",
    x = "nodesize", y = "C-Index", color = "mtry"
  ) +
  tema_tfm

ggsave(file.path(output_dir, "02_grid_search.png"),
       p_grid, width = 10, height = 5.5, dpi = 300)

# modelos con todos los datos, para la importancia de variables y el forest plot
log_msg("    Entrenando modelos finales para VIMP/forest")
modelos_full <- list()
for (nombre in names(modelos)) {
  predictores <- modelos[[nombre]]
  formula_i <- as.formula(
    paste("Surv(tiempo_mortalidad, evento_mortalidad) ~",
          paste(predictores, collapse = " + "))
  )
  cox_f <- coxph(formula_i, data = df_final, x = TRUE)
  rsf_f <- rfsrc(formula_i, data = df_final,
                 ntree = best$ntree,
                 nodesize = best$nodesize,
                 mtry = best$mtry,
                 importance = TRUE)
  modelos_full[[nombre]] <- list(cox = cox_f, rsf = rsf_f)
  rm(cox_f, rsf_f); gc()
}

# VIMP
df_vimp <- bind_rows(lapply(names(modelos_full), function(n) {
  v <- modelos_full[[n]]$rsf$importance
  tibble(Modelo = modelos_label[n], Variable = names(v), VIMP = as.numeric(v))
})) %>%
  group_by(Modelo) %>%
  mutate(Variable = reorder(Variable, VIMP)) %>%
  ungroup()

p_vimp <- ggplot(df_vimp, aes(x = VIMP, y = Variable, fill = VIMP > 0)) +
  geom_col() +
  facet_wrap(~ Modelo, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("TRUE" = "#2CA02C", "FALSE" = "#D62728"),
                    guide = "none") +
  labs(title = "VIMP, RSF tuned (mortalidad intrahospitalaria)",
       subtitle = "Verde: contribuye. Rojo: no aporta o empeora",
       x = "VIMP", y = NULL) +
  tema_tfm

ggsave(file.path(output_dir, "03_vimp_rsf.png"),
       p_vimp, width = 12, height = 9, dpi = 300)

# forest plot de Cox
df_cox <- bind_rows(lapply(names(modelos_full), function(n) {
  broom::tidy(modelos_full[[n]]$cox, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(Modelo = modelos_label[n])
}))

p_forest <- ggplot(df_cox, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high,
                      color = p.value < 0.05), size = 0.4) +
  facet_wrap(~ Modelo, scales = "free_y", ncol = 2) +
  scale_x_log10() +
  scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "grey50"),
                     labels = c("p ≥ 0.05", "p < 0.05"), name = NULL) +
  labs(title = "Hazard ratios, Cox (mortalidad intrahospitalaria)",
       subtitle = "Escala logarítmica, línea en HR = 1",
       x = "Hazard ratio (IC 95%)", y = NULL) +
  tema_tfm

ggsave(file.path(output_dir, "04_forest_cox.png"),
       p_forest, width = 12, height = 10, dpi = 300)

log_msg("Ejecución completa")
log_msg(sprintf("Resultados en: %s", normalizePath(output_dir)))
