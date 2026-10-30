# Mortalidad intrahospitalaria y perfiles de riesgo de la valoración enfermera

Código de mi Trabajo Fin de Máster (Máster Universitario en Sistemas Inteligentes, Universitat Jaume I, 2026), desarrollado en el grupo de investigación GIANT con datos del Hospital Universitario de La Plana.

**Pregunta:** ¿el perfil de riesgo que se asigna en la primera valoración enfermera (instrumento VALENF) permite predecir la mortalidad durante el ingreso?

## Resumen

- **Datos:** 8.673 episodios de hospitalización de adultos, de marzo de 2023 a febrero de 2024. En 426 (4,9 %) el paciente falleció; el alta con vida se trata como censura.
- **Estudio preliminar:** la regresión de Cox (C-index 0,852) discriminó algo mejor que el Random Survival Forest (0,842), así que el trabajo siguió con modelos estadísticos.
- **Modelo final:** un modelo loglineal de Weibull con el perfil, cuatro constantes vitales (frecuencia cardiaca, presión arterial sistólica y diastólica y saturación de oxígeno) y el tipo de ingreso y de proceso.
- **Resultados:** el perfil fue la variable que más aportó. El modelo obtuvo un C-index de 0,832 con validación cruzada y una razón de muertes observadas y esperadas (O/E) de 1,00 a 7 días.

## Scripts

| Script | Descripción |
|---|---|
| `descriptivo_datos.R` | Filtros y análisis descriptivo de los datos. |
| `estudio_preliminar_cox_rsf.R` | Regresión de Cox frente a RSF para la mortalidad. |
| `estudio_preliminar_distribuciones.R` | Comparación de distribuciones del estudio preliminar. |
| `eleccion_modelo.R` | Riesgos proporcionales, criterio de la recta y otras distribuciones. |
| `modelo_perfiles.R` | Modelo de Weibull final y su validación. |
| `comparacion_weibull_rsf.R` | Modelo de Weibull frente al RSF. |
| `comprobaciones_discusion.R` | Comprobaciones que se citan en la discusión. |
| `cox_julia.jl` | Regresión de Cox del estudio preliminar en Julia (cuaderno de Pluto). |

## Requisitos

- **R**, con los paquetes `survival`, `randomForestSRC`, `dplyr`, `tidyr`, `lubridate`, `ggplot2`, `broom`, `DBI` y `RPostgres`.
- **Julia**, solo para `cox_julia.jl`, con `CSV`, `DataFrames`, `CategoricalArrays`, `StatsModels` y `Survival`.

## Datos

Los datos no se incluyen porque contienen información clínica de los pacientes. Los scripts los leen de una carpeta `datos/`:

- `capacidad_funcional_1VALENF_LARGO_con_constantes.csv`: fichero de trabajo, para todos los scripts de R salvo `estudio_preliminar_cox_rsf.R`.
- `valenf.csv`: valoraciones VALENF, para `estudio_preliminar_cox_rsf.R`.
- `df_final_v6.csv`: tabla de episodios del estudio preliminar (`graficos_v6/df_final.rds` guardada como CSV), para `cox_julia.jl`.

`estudio_preliminar_cox_rsf.R` también consulta la base de datos del proyecto. Las credenciales se leen de variables de entorno, por ejemplo desde `~/.Renviron`:

```
VALENF_BD=...
VALENF_HOST=...
VALENF_USUARIO=...
VALENF_CLAVE=...
```

## Autor

Diego Solsona Ostal · [LinkedIn](https://www.linkedin.com/in/diego-solsona-dev)

Trabajo dirigido por Óscar Belmonte Fernández y Reyes Grangel Seguer.
