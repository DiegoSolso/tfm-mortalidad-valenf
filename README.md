# Mortalidad intrahospitalaria y perfiles de riesgo de la valoración enfermera

Código del Trabajo Fin de Máster de Diego Solsona Ostal (Máster Universitario en Sistemas Inteligentes, Universitat Jaume I, 2026).

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
