# TG_Enfermedad-de-Chagas
Implementación de las técnicas GWPCA (Geographically Weighted Principal Component Analysis) y GTWPCA (Geographically and Temporally Weighted Principal Component Analysis) para el análisis espacio-temporal de la enfermedad de Chagas en Colombia (2021–2023).

Este trabajo forma parte de un trabajo de grado orientado al estudio de patrones espaciales y su evolución temporal a partir de variables socioambientales y epidemiológicas.

## Contenido del repositorio

- **`analisis_nacional_GWPCA.R`** — Implementación nacional del método **GWPCA** utilizando las funciones definidas en `utils.R`.  
- **`analisis_nacional_GTWPCA.R`** — Implementación nacional del método **GTWPCA** utilizando las funciones definidas en `utils.R`.  
- **`analisis_municipal.R`** — Análisis municipal con **GWPCA** (no se aplica GTWPCA por falta de datos).  
- **`utils.R`** — Definición de las funciones principales para la implementación de los métodos **GWPCA** y **GTWPCA**.

---

## Requisitos

Este proyecto fue desarrollado en **R** y utiliza principalmente las siguientes librerías:
`dplyr`, `sf`, `ggplot2`, `rnaturalearth`, `GWmodel` y `gtwr`.

---

## Autora

**Liseth Lozano**  
Estudiante de Ciencia de Datos e Ingeniería de Sistemas — Colombia
