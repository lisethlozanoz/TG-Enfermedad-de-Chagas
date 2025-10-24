library(GWmodel)
library(dplyr)
library(readxl)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(lubridate)
library(ggplot2)
library(patchwork)
library(viridis)

# ==============================
# 1. Cargar y preparar datos
# ==============================
datos_norm <- read_excel("C:/Users/Natalia.LAPTOP-SRNRNE03/Downloads/Archivos TG/Datos_unificados_normalizados_VF.xlsx")

datos_epi <- datos_norm %>%
  mutate(
    FEC_NOT = ymd(FEC_NOT),
    FEC_CON = ymd(FEC_CON),
    INI_SIN = ymd(INI_SIN),
    tiempo_incubacion = as.numeric(FEC_CON - INI_SIN),
    tiempo_notificacion = as.numeric(FEC_NOT - FEC_CON),
    caso_importado = Departamento_ocurrencia != Departamento_Notificacion,
    mes_ocurrencia = month(FEC_CON)
  )

# Usar directamente la columna ANO del archivo
datos_completos <- datos_epi %>%
  group_by(Departamento_ocurrencia, ANO) %>%
  reframe(
    n_casos = n(),
    edad_promedio = mean(as.numeric(EDAD), na.rm = TRUE),
    n_mujeres = sum(SEXO == "F", na.rm = TRUE),
    mediana_incubacion = median(tiempo_incubacion, na.rm = TRUE),
    mediana_notificacion = median(tiempo_notificacion, na.rm = TRUE),
    tiempo_total_notificacion = median(as.numeric(FEC_NOT - INI_SIN), na.rm = TRUE),
    n_casos_importados = sum(caso_importado, na.rm = TRUE),
    area_predominante = as.numeric(names(which.max(table(AREA)))),
    n_area_cabecera = sum(AREA == "1", na.rm = TRUE),
    estrato_moda = as.numeric(names(which.max(table(estrato)))),
    n_regimen_contributivo = sum(TIP_SS == "C", na.rm = TRUE),
    n_regimen_excepcion = sum(TIP_SS == "P", na.rm = TRUE)
  )

# ==============================
# 2. Coordenadas departamentales
# ==============================
colombia <- ne_states(country = "Colombia", returnclass = "sf") %>%
  mutate(depto_norm = toupper(iconv(name, to = "ASCII//TRANSLIT"))) %>%
  mutate(depto_norm = gsub("´|'|\"", "", depto_norm)) %>%
  mutate(depto_norm = case_when(
    depto_norm == "GUAJIRA" ~ "LA GUAJIRA",
    depto_norm == "NORTE SANTANDER" ~ "NORTE DE SANTANDER",
    depto_norm == "BOGOTA" ~ "BOGOTÁ D.C.",
    TRUE ~ depto_norm
  )) %>%
  filter(depto_norm %in% datos_completos$Departamento_ocurrencia)

centroides <- st_centroid(colombia) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(Departamento_ocurrencia = colombia$depto_norm)

datos_completos <- datos_completos %>%
  left_join(centroides, by = "Departamento_ocurrencia") %>%
  filter(complete.cases(.))

# ==============================
# 3. Preparar matriz GTWPCA
# ==============================
vars_significativas <- c(
  "n_casos", "edad_promedio", "n_mujeres",
  "n_casos_importados", "n_area_cabecera",
  "estrato_moda", "n_regimen_contributivo", "n_regimen_excepcion"
)

X <- datos_completos %>% select(all_of(vars_significativas))
X <- X[, apply(X, 2, function(col) sd(col, na.rm = TRUE) > 0)]
X <- scale(X) %>% as.matrix()
coords <- datos_completos[, c("X", "Y")] %>% as.matrix()
times <- datos_completos$ANO
times <- as.numeric(as.character(times))

# ==============================
# 4. Validación cruzada para parámetros óptimos
# ==============================
alpha_vals <- c(0.5, 1, 2)
beta_vals  <- c(0.5, 1, 2)
bandwidth_vals <- seq(3, 10, 1)

set.seed(123)
cv_result <- cv_gtwpca(X, coords, times, alpha_vals, beta_vals, bandwidth_vals)
alpha_opt <- cv_result$best_params$alpha
beta_opt  <- cv_result$best_params$beta
bw_opt    <- cv_result$best_params$bandwidth

cat("Parámetros óptimos:\n")
print(cv_result$best_params)
cat("Error mínimo:", cv_result$best_error, "\n")

# ==============================
# 6. Ejecutar GTWPCA con parámetros óptimos
# ==============================
times <- as.numeric(datos_completos$ANO)
resultados_gtwpca <- gtwpca(X, coords, times, alpha_opt, beta_opt, bw_opt)
