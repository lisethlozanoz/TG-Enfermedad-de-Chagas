library(readxl)
library(dplyr)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(spdep)
library(ggplot2)
library(scales)
library(GWmodel)
library(zoo)
library(tidyr)

source("C:/Users/Natalia.LAPTOP-SRNRNE03/Downloads/Archivos TG/Reorganizados/utils_graficos.R")
source("C:/Users/Natalia.LAPTOP-SRNRNE03/Downloads/Archivos TG/Reorganizados/utils.R")

# ======================================================
# 2. Cargar y preparar datos epidemiológicos
# ======================================================
# Cargar los datos originales
datos_norm <- read_excel("C:/Users/Natalia.LAPTOP-SRNRNE03/Downloads/Archivos TG/Datos_unificados_normalizados_VF.xlsx")
str(datos_norm)

# Crear nuevas columnas de datos en el df
datos_epi <- datos_norm %>%
  mutate(
    FEC_NOT = ymd(FEC_NOT), 
    FEC_CON = ymd(FEC_CON),
    INI_SIN = ymd(INI_SIN),
    tiempo_incubacion = as.numeric(FEC_CON - INI_SIN), # Permite encontrar patrones de progresión de la enfermedad
    tiempo_notificacion = as.numeric(FEC_NOT - FEC_CON), # Permite encontrar demoras en el sistema de salud y posibles retrasos de diagnóstico
    caso_importado = Departamento_ocurrencia != Departamento_Notificacion, # Identifica si el contagio ocurrió fuera del departamento de notificación
    # Revela movilidad de la enfermedad y focos de transmisión
    mes_ocurrencia = month(FEC_CON) # Extraer el mes numérico (1-12) de la fecha de contagio (FEC_CON) para analizar estacionalidad
  )


datos_agregados <- datos_epi %>%
  group_by(Departamento_ocurrencia) %>%
  summarise(
    # Carga de enfermedad
    n_casos = n(),
    
    # Características demográficas
    edad_promedio = mean(as.numeric(EDAD), na.rm = TRUE),
    n_mujeres = sum(SEXO == "F", na.rm = TRUE),   
    
    # Indicadores temporales
    mediana_incubacion = median(tiempo_incubacion, na.rm = TRUE),
    mediana_notificacion = median(tiempo_notificacion, na.rm = TRUE),
    tiempo_total_notificacion = median(as.numeric(FEC_NOT - INI_SIN), na.rm = TRUE),
    
    # Patrones espaciales
    n_casos_importados = sum(caso_importado, na.rm = TRUE),  
    
    # Distribución por tipo de área
    # usamos el factor predominante (moda) y además los conteos
    area_predominante = as.numeric(names(which.max(table(AREA)))),
    n_area_cabecera = sum(AREA == "1", na.rm = TRUE),
    n_area_centro_poblado = sum(AREA == "2", na.rm = TRUE),
    n_area_rural_disperso = sum(AREA == "3", na.rm = TRUE),
    
    # Condiciones socioeconómicas
    estrato_moda = as.numeric(names(which.max(table(estrato)))), 
    # quitamos estrato_promedio por redundancia
    
    # Cobertura en salud (factores en vez de proporciones)
    n_regimen_contributivo = sum(TIP_SS == "C", na.rm = TRUE),
    n_regimen_subsidiado = sum(TIP_SS == "S", na.rm = TRUE),
    n_regimen_excepcion = sum(TIP_SS == "P", na.rm = TRUE),
    n_regimen_especial = sum(TIP_SS == "E", na.rm = TRUE),
    n_no_asegurado = sum(TIP_SS == "N", na.rm = TRUE),
    n_regimen_indeterminado = sum(TIP_SS == "I", na.rm = TRUE)
  ) %>%
  ungroup()

# ======================================================
# 4. Obtener coordenadas departamentales
# ======================================================
colombia <- ne_states(country = "colombia", returnclass = "sf") %>% # Carga los límites geográficos de los departamentos de Colombia en formato sf
  mutate(depto_norm = toupper(iconv(name, to = "ASCII//TRANSLIT"))) %>% # En depto_norm quedan los departamentos normalizados para que coincidan con mis datos
  mutate(depto_norm = gsub("´|'|\"", "", depto_norm)) %>%
  mutate(depto_norm = case_when(
    depto_norm == "GUAJIRA" ~ "LA GUAJIRA",
    depto_norm == "NORTE SANTANDER" ~ "NORTE DE SANTANDER",
    depto_norm == "BOGOTA" ~ "BOGOTÁ D.C.",
    TRUE ~ depto_norm
  ))

# Filtrar del sf Colombia solo departamentos que estén en mis datos
colombia <- colombia %>% 
  filter(depto_norm %in% datos_agregados$Departamento_ocurrencia)

# Calcular el centro geográfico (centroide) de cada departamento para tener como referencia en las técnicas
centroides <- st_centroid(colombia) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(Departamento_ocurrencia = colombia$depto_norm) # Df resultante con longitud,latitud

# ======================================================
# 5. Unir datos y preparar matriz
# ======================================================
# Agregar las coordenadas de los centroides a mis datos
datos_completos <- datos_agregados %>%
  left_join(centroides, by = "Departamento_ocurrencia") %>%
  filter(complete.cases(.))

variables_analisis <- c(
  "n_casos", 
  "edad_promedio", 
  "mediana_incubacion", 
  "mediana_notificacion", 
  "tiempo_total_notificacion",
  "n_mujeres", 
  "n_casos_importados", 
  "n_area_cabecera", 
  "n_area_centro_poblado", 
  "n_area_rural_disperso",
  "area_predominante",   
  "estrato_moda",        
  "n_regimen_contributivo",
  "n_regimen_subsidiado",
  "n_regimen_excepcion",
  "n_regimen_especial",
  "n_no_asegurado",
  "n_regimen_indeterminado"
)


# Preparar matriz de datos para GWpca
X <- datos_completos %>%
  select(all_of(variables_analisis))


X <- X[, apply(X, 2, function(col) sd(col, na.rm = TRUE) > 0)] 
X <- scale(X) %>% as.matrix() # Estandarizar

# ======================================================
# 6. Índice de Moran 
# ======================================================
# Crear lista de pesos espaciales usando centroides y k vecinos más cercanos
coords <- datos_completos[, c("X", "Y")] %>% as.matrix() # Matriz de coordenadas (longitud, latitud)
nb <- knearneigh(coords, k = 4) |> knn2nb() # Identifica los k vecinos más cercanos para cada departamento basado en sus coordenadas
lw <- nb2listw(nb, style = "W", zero.policy = TRUE) # Zero policy es TRUE para admitir que tenga departamentos sin vecinos -> Ej San Andres

# Calcular para todas las variables en X
resultados_moran <- data.frame(
  Variable = colnames(X),
  Moran_I = NA,
  p_value = NA
)

for (i in seq_along(colnames(X))) {
  var_name <- colnames(X)[i]
  test <- moran.test(datos_completos[[var_name]], lw)
  resultados_moran$Moran_I[i] <- test$estimate[["Moran I statistic"]]
  resultados_moran$p_value[i] <- test$p.value
}

print(resultados_moran)

vars_significativas <- resultados_moran %>%
  filter(p_value < 0.05) %>%
  pull(Variable)

# Crear nueva matriz solo con las variables espaciales significativas
X_sig <- X[, vars_significativas, drop = FALSE]

# ======================================================
# ACP global solo con variables espaciales significativas
# ======================================================
pca_global_sig <- prcomp(X_sig, center = TRUE, scale. = TRUE)

# Resumen del ACP
summary(pca_global_sig)
cargas_acp <- pca_global_sig$rotation
print("Cargas ACP Global:")
print(cargas_acp)

#---------------------------------------------------
# Selección de kernel y ancho de banda (fijo y adaptativo)
#---------------------------------------------------

dist_mat <- sp::spDists(coords, longlat = TRUE)
X_scaled <- scale(X_sig)
X_clean <- X_scaled[, colSums(is.finite(X_scaled)) == nrow(X_scaled), drop = FALSE]
X_clean <- X_clean[apply(X_clean, 1, function(r) all(is.finite(r))), , drop = FALSE]

eigen_global <- eigen(cov(X_clean), symmetric = TRUE)$values

bw_candidates <- seq(from = min(dist_mat[dist_mat > 0]) * 1.05, to = max(dist_mat), length.out = 8)
kernels_try <- c("bisquare", "gaussian", "exponential", "boxcar")

best_table <- tune_LOO_robusta(X_clean, coords, dist_mat, kernels = kernels_try, bw_candidates = bw_candidates, eigen_global = eigen_global)
best_row <- best_table[1, ]
best_kernel <- as.character(best_row$kernel)
best_bw <- as.numeric(best_row$bw)
cat("Mejor kernel:", best_kernel, "\n")
cat("Mejor ancho de banda:", best_bw, "\n")

#---------------------------------------------------
# Ejecución
#---------------------------------------------------

# PCA global
pca_global <- prcomp(X, scale. = TRUE)
eigen_global <- pca_global$sdev^2

# Distancias
dist_mat <- spDists(coords, longlat = TRUE)

# Búsqueda coarse
bw_coarse <- seq(min(dist_mat[dist_mat > 0]), max(dist_mat), length.out = 8)
kernels_try <- c("bisquare","gaussian")
res_coarse <- tune_LOO(X, coords, dist_mat, kernels_try, bw_coarse)

# Refinar
best0 <- res_coarse[1, ]
bw_best <- as.numeric(best0$bw)
bw_fine <- seq(bw_best*0.8, bw_best*1.2, length.out = 8)
res_fine <- tune_LOO(X, coords, dist_mat, as.character(best0$kernel), bw_fine)

# Resultado final
res_best <- rbind(res_coarse, res_fine) %>% arrange(GOF)
best_final <- res_best[1, ]
print(best_final)

# ======================================================
# 9. Ejecutar GWPCA
# ======================================================
# Mostrar resumen
resumen_gwpca <- crear_resumen_gwpca(resultados_gwpca)
print("Resumen GWPCA - Varianza explicada promedio por componente y rangos:")
print(resumen_gwpca)

# Extraer cargas (vectores propios) de cada localización
cargas_gwpca <- lapply(resultados_gwpca$autovalores, function(x) x$vectores)

# Promedio de cargas por componente en todas las localizaciones
cargas_promedio_gwpca <- Reduce("+", cargas_gwpca) / length(cargas_gwpca)

print("Cargas promedio del GWPCA:")
print(cargas_promedio_gwpca)

