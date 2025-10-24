# ANALISIS_MUNICIPAL: Implementación de Indice de Moran / ACP / GWPCA (a nivel municipal) para Cundinamarca

# ======================================================
# 1. Librerías
# ======================================================
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
# 1. Cargar datos originales y shapefile
# ======================================================
datos_norm <- read_excel("C:/Users/Natalia.LAPTOP-SRNRNE03/Downloads/Archivos TG/Datos_unificados_normalizados_VF.xlsx")

ruta_shp <- "C:/Users/Natalia.LAPTOP-SRNRNE03/Downloads/Archivos TG Reformulado 2/MGN_ADM_MPIO_GRAFICO.shp"
sf_municipal <- st_read(ruta_shp, quiet = FALSE)

# Filtrar solo Cundinamarca por código (25)
sf_muni_cund_completo <- sf_municipal %>%
  filter(dpto_ccdgo == "25")


# ======================================================
# 2. Generar variables epidemiológicas 
# ======================================================
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


# ======================================================
# 3. Agregación municipal + ANO (GTWPCA)
# ======================================================
datos_municipales_gtwpca <- datos_epi %>%
  filter(toupper(Departamento_ocurrencia) == "CUNDINAMARCA") %>%
  group_by(Departamento_ocurrencia, Municipio_ocurrencia, ANO) %>%
  summarise(
    n_casos = n(),
    edad_promedio = mean(as.numeric(EDAD), na.rm = TRUE),
    n_mujeres = sum(SEXO == "F", na.rm = TRUE),
    mediana_incubacion = median(tiempo_incubacion, na.rm = TRUE),
    mediana_notificacion = median(tiempo_notificacion, na.rm = TRUE),
    tiempo_total_notificacion = median(as.numeric(FEC_NOT - INI_SIN), na.rm = TRUE),
    n_casos_importados = sum(caso_importado, na.rm = TRUE),
    area_predominante = get_moda(AREA),
    estrato_moda = get_moda(estrato),
    n_area_cabecera = sum(AREA == "1", na.rm = TRUE),
    n_area_centro_poblado = sum(AREA == "2", na.rm = TRUE),
    n_area_rural_disperso = sum(AREA == "3", na.rm = TRUE),
    n_regimen_contributivo = sum(TIP_SS == "C", na.rm = TRUE),
    n_regimen_subsidiado = sum(TIP_SS == "S", na.rm = TRUE),
    n_regimen_excepcion = sum(TIP_SS == "P", na.rm = TRUE),
    n_regimen_especial = sum(TIP_SS == "E", na.rm = TRUE),
    n_no_asegurado = sum(TIP_SS == "N", na.rm = TRUE),
    n_regimen_indeterminado = sum(TIP_SS == "I", na.rm = TRUE),
    .groups = "drop"
  )


# ======================================================
# 4. Agregación municipal sin ANO (GWPCA)
# ======================================================
datos_municipales_gwpca <- datos_municipales_gtwpca %>%
  select(-ANO) %>%
  group_by(Departamento_ocurrencia, Municipio_ocurrencia) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# ======================================================
# 5. Unión con shapefile municipal (solo Cundinamarca)
# ======================================================
sf_municipal <- sf_municipal %>%
  mutate(
    dpto_cnmbr = toupper(iconv(dpto_cnmbr, to = "ASCII//TRANSLIT")),
    mpio_cnmbr = toupper(iconv(mpio_cnmbr, to = "ASCII//TRANSLIT"))
  )


datos_municipales_gwpca <- datos_municipales_gwpca %>%
  mutate(
    Departamento_ocurrencia = toupper(iconv(Departamento_ocurrencia, to = "ASCII//TRANSLIT")),
    Municipio_ocurrencia = toupper(iconv(Municipio_ocurrencia, to = "ASCII//TRANSLIT"))
  )

sf_muni_cund <- sf_municipal %>%
  filter(dpto_cnmbr == "CUNDINAMARCA") %>%
  inner_join(
    datos_municipales_gwpca,
    by = c("dpto_cnmbr" = "Departamento_ocurrencia",
           "mpio_cnmbr" = "Municipio_ocurrencia")
  )

# ======================================================
# 5.1. Matriz numérica y estandarizada municipal
# ======================================================
X_muni <- sf_muni_cund %>%
  st_drop_geometry() %>%
  dplyr::select(-dpto_cnmbr, -mpio_cnmbr)   # quitamos identificadores

# Seleccionar solo columnas numéricas
X_num_muni <- X_muni %>%
  select(where(is.numeric))

# Escalar matriz numérica
X_scaled_muni <- X_num_muni %>%
  as.matrix() %>%
  scale()

# Coordenadas para kernels / Moran
coords_muni <- st_coordinates(st_centroid(sf_muni_cund)) %>%
  as.data.frame()


# ======================================================
# 6. Índice de Moran Global municipal para Cundinamarca
# ======================================================

# Centroides
coords <- st_coordinates(st_centroid(sf_muni_cund))

# Vecinos contiguos
nb_contig <- poly2nb(sf_muni_cund, queen = TRUE, snap = 1e-6)

# Vecinos KNN (para nodos aislados)
nb_knn <- knn2nb(knearneigh(coords, k = 2))

# Combinar: si no tiene vecinos contiguos, usar KNN
nb_final <- nb_contig
for (i in seq_along(nb_final)) {
  if (length(nb_final[[i]]) == 0) {
    nb_final[[i]] <- nb_knn[[i]]
  }
}

# Lista de pesos espaciales
lw <- nb2listw(nb_final, style = "W", zero.policy = TRUE)

# Variables que no queremos usar en ACP
vars_geom <- c("shape_area", "shape_Leng", "shape_Area","mpio_narea")

# Aplicar Moran
X_sig_muni <- calcular_moran(X_num_muni, lw, alpha = 0.05)

# Quitar las geométricas
X_sig_muni <- X_sig_muni %>%
  dplyr::select(-any_of(vars_geom))

# ======================================================
# 7. ACP global solo con variables espaciales significativas
# ======================================================
pca_global_sig <- prcomp(X_sig_muni, center = TRUE, scale. = TRUE)

# Resumen del ACP
summary(pca_global_sig)
cargas_acp <- pca_global_sig$rotation
print("Cargas ACP Global:")
print(cargas_acp)

# ======================================================
# 8. Selección de kernel y ancho de banda (fijo y adaptativo)
# ======================================================
dist_mat_muni <- as.matrix(dist(coords_muni))
X_clean <- X_scaled_muni[, apply(X_scaled_muni, 2, function(x) all(is.finite(x)))]

X_clean <- X_clean[apply(X_clean, 1, function(x) all(is.finite(x))), ]

eigen_global <- eigen(cov(X_clean), symmetric = TRUE)$values

best_muni <- tune_LOO_robusta(X_clean, coords_muni, dist_mat_muni, eigen_global = eigen_global)

# ======================================================
# 9. GWPCA
# ======================================================
coords_muni_mat <- as.matrix(coords_muni)
stopifnot(all(is.finite(coords_muni_mat)))

# GWPCA usando la matriz limpia (sin argumento k)
resultados_gwpca_muni <- gwpca(
  X = X_clean,                   
  coordenadas = coords_muni_mat, 
  ancho_banda = 1.62845178,
  kernel = "boxcar"
)

# ======================================================
# 10. Indicador GWPCA
# ======================================================

# Extraer varianza local (2 primeras componentes)
var_local <- do.call(rbind, resultados_gwpca_muni$varianza_explicada)[, 1:2, drop = FALSE]
ponderaciones <- var_local / rowSums(var_local, na.rm = TRUE)

# Extraer y combinar las puntuaciones locales
scores_locales_list <- resultados_gwpca_muni$puntuaciones
scores_locales <- do.call(rbind, scores_locales_list)

# Verificar que tengan misma longitud
n_obs <- min(nrow(scores_locales), nrow(ponderaciones))
scores_locales <- scores_locales[1:n_obs, , drop = FALSE]
ponderaciones <- ponderaciones[1:n_obs, , drop = FALSE]

# Calcular el indicador ponderado (2 primeras componentes)
indicador_gwpca <- rowSums(scores_locales[, 1:2] * ponderaciones, na.rm = TRUE)

# Escalar de 0 a 100 (más interpretativo)
indicador_gwpca_esc <- scales::rescale(indicador_gwpca, to = c(0, 100))

# Agregar al objeto espacial
sf_muni_cund_completo$indicador_gwpca <- NA
sf_muni_cund_completo$indicador_gwpca[1:length(indicador_gwpca_esc)] <- indicador_gwpca_esc

# ======================================================
# 11. Mapa del indicador
# ======================================================
ggplot(sf_muni_cund_completo) +
  geom_sf(aes(fill = indicador_gwpca)) +
  scale_fill_viridis_c(option = "plasma", name = "Indicador GWPCA (%)") +
  theme_minimal() +
  labs(
    title = "Indicador Espacial GWPCA (2 primeras componentes)",
    subtitle = "Municipios de Cundinamarca (escala 0–100)"
  )

# ======================================================
# 12. Mapa PVT (% varianza explicada local)
# ======================================================
mapa_pvt_muni <- function(resultados_gwpca, sf_obj, id_col = "mpio_cnmbr", n_comp = 2) {
  var_list <- resultados_gwpca$varianza_explicada
  var_exp <- do.call(rbind, var_list)
  var_exp <- apply(var_exp, 2, as.numeric)
  
  n_comp_real <- min(n_comp, ncol(var_exp))
  pvt_loc <- rowSums(var_exp[, 1:n_comp_real, drop = FALSE], na.rm = TRUE) * 100
  
  # Municipios analizados
  ids_analizados <- sf_obj[[id_col]][1:length(pvt_loc)]
  df_pvt <- data.frame(
    pvt_pct = pvt_loc
  )
  df_pvt[[id_col]] <- ids_analizados
  
  # Unir al shapefile
  sf_pvt <- sf_obj %>%
    left_join(df_pvt, by = id_col)
  
  ggplot(sf_pvt) +
    geom_sf(aes(fill = pvt_pct), color = "white", size = 0.2) +
    scale_fill_viridis_c(
      option = "cividis",
      name = "PVT (%)",
      labels = scales::label_percent(accuracy = 0.1, scale = 1),
      na.value = "grey90"
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top",       
        title.hjust = 0.5,            
        title.theme = element_text(size = 14),  # ← tamaño del título más pequeño
        label.theme = element_text(size = 13)                   # ← tamaño de los porcentajes
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
}

# Ejecutar
pvt_plot_muni <- mapa_pvt_muni(resultados_gwpca_muni, sf_muni_cund_completo, n_comp = 2)

ggsave(
  filename = "PVT_GWPCA_Cundinamarca.png", 
  plot = pvt_plot_muni, width = 14, height = 6, dpi = 300
)



