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

# ======================================================
# MAPAS PVT 
# ======================================================

library(ggplot2)
library(dplyr)
library(sf)
library(patchwork)

# 1. Años en orden
anos_unicos <- sort(unique(datos_completos$ANO))  # Esto asegura 2021, 2022, 2023 en orden
pvt_todos <- list()

for (ano_map in anos_unicos) {
  idx_ano <- which(datos_completos$ANO == ano_map)
  var_locales <- t(sapply(resultados_gtwpca$pvt[idx_ano], function(x) x[1:3]))
  var_acum <- rowSums(var_locales) * 100
  
  pvt_todos[[as.character(ano_map)]] <- var_acum
}

# 2. Calcular min y max global
pvt_global <- unlist(pvt_todos)
min_pvt <- min(pvt_global, na.rm = TRUE)
max_pvt <- max(pvt_global, na.rm = TRUE)

# 3. Crear mapas con la misma escala de colores
mapas_pvt <- list()
colombia_base <- st_as_sf(colombia)

for (ano_map in anos_unicos) {
  idx_ano <- which(datos_completos$ANO == ano_map)
  var_acum <- pvt_todos[[as.character(ano_map)]]
  
  colombia_datos <- colombia_base %>%
    filter(depto_norm %in% datos_completos$Departamento_ocurrencia[idx_ano]) %>%
    mutate(PVT_acumulada = var_acum)
  
  mapa <- ggplot() +
    geom_sf(data = colombia_base, fill = "grey90", color = "white") +
    geom_sf(data = colombia_datos, aes(fill = PVT_acumulada), color = "white") +
    scale_fill_viridis_c(option = "cividis", name = "PVT (%)", limits = c(min_pvt, max_pvt)) +
    labs(title = paste(ano_map)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),      
      axis.text = element_blank(),       
      axis.ticks = element_blank(),      
      axis.title = element_blank(),      
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  
  mapas_pvt[[as.character(ano_map)]] <- mapa
}

# 4. Unir todos los mapas en una sola imagen con orden correcto
mapa_combinado <- wrap_plots(mapas_pvt, ncol = length(mapas_pvt))

# 5. Mostrar en pantalla
mapa_combinado

# 6. Guardar la imagen
ggsave("mapa_PVT_GTWPCA.png", mapa_combinado, width = 12, height = 4, dpi = 300)

# ==============================
# Mapa categórico variable ganadora - primera componente
# ==============================
primer_componente <- sapply(resultados_gtwpca$eigens, function(x) x$vectors[,1])

var_ganadora <- apply(primer_componente, 2, function(col) {
  names(col)[which.max(abs(col))]
})

colombia_sf <- st_as_sf(colombia) %>%
  mutate(
    Variable_ganadora = factor(var_ganadora)
  )

n_vars <- length(unique(colombia_sf$Variable_ganadora))
colores <- viridis::viridis(n_vars, option = "cividis")

mapa_ganadora_categorica <- ggplot() +
  geom_sf(data = colombia_sf, fill = "grey90", color = "white") +
  geom_sf(data = colombia_sf %>% filter(!is.na(Variable_ganadora)),
          aes(fill = Variable_ganadora), color = "white") +
  scale_fill_manual(values = colores, name = "Variable ganadora") +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "gray40")
  ) +
  labs(x = "Longitud", y = "Latitud", title = "Variable ganadora - primera componente")

ggsave("Variable_ganadora_CP1_categorica.png", mapa_ganadora_categorica, width = 14, height = 6, dpi = 300)

cat("✅ Mapas PVT y variable ganadora generados correctamente.\n")

# Calcular PVT acumulado para las 3 primeras componentes por año y departamento
anos_unicos <- unique(datos_completos$ANO)

# Lista para guardar resultados
pvt_acumulado <- list()

for (ano_map in anos_unicos) {
  idx_ano <- which(datos_completos$ANO == ano_map)
  
  # Extraer PVT para las 3 primeras componentes
  var_locales <- t(sapply(resultados_gtwpca$pvt[idx_ano], function(x) x[1:3]))
  var_acum <- rowSums(var_locales) * 100  # convertir a porcentaje
  
  # Crear dataframe con resultados
  df_res <- data.frame(
    Departamento = datos_completos$Departamento_ocurrencia[idx_ano],
    Ano = ano_map,
    PVT_3Comp = var_acum,
    Llega_100 = ifelse(var_acum >= 100, "Sí", "No")
  )
  
  pvt_acumulado[[as.character(ano_map)]] <- df_res
}

# Unir todos los años en un solo data.frame
pvt_acumulado_df <- do.call(rbind, pvt_acumulado)

# Ver primeros resultados
head(pvt_acumulado_df)

# ======================================================
# MAPAS DEL INDICADOR GTWPCA (PVT acumulada)
# ======================================================
anos_unicos <- sort(unique(datos_epi$ANO))
pvt_todos <- list()

# 2. Calcular PVT acumulada (3 primeras comp.) por año
for (ano_map in anos_unicos) {
  idx_ano <- which(datos_completos$ANO == ano_map)
  var_locales <- t(sapply(resultados_gtwpca$pvt[idx_ano], function(x) x[1:3]))
  pvt_todos[[as.character(ano_map)]] <- rowSums(var_locales)
}

pvt_global <- unlist(pvt_todos)
min_pvt <- min(pvt_global, na.rm = TRUE)
max_pvt <- max(pvt_global, na.rm = TRUE)

mapas_pvt <- list()
colombia_base <- st_as_sf(colombia)

for (ano_map in anos_unicos) {
  idx_ano <- which(datos_completos$ANO == ano_map)
  var_acum <- pvt_todos[[as.character(ano_map)]]
  
  # Asegurarse de no tener NAs
  var_acum[is.na(var_acum)] <- 0
  
  colombia_datos <- colombia_base %>%
    filter(depto_norm %in% datos_completos$Departamento_ocurrencia[idx_ano]) %>%
    mutate(PVT_acumulada = var_acum)
  
  mapa <- ggplot() +
    geom_sf(data = colombia_base, fill = "grey90", color = "white", size = 0.2) +
    geom_sf(data = colombia_datos, aes(fill = PVT_acumulada), color = "white", size = 0.2) +
    scale_fill_viridis(
      option = "cividis",
      direction = 1,
      name = "PVT (%)",
      labels = scales::percent_format(accuracy = 1)
    )+
    labs(title = as.character(ano_map)) +
    theme_void() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  mapas_pvt[[as.character(ano_map)]] <- mapa
}

mapa_combinado_gtwpca <- wrap_plots(mapas_pvt, ncol = length(mapas_pvt)) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

print(mapa_combinado_gtwpca)
ggsave("mapa_GTWPCA_PVT.png", mapa_combinado_gtwpca, width = 12, height = 4, dpi = 300)

# ==============================
# 7. Indicador espacio-temporal (ICP)
# ==============================
# Inicializar vector para guardar el indicador
ICP <- numeric(nrow(X))

for (i in 1:nrow(X)) {
  # Extraer autovalores y proporción de varianza por componente
  eig_i <- resultados_gtwpca$eigens[[i]]
  lambda_i <- eig_i$values              # autovalores locales
  alpha_jk <- eig_i$vectors^2           # carga al cuadrado como ponderación
  # Normalizar alpha_jk según la fórmula
  numerador <- colSums(alpha_jk * lambda_i)
  denominador <- sum(alpha_jk * lambda_i)
  alpha_j <- numerador / denominador    # ponderación α_j(u_i, v_i, t_i)
  
  # Calcular I(u_i, v_i, t_i)
  ICP[i] <- sum(alpha_j * X[i, ])
}

# Normalizar el indicador 0–100
ICP <- (ICP - min(ICP, na.rm = TRUE)) / (max(ICP, na.rm = TRUE) - min(ICP, na.rm = TRUE)) * 100

# Guardar en el data.frame completo
datos_completos$ICP_GTWPCA <- ICP

# Verificar resultados
head(datos_completos[, c("Departamento_ocurrencia", "ANO", "ICP_GTWPCA")])

# ==============================
# 8. Mapas ICP espacio-temporal
# ==============================
mapas_ICP <- list()
anos_unicos <- sort(unique(datos_completos$ANO))
colombia_base <- st_as_sf(colombia)

# Rango global del ICP para escala uniforme
min_ICP <- min(datos_completos$ICP_GTWPCA, na.rm = TRUE)
max_ICP <- max(datos_completos$ICP_GTWPCA, na.rm = TRUE)

for (ano_map in anos_unicos) {
  idx_ano <- which(datos_completos$ANO == ano_map)
  
  colombia_datos <- colombia_base %>%
    filter(depto_norm %in% datos_completos$Departamento_ocurrencia[idx_ano]) %>%
    mutate(ICP = datos_completos$ICP_GTWPCA[idx_ano])
  
  mapa <- ggplot() +
    geom_sf(data = colombia_base, fill = "grey90", color = "white", size = 0.2) +
    geom_sf(data = colombia_datos, aes(fill = ICP), color = "white", size = 0.2) +
    scale_fill_viridis_c(
      option = "cividis",
      name = "Nivel del indicador (%)",
      limits = c(min_ICP, max_ICP)
    ) +
    labs(title = as.character(ano_map)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  mapas_ICP[[as.character(ano_map)]] <- mapa
}

# Combinar todos los mapas
mapa_combinado_ICP <- wrap_plots(mapas_ICP, ncol = length(mapas_ICP)) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

# Mostrar y guardar
print(mapa_combinado_ICP)
ggsave("mapa_ICP_GTWPCA.png", mapa_combinado_ICP, width = 12, height = 4, dpi = 300)

library(ggplot2)
library(dplyr)
library(sf)
library(patchwork)
library(viridis)

# ==============================
# KDE por año con etiquetas y estilo limpio
# ==============================

coords_df <- datos_epi %>%
  mutate(ANO = as.numeric(ANO)) %>%
  left_join(centroides, by = c("Departamento_ocurrencia")) %>%
  filter(!is.na(X) & !is.na(Y))

# Asegurar CRS compatible
colombia_ll <- st_transform(colombia, crs = 4326)

anos_unicos <- sort(unique(coords_df$ANO))
mapas_kde <- list()

for (ano_map in anos_unicos) {
  df_ano <- coords_df %>% filter(ANO == ano_map)
  n_casos <- nrow(df_ano)
  
  mapa <- ggplot() +
    # Mapa base
    geom_sf(data = colombia_ll, fill = "grey90", color = "white", size = 0.2) +
    
    # KDE continuo (hotspots)
    stat_density_2d(
      data = df_ano,
      aes(x = X, y = Y, fill = after_stat(level)),
      geom = "polygon",
      alpha = 0.9,
      contour = TRUE
    ) +
    
    # Escala de color continua tipo cividis
    scale_fill_viridis_c(option = "cividis", name = "Densidad") +
    
    # Puntos de ocurrencia
    geom_point(data = df_ano, aes(x = X, y = Y),
               color = "red", size = 1.2, alpha = 0.8) +
    
    # Etiquetas
    labs(
      title = paste0("Chagas ", ano_map),
      subtitle = paste0("Casos: ", n_casos)
    ) +
    
    coord_sf() +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold", hjust = 0, vjust = -1),
      plot.subtitle = element_text(size = 10, hjust = 0, vjust = -2),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  mapas_kde[[as.character(ano_map)]] <- mapa
}

# Mostrar mapas uno al lado del otro
mapa_kde <-wrap_plots(mapas_kde, ncol = length(mapas_kde))

ggsave(
  filename = "KDE_Chagas.png",
  plot = mapa_kde,
  width = 14, height = 6, dpi = 300, bg = "white"
)
