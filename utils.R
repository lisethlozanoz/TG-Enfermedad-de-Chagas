library(GWmodel)
library(sp)
library(spdep)

# ======================================================
# 1. GWPCA
# ======================================================

# Pesos espaciales y covarianza
calcular_pesos_espaciales <- function(u, v, coordenadas, ancho_banda, kernel = "bisquare") {
  coordenada_objetivo <- matrix(c(u, v), nrow = 1)
  distancias <- gw.dist(dp.locat = coordenadas, rp.locat = coordenada_objetivo, longlat = FALSE)
  pesos <- gw.weight(distancias, ancho_banda, kernel = kernel, adaptive = FALSE)
  diag(as.vector(pesos))
}

calcular_covarianza_local <- function(X, W) {
  t(X) %*% W %*% X
}

calcular_todas_covarianzas_locales <- function(X, coordenadas, ancho_banda, kernel = "bisquare") {
  n <- nrow(X)
  sigmas_locales <- vector("list", n)
  for (i in 1:n) {
    u <- coordenadas[i, 1]
    v <- coordenadas[i, 2]
    W <- calcular_pesos_espaciales(u, v, coordenadas, ancho_banda, kernel)
    sigmas_locales[[i]] <- calcular_covarianza_local(X, W)
  }
  sigmas_locales
}

# Autovalores, varianza y puntuaciones
calcular_autovalores_locales <- function(sigmas_locales) {
  n <- length(sigmas_locales)
  autovalores_locales <- vector("list", n)
  for (i in 1:n) {
    autoval <- tryCatch({
      eigen(sigmas_locales[[i]])
    }, error = function(e) {
      warning(paste("Error en la descomposición en i =", i))
      NULL
    })
    autovalores_locales[[i]] <- if (!is.null(autoval)) {
      list(valores = autoval$values, vectores = autoval$vectors)
    } else {
      list(valores = NA, vectores = NA)
    }
  }
  autovalores_locales
}

calcular_varianza_explicada_local <- function(sigmas_locales) {
  n <- length(sigmas_locales)
  varianza_explicada <- vector("list", n)
  for (i in 1:n) {
    if (!is.null(sigmas_locales[[i]])) {
      traza <- sum(diag(sigmas_locales[[i]]))
      autovalores <- eigen(sigmas_locales[[i]])$values
      varianza_explicada[[i]] <- autovalores / traza
    } else {
      varianza_explicada[[i]] <- rep(NA, ncol(sigmas_locales[[1]]))
    }
  }
  varianza_explicada
}

calcular_puntuaciones_locales <- function(X, autovectores_locales) {
  n <- length(autovectores_locales)
  puntuaciones_locales <- vector("list", n)
  for (i in 1:n) {
    L <- autovectores_locales[[i]]$vectores
    puntuaciones_locales[[i]] <- if (is.matrix(L)) X %*% L else NA
  }
  puntuaciones_locales
}

# Función principal GWPCA
gwpca <- function(X, coordenadas, ancho_banda, kernel = "bisquare") {
  sigmas_locales <- calcular_todas_covarianzas_locales(X, coordenadas, ancho_banda, kernel)
  autovalores_locales <- calcular_autovalores_locales(sigmas_locales)
  varianza_explicada <- calcular_varianza_explicada_local(sigmas_locales)
  puntuaciones_locales <- calcular_puntuaciones_locales(X, autovalores_locales)
  
  list(
    sigmas = sigmas_locales,
    autovalores = autovalores_locales,
    varianza_explicada = varianza_explicada,
    puntuaciones = puntuaciones_locales
  )
}

# Validación cruzada
calc_gof <- function(local, global) {
  ss_res <- sum((local - global)^2)
  ss_tot <- sum(global^2)
  ss_res / ss_tot
}

gw_weights_fixed <- function(dists_row, bw, kernel) {
  GWmodel::gw.weight(dists_row, bw = bw, kernel = kernel, adaptive = FALSE)
}

tune_LOO <- function(X, coords, dist_mat, kernels, bw_candidates) {
  n <- nrow(coords)
  out <- data.frame(kernel = character(), bw = numeric(), GOF = numeric(), stringsAsFactors = FALSE)
  for (ker in kernels) {
    for (bw in bw_candidates) {
      eigenvals_all <- vector("list", n)
      for (i in seq_len(n)) {
        wvec <- gw_weights_fixed(dist_mat[i, ], bw, ker)
        wvec[i] <- 0
        sW <- sum(wvec)
        if (sW <= 0) {
          eigenvals_all[[i]] <- rep(NA, ncol(X))
        } else {
          W <- diag(wvec)
          cov_local <- t(X) %*% W %*% X / sW
          eigenvals_all[[i]] <- eigen(cov_local, symmetric = TRUE, only.values = TRUE)$values
        }
      }
      mean_eig <- Reduce("+", lapply(eigenvals_all, function(z) {
        if (all(is.na(z))) rep(0, length(eigenvals_all[[which(!sapply(eigenvals_all, function(x) all(is.na(x))))[1]]])) else z
      })) / sum(!sapply(eigenvals_all, function(x) all(is.na(x))))
      gof_val <- calc_gof(mean_eig, eigen_global)
      out <- rbind(out, data.frame(kernel = ker, bw = bw, GOF = gof_val, stringsAsFactors = FALSE))
    }
  }
  out[order(out$GOF), ]
}

tune_LOO_robusta <- function(X_scaled, coords, dist_mat, kernels = c("bisquare", "tricube", "gaussian", "exponential", "boxcar"),
                             bw_candidates = NULL, eigen_global) {
  n <- nrow(coords)
  if (is.null(bw_candidates)) {
    min_dist <- min(dist_mat[dist_mat > 0])
    bw_candidates <- seq(from = min_dist * 1.1, to = max(dist_mat), length.out = 10)
  }
  resultados <- data.frame(kernel = character(), bw = numeric(), GOF = numeric(), stringsAsFactors = FALSE)
  for (ker in kernels) {
    for (bw in bw_candidates) {
      eigenvals_all <- vector("list", n)
      for (i in seq_len(n)) {
        wvec <- gw_weights_fixed(dist_mat[i, ], bw, ker)
        wvec[i] <- 0
        if (sum(wvec) <= 0) {
          eigenvals_all[[i]] <- rep(NA, ncol(X_scaled))
        } else {
          W <- diag(wvec)
          cov_local <- t(X_scaled) %*% W %*% X_scaled / sum(wvec)
          eigenvals_all[[i]] <- eigen(cov_local, symmetric = TRUE, only.values = TRUE)$values
        }
      }
      mean_eig <- Reduce("+", lapply(eigenvals_all, function(z) if (all(is.na(z))) rep(0, ncol(X_scaled)) else z)) /
        sum(!sapply(eigenvals_all, function(x) all(is.na(x))))
      gof_val <- calc_gof(mean_eig, eigen_global)
      resultados <- rbind(resultados, data.frame(kernel = ker, bw = bw, GOF = gof_val))
    }
  }
  resultados[order(resultados$GOF), ]
}

# Resumen GWPCA
crear_resumen_gwpca <- function(resultados_gwpca, n_comp = 3) {
  n_vars <- ncol(resultados_gwpca$sigmas[[1]])
  var_locales <- sapply(resultados_gwpca$varianza_explicada, function(x) x[1:n_comp])
  var_locales <- t(var_locales)
  resumen <- data.frame(
    Componente = paste0("CP", 1:n_comp),
    Desviación.Estándar = sapply(1:n_comp, function(i) {
      mean(sapply(resultados_gwpca$autovalores, function(x) sqrt(x$valores[i])), na.rm = TRUE)
    }),
    Varianza.Explicada = colMeans(var_locales, na.rm = TRUE),
    Varianza.Acumulada = cumsum(colMeans(var_locales, na.rm = TRUE))
  )
  var_acum_local <- t(apply(var_locales, 1, cumsum))
  var_acum_min <- apply(var_acum_local, 2, min)
  var_acum_max <- apply(var_acum_local, 2, max)
  resumen$Varianza.Acumulada.Min <- var_acum_min
  resumen$Varianza.Acumulada.Max <- var_acum_max
  return(resumen)
}

crear_resumen_gwpca_generico <- function(resultados_gwpca, ids, n_comp = 2) {
  var_local <- resultados_gwpca$varianza_explicada
  var_local_mat <- do.call(rbind, lapply(var_local, function(x) x[1:n_comp]))
  var_prom <- colMeans(var_local_mat, na.rm = TRUE)
  var_acum <- cumsum(var_prom)
  resumen_global <- data.frame(
    Componente = paste0("CP", 1:n_comp),
    Varianza.Explicada.Promedio = var_prom,
    Varianza.Acumulada = var_acum
  )
  resumen_local <- as.data.frame(var_local_mat)
  colnames(resumen_local) <- paste0("CP", 1:n_comp)
  resumen_local$id <- ids
  list(
    resumen_global = resumen_global,
    resumen_local = resumen_local,
    varianza_local = var_local_mat
  )
}

# ======================================================
# 2. GTWPCA
# ======================================================

calculate_spatiotemporal_distances <- function(coords, times, alpha, beta) {
  n <- nrow(coords)
  dST <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) {
    spatial_sq <- sum((coords[i, ] - coords[j, ])^2)
    temporal_sq <- (times[i] - times[j])^2
    dST[i, j] <- alpha * spatial_sq + beta * temporal_sq
  }
  sqrt(dST)
}

bisquare_kernel_st <- function(dists, bandwidth) {
  ifelse(dists < bandwidth, (1 - (dists / bandwidth)^2)^2, 0)
}

calculate_spatiotemporal_weights <- function(i, dST, bandwidth) {
  diag(bisquare_kernel_st(dST[i, ], bandwidth))
}

calculate_local_covariance <- function(X, W) {
  t(X) %*% W %*% X
}

calculate_local_eigen <- function(local_sigmas) {
  lapply(local_sigmas, function(Sigma) {
    eig <- eigen(Sigma)
    list(values = eig$values, vectors = eig$vectors)
  })
}

calculate_local_scores <- function(X, local_eigenvectors) {
  lapply(local_eigenvectors, function(eig) X %*% eig$vectors)
}

calculate_local_pvt <- function(local_eigen) {
  lapply(local_eigen, function(eig) eig$values / sum(eig$values))
}

gtwpca <- function(X, coords, times, alpha, beta, bandwidth) {
  dST <- calculate_spatiotemporal_distances(coords, times, alpha, beta)
  n <- nrow(X)
  local_sigmas <- vector("list", n)
  for (i in 1:n) {
    W <- calculate_spatiotemporal_weights(i, dST, bandwidth)
    local_sigmas[[i]] <- calculate_local_covariance(X, W)
  }
  local_eigens <- calculate_local_eigen(local_sigmas)
  local_scores <- calculate_local_scores(X, local_eigens)
  local_pvt <- calculate_local_pvt(local_eigens)
  list(
    sigmas = local_sigmas,
    autovalores = local_eigens,
    puntuaciones = local_scores,
    varianza_explicada = local_pvt
  )
}

cv_gtwpca <- function(X, coords, times, alpha_vals, beta_vals, bandwidth_vals) {
  best_error <- Inf
  best_params <- list()
  n <- nrow(X)
  p <- ncol(X)
  for (alpha in alpha_vals) {
    for (beta in beta_vals) {
      for (bw in bandwidth_vals) {
        resultados <- gtwpca(X, coords, times, alpha, beta, bw)
        X_hat <- matrix(NA, n, p)
        for (i in 1:n) {
          loadings_i <- resultados$eigens[[i]]$vectors
          scores_i <- resultados$scores[[i]][i, , drop = FALSE]
          X_hat[i, ] <- scores_i %*% t(loadings_i)
        }
        error <- mean((X - X_hat)^2, na.rm = TRUE)
        if (error < best_error) {
          best_error <- error
          best_params <- list(alpha = alpha, beta = beta, bandwidth = bw)
        }
      }
    }
  }
  list(best_params = best_params, best_error = best_error)
}

# ======================================================
# 3. FUNCIONES GENERALES
# ======================================================

calcular_moran <- function(X, lw, alpha = 0.05) {
  resultados_moran <- data.frame(
    Variable = character(),
    Moran_I = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  for (var_name in colnames(X)) {
    x <- X[[var_name]]
    if (!is.numeric(x) || length(unique(na.omit(x))) < 3) next
    test <- tryCatch(
      { moran.test(x, lw, zero.policy = TRUE, na.action = na.exclude) },
      error = function(e) NULL
    )
    if (!is.null(test)) {
      resultados_moran <- rbind(resultados_moran, data.frame(
        Variable = var_name,
        Moran_I = test$estimate[["Moran I statistic"]],
        p_value = test$p.value
      ))
    }
  }
  print(resultados_moran)
  vars_significativas <- resultados_moran %>%
    dplyr::filter(p_value < alpha) %>%
    dplyr::pull(Variable)
  X[, vars_significativas, drop = FALSE]
}

get_moda <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

