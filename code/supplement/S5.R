# code/supplement/S5.R
#
# Reproduces Table S5 (observation noise vs. observation+process noise) for
# the STGP model.
#
#   - "STGP (observation noise included)" is the main paper's STGP
#     implementation -- NOT rerun here, pulled directly from
#     results/final_results.csv (method "STGP (ours)", Table 2 / Table 3).
#   - "STGP (both observation and process noise)" is implemented here,
#     adapted from Tablle_2_3_delta.R: the first-stage process variance
#     sigma_star_2 = s_2 - nug (thinned-twinGP variance minus the nugget) is
#     folded into the second-stage covariance via a scaling parameter
#     lambda/delta, applied through a first-order Neumann approximation of
#     K^{-1} (see fit_stage2_given_delta / apply_Sinv_1st below). Length-
#     scales theta are fit once per turbine under the noise-free model and
#     held fixed while delta is tuned by minimizing RMSE on a validation set
#     of 1,000 samples from 4 held-out turbines, over the grid
#     10^seq(-5, 2, by = 1) (as in the original delta script).
#
# Output: results/supplement/table_s5.csv
#   (Method, RMSE_2017, NLPD_2017, RMSE_2018, NLPD_2018)

testset  <- c(1:66)
testset1 <- c(1:46, 48:50, 52, 54:60, 62:66)

data_dir  <- file.path("data")
input_folder <- file.path("data", "processed data")
terrain_data_path <- file.path("data", "weightedTerrainData.csv")
results_dir <- file.path("results", "supplement")
final_results_path <- file.path("results", "final_results.csv")

stopifnot(dir.exists(input_folder))
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
table_s5_path <- file.path(results_dir, "table_s5.csv")

# -------------------------
# Inline table_s5.csv writer (no separate helper file)
# -------------------------
update_table_s5 <- function(method, rmse_2017 = NA_real_, nlpd_2017 = NA_real_,
                             rmse_2018 = NA_real_, nlpd_2018 = NA_real_) {
  cols <- c("Method", "RMSE_2017", "NLPD_2017", "RMSE_2018", "NLPD_2018")
  if (file.exists(table_s5_path)) {
    df <- read.csv(table_s5_path, stringsAsFactors = FALSE)
  } else {
    df <- data.frame(matrix(ncol = length(cols), nrow = 0)); colnames(df) <- cols
  }
  match_idx <- which(df$Method == method)
  new_row <- data.frame(Method = method, RMSE_2017 = rmse_2017, NLPD_2017 = nlpd_2017,
                         RMSE_2018 = rmse_2018, NLPD_2018 = nlpd_2018, stringsAsFactors = FALSE)
  if (length(match_idx) > 0) df[match_idx[1], cols] <- new_row[cols] else df <- rbind(df, new_row[cols])
  write.csv(df, table_s5_path, row.names = FALSE)
  cat("[table_s5] updated", table_s5_path, "->", method, "\n")
}

# -------------------------
# Pull "STGP (observation noise included)" row from results/final_results.csv
# -------------------------
if (file.exists(final_results_path)) {
  fr <- read.csv(final_results_path, stringsAsFactors = FALSE)
  method_col <- if ("Method" %in% names(fr)) "Method" else names(fr)[1]
  table_col  <- if ("Table" %in% names(fr)) "Table" else if ("table_id" %in% names(fr)) "table_id" else NA
  rmse_col   <- if ("RMSE" %in% names(fr)) "RMSE" else NA
  nlpd_col   <- if ("NLPD" %in% names(fr)) "NLPD" else NA

  get_val <- function(table_id, col) {
    if (is.na(table_col) || is.na(col)) return(NA_real_)
    v <- fr[fr[[method_col]] == "STGP (ours)" & fr[[table_col]] == table_id, col]
    if (length(v) > 0) v[1] else NA_real_
  }
  obs_rmse_2017 <- get_val("Table 2", rmse_col); obs_nlpd_2017 <- get_val("Table 2", nlpd_col)
  obs_rmse_2018 <- get_val("Table 3", rmse_col); obs_nlpd_2018 <- get_val("Table 3", nlpd_col)
} else {
  cat("[S5] warning:", final_results_path, "not found; observation-noise row left NA\n")
  obs_rmse_2017 <- NA_real_; obs_nlpd_2017 <- NA_real_
  obs_rmse_2018 <- NA_real_; obs_nlpd_2018 <- NA_real_
}

update_table_s5("STGP (observation noise included)",
                 rmse_2017 = obs_rmse_2017, nlpd_2017 = obs_nlpd_2017,
                 rmse_2018 = obs_rmse_2018, nlpd_2018 = obs_nlpd_2018)

# -------------------------
# Data loading (same as Tablle_2_3_delta.R)
# -------------------------
read_vec_csv <- function(path, colname = NULL) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(colname)) {
    if (!colname %in% names(df)) stop("Column '", colname, "' not found in: ", path)
    v <- df[[colname]]
  } else {
    if (ncol(df) != 1) stop("Expected 1 column in: ", path, " but found ", ncol(df))
    v <- df[[1]]
  }
  if (is.character(v)) suppressWarnings(vn <- as.numeric(v)) else vn <- v
  if (is.numeric(vn) && !all(is.na(vn))) return(vn)
  return(v)
}

read_mat_csv <- function(path) {
  df_try <- try(read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
  if (!inherits(df_try, "try-error")) {
    df <- df_try
    rn <- rownames(df)
    if (any(is.na(rn)) || any(rn == "") || any(duplicated(rn))) {
      df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    }
  } else {
    df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  df_num <- df
  for (j in seq_along(df_num)) suppressWarnings(df_num[[j]] <- as.numeric(df_num[[j]]))
  na_ratio <- function(x) mean(is.na(x))
  before <- mean(vapply(df, na_ratio, numeric(1)))
  after  <- mean(vapply(df_num, na_ratio, numeric(1)))
  if (after <= before + 1e-12) df <- df_num
  as.matrix(df)
}

temp_vector   <- read_vec_csv(file.path(input_folder, "temp_vector.csv"), colname = "temp")
speed_vector  <- read_vec_csv(file.path(input_folder, "speed_vector.csv"), colname = "speed")
power_matrix  <- read_mat_csv(file.path(input_folder, "power_matrix.csv"))
sd_matrix     <- read_mat_csv(file.path(input_folder, "sd_matrix.csv"))
sigma2_matrix <- read_mat_csv(file.path(input_folder, "sigma2_matrix.csv"))

library(dplyr)

turbine_files <- list.files(data_dir, pattern = "Turbine[1-9]{1}_2017.csv|Turbine[1-6][0-9]_2017.csv", full.names = TRUE)
set.seed(15)
turbine_data_list <- list()
for (i in seq_along(turbine_files)) turbine_data_list[[i]] <- read.csv(turbine_files[i])

all_speeds <- c()
for (i in seq_along(turbine_data_list)) all_speeds <- c(all_speeds, turbine_data_list[[i]]$wind_speed)
all_temp <- c()
for (i in seq_along(turbine_data_list)) all_temp <- c(all_temp, turbine_data_list[[i]]$temperature)

terrain_data <- read.csv(terrain_data_path)
scale_01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
scaled_terrain_data <- as.data.frame(lapply(terrain_data, scale_01))

scale_01_speed <- function(x) (x - min(all_speeds)) / (max(all_speeds) - min(all_speeds))
scale_01_temp  <- function(x) (x - min(all_temp))   / (max(all_temp)   - min(all_temp))

eps_jit <- 1e-8
eps_inv <- 1e-8

z_raw  <- scale_01_speed(speed_vector)
z1_raw <- scale_01_temp(temp_vector)
ord_global <- order(z_raw)
z_sorted  <- z_raw[ord_global]
z1_sorted <- z1_raw[ord_global]
Ez_global  <- as.matrix(dist(z_sorted,  diag = TRUE, upper = TRUE))
Ez1_global <- as.matrix(dist(z1_sorted, diag = TRUE, upper = TRUE))
m_global <- length(z_sorted)

cap_outliers_row <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  upper <- q[2] + 1.5 * iqr
  x[x > upper] <- upper
  x
}

# =========================================================
# Stage-2 fit given a fixed delta (observation + process noise via
# first-order Neumann approximation of K^{-1})
# =========================================================
fit_stage2_given_delta <- function(vec.y, n, m, Rx, Rz, sigma_star_2, delta,
                                    eps_inv = 0, cap_fun = NULL) {
  N <- n * m
  oneN <- rep(1, N)

  sigma_star_2 <- pmax(sigma_star_2, 0)
  if (!is.null(cap_fun)) sigma_star_2 <- t(apply(sigma_star_2, 1, cap_fun))

  Rxinv <- solve(Rx + eps_inv * diag(n))
  Rzinv <- solve(Rz + eps_inv * diag(m))

  apply_Rinv <- function(v) {
    V <- matrix(v, nrow = n, ncol = m, byrow = TRUE)
    out <- Rxinv %*% V %*% Rzinv
    c(t(out))
  }

  a0  <- apply_Rinv(oneN)
  b0  <- apply_Rinv(vec.y)
  mu0 <- sum(b0) / sum(a0)

  z0     <- vec.y - mu0 * oneN
  u0     <- apply_Rinv(z0)
  tau2_0 <- max((1 / N) * sum(z0 * u0), 1e-12)

  sigma_star_vec <- c(t(sigma_star_2))
  lambda_eff     <- delta * (sigma_star_vec / tau2_0)

  apply_Sinv_1st <- function(b) {
    u  <- apply_Rinv(b)
    u2 <- apply_Rinv(lambda_eff * u)
    u - u2
  }

  a  <- apply_Sinv_1st(oneN)
  b  <- apply_Sinv_1st(vec.y)
  mu <- sum(b) / sum(a)

  z_res    <- vec.y - mu * oneN
  Sz       <- apply_Sinv_1st(z_res)
  tau2_hat <- max((1 / N) * sum(z_res * Sz), 1e-12)

  A <- matrix(a, nrow = n, ncol = m, byrow = TRUE)
  B <- matrix(b, nrow = n, ncol = m, byrow = TRUE)

  list(A = A, B = B, mu = mu, sigma2 = tau2_hat, Rxinv = Rxinv, Rzinv = Rzinv,
       apply_Sinv_1st = apply_Sinv_1st, vec_y = vec.y, lambda_eff = lambda_eff)
}

# Lightweight LK mean predictor for delta-tuning validation RMSE
predict_LK_mean_light <- function(u, x, z, z1, theta, A, B, z.ev, z.ev1) {
  Auu <- t(t(x) - u)
  rx <- 1 / (1 + (Auu[,1]/theta[1])^2 + (Auu[,2]/theta[2])^2 + (Auu[,3]/theta[3])^2)
  Ax <- as.numeric(crossprod(rx, A))
  Bx <- as.numeric(crossprod(rx, B))
  Dz  <- abs(outer(z,  z.ev,  "-"))
  Dz1 <- abs(outer(z1, z.ev1, "-"))
  RZ  <- exp(-Dz/theta[4] - Dz1/theta[5])
  denom <- as.numeric(crossprod(Ax, RZ))
  num   <- as.numeric(crossprod(Bx, RZ))
  num / denom
}

# Delta grid search: minimizes RMSE on a validation set of val_P samples
# from 4 held-out turbines
tune_delta_rmse_light <- function(delta_grid, vec.y, n, m, Rx, Rz, sigma_star_2,
                                   x_all, x_train, z, z1, theta,
                                   val_turbs, val_P = 1000, data_dir, eps_inv = 0,
                                   cap_fun = NULL, seed = 123, stability_check = TRUE) {
  set.seed(seed)

  val_cache <- vector("list", length(val_turbs))
  names(val_cache) <- as.character(val_turbs)
  for (k in seq_along(val_turbs)) {
    tv <- val_turbs[k]
    tv_data <- read.csv(file.path(data_dir, paste0("Turbine", tv, "_2017.csv")))
    tv_data$scaled_wind_speed  <- scale_01_speed(tv_data$wind_speed)
    tv_data$scaled_temperature <- scale_01_temp(tv_data$temperature)
    idx <- sample.int(nrow(tv_data), min(val_P, nrow(tv_data)))
    val_cache[[k]] <- list(u = as.numeric(x_all[tv, ]),
                            z = tv_data$scaled_wind_speed[idx],
                            z1 = tv_data$scaled_temperature[idx],
                            y = tv_data$power[idx])
  }

  rmse_delta <- rep(NA_real_, length(delta_grid))
  for (dd in seq_along(delta_grid)) {
    delta <- delta_grid[dd]
    fit2 <- fit_stage2_given_delta(vec.y, n, m, Rx, Rz, sigma_star_2, delta,
                                    eps_inv = eps_inv, cap_fun = cap_fun)
    se_all <- numeric(0); ok <- TRUE
    for (k in seq_along(val_cache)) {
      item <- val_cache[[k]]
      yhat <- predict_LK_mean_light(item$u, x_train, z, z1, theta, fit2$A, fit2$B, item$z, item$z1)
      if (stability_check) {
        if (any(!is.finite(yhat))) { ok <- FALSE; break }
        if (max(abs(yhat), na.rm = TRUE) > 1e6) { ok <- FALSE; break }
      }
      se_all <- c(se_all, (item$y - yhat)^2)
    }
    rmse_delta[dd] <- if (ok) sqrt(mean(se_all, na.rm = TRUE)) else Inf
  }
  best_id <- which.min(rmse_delta)
  list(delta_best = delta_grid[best_id], rmse_delta = rmse_delta)
}

# =========================================================
# Main LOO loop: fit theta (noise-free), tune delta, predict, score
# =========================================================
rmse_list <- c(); rmse_list_18 <- c()
nlpd_list <- c(); nlpd_list_18 <- c()

for (i in testset) {
  cat("=== S5 Turbine", i, "===\n")
  set.seed(i)

  test_data <- read.csv(file.path(data_dir, paste0("Turbine", i, "_2017.csv")))
  test_data$scaled_wind_speed  <- scale_01_speed(test_data$wind_speed)
  test_data$scaled_temperature <- scale_01_temp(test_data$temperature)
  z.ev  <- test_data$scaled_wind_speed
  z.ev1 <- test_data$scaled_temperature
  f.ev  <- test_data$power

  do_2018 <- i %in% testset1
  if (do_2018) {
    test_data1 <- read.csv(file.path(data_dir, paste0("Turbine", i, "_2018.csv")))
    test_data1$scaled_wind_speed  <- scale_01_speed(test_data1$wind_speed)
    test_data1$scaled_temperature <- scale_01_temp(test_data1$temperature)
    z.ev_18  <- test_data1$scaled_wind_speed
    z.ev1_18 <- test_data1$scaled_temperature
    f.ev_18  <- test_data1$power
  }

  x_all <- as.matrix(scaled_terrain_data[, 2:4])
  u <- as.numeric(x_all[i, ])
  x <- x_all[-i, , drop = FALSE]
  n <- nrow(x)

  z <- z_sorted; z1 <- z1_sorted
  Ez <- Ez_global; Ez1 <- Ez1_global
  m <- m_global; N <- n * m

  y <- power_matrix[, ord_global, drop = FALSE]
  y <- y[-i, , drop = FALSE]
  vec.y <- c(t(y))

  nug <- sigma2_matrix[-i, ord_global]
  s_2 <- (sd_matrix[-i, ord_global])^2

  Ex1 <- as.matrix(dist(x[, 1], diag = TRUE, upper = TRUE))
  Ex2 <- as.matrix(dist(x[, 2], diag = TRUE, upper = TRUE))
  Ex3 <- as.matrix(dist(x[, 3], diag = TRUE, upper = TRUE))

  basis.x <- function(h, th123) 1 / (1 + sum((h / th123)^2))
  basis.z <- function(h, h1, th45) exp(-abs(h) / th45[1] - abs(h1) / th45[2])
  r.x <- function(u0, th) {
    Auu <- t(t(x) - u0)
    apply(Auu, 1, function(row) basis.x(row, th[1:3]))
  }
  r.z <- function(v, v1, th) basis.z(v - z, v1 - z1, th[4:5])

  # ---- fit length-scales theta under the noise-free (Exponential/RQ) model ----
  ML_theta <- function(theta) {
    if (any(!is.finite(theta)) || any(theta <= 0)) return(Inf)
    Rx <- 1/(1 + (Ex1/theta[1])^2 + (Ex2/theta[2])^2 + (Ex3/theta[3])^2) + eps_jit*diag(n)
    Rz <- exp(-Ez/theta[4] - Ez1/theta[5]) + eps_jit*diag(m)
    Rx_ch <- tryCatch(chol(Rx), error = function(e) NULL)
    Rz_ch <- tryCatch(chol(Rz), error = function(e) NULL)
    if (is.null(Rx_ch) || is.null(Rz_ch)) return(Inf)
    logdetRx <- 2 * sum(log(diag(Rx_ch)))
    logdetRz <- 2 * sum(log(diag(Rz_ch)))
    logdetR  <- m * logdetRx + n * logdetRz
    Rxinv <- solve(Rx); Rzinv <- solve(Rz)
    apply_Rinv <- function(v) {
      V <- matrix(v, nrow = n, ncol = m, byrow = TRUE)
      c(t(Rxinv %*% V %*% Rzinv))
    }
    a0 <- apply_Rinv(rep(1, N)); b0 <- apply_Rinv(vec.y)
    mu0 <- sum(b0) / sum(a0)
    z0 <- vec.y - mu0 * rep(1, N)
    u0 <- apply_Rinv(z0)
    s2 <- sum(z0 * u0) / N
    if (!is.finite(s2) || s2 <= 0) return(Inf)
    N * log(s2) + logdetR
  }

  fit_theta <- nlminb(start = log(c(0.1, 0.1, 0.1, 1, 1)),
                       objective = function(phi) ML_theta(exp(phi)),
                       control = list(iter.max = 250, rel.tol = 1e-6))
  theta <- exp(fit_theta$par)

  # ---- process variance + delta tuning ----
  sigma_star_2 <- pmax(s_2 - nug, 0)
  sigma_star_2 <- t(apply(sigma_star_2, 1, cap_outliers_row))

  set.seed(123)
  cand_turbs <- setdiff(seq_len(nrow(x_all)), i)
  val_turbs  <- sample(cand_turbs, 4)
  delta_grid <- 10^seq(-5, 2, by = 1)

  Rx2 <- 1/(1 + (Ex1/theta[1])^2 + (Ex2/theta[2])^2 + (Ex3/theta[3])^2)
  Rz2 <- exp(-Ez/theta[4] - Ez1/theta[5])

  delta_fit <- tune_delta_rmse_light(
    delta_grid = delta_grid, vec.y = vec.y, n = n, m = m, Rx = Rx2, Rz = Rz2,
    sigma_star_2 = sigma_star_2, x_all = x_all, x_train = x, z = z, z1 = z1,
    theta = theta, val_turbs = val_turbs, val_P = 1000, data_dir = data_dir,
    eps_inv = eps_inv, cap_fun = cap_outliers_row, seed = 123, stability_check = TRUE
  )
  delta_best <- delta_fit$delta_best
  cat("  Turbine", i, "selected delta:", format(delta_best, scientific = TRUE), "\n")

  fit2 <- fit_stage2_given_delta(vec.y, n, m, Rx2, Rz2, sigma_star_2, delta_best,
                                  eps_inv = eps_inv, cap_fun = cap_outliers_row)
  A <- fit2$A; B <- fit2$B; mu <- fit2$mu; sigma2 <- fit2$sigma2
  Rxinv <- fit2$Rxinv; Rzinv <- fit2$Rzinv
  tau2 <- 0  # predictive-variance correction beyond mspe_latent, as in the original delta script

  rx_vec <- r.x(u, theta)
  nearest_col <- vapply(z.ev, function(zz) which.min(abs(z - zz)), integer(1))
  nug_scalar  <- vapply(nearest_col, function(k) median(nug[, k], na.rm = TRUE), numeric(1))

  yhat <- numeric(length(z.ev)); yhat_sd <- numeric(length(z.ev))
  for (jj in seq_along(z.ev)) {
    rz_vec <- r.z(z.ev[jj], z.ev1[jj], theta)
    denom <- as.numeric(t(rx_vec) %*% A %*% rz_vec)
    num   <- as.numeric(t(rx_vec) %*% B %*% rz_vec)
    yhat[jj] <- num / denom
    alpha <- as.numeric(t(rx_vec) %*% Rxinv %*% rx_vec)
    beta  <- as.numeric(t(rz_vec) %*% Rzinv %*% rz_vec)
    ab <- alpha * beta
    mspe_latent <- sigma2 * (1 - ab + ab * (1 - denom)^2 / denom^2)
    yhat_sd[jj] <- sqrt(pmax(mspe_latent + tau2 + nug_scalar[jj], 1e-12))
  }

  nlpds <- 0.5*log(2*pi*yhat_sd^2) + 0.5*((f.ev - yhat)^2)/(yhat_sd^2)
  rmse_list <- c(rmse_list, sqrt(mean((f.ev - yhat)^2, na.rm = TRUE)))
  nlpd_list <- c(nlpd_list, mean(nlpds, na.rm = TRUE))

  if (do_2018) {
    nearest_col_18 <- vapply(z.ev_18, function(zz) which.min(abs(z - zz)), integer(1))
    nug_scalar_18  <- vapply(nearest_col_18, function(k) median(nug[, k], na.rm = TRUE), numeric(1))

    yhat_18 <- numeric(length(z.ev_18)); yhat_sd_18 <- numeric(length(z.ev_18))
    for (jj in seq_along(z.ev_18)) {
      rz_vec_18 <- r.z(z.ev_18[jj], z.ev1_18[jj], theta)
      denom <- as.numeric(t(rx_vec) %*% A %*% rz_vec_18)
      num   <- as.numeric(t(rx_vec) %*% B %*% rz_vec_18)
      yhat_18[jj] <- num / denom
      alpha <- as.numeric(t(rx_vec) %*% Rxinv %*% rx_vec)
      beta  <- as.numeric(t(rz_vec_18) %*% Rzinv %*% rz_vec_18)
      ab <- alpha * beta
      mspe_latent_18 <- sigma2 * (1 - ab + ab * (1 - denom)^2 / denom^2)
      yhat_sd_18[jj] <- sqrt(pmax(mspe_latent_18 + tau2 + nug_scalar_18[jj], 1e-12))
    }

    nlpds18 <- 0.5*log(2*pi*yhat_sd_18^2) + 0.5*((f.ev_18 - yhat_18)^2)/(yhat_sd_18^2)
    rmse_list_18 <- c(rmse_list_18, sqrt(mean((f.ev_18 - yhat_18)^2, na.rm = TRUE)))
    nlpd_list_18 <- c(nlpd_list_18, mean(nlpds18, na.rm = TRUE))
  }
}

# =========================================================
# Aggregate + write results/supplement/table_s5.csv
# =========================================================
update_table_s5(
  "STGP (both observation and process noise)",
  rmse_2017 = mean(rmse_list, na.rm = TRUE),
  nlpd_2017 = mean(nlpd_list, na.rm = TRUE),
  rmse_2018 = mean(rmse_list_18, na.rm = TRUE),
  nlpd_2018 = mean(nlpd_list_18, na.rm = TRUE)
)
