# code/supplement/S4.R
#
# Reproduces Table S4 (kernel sensitivity) for the STGP model, looping over
# Temporal x Spatial kernel combinations. Reuses the exact LOO / distance-
# matrix setup from code/Table2-Table3(STGP).R; only the covariance
# functions (Rx, Rz), their gradients, and the prediction basis functions
# change per combination. Reports RMSE only (LK prediction, matching how
# the main results.csv "STGP (ours)" row is computed), not NLPD.
#
# The "Exponential | RQ" row is the main STGP model already computed in
# Table 2 / Table 3 -- it is NOT refit here, just pulled from
# results/final_results.csv. The other 6 combinations are refit from scratch:
#   Exponential | Exponential
#   Exponential | Matern(1.5)
#   Exponential | Matern(2.5)
#   Matern(1.5) | RQ
#   Matern(2.5) | RQ
#   RQ          | RQ
#
# NOTE: each combination is its own 66-turbine LOO GP fit, same cost as the
# main STGP run (~4-5 hours per combo per the main README) -- expect this
# script to take on the order of a day to complete all 6.
#
# Output: results/supplement/table_s4.csv
#   (Temporal_Kernel, Spatial_Kernel, RMSE_2017, RMSE_2018)

testset  <- c(1:66)
testset1 <- c(1:46, 48:50, 52, 54:60, 62:66)

data_dir <- file.path("data")
data_path <- data_dir
input_folder <- file.path("data", "processed data")
terrain_data_path <- file.path("data", "weightedTerrainData.csv")
results_dir <- file.path("results", "supplement")
final_results_path <- file.path("results", "final_results.csv")

stopifnot(dir.exists(input_folder))
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
table_s4_path <- file.path(results_dir, "table_s4.csv")

# -------------------------
# Inline table_s4.csv writer (no separate helper file)
# -------------------------
update_table_s4 <- function(temporal_kernel, spatial_kernel, rmse_2017, rmse_2018) {
  cols <- c("Temporal_Kernel", "Spatial_Kernel", "RMSE_2017", "RMSE_2018")
  if (file.exists(table_s4_path)) {
    df <- read.csv(table_s4_path, stringsAsFactors = FALSE)
  } else {
    df <- data.frame(matrix(ncol = length(cols), nrow = 0))
    colnames(df) <- cols
  }
  match_idx <- which(df$Temporal_Kernel == temporal_kernel & df$Spatial_Kernel == spatial_kernel)
  new_row <- data.frame(Temporal_Kernel = temporal_kernel, Spatial_Kernel = spatial_kernel,
                         RMSE_2017 = rmse_2017, RMSE_2018 = rmse_2018, stringsAsFactors = FALSE)
  if (length(match_idx) > 0) {
    df[match_idx[1], cols] <- new_row[cols]
  } else {
    df <- rbind(df, new_row[cols])
  }
  write.csv(df, table_s4_path, row.names = FALSE)
  cat("[table_s4] updated", table_s4_path, "->", temporal_kernel, "|", spatial_kernel, "\n")
}

# -------------------------
# Pull the main (Exponential | RQ) row from results/final_results.csv,
# instead of refitting the already-computed main model.
# -------------------------
kernel_label <- function(name) {
  switch(name,
    "RQ" = "Rational Quadratic (alpha=1)",
    "Exponential" = "Exponential",
    "Matern1.5" = "Matern (nu=1.5)",
    "Matern2.5" = "Matern (nu=2.5)"
  )
}

if (file.exists(final_results_path)) {
  fr <- read.csv(final_results_path, stringsAsFactors = FALSE)
  method_col <- if ("Method" %in% names(fr)) "Method" else names(fr)[1]
  table_col  <- if ("Table" %in% names(fr)) "Table" else if ("table_id" %in% names(fr)) "table_id" else NA
  rmse_col   <- if ("RMSE" %in% names(fr)) "RMSE" else NA

  if (!is.na(table_col) && !is.na(rmse_col)) {
    r17 <- fr[fr[[method_col]] == "STGP (ours)" & fr[[table_col]] == "Table 2", rmse_col]
    r18 <- fr[fr[[method_col]] == "STGP (ours)" & fr[[table_col]] == "Table 3", rmse_col]
    r17 <- if (length(r17) > 0) r17[1] else NA_real_
    r18 <- if (length(r18) > 0) r18[1] else NA_real_
  } else {
    cat("[S4] warning: unexpected schema in", final_results_path, "; main row left NA\n")
    r17 <- NA_real_; r18 <- NA_real_
  }
} else {
  cat("[S4] warning:", final_results_path, "not found; main row left NA\n")
  r17 <- NA_real_; r18 <- NA_real_
}

update_table_s4(kernel_label("Exponential"), kernel_label("RQ"), r17, r18)

# -------------------------
# Data loading (verbatim from Table2-Table3(STGP).R)
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

temp_vector  <- read_vec_csv(file.path(input_folder, "temp_vector.csv"), colname = "temp")
speed_vector <- read_vec_csv(file.path(input_folder, "speed_vector.csv"), colname = "speed")

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

power_matrix <- read_mat_csv(file.path(input_folder, "power_matrix.csv"))

library(dplyr)

turbine_files <- list.files(data_path, pattern = "Turbine[1-9]{1}_2017.csv|Turbine[1-6][0-9]_2017.csv", full.names = TRUE)
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

# =========================================================
# Kernel function factory
# =========================================================
get_spatial_kernel <- function(name) {
  if (name == "RQ") {
    list(
      Rx = function(Ex1, Ex2, Ex3, theta) 1/(1 + (Ex1/theta[1])^2 + (Ex2/theta[2])^2 + (Ex3/theta[3])^2),
      dRx = function(k, Exk, theta, Rx) 2 * Rx^2 * (Exk^2) / theta[k]^3,
      basis = function(h, theta) 1/(1 + sum((h/theta[1:3])^2))
    )
  } else if (name == "Exponential") {
    list(
      Rx = function(Ex1, Ex2, Ex3, theta) exp(-Ex1/theta[1] - Ex2/theta[2] - Ex3/theta[3]),
      dRx = function(k, Exk, theta, Rx) Rx * (Exk / theta[k]^2),
      basis = function(h, theta) exp(-sum(abs(h)/theta[1:3]))
    )
  } else if (name == "Matern1.5") {
    list(
      Rx = function(Ex1, Ex2, Ex3, theta) {
        ((1 + sqrt(3)*Ex1/theta[1]) * exp(-sqrt(3)*Ex1/theta[1])) *
        ((1 + sqrt(3)*Ex2/theta[2]) * exp(-sqrt(3)*Ex2/theta[2])) *
        ((1 + sqrt(3)*Ex3/theta[3]) * exp(-sqrt(3)*Ex3/theta[3]))
      },
      dRx = function(k, Exk, theta, Rx) {
        a_k <- sqrt(3) * Exk / theta[k]
        Rx * (a_k^2) / (theta[k] * (1 + a_k))
      },
      basis = function(h, theta) {
        a <- sqrt(3) * abs(h) / theta[1:3]
        prod((1 + a) * exp(-a))
      }
    )
  } else if (name == "Matern2.5") {
    list(
      Rx = function(Ex1, Ex2, Ex3, theta) {
        ((1 + sqrt(5)*Ex1/theta[1] + (5/3)*(Ex1^2)/theta[1]^2) * exp(-sqrt(5)*Ex1/theta[1])) *
        ((1 + sqrt(5)*Ex2/theta[2] + (5/3)*(Ex2^2)/theta[2]^2) * exp(-sqrt(5)*Ex2/theta[2])) *
        ((1 + sqrt(5)*Ex3/theta[3] + (5/3)*(Ex3^2)/theta[3]^2) * exp(-sqrt(5)*Ex3/theta[3]))
      },
      dRx = function(k, Exk, theta, Rx) {
        a_k <- sqrt(5) * Exk / theta[k]
        Rx * (a_k^2 * (1 + a_k)) / (3 * theta[k] * (1 + a_k + (a_k^2)/3))
      },
      basis = function(h, theta) {
        a <- sqrt(5) * abs(h) / theta[1:3]
        prod((1 + a + (a^2)/3) * exp(-a))
      }
    )
  } else stop("Unknown spatial kernel: ", name)
}

get_temporal_kernel <- function(name) {
  if (name == "Exponential") {
    list(
      Rz = function(Ez, Ez1, theta) exp(-Ez/theta[4] - Ez1/theta[5]),
      dRz = function(k, Ez_k, theta, Rz) Rz * Ez_k / theta[k]^2,
      basis = function(h, h1, theta) exp(-abs(h)/theta[4] - abs(h1)/theta[5])
    )
  } else if (name == "RQ") {
    list(
      Rz = function(Ez, Ez1, theta) 1/(1 + (Ez/theta[4])^2 + (Ez1/theta[5])^2),
      dRz = function(k, Ez_k, theta, Rz) Rz^2 * (2 * (Ez_k^2) / theta[k]^3),
      basis = function(h, h1, theta) 1/(1 + (h/theta[4])^2 + (h1/theta[5])^2)
    )
  } else if (name == "Matern1.5") {
    list(
      Rz = function(Ez, Ez1, theta) {
        ((1 + sqrt(3)*Ez/theta[4]) * exp(-sqrt(3)*Ez/theta[4])) *
        ((1 + sqrt(3)*Ez1/theta[5]) * exp(-sqrt(3)*Ez1/theta[5]))
      },
      dRz = function(k, Ez_k, theta, Rz) {
        c_k <- sqrt(3) * Ez_k
        a_k <- c_k / theta[k]
        Rz * (c_k^2) / (theta[k]^3 * (1 + a_k))
      },
      basis = function(h, h1, theta) {
        a  <- sqrt(3) * abs(h)  / theta[4]
        a1 <- sqrt(3) * abs(h1) / theta[5]
        (1 + a) * exp(-a) * (1 + a1) * exp(-a1)
      }
    )
  } else if (name == "Matern2.5") {
    list(
      Rz = function(Ez, Ez1, theta) {
        ((1 + sqrt(5)*Ez/theta[4] + (5/3)*(Ez^2)/theta[4]^2) * exp(-sqrt(5)*Ez/theta[4])) *
        ((1 + sqrt(5)*Ez1/theta[5] + (5/3)*(Ez1^2)/theta[5]^2) * exp(-sqrt(5)*Ez1/theta[5]))
      },
      dRz = function(k, Ez_k, theta, Rz) {
        a_k <- sqrt(5) * Ez_k / theta[k]
        Rz * (a_k^2 * (1 + a_k)) / (3 * theta[k] * (1 + a_k + (a_k^2)/3))
      },
      basis = function(h, h1, theta) {
        a  <- sqrt(5) * abs(h)  / theta[4]
        a1 <- sqrt(5) * abs(h1) / theta[5]
        (1 + a + (a^2)/3) * exp(-a) * (1 + a1 + (a1^2)/3) * exp(-a1)
      }
    )
  } else stop("Unknown temporal kernel: ", name)
}

# 6 combos to refit (Exponential | RQ excluded -- pulled from final_results.csv above)
kernel_combos <- list(
  list(temporal = "Exponential", spatial = "Exponential"),
  list(temporal = "Exponential", spatial = "Matern1.5"),
  list(temporal = "Exponential", spatial = "Matern2.5"),
  list(temporal = "Matern1.5",   spatial = "RQ"),
  list(temporal = "Matern2.5",   spatial = "RQ"),
  list(temporal = "RQ",          spatial = "RQ")
)

# =========================================================
# Generic ML() builder for one kernel combo
# =========================================================
make_ML <- function(sp, tp) {
  function(theta, Ex1, Ex2, Ex3, Ez, Ez1, n, m, y, vec.y) {
    eps <- 1e-10

    Rx <- sp$Rx(Ex1, Ex2, Ex3, theta) + eps * diag(n)
    Rxinv <- solve(Rx)
    Rz <- tp$Rz(Ez, Ez1, theta) + eps * diag(m)
    Rzinv <- solve(Rz)

    a_mat <- Rxinv %*% rep(1, n) %*% t(rep(1, m)) %*% Rzinv
    b_mat <- Rxinv %*% as.matrix(y) %*% Rzinv
    a <- c(t(a_mat)); b <- c(t(b_mat))
    mu <- sum(b) / sum(a)
    sigma2 <- mean((vec.y - mu * a) * (b - mu * a))

    logdetRx <- determinant(Rx, logarithm = TRUE)$modulus[1]
    logdetRz <- determinant(Rz, logarithm = TRUE)$modulus[1]
    val <- m * n * log(sigma2) + m * logdetRx + n * logdetRz

    grad <- numeric(5)
    denom_a <- sum(a); denom_b <- sum(b)
    dmu <- function(da, db) (sum(db) * denom_a - denom_b * sum(da)) / denom_a^2

    for (k in 1:3) {
      Exk <- switch(k, Ex1, Ex2, Ex3)
      dRx <- sp$dRx(k, Exk, theta, Rx)
      dKx <- -Rxinv %*% dRx %*% Rxinv
      da_mat <- dKx %*% rep(1, n) %*% t(rep(1, m)) %*% Rzinv
      db_mat <- dKx %*% y %*% Rzinv
      da <- c(t(da_mat)); db <- c(t(db_mat))
      d_mu <- dmu(da, db)
      d_sigma2 <- mean(-d_mu * a * (b - mu * a) - mu * da * (b - mu * a) +
                          (vec.y - mu * a) * (db - d_mu * a - mu * da))
      grad[k] <- (m * n / sigma2) * d_sigma2 + m * sum(Rxinv * dRx)
    }

    for (k in 4:5) {
      Ez_k <- if (k == 4) Ez else Ez1
      dRz <- tp$dRz(k, Ez_k, theta, Rz)
      dKz <- -Rzinv %*% dRz %*% Rzinv
      da_mat <- Rxinv %*% rep(1, n) %*% t(rep(1, m)) %*% dKz
      db_mat <- Rxinv %*% y %*% dKz
      da <- c(t(da_mat)); db <- c(t(db_mat))
      d_mu <- dmu(da, db)
      d_sigma2 <- mean(-d_mu * a * (b - mu * a) - mu * da * (b - mu * a) +
                          (vec.y - mu * a) * (db - d_mu * a - mu * da))
      grad[k] <- (m * n / sigma2) * d_sigma2 + n * sum(Rzinv * dRz)
    }

    attr(val, "gradient") <- grad
    val
  }
}

# =========================================================
# Main loop: turbines outer (distance matrices computed once),
# kernel combos inner (reuse distance matrices, refit per combo)
# =========================================================
combo_key <- function(kc) paste(kc$temporal, kc$spatial, sep = "__")
rmse17_store <- setNames(vector("list", length(kernel_combos)), sapply(kernel_combos, combo_key))
rmse18_store <- setNames(vector("list", length(kernel_combos)), sapply(kernel_combos, combo_key))
for (nm in names(rmse17_store)) { rmse17_store[[nm]] <- c(); rmse18_store[[nm]] <- c() }

for (i in testset) {
  cat("=== S4 Turbine", i, "===\n")
  set.seed(i)

  scale_01_test <- function(x) (x - min(all_speeds)) / (max(all_speeds) - min(all_speeds))
  scale_01_test_temperature <- function(x) (x - min(all_temp)) / (max(all_temp) - min(all_temp))

  test_data <- read.csv(file.path(data_path, paste0("Turbine", i, "_2017.csv")))
  test_data$scaled_wind_speed  <- scale_01_test(test_data$wind_speed)
  test_data$scaled_temperature <- scale_01_test_temperature(test_data$temperature)
  z.ev  <- test_data$scaled_wind_speed
  z.ev1 <- test_data$scaled_temperature
  f.ev  <- test_data$power

  do_2018 <- i %in% testset1
  if (do_2018) {
    test_data1 <- read.csv(file.path(data_path, paste0("Turbine", i, "_2018.csv")))
    test_data1$scaled_wind_speed  <- scale_01_test(test_data1$wind_speed)
    test_data1$scaled_temperature <- scale_01_test_temperature(test_data1$temperature)
    z.ev_18  <- test_data1$scaled_wind_speed
    z.ev1_18 <- test_data1$scaled_temperature
    f.ev_18  <- test_data1$power
  }

  x <- scaled_terrain_data[, 2:4]
  x <- x[-i, , drop = FALSE]

  z  <- scale_01_test(speed_vector)
  z1 <- scale_01_test_temperature(temp_vector)
  ord <- order(z)
  z <- z[ord]; z1 <- z1[ord]

  y <- power_matrix[, ord, drop = FALSE]
  y <- y[-i, , drop = FALSE]

  n <- nrow(x); m <- length(z)
  vec.y <- c(t(y))

  Ex1 <- as.matrix(dist(x[, 1], diag = TRUE, upper = TRUE))
  Ex2 <- as.matrix(dist(x[, 2], diag = TRUE, upper = TRUE))
  Ex3 <- as.matrix(dist(x[, 3], diag = TRUE, upper = TRUE))
  Ez  <- as.matrix(dist(z,  diag = TRUE, upper = TRUE))
  Ez1 <- as.matrix(dist(z1, diag = TRUE, upper = TRUE))

  u <- c(scaled_terrain_data[i, 2], scaled_terrain_data[i, 3], scaled_terrain_data[i, 4])

  for (kc in kernel_combos) {
    key <- combo_key(kc)
    sp <- get_spatial_kernel(kc$spatial)
    tp <- get_temporal_kernel(kc$temporal)
    ML <- make_ML(sp, tp)

    ML_log <- function(phi) {
      theta <- exp(phi)
      out <- ML(theta, Ex1, Ex2, Ex3, Ez, Ez1, n, m, y, vec.y)
      g <- attr(out, "gradient") * theta
      attr(out, "gradient") <- g
      out
    }

    phi0 <- log(c(0.1, 0.1, 0.1, 1, 1))
    fit <- nlminb(
      start = phi0,
      objective = function(p) ML_log(p),
      gradient = function(p) attr(ML_log(p), "gradient"),
      control = list(iter.max = 200, rel.tol = 1e-8)
    )
    theta <- exp(fit$par)
    cat("  [", key, "] turbine", i, "convergence:", fit$convergence, "\n")

    Rx <- sp$Rx(Ex1, Ex2, Ex3, theta)
    Rxinv <- solve(Rx + 1e-8 * diag(n))
    Rz <- tp$Rz(Ez, Ez1, theta)
    Rzinv <- solve(Rz + 1e-8 * diag(m))

    a <- c(t((Rxinv %*% rep(1, n)) %*% (t(rep(1, m)) %*% Rzinv)))
    b <- c(t(Rxinv %*% as.matrix(y) %*% Rzinv))
    A <- matrix(a, nrow = n, ncol = m, byrow = TRUE)
    B <- matrix(b, nrow = n, ncol = m, byrow = TRUE)

    r.x_fun <- function(uu) {
      Auu <- t(t(x) - uu)
      apply(Auu, 1, function(h) sp$basis(h, theta))
    }
    r.z_fun <- function(v, v1) tp$basis(v - z, v1 - z1, theta)

    rx_vec <- r.x_fun(u)

    yhat_lk <- numeric(length(z.ev))
    for (j in seq_along(z.ev)) {
      rz_vec <- r.z_fun(z.ev[j], z.ev1[j])
      denom <- as.numeric(t(rx_vec) %*% A %*% rz_vec)
      num   <- as.numeric(t(rx_vec) %*% B %*% rz_vec)
      yhat_lk[j] <- num / denom
    }
    rmse17 <- sqrt(mean((f.ev - yhat_lk)^2, na.rm = TRUE))
    rmse17_store[[key]] <- c(rmse17_store[[key]], rmse17)

    if (do_2018) {
      yhat_lk_18 <- numeric(length(z.ev_18))
      for (j in seq_along(z.ev_18)) {
        rz_vec_18 <- r.z_fun(z.ev_18[j], z.ev1_18[j])
        denom <- as.numeric(t(rx_vec) %*% A %*% rz_vec_18)
        num   <- as.numeric(t(rx_vec) %*% B %*% rz_vec_18)
        yhat_lk_18[j] <- num / denom
      }
      rmse18 <- sqrt(mean((f.ev_18 - yhat_lk_18)^2, na.rm = TRUE))
      rmse18_store[[key]] <- c(rmse18_store[[key]], rmse18)
    }
  }
}

# =========================================================
# Aggregate + write results/supplement/table_s4.csv
# =========================================================
for (kc in kernel_combos) {
  key <- combo_key(kc)
  rmse_2017 <- mean(rmse17_store[[key]], na.rm = TRUE)
  rmse_2018 <- mean(rmse18_store[[key]], na.rm = TRUE)
  update_table_s4(kernel_label(kc$temporal), kernel_label(kc$spatial), rmse_2017, rmse_2018)
}
