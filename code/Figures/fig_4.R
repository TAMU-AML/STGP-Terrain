# code/figures/fig_4.R
#
# Reproduces Figure 4: comparison of Limit Kriging (LK) vs Ordinary Kriging
# (OK) predictors.
#
#   Left:  boxplot of RMSE across all 66 turbines (LOTO experiment), LK vs
#          OK. Read directly from results/intermediate/table2-table3.csv
#          (produced by code/Table2-Table3(STGP).R) -- NOT recomputed here.
#   Right: power-curve scatter for Turbine #60 only. The STGP model is
#          refit leave-one-out for turbine 60, then LK and OK predictions
#          are evaluated at that turbine's real 2017 test points and
#          plotted against the actual observed data.
#
# Output:
#   results/figures/fig4_lk_ok_rmse_boxplot.pdf
#   results/figures/fig4_turbine60_lk_ok.pdf

data_dir <- file.path("data")
input_folder <- file.path("data", "processed data")
terrain_data_path <- file.path("data", "weightedTerrainData.csv")
intermediate_path <- file.path("results", "intermediate", "table2-table3.csv")
figures_dir <- file.path("results", "figures")

stopifnot(dir.exists(input_folder))
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)

# =========================================================
# LEFT PANEL: RMSE boxplot, LK vs OK (from existing intermediate results)
# =========================================================
if (!file.exists(intermediate_path)) {
  stop("Missing ", intermediate_path,
       " -- run code/Table2-Table3(STGP).R first to produce it.")
}

res <- read.csv(intermediate_path, stringsAsFactors = FALSE)
stopifnot(all(c("RMSE", "RMSE_OK_2017") %in% names(res)))

pdf(file.path(figures_dir, "fig4_lk_ok_rmse_boxplot.pdf"), width = 4, height = 5)
par(mar = c(3, 4, 1, 1))
boxplot(res$RMSE, res$RMSE_OK_2017,
        names = c("LK", "OK"),
        col = c("forestgreen", "red"),
        ylab = "RMSE across 66 Turbines")
dev.off()

# =========================================================
# RIGHT PANEL: Turbine #60 -- LK and OK predictions vs actual data
# =========================================================
TARGET_TURBINE <- 60

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

terrain_data <- read.csv(terrain_data_path)
scale_01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
scaled_terrain_data <- as.data.frame(lapply(terrain_data, scale_01))

turbine_files <- list.files(data_dir, pattern = "Turbine[1-9]{1}_2017.csv|Turbine[1-6][0-9]_2017.csv", full.names = TRUE)
set.seed(15)
turbine_data_list <- list()
for (k in seq_along(turbine_files)) turbine_data_list[[k]] <- read.csv(turbine_files[k])

all_speeds <- c()
for (k in seq_along(turbine_data_list)) all_speeds <- c(all_speeds, turbine_data_list[[k]]$wind_speed)
all_temp <- c()
for (k in seq_along(turbine_data_list)) all_temp <- c(all_temp, turbine_data_list[[k]]$temperature)

scale_01_test <- function(x) (x - min(all_speeds)) / (max(all_speeds) - min(all_speeds))
scale_01_test_temperature <- function(x) (x - min(all_temp)) / (max(all_temp) - min(all_temp))

i <- TARGET_TURBINE
cat("=== Fig4 refit: Turbine", i, "===\n")

test_data <- read.csv(file.path(data_dir, paste0("Turbine", i, "_2017.csv")))
test_data$scaled_wind_speed  <- scale_01_test(test_data$wind_speed)
test_data$scaled_temperature <- scale_01_test_temperature(test_data$temperature)
z.ev  <- test_data$scaled_wind_speed
z.ev1 <- test_data$scaled_temperature
f.ev  <- test_data$power
wind_speed_ev <- test_data$wind_speed  # unscaled, for plotting x-axis

x <- scaled_terrain_data[, 2:4]
x <- x[-i, , drop = FALSE]
y <- power_matrix

z  <- scale_01_test(speed_vector)
z1 <- scale_01_test_temperature(temp_vector)
ord <- order(z)
z  <- z[ord]
z1 <- z1[ord]
y  <- power_matrix[, ord, drop = FALSE]
y  <- y[-i, , drop = FALSE]

n <- nrow(x)
m <- length(z)
vec.y <- c(t(y))

Ex1 <- as.matrix(dist(x[, 1], diag = TRUE, upper = TRUE))
Ex2 <- as.matrix(dist(x[, 2], diag = TRUE, upper = TRUE))
Ex3 <- as.matrix(dist(x[, 3], diag = TRUE, upper = TRUE))
Ez  <- as.matrix(dist(z,  diag = TRUE, upper = TRUE))
Ez1 <- as.matrix(dist(z1, diag = TRUE, upper = TRUE))

ML <- function(para) {
  theta <- para
  eps <- 1e-10

  Rx <- 1 / (1 + (Ex1/theta[1])^2 + (Ex2/theta[2])^2 + (Ex3/theta[3])^2)
  Rx <- Rx + eps * diag(n)
  Rxinv <- solve(Rx)

  Rz <- exp(-Ez/theta[4] - Ez1/theta[5])
  Rz <- Rz + eps * diag(m)
  Rzinv <- solve(Rz)

  a_mat <- Rxinv %*% rep(1, n) %*% t(rep(1, m)) %*% Rzinv
  b_mat <- Rxinv %*% as.matrix(y) %*% Rzinv
  a <- c(t(a_mat)); b <- c(t(b_mat))
  mu <- sum(b) / sum(a)
  sigma2 <- (1/(m*n)) * sum((vec.y - mu) * (b - mu * a))

  logdetRx <- determinant(Rx, logarithm = TRUE)$modulus[1]
  logdetRz <- determinant(Rz, logarithm = TRUE)$modulus[1]
  val <- m * n * log(sigma2) + m * logdetRx + n * logdetRz

  grad <- numeric(5)
  denom_a <- sum(a); denom_b <- sum(b)
  dmu <- function(da, db) (sum(db) * denom_a - denom_b * sum(da)) / denom_a^2

  for (k in 1:3) {
    Exk <- switch(k, Ex1, Ex2, Ex3)
    dRx  <- 2 * Rx^2 * (Exk^2) / theta[k]^3
    dKx  <- -Rxinv %*% dRx %*% Rxinv
    da_mat <- dKx %*% rep(1, n) %*% t(rep(1, m)) %*% Rzinv
    db_mat <- dKx %*% y %*% Rzinv
    da <- c(t(da_mat)); db <- c(t(db_mat))
    d_mu <- dmu(da, db)
    d_sigma2 <- mean(-d_mu*a*(b-mu*a) - mu*da*(b-mu*a) + (vec.y-mu*a)*(db-d_mu*a-mu*da))
    grad[k] <- (m * n / sigma2) * d_sigma2 + m * sum(Rxinv * dRx)
  }

  for (k in 4:5) {
    Ez_k <- if (k == 4) Ez else Ez1
    dRz  <- Rz * Ez_k / theta[k]^2
    dKz  <- -Rzinv %*% dRz %*% Rzinv
    da_mat <- Rxinv %*% rep(1, n) %*% t(rep(1, m)) %*% dKz
    db_mat <- Rxinv %*% y %*% dKz
    da <- c(t(da_mat)); db <- c(t(db_mat))
    d_mu <- dmu(da, db)
    d_sigma2 <- mean(-d_mu*a*(b-mu*a) - mu*da*(b-mu*a) + (vec.y-mu*a)*(db-d_mu*a-mu*da))
    grad[k] <- (m * n / sigma2) * d_sigma2 + n * sum(Rzinv * dRz)
  }

  attr(val, "gradient") <- grad
  val
}

ML_log <- function(phi) {
  theta <- exp(phi)
  out   <- ML(theta)
  g     <- attr(out, "gradient") * theta
  attr(out, "gradient") <- g
  out
}

fit <- nlminb(
  start     = log(c(0.1, 0.1, 0.1, 1, 1)),
  objective = function(p) ML_log(p),
  gradient  = function(p) attr(ML_log(p), "gradient"),
  control   = list(iter.max = 200, rel.tol = 1e-8)
)
theta <- exp(fit$par)
cat("Turbine", i, "nlminb convergence code:", fit$convergence, "\n")

Rx    <- 1/(1 + (Ex1/theta[1])^2 + (Ex2/theta[2])^2 + (Ex3/theta[3])^2)
Rxinv <- solve(Rx + 1e-8 * diag(n))
Rz    <- exp(-Ez/theta[4] - Ez1/theta[5])
Rzinv <- solve(Rz + 1e-8 * diag(m))

a <- c(t((Rxinv %*% rep(1, n)) %*% (t(rep(1, m)) %*% Rzinv)))
b <- c(t(Rxinv %*% as.matrix(y) %*% Rzinv))
mu <- sum(b) / sum(a)

A <- matrix(a, nrow = n, ncol = m, byrow = TRUE)
B <- matrix(b, nrow = n, ncol = m, byrow = TRUE)

# OK coefficient (same form as code/Table2-Table3(STGP).R)
COEF <- Rxinv %*% (y - mu) %*% Rzinv

basis.x <- function(h) 1/(1 + sum((h/theta[1:3])^2))
basis.z <- function(h, h1) exp(-abs(h)/theta[4] - abs(h1)/theta[5])
r.x <- function(u) {
  Auu <- t(t(x) - u)
  apply(Auu, 1, basis.x)
}
r.z <- function(v, v1) basis.z(v - z, v1 - z1)

u <- c(scaled_terrain_data[i, 2], scaled_terrain_data[i, 3], scaled_terrain_data[i, 4])
rx_vec <- r.x(u)

yhat_lk <- numeric(length(z.ev))
yhat_ok <- numeric(length(z.ev))

for (j in seq_along(z.ev)) {
  rz_vec <- r.z(z.ev[j], z.ev1[j])

  denom <- as.numeric(t(rx_vec) %*% A %*% rz_vec)
  num   <- as.numeric(t(rx_vec) %*% B %*% rz_vec)
  yhat_lk[j] <- num / denom

  yhat_ok[j] <- as.numeric(mu + t(rx_vec %*% COEF %*% rz_vec))
}

rmse_lk <- sqrt(mean((f.ev - yhat_lk)^2, na.rm = TRUE))
rmse_ok <- sqrt(mean((f.ev - yhat_ok)^2, na.rm = TRUE))
cat("Turbine", i, "RMSE -- LK:", rmse_lk, " OK:", rmse_ok, "\n")

write.csv(
  data.frame(wind_speed = wind_speed_ev, actual_power = f.ev,
             pred_lk = yhat_lk, pred_ok = yhat_ok),
  file.path(figures_dir, paste0("fig4_turbine", TARGET_TURBINE, "_predictions.csv")),
  row.names = FALSE
)

# -------------------------------------------------
# Plot: actual data (black) with LK (green) and OK (red) predictions,
# one sub-panel each, matching Figure 4's right-hand layout
# -------------------------------------------------
pdf(file.path(figures_dir, "fig4_turbine60_lk_ok.pdf"), width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(5, 5, 2, 1), cex.lab = 1.2, cex.axis = 1.0)

plot(wind_speed_ev, f.ev, pch = 16, cex = 0.4, col = "black",
     xlab = "Wind Speed", ylab = "Scaled Power", main = "LK")
points(wind_speed_ev, yhat_lk, pch = 16, cex = 0.4, col = "forestgreen")

plot(wind_speed_ev, f.ev, pch = 16, cex = 0.4, col = "black",
     xlab = "Wind Speed", ylab = "Scaled Power", main = "OK")
points(wind_speed_ev, yhat_ok, pch = 16, cex = 0.4, col = "red")

dev.off()
