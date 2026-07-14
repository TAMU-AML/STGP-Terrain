# code/figures/fig_5.R
#
# Reproduces Figure 5: Accumulated Local Effects (ALE) of Power with respect
# to each covariate -- temporal (wind speed, temperature) and spatial
# (slope, RIX, ridge height) -- from the STGP model fit on the FULL dataset
# (all 66 turbines, no leave-one-out).
#
# For the temporal ALE curves, the model is evaluated once per turbine
# (holding that turbine's own terrain fixed) and the resulting curves are
# averaged. For the spatial ALE curves, the model is evaluated once per
# time/wind grid point (holding speed/temperature fixed) and averaged.
#
# Requires the ALEPlot package: install.packages("ALEPlot")
#
# Output:
#   results/figures/fig5_temporal_ale.pdf  (speed, temperature -- 2 panels)
#   results/figures/fig5_spatial_ale.pdf   (slope, RIX, ridge -- 3 panels)
#   results/figures/fig5_ale_ranges.csv    (relative ALE range per covariate)

data_dir <- file.path("data")
input_folder <- file.path("data", "processed data")
terrain_data_path <- file.path("data", "weightedTerrainData.csv")
figures_dir <- file.path("results", "figures")

stopifnot(dir.exists(input_folder))
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(ALEPlot)

# -------------------------------------------------
# Data loading (same convention as Table2-Table3(STGP).R)
# -------------------------------------------------
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

temp_vector  <- read_vec_csv(file.path(input_folder, "temp_vector.csv"), colname = "temp")
speed_vector <- read_vec_csv(file.path(input_folder, "speed_vector.csv"), colname = "speed")
power_matrix <- read_mat_csv(file.path(input_folder, "power_matrix.csv"))

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

# -------------------------------------------------
# Fit STGP on the FULL dataset (all 66 turbines, no leave-one-out)
# -------------------------------------------------
x <- scaled_terrain_data[, 2:4]
y <- power_matrix
z  <- scale_01_test(speed_vector)
z1 <- scale_01_test_temperature(temp_vector)

sorted_indices <- order(z)
z  <- z[sorted_indices]
z1 <- z1[sorted_indices]
y  <- power_matrix[, sorted_indices]

n <- dim(x)[1]
m <- length(z)
vec.y <- c(t(y))

Ex1 <- as.matrix(dist(x[, 1], diag = TRUE, upper = TRUE))
Ex2 <- as.matrix(dist(x[, 2], diag = TRUE, upper = TRUE))
Ex3 <- as.matrix(dist(x[, 3], diag = TRUE, upper = TRUE))
Ez  <- as.matrix(dist(z,  diag = TRUE, upper = TRUE))
Ez1 <- as.matrix(dist(z1, diag = TRUE, upper = TRUE))

ML <- function(para) {
  Rx <- 1 / (1 + (Ex1/para[1])^2 + (Ex2/para[2])^2 + (Ex3/para[3])^2)
  Rxinv <- solve(Rx + 1e-8 * diag(n))
  Rz <- exp(-Ez/para[4] - Ez1/para[5])
  Rzinv <- solve(Rz + 1e-8 * diag(m))
  a <- c(t((Rxinv %*% rep(1, n)) %*% (t(rep(1, m)) %*% Rzinv)))
  b <- c(t(Rxinv %*% as.matrix(y) %*% Rzinv))
  mu <- sum(b) / sum(a)
  sigma2 <- (1/(m*n)) * sum((vec.y - mu) * (b - mu * a))
  m*n*log(sigma2) + m*determinant(Rx, logarithm = TRUE)$mod[1] + n*determinant(Rz, logarithm = TRUE)$mod[1]
}

ini <- c(0.1, 0.1, 0.1, 1, 1)
a.opt <- optim(ini, ML, lower = ini/100, upper = ini*100, method = "L-BFGS-B")
theta <- a.opt$par
cat("Fitted theta:", theta, "\n")

Rx <- 1/(1 + (Ex1/theta[1])^2 + (Ex2/theta[2])^2 + (Ex3/theta[3])^2)
Rxinv <- solve(Rx + 1e-8 * diag(n))
Rz <- exp(-Ez/theta[4] - Ez1/theta[5])
Rzinv <- solve(Rz + 1e-8 * diag(m))
a <- c(t((Rxinv %*% rep(1, n)) %*% (t(rep(1, m)) %*% Rzinv)))
b <- c(t(Rxinv %*% as.matrix(y) %*% Rzinv))

# -------------------------------------------------
# ALE basis functions
# -------------------------------------------------
basis.x <- function(h, theta) 1 / (1 + sum((h / theta[1:3])^2))
basis.z <- function(h, h1, theta) exp(-abs(h) / theta[4] - abs(h1) / theta[5])

r.x <- function(u, x, theta) {
  t(t(x) - u) |> apply(1, function(h) basis.x(h, theta))
}
r.z <- function(v, v1, z, z1, theta) basis.z(v - z, v1 - z1, theta)

# -------------------------------------------------
# Temporal ALE (speed, temperature): one model call per turbine, averaged
# -------------------------------------------------
K <- 200
num_turbines <- nrow(x)
X_temporal <- data.frame(speed = z, temperature = z1)

f_speed_all <- matrix(0, nrow = K + 1, ncol = num_turbines)
f_temp_all  <- matrix(0, nrow = K + 1, ncol = num_turbines)

for (i in 1:num_turbines) {
  u_fixed <- as.numeric(x[i, ])
  gp_model <- list(
    A = matrix(a, nrow = n, ncol = m, byrow = TRUE),
    B = matrix(b, nrow = n, ncol = m, byrow = TRUE),
    x = x, z = z, z1 = z1, theta = theta, u = u_fixed
  )
  predict.gp <- function(X.model, newdata) {
    apply(newdata, 1, function(row) {
      rx <- r.x(X.model$u, X.model$x, X.model$theta)
      rz <- r.z(row["speed"], row["temperature"], X.model$z, X.model$z1, X.model$theta)
      as.numeric((t(rx) %*% X.model$B %*% rz) / (t(rx) %*% X.model$A %*% rz))
    })
  }
  ALE.speed <- ALEPlot(X_temporal, gp_model, predict.gp, J = 1, K = K)
  ALE.temp  <- ALEPlot(X_temporal, gp_model, predict.gp, J = 2, K = K)
  f_speed_all[, i] <- ALE.speed$f.values
  f_temp_all[, i]  <- ALE.temp$f.values
}

mean_speed <- rowMeans(f_speed_all)
mean_temp  <- rowMeans(f_temp_all)
x_vals_speed <- ALE.speed$x.values
x_vals_temp  <- ALE.temp$x.values

# -------------------------------------------------
# Spatial ALE (slope, RIX, ridge): one model call per time point, averaged
# -------------------------------------------------
K_spatial <- nrow(x) - 1  # matches original script's K = 65 for n = 66 terrain rows
num_time <- length(z)
x_spatial <- data.frame(slope = x[, 1], rix = x[, 2], ridge = x[, 3])

f_slope_all <- matrix(0, nrow = K_spatial + 1, ncol = num_time)
f_rix_all   <- matrix(0, nrow = K_spatial + 1, ncol = num_time)
f_ridge_all <- matrix(0, nrow = K_spatial + 1, ncol = num_time)

for (j in 1:num_time) {
  gp_model_spatial <- list(
    A = matrix(a, nrow = n, ncol = m, byrow = TRUE),
    B = matrix(b, nrow = n, ncol = m, byrow = TRUE),
    x = x, z = z, z1 = z1, theta = theta, v = z[j], v1 = z1[j]
  )
  predict.gp.spatial <- function(X.model, newdata) {
    apply(newdata, 1, function(row) {
      rx <- r.x(row, X.model$x, X.model$theta)
      rz <- r.z(X.model$v, X.model$v1, X.model$z, X.model$z1, X.model$theta)
      as.numeric((t(rx) %*% X.model$B %*% rz) / (t(rx) %*% X.model$A %*% rz))
    })
  }
  ALE.slope <- ALEPlot(x_spatial, gp_model_spatial, predict.gp.spatial, J = 1, K = K_spatial)
  ALE.rix   <- ALEPlot(x_spatial, gp_model_spatial, predict.gp.spatial, J = 2, K = K_spatial)
  ALE.ridge <- ALEPlot(x_spatial, gp_model_spatial, predict.gp.spatial, J = 3, K = K_spatial)
  f_slope_all[, j] <- ALE.slope$f.values
  f_rix_all[, j]   <- ALE.rix$f.values
  f_ridge_all[, j] <- ALE.ridge$f.values
}

mean_slope <- rowMeans(f_slope_all)
mean_rix   <- rowMeans(f_rix_all)
mean_ridge <- rowMeans(f_ridge_all)
x_vals_spatial <- ALE.slope$x.values

# -------------------------------------------------
# Plots
# -------------------------------------------------
pdf(file.path(figures_dir, "fig5_spatial_ale.pdf"), width = 11, height = 4)
par(mfrow = c(1, 3), mar = c(5, 5, 2, 1), cex.lab = 1.5, cex.axis = 1.3)

plot(x_vals_spatial, f_slope_all[, 1], type = "l", col = "gray", lwd = 0.5,
     ylim = range(f_slope_all), xlab = "Slope", ylab = "ALE of Power")
for (j in 2:num_time) lines(x_vals_spatial, f_slope_all[, j], col = "gray", lwd = 0.5)
lines(x_vals_spatial, mean_slope, col = "black", lwd = 2)

plot(x_vals_spatial, f_rix_all[, 1], type = "l", col = "gray", lwd = 0.5,
     ylim = range(f_rix_all), xlab = "RIX", ylab = "ALE of Power")
for (j in 2:num_time) lines(x_vals_spatial, f_rix_all[, j], col = "gray", lwd = 0.5)
lines(x_vals_spatial, mean_rix, col = "black", lwd = 2)

plot(x_vals_spatial, f_ridge_all[, 1], type = "l", col = "gray", lwd = 0.5,
     ylim = range(f_ridge_all), xlab = "Ridge", ylab = "ALE of Power")
for (j in 2:num_time) lines(x_vals_spatial, f_ridge_all[, j], col = "gray", lwd = 0.5)
lines(x_vals_spatial, mean_ridge, col = "black", lwd = 2)

dev.off()

pdf(file.path(figures_dir, "fig5_temporal_ale.pdf"), width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(5, 5, 2, 1), cex.lab = 1.2, cex.axis = 1.0)

plot(x_vals_speed, f_speed_all[, 1], type = "l", col = "gray", lwd = 0.5,
     ylim = range(f_speed_all), xlab = "Speed", ylab = "ALE of Power")
for (i in 2:num_turbines) lines(x_vals_speed, f_speed_all[, i], col = "gray", lwd = 0.5)
lines(x_vals_speed, mean_speed, col = "black", lwd = 2)

plot(x_vals_temp, f_temp_all[, 1], type = "l", col = "gray", lwd = 0.5,
     ylim = range(f_temp_all), xlab = "Temperature", ylab = "ALE of Power")
for (i in 2:num_turbines) lines(x_vals_temp, f_temp_all[, i], col = "gray", lwd = 0.5)
lines(x_vals_temp, mean_temp, col = "black", lwd = 2)

dev.off()

# -------------------------------------------------
# Feature-importance ranges (relative ALE range -> weights used in d_WD)
# -------------------------------------------------
range_speed <- max(mean_speed) - min(mean_speed)
range_temp  <- max(mean_temp)  - min(mean_temp)
range_slope <- max(mean_slope) - min(mean_slope)
range_rix   <- max(mean_rix)   - min(mean_rix)
range_ridge <- max(mean_ridge) - min(mean_ridge)

cat("ALE ranges (speed, temp, slope, rix, ridge):\n")
print(c(speed = range_speed, temp = range_temp, slope = range_slope, rix = range_rix, ridge = range_ridge))

write.csv(
  data.frame(covariate = c("speed", "temperature", "slope", "rix", "ridge"),
             ale_range = c(range_speed, range_temp, range_slope, range_rix, range_ridge)),
  file.path(figures_dir, "fig5_ale_ranges.csv"), row.names = FALSE
)
