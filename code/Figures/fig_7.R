# code/figures/fig_7.R
#
# Reproduces Figure 7: naive vs. terrain-adjusted turbine underperformance.
#
# For each turbine i:
#   P1 = turbine i's own total power over the support-point grid (actual)
#   P2 = average of that total across all OTHER turbines (fleet baseline)
#   P3 = STGP-predicted total power over the grid, using turbine i's own
#        terrain, with the STGP model fit leave-one-out (turbine i excluded)
#   P4 = same prediction, but using the AVERAGE terrain of all other
#        turbines instead of turbine i's own terrain
#
#   naive_underperf     = (P2 - P1) / P2 * 100   -- deviation from fleet avg
#   terrain_effect       = (P4 - P3) / P2 * 100   -- how much of that the
#                                                     model attributes to
#                                                     turbine i's terrain
#   adjusted_underperf  = naive_underperf - terrain_effect
#
# Turbines where the sign of naive vs. adjusted underperformance disagrees
# are flagged as a false alarm (naive flags underperformance, terrain-
# adjusted does not) or a missed case (the reverse).
#
# NOTE: this refits the STGP model leave-one-out for all 66 turbines, so
# it's on the same order of runtime as the main STGP results (~4-5 hours).
#
# Output:
#   results/figures/fig7_underperformance.pdf
#   results/figures/fig7_underperformance_data.csv

data_dir <- file.path("data")
input_folder <- file.path("data", "processed data")
terrain_data_path <- file.path("data", "weightedTerrainData.csv")
figures_dir <- file.path("results", "figures")

stopifnot(dir.exists(input_folder))
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)

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
# P1, P2: naive fleet comparison (no model needed)
# -------------------------------------------------
power_sums <- rowSums(power_matrix)  # total power per turbine over the grid

P1 <- power_sums
P2 <- vapply(seq_along(power_sums), function(i) mean(power_sums[-i]), numeric(1))

# -------------------------------------------------
# P3, P4: STGP-predicted totals (LOO fit per turbine), own vs. average terrain
# -------------------------------------------------
P3 <- numeric(66)
P4 <- numeric(66)

for (i in 1:66) {
  cat("=== Fig7 Turbine", i, "===\n")

  x <- scaled_terrain_data[, 2:4]
  x <- x[-i, ]
  y <- power_matrix
  z  <- scale_01_test(speed_vector)
  z1 <- scale_01_test_temperature(temp_vector)

  sorted_indices <- order(z)
  z  <- z[sorted_indices]
  z1 <- z1[sorted_indices]
  y  <- power_matrix[, sorted_indices]
  y  <- y[-i, ]

  n <- dim(x)[1]
  m <- length(z)
  vec.y <- c(t(y))

  Ex1 <- as.matrix(dist(x[,1], diag = TRUE, upper = TRUE))
  Ex2 <- as.matrix(dist(x[,2], diag = TRUE, upper = TRUE))
  Ex3 <- as.matrix(dist(x[,3], diag = TRUE, upper = TRUE))
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
    y_vec <- vec.y
    sigma2 <- mean((y_vec - mu*a) * (b - mu*a))

    logdetRx <- determinant(Rx, logarithm = TRUE)$modulus[1]
    logdetRz <- determinant(Rz, logarithm = TRUE)$modulus[1]
    val <- m*n*log(sigma2) + m*logdetRx + n*logdetRz

    grad <- numeric(5)
    denom_a <- sum(a); denom_b <- sum(b)
    dmu <- function(da, db) (sum(db)*denom_a - denom_b*sum(da)) / denom_a^2

    for (k in 1:3) {
      Exk <- switch(k, Ex1, Ex2, Ex3)
      dRx <- 2 * Rx^2 * (Exk^2) / theta[k]^3
      dKx <- -Rxinv %*% dRx %*% Rxinv
      da_mat <- dKx %*% rep(1, n) %*% t(rep(1, m)) %*% Rzinv
      db_mat <- dKx %*% y %*% Rzinv
      da <- c(t(da_mat)); db <- c(t(db_mat))
      d_mu <- dmu(da, db)
      d_sigma2 <- mean(-d_mu*a*(b-mu*a) - mu*da*(b-mu*a) + (y_vec-mu*a)*(db-d_mu*a-mu*da))
      grad[k] <- (m*n/sigma2)*d_sigma2 + m*sum(Rxinv*dRx)
    }

    for (k in 4:5) {
      Ez_k <- if (k == 4) Ez else Ez1
      dRz <- Rz * Ez_k / theta[k]^2
      dKz <- -Rzinv %*% dRz %*% Rzinv
      da_mat <- Rxinv %*% rep(1, n) %*% t(rep(1, m)) %*% dKz
      db_mat <- Rxinv %*% y %*% dKz
      da <- c(t(da_mat)); db <- c(t(db_mat))
      d_mu <- dmu(da, db)
      d_sigma2 <- mean(-d_mu*a*(b-mu*a) - mu*da*(b-mu*a) + (y_vec-mu*a)*(db-d_mu*a-mu*da))
      grad[k] <- (m*n/sigma2)*d_sigma2 + n*sum(Rzinv*dRz)
    }

    attr(val, "gradient") <- grad
    val
  }

  ML_log <- function(phi) {
    theta <- exp(phi)
    out <- ML(theta)
    g <- attr(out, "gradient") * theta
    attr(out, "gradient") <- g
    out
  }

  fit <- nlminb(start = log(c(0.1,0.1,0.1,1,1)),
                objective = function(p) ML_log(p),
                gradient = function(p) attr(ML_log(p), "gradient"),
                control = list(iter.max = 200, rel.tol = 1e-8))
  theta <- exp(fit$par)

  Rx <- 1/(1+(Ex1/theta[1])^2+(Ex2/theta[2])^2+(Ex3/theta[3])^2)
  Rxinv <- solve(Rx + 1e-8*diag(n))
  Rz <- exp(-Ez/theta[4] - Ez1/theta[5])
  Rzinv <- solve(Rz + 1e-8*diag(m))
  a <- c(t((Rxinv %*% rep(1,n)) %*% (t(rep(1,m)) %*% Rzinv)))
  b <- c(t(Rxinv %*% as.matrix(y) %*% Rzinv))

  basis.x <- function(h) 1/(1 + sum((h/theta[1:3])^2))
  basis.z <- function(h, h1) exp(-abs(h)/theta[4] - abs(h1)/theta[5])
  r.x <- function(u) apply(t(t(x) - u), 1, basis.x)
  r.z <- function(v, v1) basis.z(v - z, v1 - z1)

  A <- matrix(a, nrow = n, ncol = m, byrow = TRUE)
  B <- matrix(b, nrow = n, ncol = m, byrow = TRUE)

  predict_total <- function(u) {
    rx_vec <- r.x(u)
    yhat <- numeric(length(z))
    for (j in seq_along(z)) {
      rz_vec <- r.z(z[j], z1[j])
      denom <- as.numeric(t(rx_vec) %*% A %*% rz_vec)
      num   <- as.numeric(t(rx_vec) %*% B %*% rz_vec)
      yhat[j] <- num / denom
    }
    sum(yhat)
  }

  u_own <- c(scaled_terrain_data[i, 2], scaled_terrain_data[i, 3], scaled_terrain_data[i, 4])
  P3[i] <- predict_total(u_own)

  u_avg <- colMeans(scaled_terrain_data[-i, 2:4])
  P4[i] <- predict_total(u_avg)
}

# -------------------------------------------------
# Underperformance rates
# -------------------------------------------------
naive_underperf    <- (P2 - P1) / P2 * 100
terrain_effect     <- (P4 - P3) / P2 * 100
adjusted_underperf <- naive_underperf - terrain_effect

false_alarm <- which(naive_underperf > 0 & adjusted_underperf < 0)
missed_case <- which(naive_underperf < 0 & adjusted_underperf > 0)

write.csv(
  data.frame(Turbine_ID = 1:66, P1 = P1, P2 = P2, P3 = P3, P4 = P4,
             Naive_Underperf = naive_underperf, Terrain_Effect = terrain_effect,
             Adjusted_Underperf = adjusted_underperf),
  file.path(figures_dir, "fig7_underperformance_data.csv"), row.names = FALSE
)

# -------------------------------------------------
# Plot
# -------------------------------------------------
pdf(file.path(figures_dir, "fig7_underperformance.pdf"), width = 8, height = 5)
par(mar = c(5.1, 4.1, 2, 2.1))

plot(naive_underperf, type = "l", col = "darkred", lwd = 2,
     xlab = "Turbine Index", ylab = "Underperformance Rate (%)", main = "")
lines(adjusted_underperf, col = "darkblue", lwd = 2)
abline(h = 0, col = "gray60", lty = 2)

points(false_alarm, adjusted_underperf[false_alarm], col = "forestgreen", pch = 19)
points(missed_case, adjusted_underperf[missed_case], col = "purple", pch = 19)

text(false_alarm, adjusted_underperf[false_alarm] - 3,
     labels = false_alarm, col = "forestgreen", cex = 0.9)
text(missed_case, adjusted_underperf[missed_case] + 4,
     labels = missed_case, col = "purple", cex = 0.9)

legend("topleft",
       legend = c("Naive", "Adjusted", "+Naive / -Adjusted", "-Naive / +Adjusted"),
       col = c("darkred", "darkblue", "forestgreen", "purple"),
       lwd = c(2, 2, NA, NA), pch = c(NA, NA, 19, 19),
       bty = "n", cex = 0.8)

dev.off()
