# code/supplement/S1.R
#
# Reproduces the X+S ("terrain as extra covariates") TwinGP row of Table S1.
# TwinGP(x) is NOT recomputed here -- it's identical to the "TwinGP" row in
# Table 2 / Table 3 of the main results.
#
# Reuses the data loading, LOO slicing, and the (previously commented-out)
# TwinGP(x+s) block from code/Table2-Table3(twinGP+Binning).R.
# Writes results/supplement/table_s1.csv directly -- no separate helper file.

data_path <- file.path("data")
terrain_path <- file.path("data", "weightedTerrainData.csv")
results_dir <- file.path("results", "supplement")

library(data.table)
library(twingp)
library(dplyr)
library(readr)
library(tidyr)

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
table_s1_path <- file.path(results_dir, "table_s1.csv")

# -------------------------
# Inline table_s1.csv writer (no separate helper file)
# -------------------------
update_table_s1 <- function(method, version,
                             rmse_2017 = NA_real_, rmse_2018 = NA_real_,
                             nlpd_2017 = NA_real_, nlpd_2018 = NA_real_) {
  cols <- c("Method", "Version", "RMSE_2017", "RMSE_2018", "NLPD_2017", "NLPD_2018")

  if (file.exists(table_s1_path)) {
    df <- read.csv(table_s1_path, stringsAsFactors = FALSE)
  } else {
    df <- data.frame(matrix(ncol = length(cols), nrow = 0))
    colnames(df) <- cols
  }

  match_idx <- which(df$Method == method & df$Version == version)
  new_row <- data.frame(
    Method = method, Version = version,
    RMSE_2017 = rmse_2017, RMSE_2018 = rmse_2018,
    NLPD_2017 = nlpd_2017, NLPD_2018 = nlpd_2018,
    stringsAsFactors = FALSE
  )

  if (length(match_idx) > 0) {
    df[match_idx[1], cols] <- new_row[cols]
  } else {
    df <- rbind(df, new_row[cols])
  }

  write.csv(df, table_s1_path, row.names = FALSE)
  cat("[table_s1] updated", table_s1_path, "->", method, "(", version, ")\n")
}

turbine_ids <- 1:66
testset_2018 <- c(1:46, 48:50, 52, 54:60, 62:66)

terrain_data <- read.csv(terrain_path)
scale_01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
terrain_data[, 2:4] <- lapply(terrain_data[, 2:4, drop = FALSE], scale_01)

terrain_mat <- as.matrix(terrain_data[1:66, 2:4])
storage.mode(terrain_mat) <- "double"

# -------------------------
# Load all turbine-year CSVs ONCE (cache)
# -------------------------
cache_key <- function(id, year) sprintf("T%02d_%d", id, year)
data_cache <- vector("list", length = length(turbine_ids) * 2)
names(data_cache) <- as.vector(outer(sprintf("T%02d", turbine_ids), c("2017", "2018"), paste, sep = "_"))

for (id in turbine_ids) {
  for (yr in c(2017, 2018)) {
    f <- sprintf("%s/Turbine%d_%d.csv", data_path, id, yr)
    d <- tryCatch(fread(file = f, showProgress = FALSE), error = function(e) NULL)
    if (!is.null(d)) {
      d[, wind_speed := as.numeric(wind_speed)]
      d[, temperature := as.numeric(temperature)]
      d[, power := as.numeric(power)]
    }
    data_cache[[cache_key(id, yr)]] <- d
  }
}
get_data <- function(id, year) data_cache[[cache_key(id, year)]]

# -------------------------
# Full 2017 training pool (built once)
# -------------------------
X2017_list <- vector("list", 66); y2017_list <- vector("list", 66); tid2017_list <- vector("list", 66)
for (j in turbine_ids) {
  d <- get_data(j, 2017)
  if (!is.null(d)) {
    X2017_list[[j]] <- cbind(d$wind_speed, d$temperature)
    y2017_list[[j]] <- d$power
    tid2017_list[[j]] <- rep(j, nrow(d))
  } else {
    X2017_list[[j]] <- matrix(numeric(0), ncol = 2)
    y2017_list[[j]] <- numeric(0)
    tid2017_list[[j]] <- integer(0)
  }
}
X2017  <- do.call(rbind, X2017_list)
y2017  <- as.numeric(unlist(y2017_list, use.names = FALSE))
tid2017 <- as.integer(unlist(tid2017_list, use.names = FALSE))
S2017 <- terrain_mat[tid2017, , drop = FALSE]

S_test_for <- function(i, n) matrix(rep(terrain_mat[i, ], each = n), nrow = n, ncol = 3)

get_train_loo_2017 <- function(leave_out_id) {
  idx <- tid2017 != leave_out_id
  list(X = X2017[idx, , drop = FALSE], y = y2017[idx],
       tid = tid2017[idx], S = S2017[idx, , drop = FALSE])
}

# -------------------------
# TwinGP(x+s) LOO
# -------------------------
results_long <- data.frame()

for (year in c(2017, 2018)) {
  test_ids <- if (year == 2017) turbine_ids else testset_2018

  for (i in test_ids) {
    cat("[S1 TwinGP X+S] Turbine", i, "Year", year, "\n")

    test_data <- get_data(i, year)
    if (is.null(test_data)) next

    test_speed <- test_data$wind_speed
    test_temp  <- test_data$temperature
    test_power <- test_data$power

    X_test  <- cbind(test_speed, test_temp)
    Xs_test <- cbind(X_test, S_test_for(i, nrow(X_test)))

    tr <- get_train_loo_2017(i)
    Xs_train <- cbind(tr$X, tr$S)
    y_train  <- tr$y

    set.seed(i)
    t1 <- Sys.time()
    twin_out_s <- twingp(x = Xs_train, y = y_train, x_test = Xs_test)
    t2 <- Sys.time()

    pred_s <- as.numeric(twin_out_s$mu)
    pred_sd_s <- as.numeric(twin_out_s$sigma)

    rmse_s <- sqrt(mean((pred_s - test_power)^2, na.rm = TRUE))
    nlpd_s <- mean(0.5 * log(2 * pi * pred_sd_s^2) +
                     0.5 * ((test_power - pred_s)^2) / (pred_sd_s^2),
                   na.rm = TRUE)

    results_long <- rbind(results_long, data.frame(
      Method = "TwinGP", Version = "X+S", Turbine = i, Year = year,
      RMSE = rmse_s, NLPD = nlpd_s,
      Runtime = round(as.numeric(difftime(t2, t1, units = "secs")), 4)
    ))
  }
}

write_csv(results_long, file.path(results_dir, "s1_twingp_long.csv"))

summary_table <- results_long %>%
  group_by(Year) %>%
  summarise(RMSE = mean(RMSE, na.rm = TRUE), NLPD = mean(NLPD, na.rm = TRUE), .groups = "drop")

rmse_2017 <- summary_table$RMSE[summary_table$Year == 2017]
rmse_2018 <- summary_table$RMSE[summary_table$Year == 2018]
nlpd_2017 <- summary_table$NLPD[summary_table$Year == 2017]
nlpd_2018 <- summary_table$NLPD[summary_table$Year == 2018]

update_table_s1("TwinGP", "X+S",
                 rmse_2017 = rmse_2017, rmse_2018 = rmse_2018,
                 nlpd_2017 = nlpd_2017, nlpd_2018 = nlpd_2018)
