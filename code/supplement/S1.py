# -*- coding: utf-8 -*-
"""
code/supplement/S1.py

Reproduces the Table S1 rows NOT already covered by the main results:
  - XGBoost:        X+S   (X-only already in Table 2/3)
  - Multi-layer NN: X+S   (X-only already in Table 2/3)
  - Bayesian NN:    X     (X+S already in Table 2/3 -- its code always
                           concatenates terrain, so the existing "BNN" row
                           in the main results IS the X+S version)
TwinGP(X+S) is handled separately in S1.R.

Terrain scaling follows each method's own existing convention:
  - XGBoost: raw terrain, no scaling (matches Table2-Table3(XGBoost).ipynb)
  - NN / BNN: terrain standardized once globally via StandardScaler
    (matches Table2-Table3(BNN).ipynb; applied to NN by analogy since NN
    also standardizes its dynamic features, unlike XGBoost)

Output: results/supplement/table_s1.csv (Method, Version, RMSE_2017, RMSE_2018)
Written directly by this script -- no separate helper module.
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from sklearn.preprocessing import StandardScaler

# -------------------------------------------------
# Paths
# -------------------------------------------------
PROJECT_ROOT = Path.cwd()
if PROJECT_ROOT.name in ("code", "supplement"):
    PROJECT_ROOT = PROJECT_ROOT.parent
    if PROJECT_ROOT.name == "code":
        PROJECT_ROOT = PROJECT_ROOT.parent

DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "supplement"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

TERRAIN_PATH = DATA_DIR / "weightedTerrainData.csv"
OUT_LONG = RESULTS_DIR / "s1_python_long.csv"
TABLE_S1_PATH = RESULTS_DIR / "table_s1.csv"

if not TERRAIN_PATH.exists():
    raise FileNotFoundError(f"Terrain file not found: {TERRAIN_PATH}")


# -------------------------------------------------
# Inline table_s1.csv writer (no separate helper file)
# -------------------------------------------------
def update_table_s1(method, version, rmse_2017=np.nan, rmse_2018=np.nan,
                     nlpd_2017=np.nan, nlpd_2018=np.nan):
    cols = ["Method", "Version", "RMSE_2017", "RMSE_2018", "NLPD_2017", "NLPD_2018"]
    df = pd.read_csv(TABLE_S1_PATH) if TABLE_S1_PATH.exists() else pd.DataFrame(columns=cols)

    mask = (df["Method"] == method) & (df["Version"] == version)
    row = {"Method": method, "Version": version,
           "RMSE_2017": rmse_2017, "RMSE_2018": rmse_2018,
           "NLPD_2017": nlpd_2017, "NLPD_2018": nlpd_2018}

    if mask.any():
        for k, v in row.items():
            df.loc[mask, k] = v
    else:
        df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)

    df[cols].to_csv(TABLE_S1_PATH, index=False)
    print(f"[table_s1] updated {TABLE_S1_PATH} -> {method} ({version})")


# -------------------------------------------------
# IDs (same split as Table2-3)
# -------------------------------------------------
TURBINE_IDS = list(range(1, 67))
TESTSET_2018 = list(range(1, 47)) + [48, 49, 50, 52] + list(range(54, 61)) + list(range(62, 67))

FEATURES = ["wind_speed", "temperature"]
TARGET = "power"

# -------------------------------------------------
# Terrain: two versions, one per convention
# -------------------------------------------------
terrain_df = pd.read_csv(TERRAIN_PATH)
terrain_cols = terrain_df.columns[1:4]
terrain_mat_raw = terrain_df.loc[:65, terrain_cols].values          # XGBoost: unscaled
terrain_mat_scaled = StandardScaler().fit_transform(
    terrain_mat_raw.astype(np.float32)
)                                                                     # NN / BNN: standardized

# -------------------------------------------------
# Helpers
# -------------------------------------------------
def load_turbine_csv(tid, year):
    f = DATA_DIR / f"Turbine{tid}_{year}.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f)
    for c in ["wind_speed", "temperature", "power"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def rmse(yhat, y):
    yhat = np.asarray(yhat, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(yhat) & np.isfinite(y)
    if not np.any(m):
        return np.nan
    return float(np.sqrt(np.mean((yhat[m] - y[m]) ** 2)))


data_cache = {}
for tid in TURBINE_IDS:
    for yr in (2017, 2018):
        data_cache[(tid, yr)] = load_turbine_csv(tid, yr)

# Full 2017 pool, built once
X2017_list, y2017_list, tid2017_list = [], [], []
for tid in TURBINE_IDS:
    df = data_cache[(tid, 2017)]
    if df is None:
        continue
    X2017_list.append(df[FEATURES].values)
    y2017_list.append(df[TARGET].values)
    tid2017_list.append(np.full(len(df), tid, dtype=int))

X2017 = np.vstack(X2017_list)
y2017 = np.concatenate(y2017_list)
tid2017 = np.concatenate(tid2017_list)


def run_loo(method_name, fit_predict_fn, terrain_mat):
    """
    fit_predict_fn(Xs_train, y_train, Xs_test, seed) -> yhat
    LOO across all 66 turbines (2017) and the Table3 subset (2018),
    with `terrain_mat` concatenated to [wind_speed, temperature].
    """
    S2017 = terrain_mat[tid2017 - 1]
    results = []
    for year in (2017, 2018):
        test_ids = TURBINE_IDS if year == 2017 else TESTSET_2018
        for i in test_ids:
            df_test = data_cache[(i, year)]
            if df_test is None:
                continue
            print(f"[{method_name} X+S] Turbine {i}, Year {year}")

            X_test = df_test[FEATURES].values
            y_test = df_test[TARGET].values
            S_test = np.tile(terrain_mat[i - 1], (len(X_test), 1))
            Xs_test = np.hstack([X_test, S_test])

            mask = tid2017 != i
            Xs_train = np.hstack([X2017[mask], S2017[mask]])
            y_train = y2017[mask]

            t0 = time.time()
            pred = fit_predict_fn(Xs_train, y_train, Xs_test, seed=i)
            t1 = time.time()

            results.append({
                "Method": method_name, "Version": "X+S", "Turbine": i, "Year": year,
                "RMSE": rmse(pred, y_test), "Runtime": round(t1 - t0, 4),
            })
    return pd.DataFrame(results)


# -------------------------------------------------
# XGBoost (LightGBM) X+S -- raw terrain, no scaling
# -------------------------------------------------
def fit_xgb(X_train, y_train, X_test, seed):
    model = LGBMRegressor(
        objective="regression", n_estimators=200, learning_rate=0.1,
        max_depth=8, subsample=0.8, colsample_bytree=0.8,
        random_state=seed, n_jobs=-1,
    )
    model.fit(X_train, y_train)
    return model.predict(X_test)


# -------------------------------------------------
# Multi-layer NN X+S (8-16-8 architecture, matches Table2-3)
# Dynamic features [wind_speed, temperature] get a per-fold StandardScaler,
# exactly as Table2-Table3(NN).ipynb does; terrain columns are already
# globally standardized (terrain_mat_scaled) and pass through unchanged --
# same split responsibility as Table2-Table3(BNN).ipynb.
# -------------------------------------------------
def fit_nn(X_train, y_train, X_test, seed):
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers, Sequential

    np.random.seed(seed)
    tf.random.set_seed(seed)

    n_dyn = len(FEATURES)
    dyn_scaler = StandardScaler()
    X_train_s = X_train.copy()
    X_test_s = X_test.copy()
    X_train_s[:, :n_dyn] = dyn_scaler.fit_transform(X_train[:, :n_dyn])
    X_test_s[:, :n_dyn] = dyn_scaler.transform(X_test[:, :n_dyn])

    model = Sequential([
        layers.Dense(8, activation="relu"),
        layers.Dense(16, activation="relu"),
        layers.Dense(8, activation="relu"),
        layers.Dense(1),
    ])
    model.compile(
        loss=keras.losses.MeanSquaredError(),
        optimizer=keras.optimizers.Adam(learning_rate=0.001),
        metrics=["mae", "mse"],
    )
    model.fit(X_train_s, y_train, epochs=100, batch_size=2048, verbose=0, shuffle=True)
    pred = model.predict(X_test_s, batch_size=2048, verbose=0).reshape(-1)
    return np.clip(pred, 0.0, None)


# -------------------------------------------------
# Bayesian NN -- X ONLY.
#
# Table2-Table3(BNN).ipynb always concatenates terrain (build_xy_with_terrain
# has no X-only branch), so the existing "BNN" row in the main results
# (Table 2/3) is already the X+S version. Here we build the X-only variant
# by dropping terrain from that same procedure.
# -------------------------------------------------
def run_bnn_loo_xonly():
    import random
    import torch
    import torch.nn as nn
    import torch.optim as optim
    import torchbnn as bnn
    from torch.utils.data import TensorDataset, DataLoader

    SEED = 15
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    BATCH_SIZE, EPOCHS, LR = 4096, 6, 1e-3
    H, PRIOR_MU, PRIOR_SIGMA, KL_WEIGHT, MC_DRAWS = 128, 0.0, 0.1, 0.01, 40

    def build_xy(df):
        d = df[FEATURES + [TARGET]].dropna().copy()
        X = d[FEATURES].to_numpy(np.float32)
        y = d[TARGET].to_numpy(np.float32)
        return X, y

    def seed_all(seed=SEED):
        random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)

    def nlpd_gaussian(y, mu, sd, sd_floor=1e-6):
        sd = np.maximum(sd, sd_floor)
        return float(np.mean(0.5 * np.log(2.0 * np.pi * sd**2) + 0.5 * ((y - mu) ** 2) / (sd ** 2)))

    def make_model(d_in):
        return nn.Sequential(
            bnn.BayesLinear(prior_mu=PRIOR_MU, prior_sigma=PRIOR_SIGMA, in_features=d_in, out_features=H),
            nn.ReLU(),
            bnn.BayesLinear(prior_mu=PRIOR_MU, prior_sigma=PRIOR_SIGMA, in_features=H, out_features=H),
            nn.ReLU(),
            bnn.BayesLinear(prior_mu=PRIOR_MU, prior_sigma=PRIOR_SIGMA, in_features=H, out_features=1),
        ).to(DEVICE)

    mse_loss = nn.MSELoss()
    kl_loss = bnn.BKLLoss(reduction="mean", last_layer_only=False)

    def train_bnn(model, X_train, y_train):
        X_tensor = torch.tensor(X_train, dtype=torch.float32)
        y_tensor = torch.tensor(y_train, dtype=torch.float32).view(-1, 1)
        loader = DataLoader(TensorDataset(X_tensor, y_tensor), batch_size=BATCH_SIZE,
                             shuffle=True, num_workers=0, drop_last=False)
        optimizer = optim.Adam(model.parameters(), lr=LR)
        model.train()
        for _ in range(EPOCHS):
            for xb, yb in loader:
                xb, yb = xb.to(DEVICE), yb.to(DEVICE)
                pred = model(xb)
                loss = mse_loss(pred, yb) + KL_WEIGHT * kl_loss(model)
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()
        return model

    @torch.no_grad()
    def predict_mc(model, X_test, mc_draws=MC_DRAWS, pred_batch_size=8192):
        model.eval()
        X_tensor = torch.tensor(X_test, dtype=torch.float32)
        n = X_tensor.shape[0]
        draws = []
        for _ in range(mc_draws):
            preds = []
            for start in range(0, n, pred_batch_size):
                xb = X_tensor[start:start + pred_batch_size].to(DEVICE)
                preds.append(model(xb).cpu().numpy().reshape(-1))
            draws.append(np.concatenate(preds))
        draws = np.stack(draws, axis=0)
        return draws.mean(axis=0), draws.std(axis=0)

    seed_all(SEED)
    rows = []

    for i in TURBINE_IDS:
        print(f"[BNN X-only] Turbine {i}")
        test17 = data_cache[(i, 2017)]
        test18 = data_cache[(i, 2018)] if i in TESTSET_2018 else None
        if test17 is None or not set(FEATURES + [TARGET]).issubset(test17.columns):
            continue

        X_train_parts, y_train_parts = [], []
        for j in TURBINE_IDS:
            if j == i:
                continue
            tr = data_cache[(j, 2017)]
            if tr is None or not set(FEATURES + [TARGET]).issubset(tr.columns):
                continue
            Xj, yj = build_xy(tr)
            if len(yj) == 0:
                continue
            X_train_parts.append(Xj)
            y_train_parts.append(yj)
        if not X_train_parts:
            continue

        X_train = np.vstack(X_train_parts).astype(np.float32)
        y_train = np.concatenate(y_train_parts).astype(np.float32)

        x_scaler = StandardScaler()
        X_train_s = x_scaler.fit_transform(X_train).astype(np.float32)

        y_mean, y_std = float(np.mean(y_train)), float(np.std(y_train))
        if not np.isfinite(y_std) or y_std < 1e-8:
            y_std = 1.0
        y_train_s = ((y_train - y_mean) / y_std).astype(np.float32)

        model = make_model(d_in=X_train_s.shape[1])
        model = train_bnn(model, X_train_s, y_train_s)

        X_te17, y_te17 = build_xy(test17)
        X_te17_s = x_scaler.transform(X_te17).astype(np.float32)
        pred17_s, pred17_sd_s = predict_mc(model, X_te17_s)
        pred17 = np.clip(pred17_s * y_std + y_mean, 0.0, None)
        pred17_sd = pred17_sd_s * y_std
        rmse17 = rmse(pred17, y_te17)
        nlpd17 = nlpd_gaussian(y_te17, pred17, pred17_sd)

        rmse18, nlpd18 = np.nan, np.nan
        if test18 is not None and set(FEATURES + [TARGET]).issubset(test18.columns):
            X_te18, y_te18 = build_xy(test18)
            if len(y_te18) > 0:
                X_te18_s = x_scaler.transform(X_te18).astype(np.float32)
                pred18_s, pred18_sd_s = predict_mc(model, X_te18_s)
                pred18 = np.clip(pred18_s * y_std + y_mean, 0.0, None)
                pred18_sd = pred18_sd_s * y_std
                rmse18 = rmse(pred18, y_te18)
                nlpd18 = nlpd_gaussian(y_te18, pred18, pred18_sd)

        rows.append({"Method": "Bayesian NN", "Version": "X", "Turbine": i,
                      "RMSE_2017": rmse17, "NLPD_2017": nlpd17,
                      "RMSE_2018": rmse18, "NLPD_2018": nlpd18})

        del model
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    return pd.DataFrame(rows)


# -------------------------------------------------
# Run everything, save long results, update table_s1.csv directly
# -------------------------------------------------
if __name__ == "__main__":
    all_long = []

    for method_name, fn, terrain_mat in [
        ("XGBoost", fit_xgb, terrain_mat_raw),
        ("Multi-layer NN", fit_nn, terrain_mat_scaled),
    ]:
        long_df = run_loo(method_name, fn, terrain_mat)
        all_long.append(long_df)
        update_table_s1(
            method_name, "X+S",
            rmse_2017=long_df.loc[long_df["Year"] == 2017, "RMSE"].mean(),
            rmse_2018=long_df.loc[long_df["Year"] == 2018, "RMSE"].mean(),
        )

    pd.concat(all_long, ignore_index=True).to_csv(OUT_LONG, index=False)
    print("Saved:", OUT_LONG)

    # Bayesian NN: X-only (X+S already lives in results/final_results.csv)
    bnn_df = run_bnn_loo_xonly()
    bnn_df.to_csv(RESULTS_DIR / "s1_bnn_X_long.csv", index=False)
    update_table_s1(
        "Bayesian NN", "X",
        rmse_2017=bnn_df["RMSE_2017"].mean(),
        rmse_2018=bnn_df["RMSE_2018"].mean(),
        nlpd_2017=bnn_df["NLPD_2017"].mean(),
        nlpd_2018=bnn_df["NLPD_2018"].mean(),
    )

    print("Done.")
