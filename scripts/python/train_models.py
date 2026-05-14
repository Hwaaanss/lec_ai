#!/usr/bin/env python
from __future__ import annotations

import argparse
import itertools
import json
import os
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MPLCONFIG_DIR = PROJECT_ROOT / ".cache" / "matplotlib"
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap
import yaml
from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    roc_auc_score,
)
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.svm import SVC
from sklearn.exceptions import ConvergenceWarning
from xgboost import XGBClassifier

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="Tight layout not applied.*", category=UserWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    return parser.parse_args()


def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def build_preprocessor(df: pd.DataFrame, drop_cols: Iterable[str]) -> Tuple[ColumnTransformer, List[str], List[str]]:
    feature_cols = [c for c in df.columns if c not in set(drop_cols)]
    numeric_cols = [c for c in feature_cols if pd.api.types.is_numeric_dtype(df[c])]
    categorical_cols = [c for c in feature_cols if c not in numeric_cols]

    numeric_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )
    categorical_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", OneHotEncoder(handle_unknown="ignore")),
        ]
    )

    preprocessor = ColumnTransformer(
        transformers=[
            ("num", numeric_pipe, numeric_cols),
            ("cat", categorical_pipe, categorical_cols),
        ]
    )
    return preprocessor, numeric_cols, categorical_cols


def model_specs(random_state: int) -> Dict[str, Tuple[object, List[dict]]]:
    return {
        "ridge_logreg": (
            LogisticRegression(
                penalty="l2",
                solver="liblinear",
                class_weight="balanced",
                max_iter=5000,
                random_state=random_state,
            ),
            [{"model__C": [0.01, 0.1, 1.0, 10.0, 50.0]}],
        ),
        "elastic_net_logreg": (
            LogisticRegression(
                penalty="elasticnet",
                solver="saga",
                class_weight="balanced",
                max_iter=20000,
                tol=1e-3,
                random_state=random_state,
            ),
            [
                {
                    "model__C": [0.01, 0.1, 1.0, 10.0],
                    "model__l1_ratio": [0.1, 0.3, 0.5, 0.7, 0.9],
                }
            ],
        ),
        "linear_svm": (
            SVC(
                kernel="linear",
                probability=True,
                class_weight="balanced",
                random_state=random_state,
            ),
            [{"model__C": [0.01, 0.1, 1.0, 10.0, 50.0]}],
        ),
        "xgboost": (
            XGBClassifier(
                objective="binary:logistic",
                eval_metric="logloss",
                random_state=random_state,
                scale_pos_weight=1.0,
                tree_method="hist",
                n_jobs=1,
            ),
            [
                {
                    "model__max_depth": [2, 4],
                    "model__learning_rate": [0.05],
                    "model__n_estimators": [120],
                    "model__subsample": [0.8],
                    "model__colsample_bytree": [0.8],
                }
            ],
        ),
    }


def iter_param_grid(grid_list: List[dict]) -> Iterable[dict]:
    from itertools import product

    for grid in grid_list:
        keys = list(grid.keys())
        vals = [grid[k] for k in keys]
        for combo in product(*vals):
            yield dict(zip(keys, combo))


def optimal_threshold(y_true: np.ndarray, y_prob: np.ndarray) -> float:
    thresholds = np.unique(np.clip(y_prob, 0.0, 1.0))
    best_threshold = 0.5
    best_score = -np.inf
    for threshold in thresholds:
        pred = (y_prob >= threshold).astype(int)
        score = balanced_accuracy_score(y_true, pred)
        if score > best_score:
            best_score = score
            best_threshold = float(threshold)
    return best_threshold


def evaluate_probs(y_true: np.ndarray, y_prob: np.ndarray, threshold: float) -> dict:
    y_pred = (y_prob >= threshold).astype(int)
    return {
        "auroc": float(roc_auc_score(y_true, y_prob)),
        "auprc": float(average_precision_score(y_true, y_prob)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "threshold": float(threshold),
    }


def select_model_with_inner_cv(
    X: pd.DataFrame,
    y: np.ndarray,
    strata: np.ndarray,
    base_pipeline: Pipeline,
    param_grid_list: List[dict],
    inner_folds: int,
    inner_repeats: int,
    random_state: int,
) -> Tuple[Pipeline, dict, float]:
    cv = RepeatedStratifiedKFold(n_splits=inner_folds, n_repeats=inner_repeats, random_state=random_state)
    best_score = -np.inf
    best_params = None
    best_threshold = 0.5

    for params in iter_param_grid(param_grid_list):
        fold_scores = []
        pooled_y = []
        pooled_probs = []
        for train_idx, valid_idx in cv.split(X, strata):
            model = clone(base_pipeline)
            model.set_params(**params)
            model.fit(X.iloc[train_idx], y[train_idx])
            probs = model.predict_proba(X.iloc[valid_idx])[:, 1]
            pooled_y.append(y[valid_idx])
            pooled_probs.append(probs)
            fold_scores.append(roc_auc_score(y[valid_idx], probs))
        mean_score = float(np.mean(fold_scores))
        if mean_score > best_score:
            best_score = mean_score
            best_params = params
            best_threshold = optimal_threshold(np.concatenate(pooled_y), np.concatenate(pooled_probs))

    final_model = clone(base_pipeline)
    final_model.set_params(**best_params)
    final_model.fit(X, y)
    return final_model, best_params, best_threshold


def transformed_feature_names(fitted_pipeline: Pipeline, numeric_cols: List[str], categorical_cols: List[str]) -> List[str]:
    pre = fitted_pipeline.named_steps["pre"]
    feature_names = list(numeric_cols)
    if categorical_cols:
        ohe = pre.named_transformers_["cat"].named_steps["onehot"]
        feature_names.extend(ohe.get_feature_names_out(categorical_cols).tolist())
    return feature_names


def summarize_linear_importance(model: Pipeline, feature_names: List[str]) -> pd.DataFrame:
    coef = model.named_steps["model"].coef_
    if hasattr(coef, "toarray"):
        coef = coef.toarray()
    coef = np.asarray(coef).ravel()
    out = pd.DataFrame({"feature": feature_names, "coefficient": coef})
    out["abs_coefficient"] = out["coefficient"].abs()
    return out.sort_values("abs_coefficient", ascending=False)


def summarize_xgb_importance(model: Pipeline, feature_names: List[str], X: pd.DataFrame, output_prefix: Path) -> pd.DataFrame:
    transformed = model.named_steps["pre"].transform(X)
    booster = model.named_steps["model"]
    importance = pd.DataFrame(
        {
            "feature": feature_names,
            "gain": booster.feature_importances_,
        }
    ).sort_values("gain", ascending=False)

    explainer = shap.TreeExplainer(booster)
    shap_values = explainer.shap_values(transformed)
    plt.figure(figsize=(13, 8))
    shap.summary_plot(shap_values, transformed, feature_names=feature_names, show=False, max_display=20, plot_size=(13, 8))
    plt.gcf().subplots_adjust(left=0.34, right=0.98, top=0.95, bottom=0.10)
    plt.savefig(output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.2)
    plt.close()
    return importance


def require_torch_mps():
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    if not torch.backends.mps.is_available():
        raise RuntimeError("PASNet requires torch MPS, but torch.backends.mps.is_available() is False.")
    return torch, nn, F, torch.device("mps")


def load_pasnet_inputs(processed_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.read_csv(processed_dir / "pasnet_expression_primary.tsv", sep="\t")
    metadata = pd.read_csv(processed_dir / "metadata_primary.tsv", sep="\t")
    gene_sets = pd.read_csv(processed_dir / "pasnet_gene_sets.tsv", sep="\t")
    df = expression.merge(
        metadata[["geo_accession", "response_primary", "response_binary", "source"]],
        on="geo_accession",
        how="inner",
    )
    if df.empty:
        raise RuntimeError("PASNet input merge produced no rows.")
    return df, gene_sets


def pasnet_param_grid(round_name: str) -> List[dict]:
    base = {
        "dropout": [0.1, 0.3, 0.5],
        "hidden_dim": [16, 32, 64],
        "lr": [1e-3, 3e-4],
        "weight_decay": [1e-4, 1e-3],
        "top_genes": [5000],
        "min_genes": [5],
        "pathway_mode": ["trainable", "mean"],
    }
    if round_name == "extended":
        base["top_genes"] = [3000, 5000, None]
        base["min_genes"] = [5, 10]
        base["dropout"] = [0.3, 0.5, 0.6]
        base["weight_decay"] = [1e-3, 3e-3]
    keys = list(base)
    return [dict(zip(keys, values)) for values in itertools.product(*(base[k] for k in keys))]


def prepare_pasnet_arrays(
    train_df: pd.DataFrame,
    apply_df: pd.DataFrame,
    gene_sets: pd.DataFrame,
    top_genes: int | None,
    min_genes: int,
) -> Tuple[np.ndarray, np.ndarray, List[str], List[str]]:
    drop_cols = {"geo_accession", "response_primary", "response_binary", "source"}
    all_gene_cols = [c for c in train_df.columns if c not in drop_cols]
    variances = train_df[all_gene_cols].var(axis=0).sort_values(ascending=False)
    selected_genes = variances.index.tolist() if top_genes is None else variances.head(top_genes).index.tolist()

    selected_set = set(selected_genes)
    gs = gene_sets.loc[gene_sets["gene_symbol"].isin(selected_set)].copy()
    pathway_sizes = gs.groupby("pathway")["gene_symbol"].nunique()
    keep_pathways = pathway_sizes[pathway_sizes >= min_genes].index.tolist()
    if not keep_pathways:
        raise RuntimeError("PASNet gene-set filtering removed all pathways.")

    gs = gs.loc[gs["pathway"].isin(keep_pathways)]
    used_genes = sorted(gs["gene_symbol"].unique().tolist())
    pathways = sorted(keep_pathways)
    gene_index = {gene: idx for idx, gene in enumerate(used_genes)}
    pathway_index = {pathway: idx for idx, pathway in enumerate(pathways)}
    mask = np.zeros((len(pathways), len(used_genes)), dtype=np.float32)
    for row in gs.itertuples(index=False):
        mask[pathway_index[row.pathway], gene_index[row.gene_symbol]] = 1.0

    train_x = train_df[used_genes].to_numpy(dtype=np.float32)
    apply_x = apply_df[used_genes].to_numpy(dtype=np.float32)
    mean = train_x.mean(axis=0, keepdims=True)
    std = train_x.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    apply_x = (apply_x - mean) / std
    return apply_x.astype(np.float32), mask, used_genes, pathways


def train_pasnet_model(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_valid: np.ndarray,
    y_valid: np.ndarray,
    mask: np.ndarray,
    params: dict,
    seed: int,
) -> Tuple[object, float]:
    torch, nn, F, device = require_torch_mps()
    torch.manual_seed(seed)
    np.random.seed(seed)

    class MaskedLinear(nn.Module):
        def __init__(self, mask_array: np.ndarray, mode: str):
            super().__init__()
            self.mode = mode
            self.bias = nn.Parameter(torch.zeros(mask_array.shape[0]))
            if mode == "mean":
                denom = np.maximum(mask_array.sum(axis=1, keepdims=True), 1.0)
                normalized_mask = mask_array / denom
                self.register_buffer("mask", torch.tensor(normalized_mask, dtype=torch.float32))
                self.scale = nn.Parameter(torch.ones(mask_array.shape[0]))
            else:
                self.register_buffer("mask", torch.tensor(mask_array, dtype=torch.float32))
                self.weight = nn.Parameter(torch.empty(mask_array.shape[0], mask_array.shape[1]))
                nn.init.xavier_uniform_(self.weight)

        def forward(self, x):
            if self.mode == "mean":
                return F.linear(x, self.mask, None) * self.scale + self.bias
            return F.linear(x, self.weight * self.mask, self.bias)

    class PASNet(nn.Module):
        def __init__(self, mask_array: np.ndarray, hidden_dim: int, dropout: float):
            super().__init__()
            self.gene_to_pathway = MaskedLinear(mask_array, params.get("pathway_mode", "trainable"))
            self.pathway_bn = nn.BatchNorm1d(mask_array.shape[0])
            self.dropout = nn.Dropout(dropout)
            self.hidden = nn.Linear(mask_array.shape[0], hidden_dim)
            self.hidden_bn = nn.BatchNorm1d(hidden_dim)
            self.out = nn.Linear(hidden_dim, 1)

        def forward(self, x):
            x = self.gene_to_pathway(x)
            x = self.pathway_bn(x)
            x = F.relu(x)
            x = self.dropout(x)
            x = self.hidden(x)
            x = self.hidden_bn(x)
            x = F.relu(x)
            x = self.dropout(x)
            return self.out(x).squeeze(-1)

    model = PASNet(mask, params["hidden_dim"], params["dropout"]).to(device)
    x_train = torch.tensor(X_train, dtype=torch.float32, device=device)
    y_train_t = torch.tensor(y_train.astype(np.float32), dtype=torch.float32, device=device)
    x_valid = torch.tensor(X_valid, dtype=torch.float32, device=device)

    pos = max(float((y_train == 1).sum()), 1.0)
    neg = max(float((y_train == 0).sum()), 1.0)
    criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(neg / pos, dtype=torch.float32, device=device))
    optimizer = torch.optim.AdamW(model.parameters(), lr=params["lr"], weight_decay=params["weight_decay"])

    best_state = None
    best_score = -np.inf
    stale_epochs = 0
    max_epochs = int(params.get("epochs", 500))
    patience = int(params.get("patience", 60))

    for _ in range(max_epochs):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        loss = criterion(model(x_train), y_train_t)
        loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            valid_prob = torch.sigmoid(model(x_valid)).cpu().numpy()
        score = roc_auc_score(y_valid, valid_prob)
        if score > best_score + 1e-5:
            best_score = float(score)
            best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            stale_epochs = 0
        else:
            stale_epochs += 1
        if stale_epochs >= patience:
            break

    if best_state is not None:
        model.load_state_dict({key: value.to(device) for key, value in best_state.items()})
    return model, best_score


def predict_pasnet(model, X: np.ndarray) -> np.ndarray:
    torch, _, _, device = require_torch_mps()
    model.eval()
    with torch.no_grad():
        probs = torch.sigmoid(model(torch.tensor(X, dtype=torch.float32, device=device))).cpu().numpy()
    return probs.astype(float)


def fit_pasnet_with_inner_cv(
    train_df: pd.DataFrame,
    y_train: np.ndarray,
    strata_train: np.ndarray,
    gene_sets: pd.DataFrame,
    params_grid: List[dict],
    inner_folds: int,
    random_state: int,
) -> Tuple[object, dict, float, np.ndarray, List[str], List[str]]:
    cv = RepeatedStratifiedKFold(n_splits=inner_folds, n_repeats=1, random_state=random_state)
    best_score = -np.inf
    best_params = None
    best_threshold = 0.5

    for params in params_grid:
        fold_scores = []
        pooled_y = []
        pooled_probs = []
        for fold_seed, (inner_train_idx, inner_valid_idx) in enumerate(cv.split(train_df, strata_train), start=1):
            inner_train = train_df.iloc[inner_train_idx]
            inner_valid = train_df.iloc[inner_valid_idx]
            X_inner_train, mask, _, _ = prepare_pasnet_arrays(
                inner_train, inner_train, gene_sets, params["top_genes"], params["min_genes"]
            )
            X_inner_valid, _, _, _ = prepare_pasnet_arrays(
                inner_train, inner_valid, gene_sets, params["top_genes"], params["min_genes"]
            )
            model, score = train_pasnet_model(
                X_inner_train,
                y_train[inner_train_idx],
                X_inner_valid,
                y_train[inner_valid_idx],
                mask,
                {**params, "epochs": 500, "patience": 60},
                seed=random_state + fold_seed,
            )
            probs = predict_pasnet(model, X_inner_valid)
            pooled_y.append(y_train[inner_valid_idx])
            pooled_probs.append(probs)
            fold_scores.append(score)
        mean_score = float(np.mean(fold_scores))
        if mean_score > best_score:
            best_score = mean_score
            best_params = params
            best_threshold = optimal_threshold(np.concatenate(pooled_y), np.concatenate(pooled_probs))

    X_final, mask, used_genes, pathways = prepare_pasnet_arrays(train_df, train_df, gene_sets, best_params["top_genes"], best_params["min_genes"])
    model, _ = train_pasnet_model(
        X_final,
        y_train,
        X_final,
        y_train,
        mask,
        {**best_params, "epochs": 500, "patience": 60},
        seed=random_state,
    )
    return model, best_params, best_threshold, mask, used_genes, pathways


def summarize_pasnet_importance(model, pathways: List[str], repeat_idx: int, fold_idx: int) -> pd.DataFrame:
    pathway_weight = model.hidden.weight.detach().cpu().numpy()
    out_weight = model.out.weight.detach().cpu().numpy().ravel()
    importance = np.abs(pathway_weight.T @ out_weight)
    out = pd.DataFrame({"feature": pathways, "gain": importance})
    out["feature_set"] = "gene_pathway"
    out["model"] = "pasnet"
    out["repeat"] = repeat_idx
    out["fold"] = fold_idx
    return out.sort_values("gain", ascending=False)


def run_pasnet_nested_cv(
    df: pd.DataFrame,
    gene_sets: pd.DataFrame,
    baseline_best_auroc: float,
    outer_folds: int,
    outer_repeats: int,
    inner_folds: int,
    random_state: int,
    tables_dir: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y = df["response_binary"].astype(int).to_numpy()
    strata = (df["source"].astype(str) + "__" + df["response_primary"].astype(str)).to_numpy()
    outer_cv = RepeatedStratifiedKFold(n_splits=outer_folds, n_repeats=outer_repeats, random_state=random_state)
    fold_rows = []
    importance_frames = []

    params_grid = pasnet_param_grid("initial")
    for split_idx, (train_idx, test_idx) in enumerate(outer_cv.split(df, strata), start=1):
        repeat_idx = (split_idx - 1) // outer_folds + 1
        fold_idx = (split_idx - 1) % outer_folds + 1
        print(f"[modeling] start feature_set=gene_pathway model=pasnet repeat={repeat_idx} fold={fold_idx}", flush=True)
        train_df = df.iloc[train_idx].copy()
        test_df = df.iloc[test_idx].copy()
        y_train, y_test = y[train_idx], y[test_idx]
        model, best_params, threshold, _, _, pathways = fit_pasnet_with_inner_cv(
            train_df,
            y_train,
            strata[train_idx],
            gene_sets,
            params_grid,
            inner_folds=inner_folds,
            random_state=random_state + split_idx,
        )
        X_test, _, _, _ = prepare_pasnet_arrays(train_df, test_df, gene_sets, best_params["top_genes"], best_params["min_genes"])
        test_probs = predict_pasnet(model, X_test)
        metrics = evaluate_probs(y_test, test_probs, threshold)
        fold_rows.append(
            {
                "feature_set": "gene_pathway",
                "model": "pasnet",
                "repeat": repeat_idx,
                "fold": fold_idx,
                **metrics,
                "best_params": json.dumps(best_params, ensure_ascii=False),
            }
        )
        importance_frames.append(summarize_pasnet_importance(model, pathways, repeat_idx, fold_idx))
        print(f"[modeling] done feature_set=gene_pathway model=pasnet repeat={repeat_idx} fold={fold_idx}", flush=True)

    fold_df = pd.DataFrame(fold_rows)
    importance_df = pd.concat(importance_frames, ignore_index=True)
    mean_auroc = float(fold_df["auroc"].mean())
    target_auroc = min(0.76, baseline_best_auroc - 0.02)
    tuning_note = (
        "PASNet reached the AUROC target."
        if mean_auroc >= target_auroc
        else "PASNet did not reach the baseline AUROC target after the fixed non-leaky tuning grid."
    )
    tuning_df = pd.DataFrame(
        [
            {
                "model": "pasnet",
                "feature_set": "gene_pathway",
                "auroc_mean": mean_auroc,
                "baseline_best_auroc": baseline_best_auroc,
                "target_auroc": target_auroc,
                "note": tuning_note,
            }
        ]
    )

    fold_df.to_csv(tables_dir / "nested_cv_gene_pathway_pasnet.tsv", sep="\t", index=False)
    importance_df.to_csv(tables_dir / "feature_importance_gene_pathway_pasnet.tsv", sep="\t", index=False)
    tuning_df.to_csv(tables_dir / "pasnet_tuning_summary.tsv", sep="\t", index=False)
    return fold_df, importance_df, tuning_df


def run_pasnet_source_transfer(
    df: pd.DataFrame,
    gene_sets: pd.DataFrame,
    inner_folds: int,
    random_state: int,
) -> List[dict]:
    rows = []
    params_grid = pasnet_param_grid("initial")
    for train_source, test_source in [("MDACC", "ISPY"), ("ISPY", "MDACC")]:
        train_df = df.loc[df["source"] == train_source].copy()
        test_df = df.loc[df["source"] == test_source].copy()
        if train_df.empty or test_df.empty:
            continue
        y_train = train_df["response_binary"].astype(int).to_numpy()
        y_test = test_df["response_binary"].astype(int).to_numpy()
        strata_train = (train_df["source"].astype(str) + "__" + train_df["response_primary"].astype(str)).to_numpy()
        model, best_params, threshold, _, _, _ = fit_pasnet_with_inner_cv(
            train_df,
            y_train,
            strata_train,
            gene_sets,
            params_grid,
            inner_folds=inner_folds,
            random_state=random_state,
        )
        X_test, _, _, _ = prepare_pasnet_arrays(train_df, test_df, gene_sets, best_params["top_genes"], best_params["min_genes"])
        probs = predict_pasnet(model, X_test)
        metrics = evaluate_probs(y_test, probs, threshold)
        rows.append(
            {
                "feature_set": "gene_pathway",
                "model": "pasnet",
                "train_source": train_source,
                "test_source": test_source,
                **metrics,
                "best_params": json.dumps(best_params, ensure_ascii=False),
            }
        )
    return rows


def run_nested_cv(
    df: pd.DataFrame,
    feature_set_name: str,
    model_name: str,
    model,
    param_grid_list: List[dict],
    outer_folds: int,
    outer_repeats: int,
    inner_folds: int,
    inner_repeats: int,
    random_state: int,
    tables_dir: Path,
    models_dir: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = df.copy()
    y = df["response_binary"].astype(int).to_numpy()
    strata = (df["source"].astype(str) + "__" + df["response_primary"].astype(str)).to_numpy()
    drop_cols = ["geo_accession", "source", "response_primary", "response_binary"]
    X = df.drop(columns=drop_cols)
    preprocessor, numeric_cols, categorical_cols = build_preprocessor(df, drop_cols=drop_cols)

    outer_cv = RepeatedStratifiedKFold(n_splits=outer_folds, n_repeats=outer_repeats, random_state=random_state)
    fold_rows: List[dict] = []
    importance_frames: List[pd.DataFrame] = []

    for split_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, strata), start=1):
        repeat_idx = (split_idx - 1) // outer_folds + 1
        fold_idx = (split_idx - 1) % outer_folds + 1
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        strata_train = strata[train_idx]

        model_for_outer = clone(model)
        if model_name == "xgboost":
            pos = max(float((y_train == 1).sum()), 1.0)
            neg = max(float((y_train == 0).sum()), 1.0)
            model_for_outer.set_params(scale_pos_weight=neg / pos, verbosity=0)
        base_pipeline = Pipeline([("pre", preprocessor), ("model", model_for_outer)])

        fitted_model, best_params, threshold = select_model_with_inner_cv(
            X_train,
            y_train,
            strata_train,
            base_pipeline,
            param_grid_list,
            inner_folds,
            inner_repeats,
            random_state,
        )

        test_probs = fitted_model.predict_proba(X_test)[:, 1]
        metrics = evaluate_probs(y_test, test_probs, threshold)
        fold_rows.append(
            {
                "feature_set": feature_set_name,
                "model": model_name,
                "repeat": repeat_idx,
                "fold": fold_idx,
                **metrics,
                "best_params": json.dumps(best_params, ensure_ascii=False),
            }
        )

        feature_names = transformed_feature_names(fitted_model, numeric_cols, categorical_cols)
        if model_name in {"ridge_logreg", "elastic_net_logreg"}:
            importance = summarize_linear_importance(fitted_model, feature_names)
        elif model_name == "linear_svm":
            coef = fitted_model.named_steps["model"].coef_
            if hasattr(coef, "toarray"):
                coef = coef.toarray()
            coef = np.asarray(coef).ravel()
            importance = pd.DataFrame({"feature": feature_names, "coefficient": coef})
            importance["abs_coefficient"] = importance["coefficient"].abs()
            importance = importance.sort_values("abs_coefficient", ascending=False)
        else:
            shap_path = models_dir / f"shap_{feature_set_name}_{model_name}_repeat{repeat_idx}_fold{fold_idx}"
            importance = summarize_xgb_importance(fitted_model, feature_names, X_train, shap_path)
        importance["feature_set"] = feature_set_name
        importance["model"] = model_name
        importance["repeat"] = repeat_idx
        importance["fold"] = fold_idx
        importance_frames.append(importance)

    fold_df = pd.DataFrame(fold_rows)
    importance_df = pd.concat(importance_frames, ignore_index=True)

    fold_df.to_csv(tables_dir / f"nested_cv_{feature_set_name}_{model_name}.tsv", sep="\t", index=False)
    importance_df.to_csv(tables_dir / f"feature_importance_{feature_set_name}_{model_name}.tsv", sep="\t", index=False)

    return fold_df, importance_df


def run_source_transfer(
    df: pd.DataFrame,
    feature_set_name: str,
    model_name: str,
    model,
    param_grid_list: List[dict],
    inner_folds: int,
    inner_repeats: int,
    random_state: int,
) -> List[dict]:
    rows = []
    for train_source, test_source in [("MDACC", "ISPY"), ("ISPY", "MDACC")]:
        train_df = df.loc[df["source"] == train_source].copy()
        test_df = df.loc[df["source"] == test_source].copy()
        if train_df.empty or test_df.empty:
            continue

        y_train = train_df["response_binary"].astype(int).to_numpy()
        strata_train = (train_df["source"].astype(str) + "__" + train_df["response_primary"].astype(str)).to_numpy()
        drop_cols = ["geo_accession", "source", "response_primary", "response_binary"]
        preprocessor, _, _ = build_preprocessor(train_df, drop_cols=drop_cols)
        model_for_transfer = clone(model)
        if model_name == "xgboost":
            pos = max(float((y_train == 1).sum()), 1.0)
            neg = max(float((y_train == 0).sum()), 1.0)
            model_for_transfer.set_params(scale_pos_weight=neg / pos, verbosity=0)
        base_pipeline = Pipeline([("pre", preprocessor), ("model", model_for_transfer)])
        fitted_model, best_params, threshold = select_model_with_inner_cv(
            train_df.drop(columns=drop_cols),
            y_train,
            strata_train,
            base_pipeline,
            param_grid_list,
            inner_folds=inner_folds,
            inner_repeats=inner_repeats,
            random_state=random_state,
        )

        test_probs = fitted_model.predict_proba(test_df.drop(columns=drop_cols))[:, 1]
        metrics = evaluate_probs(test_df["response_binary"].astype(int).to_numpy(), test_probs, threshold)
        rows.append(
            {
                "feature_set": feature_set_name,
                "model": model_name,
                "train_source": train_source,
                "test_source": test_source,
                **metrics,
                "best_params": json.dumps(best_params, ensure_ascii=False),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    paths = cfg["paths"]
    modeling_cfg = cfg["modeling"]

    tables_dir = Path(paths["tables_dir"])
    plots_dir = Path(paths["plots_dir"]) / "modeling"
    models_dir = Path(paths["models_dir"])
    ensure_dir(tables_dir)
    ensure_dir(plots_dir)
    ensure_dir(models_dir)

    datasets = {
        "clinical_only": pd.read_csv(Path(paths["processed_dir"]) / "model_input_clinical_primary.tsv", sep="\t"),
        "bi_only": pd.read_csv(Path(paths["processed_dir"]) / "model_input_bi_primary.tsv", sep="\t"),
        "clinical_plus_bi": pd.read_csv(Path(paths["processed_dir"]) / "model_input_combined_primary.tsv", sep="\t"),
    }

    all_fold_results = []
    all_transfer_results = []
    top_signature_rows = []
    outer_repeats = modeling_cfg.get("outer_repeats", 1)
    inner_repeats = modeling_cfg.get("inner_repeats", 1)

    for feature_set_name, df in datasets.items():
        specs = model_specs(modeling_cfg["random_state"])
        for model_name, (model, param_grid_list) in specs.items():
            print(f"[modeling] start feature_set={feature_set_name} model={model_name}", flush=True)
            fold_df, importance_df = run_nested_cv(
                df=df,
                feature_set_name=feature_set_name,
                model_name=model_name,
                model=model,
                param_grid_list=param_grid_list,
                outer_folds=modeling_cfg["outer_folds"],
                outer_repeats=outer_repeats,
                inner_folds=modeling_cfg["inner_folds"],
                inner_repeats=inner_repeats,
                random_state=modeling_cfg["random_state"],
                tables_dir=tables_dir,
                models_dir=models_dir,
            )
            all_fold_results.append(fold_df)

            avg_imp = (
                importance_df.groupby("feature", as_index=False)
                .agg(
                    mean_abs_importance=("abs_coefficient" if "abs_coefficient" in importance_df.columns else "gain", "mean")
                )
                .sort_values("mean_abs_importance", ascending=False)
                .head(15)
            )
            avg_imp["feature_set"] = feature_set_name
            avg_imp["model"] = model_name
            top_signature_rows.append(avg_imp)

            transfer_rows = run_source_transfer(
                df=df,
                feature_set_name=feature_set_name,
                model_name=model_name,
                model=model,
                param_grid_list=param_grid_list,
                inner_folds=modeling_cfg["inner_folds"],
                inner_repeats=inner_repeats,
                random_state=modeling_cfg["random_state"],
            )
            all_transfer_results.extend(transfer_rows)
            print(f"[modeling] done feature_set={feature_set_name} model={model_name}", flush=True)

    baseline_perf = pd.concat(all_fold_results, ignore_index=True)
    baseline_best_auroc = float(
        baseline_perf.groupby(["feature_set", "model"])["auroc"].mean().max()
    )
    pasnet_df, pasnet_gene_sets = load_pasnet_inputs(Path(paths["processed_dir"]))
    pasnet_fold_df, pasnet_importance_df, _ = run_pasnet_nested_cv(
        df=pasnet_df,
        gene_sets=pasnet_gene_sets,
        baseline_best_auroc=baseline_best_auroc,
        outer_folds=modeling_cfg["outer_folds"],
        outer_repeats=outer_repeats,
        inner_folds=modeling_cfg["inner_folds"],
        random_state=modeling_cfg["random_state"],
        tables_dir=tables_dir,
    )
    all_fold_results.append(pasnet_fold_df)
    pasnet_avg_imp = (
        pasnet_importance_df.groupby("feature", as_index=False)
        .agg(mean_abs_importance=("gain", "mean"))
        .sort_values("mean_abs_importance", ascending=False)
        .head(15)
    )
    pasnet_avg_imp["feature_set"] = "gene_pathway"
    pasnet_avg_imp["model"] = "pasnet"
    top_signature_rows.append(pasnet_avg_imp)
    all_transfer_results.extend(
        run_pasnet_source_transfer(
            df=pasnet_df,
            gene_sets=pasnet_gene_sets,
            inner_folds=modeling_cfg["inner_folds"],
            random_state=modeling_cfg["random_state"],
        )
    )

    model_perf = pd.concat(all_fold_results, ignore_index=True)
    model_perf.to_csv(tables_dir / "model_performance_nested_cv.tsv", sep="\t", index=False)

    perf_summary = (
        model_perf.groupby(["feature_set", "model"], as_index=False)
        .agg(
            auroc_mean=("auroc", "mean"),
            auroc_sd=("auroc", "std"),
            auprc_mean=("auprc", "mean"),
            auprc_sd=("auprc", "std"),
            balanced_accuracy_mean=("balanced_accuracy", "mean"),
            balanced_accuracy_sd=("balanced_accuracy", "std"),
        )
        .sort_values(["auroc_mean", "auprc_mean"], ascending=False)
    )
    perf_summary.to_csv(tables_dir / "model_performance_summary.tsv", sep="\t", index=False)

    if all_transfer_results:
        transfer_df = pd.DataFrame(all_transfer_results)
        transfer_df.to_csv(tables_dir / "model_performance_source_transfer.tsv", sep="\t", index=False)

    top_signatures = pd.concat(top_signature_rows, ignore_index=True)
    top_signatures.to_csv(tables_dir / "top_predictive_signatures.tsv", sep="\t", index=False)

    plt.figure(figsize=(12, 7))
    sns.barplot(data=perf_summary, x="model", y="auroc_mean", hue="feature_set")
    plt.title("Nested CV AUROC by model and feature set")
    plt.ylabel("Mean AUROC")
    plt.xlabel("")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(plots_dir / "nested_cv_auroc_barplot.png", dpi=300)
    plt.close()

    plt.figure(figsize=(12, 7))
    sns.barplot(data=perf_summary, x="model", y="auprc_mean", hue="feature_set")
    plt.title("Nested CV AUPRC by model and feature set")
    plt.ylabel("Mean AUPRC")
    plt.xlabel("")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(plots_dir / "nested_cv_auprc_barplot.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
