"""Grouped county validation of climate predictors for vegetation cover."""
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit
from sklearn.cluster import KMeans

from inputs import RESULTS, WORK

PENALTIES = np.array([0, .001, .003, .01, .03, .1, .3, 1, 3])
SEED = 20261030


def read_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def scores(train, test):
    mean, sd = train.mean(axis=0), train.std(axis=0, ddof=1)
    z, v = (train - mean) / sd, (test - mean) / sd
    values, vectors = np.linalg.eigh(np.corrcoef(z.T))
    direction = vectors[:, -1] * np.sign(vectors[0, -1])
    return z, v, z @ direction / np.sqrt(values[-1]), v @ direction / np.sqrt(values[-1])


def linear(a, b, y):
    design = np.column_stack([np.ones(len(a)), a])
    return np.column_stack([np.ones(len(b)), b]) @ np.linalg.lstsq(design, y, rcond=None)[0]


def fractional(a, b, y):
    design = np.column_stack([np.ones(len(a)), a])
    test = np.column_stack([np.ones(len(b)), b])

    def loss(beta):
        eta = design @ beta
        return np.mean(np.logaddexp(0, eta) - y * eta)

    def gradient(beta):
        return design.T @ (expit(design @ beta) - y) / len(y)

    def hessian(beta):
        p = expit(design @ beta)
        return design.T @ ((p * (1 - p))[:, None] * design) / len(y)

    start = np.zeros(design.shape[1])
    start[0] = logit(y.mean())
    fit = minimize(loss, start, jac=gradient, hess=hessian, method="trust-exact",
                   options={"gtol": 1e-10, "maxiter": 200})
    if max(abs(gradient(fit.x))) >= 1e-8:
        raise ValueError(fit.message)
    return expit(test @ fit.x)


def ridge(a, b, y, penalty):
    centred = a - a.mean()
    slope = centred @ (y - y.mean()) / (centred @ centred + len(a) * penalty)
    return y.mean() + (b - a.mean()) * slope


def choose_penalty(x, y, codes, xy):
    groups, idx = np.unique(codes, return_inverse=True)
    centres = np.array([xy[codes == g].mean(0) for g in groups])
    labels = KMeans(n_clusters=3, n_init=30, random_state=SEED).fit_predict(centres)[idx]
    errors = np.zeros(len(PENALTIES))
    for fold in range(3):
        tr, te = labels != fold, labels == fold
        _, _, a, b = scores(x[tr], x[te])
        for j, penalty in enumerate(PENALTIES):
            errors[j] += np.sum((ridge(a, b, y[tr], penalty) - y[te]) ** 2)
    errors /= len(y)
    return PENALTIES[np.flatnonzero(errors <= errors.min() + 1e-14)[0]]


def run(zonal_table=None):
    rows = [r for r in read_rows(RESULTS / "counties.csv") if r["included"] == "1"]
    if zonal_table:
        fresh = {int(r["ZoneFID"]): r for r in read_rows(zonal_table) if r["included"] == "1"}
        if fresh.keys() != {int(r["ZoneFID"]) for r in rows}:
            raise ValueError("County join keys differ")
        for row in rows:
            new = fresh[int(row["ZoneFID"])]
            if int(new["ADCODE99"]) != int(row["ADCODE99"]):
                raise ValueError("County codes differ")
            row.update({key: new[key] for key in ("jiangshui", "zhengfa", "percent")})
    # The original regression import stored 12 significant digits.
    def column(name):
        return np.array([float(format(float(r[name]), ".12g")) for r in rows])
    x = np.column_stack([column("jiangshui"), column("zhengfa")])
    y = column("percent")
    xy = np.array([[float(r["x"]), float(r["y"])] for r in rows])
    codes = np.array([int(r["ADCODE99"]) for r in rows])
    groups, idx = np.unique(codes, return_inverse=True)
    schemes = {f"spatial{k}": np.array([int(r[f"fold{k}"]) for r in rows]) for k in (4, 5, 8)}
    for seed in range(SEED, SEED + 10):
        order = np.random.default_rng(seed).permutation(len(groups))
        assigned = np.empty(len(groups), int)
        assigned[order] = np.arange(len(groups)) % 5
        schemes[f"group_random5_seed{seed}"] = assigned[idx]

    output = WORK / "tables" / "county_predictions.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["scheme", "fold", "ZoneFID", "ADCODE99", "observed", "model", "predicted"])
        for scheme, labels in schemes.items():
            for fold in np.unique(labels):
                tr, te = labels != fold, labels == fold
                assert not set(codes[tr]) & set(codes[te])
                z, v, a, b = scores(x[tr], x[te])
                keep_one = (1 + abs(np.corrcoef(z.T)[0, 1])) / 2 >= .85
                ols = linear(z, v, y[tr])
                forecasts = {"OLS_two": ols, "PCR_85": linear(a, b, y[tr]) if keep_one else ols,
                             "P_only": linear(z[:, 0], v[:, 0], y[tr]),
                             "E_only": linear(z[:, 1], v[:, 1], y[tr]),
                             "training_mean": np.full(te.sum(), y[tr].mean()),
                             "PC1_fractional_logit": fractional(a, b, y[tr])}
                if scheme.startswith("spatial"):
                    penalty = choose_penalty(x[tr], y[tr], codes[tr], xy[tr])
                    forecasts["PC1_ridge"] = ridge(a, b, y[tr], penalty)
                for model, values in forecasts.items():
                    if scheme == "spatial5" and model in ("OLS_two", "PCR_85"):
                        field = {"OLS_two": "ols5", "PCR_85": "pcr5"}[model]
                        expected = [float(rows[i][field]) for i in np.flatnonzero(te)]
                        np.testing.assert_allclose(values, expected, rtol=0, atol=2e-8)
                    for i, value in zip(np.flatnonzero(te), values):
                        writer.writerow([scheme, fold, rows[i]["ZoneFID"], codes[i], y[i], model, value])
