"""The three latent heads: cleavage (MLP), TAP (ridge), MHC (projection + centroid).

Checkpoints are plain numpy/torch-tensor arrays plus JSON metadata -- no
pickled sklearn/torch objects -- so a checkpoint directory stays
inspectable and diffable, per the project's "prefer structured outputs"
rule.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import torch
from torch import nn

from re_agent.e2e_pls.schema import EMBEDDING_DIM, ESM3_MODEL_ID

MODEL_VERSION = "e2e_pls-heads-v1"


# --------------------------------------------------------------------------
# Shared calibration primitive: isotonic regression, stored as breakpoints
# so inference is a plain np.interp (no sklearn object in the checkpoint).
# --------------------------------------------------------------------------


@dataclass
class IsotonicCalibrator:
    x_thresholds: list[float]
    y_thresholds: list[float]

    @classmethod
    def fit(cls, x: np.ndarray, y: np.ndarray) -> IsotonicCalibrator:
        from sklearn.isotonic import IsotonicRegression

        ir = IsotonicRegression(out_of_bounds="clip")
        ir.fit(x, y)
        return cls(x_thresholds=ir.X_thresholds_.tolist(), y_thresholds=ir.y_thresholds_.tolist())

    def predict(self, x: np.ndarray) -> np.ndarray:
        return np.interp(x, self.x_thresholds, self.y_thresholds)

    def to_dict(self) -> dict:
        return {"x_thresholds": self.x_thresholds, "y_thresholds": self.y_thresholds}

    @classmethod
    def from_dict(cls, d: dict) -> IsotonicCalibrator:
        return cls(x_thresholds=d["x_thresholds"], y_thresholds=d["y_thresholds"])


# --------------------------------------------------------------------------
# Cleavage head
# --------------------------------------------------------------------------


class _CleavageMLP(nn.Module):
    def __init__(self, input_dim: int = EMBEDDING_DIM, hidden_dim: int = 128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 2),
        )

    def forward(self, n_vec: torch.Tensor, c_vec: torch.Tensor) -> torch.Tensor:
        return self.net(torch.cat([n_vec, c_vec], dim=-1))  # logits, shape (..., 2)


@dataclass
class CleavageHead:
    mlp: _CleavageMLP
    calibrator_n: IsotonicCalibrator | None = None
    calibrator_c: IsotonicCalibrator | None = None

    @classmethod
    def new(cls, hidden_dim: int = 128) -> CleavageHead:
        return cls(mlp=_CleavageMLP(hidden_dim=hidden_dim))

    def raw_logits(self, n_vecs: np.ndarray, c_vecs: np.ndarray) -> np.ndarray:
        self.mlp.eval()
        with torch.no_grad():
            logits = self.mlp(torch.from_numpy(n_vecs).float(), torch.from_numpy(c_vecs).float())
        return logits.numpy()

    def fit(
        self,
        n_vecs: np.ndarray,
        c_vecs: np.ndarray,
        cleave_n_targets: np.ndarray,
        cleave_c_targets: np.ndarray,
        epochs: int = 300,
        lr: float = 1e-2,
        calibration_frac: float = 0.2,
        seed: int = 0,
    ) -> dict:
        rng = np.random.default_rng(seed)
        n = len(n_vecs)
        idx = rng.permutation(n)
        n_cal = max(1, int(n * calibration_frac))
        cal_idx, train_idx = idx[:n_cal], idx[n_cal:]

        x_n = torch.from_numpy(n_vecs[train_idx]).float()
        x_c = torch.from_numpy(c_vecs[train_idx]).float()
        y = torch.from_numpy(
            np.stack([cleave_n_targets[train_idx], cleave_c_targets[train_idx]], axis=1)
        ).float()

        self.mlp.train()
        optimizer = torch.optim.Adam(self.mlp.parameters(), lr=lr)
        loss_fn = nn.BCEWithLogitsLoss()
        losses = []
        for _ in range(epochs):
            optimizer.zero_grad()
            logits = self.mlp(x_n, x_c)
            loss = loss_fn(logits, y)
            loss.backward()
            optimizer.step()
            losses.append(float(loss.item()))

        cal_logits = self.raw_logits(n_vecs[cal_idx], c_vecs[cal_idx])
        cal_probs = 1 / (1 + np.exp(-cal_logits))
        self.calibrator_n = IsotonicCalibrator.fit(cal_probs[:, 0], cleave_n_targets[cal_idx])
        self.calibrator_c = IsotonicCalibrator.fit(cal_probs[:, 1], cleave_c_targets[cal_idx])
        return {"final_loss": losses[-1], "n_train": len(train_idx), "n_calibration": len(cal_idx)}

    def predict(self, n_vecs: np.ndarray, c_vecs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        logits = self.raw_logits(n_vecs, c_vecs)
        probs = 1 / (1 + np.exp(-logits))
        n_prob = self.calibrator_n.predict(probs[:, 0]) if self.calibrator_n else probs[:, 0]
        c_prob = self.calibrator_c.predict(probs[:, 1]) if self.calibrator_c else probs[:, 1]
        return n_prob, c_prob

    def save(self, path: Path) -> None:
        path.mkdir(parents=True, exist_ok=True)
        torch.save(self.mlp.state_dict(), path / "mlp.pt")
        hidden_dim = self.mlp.net[0].out_features
        meta = {
            "hidden_dim": hidden_dim,
            "calibrator_n": self.calibrator_n.to_dict() if self.calibrator_n else None,
            "calibrator_c": self.calibrator_c.to_dict() if self.calibrator_c else None,
        }
        (path / "meta.json").write_text(json.dumps(meta))

    @classmethod
    def load(cls, path: Path) -> CleavageHead:
        meta = json.loads((path / "meta.json").read_text())
        mlp = _CleavageMLP(hidden_dim=meta["hidden_dim"])
        mlp.load_state_dict(torch.load(path / "mlp.pt", weights_only=True))
        return cls(
            mlp=mlp,
            calibrator_n=IsotonicCalibrator.from_dict(meta["calibrator_n"])
            if meta["calibrator_n"]
            else None,
            calibrator_c=IsotonicCalibrator.from_dict(meta["calibrator_c"])
            if meta["calibrator_c"]
            else None,
        )


# --------------------------------------------------------------------------
# TAP head: ridge regression + bootstrap uncertainty
# --------------------------------------------------------------------------


@dataclass
class TapHead:
    coef: np.ndarray  # (D,)
    intercept: float
    bootstrap_coef: np.ndarray  # (B, D)
    bootstrap_intercept: np.ndarray  # (B,)
    alpha: float

    @classmethod
    def fit(
        cls, x: np.ndarray, y: np.ndarray, alpha: float = 10.0, n_bootstrap: int = 50, seed: int = 0
    ) -> TapHead:
        from sklearn.linear_model import Ridge

        primary = Ridge(alpha=alpha).fit(x, y)
        rng = np.random.default_rng(seed)
        n = len(x)
        boot_coef, boot_intercept = [], []
        for _ in range(n_bootstrap):
            sample = rng.integers(0, n, n)
            model = Ridge(alpha=alpha).fit(x[sample], y[sample])
            boot_coef.append(model.coef_)
            boot_intercept.append(model.intercept_)
        return cls(
            coef=primary.coef_,
            intercept=float(primary.intercept_),
            bootstrap_coef=np.array(boot_coef),
            bootstrap_intercept=np.array(boot_intercept),
            alpha=alpha,
        )

    def predict(self, x: np.ndarray) -> np.ndarray:
        return x @ self.coef + self.intercept

    def predict_with_uncertainty(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        mean = self.predict(x)
        bootstrap_preds = x @ self.bootstrap_coef.T + self.bootstrap_intercept  # (n, B)
        return mean, bootstrap_preds.std(axis=1)

    def save(self, path: Path) -> None:
        path.mkdir(parents=True, exist_ok=True)
        np.savez(
            path / "weights.npz",
            coef=self.coef,
            bootstrap_coef=self.bootstrap_coef,
            bootstrap_intercept=self.bootstrap_intercept,
        )
        (path / "meta.json").write_text(
            json.dumps({"intercept": self.intercept, "alpha": self.alpha})
        )

    @classmethod
    def load(cls, path: Path) -> TapHead:
        arrays = np.load(path / "weights.npz")
        meta = json.loads((path / "meta.json").read_text())
        return cls(
            coef=arrays["coef"],
            intercept=meta["intercept"],
            bootstrap_coef=arrays["bootstrap_coef"],
            bootstrap_intercept=arrays["bootstrap_intercept"],
            alpha=meta["alpha"],
        )


# --------------------------------------------------------------------------
# MHC head: PCA projection + per-allele binder centroid + calibration
# --------------------------------------------------------------------------


@dataclass
class MhcHead:
    projection_weight: np.ndarray  # (output_dim, EMBEDDING_DIM)
    projection_bias: np.ndarray  # (output_dim,)
    centroids: dict[str, np.ndarray] = field(
        default_factory=dict
    )  # allele -> (output_dim,) unit vector
    calibrators: dict[str, IsotonicCalibrator] = field(default_factory=dict)

    @property
    def output_dim(self) -> int:
        return self.projection_weight.shape[0]

    def project(self, pooled: np.ndarray) -> np.ndarray:
        """pooled: (N, EMBEDDING_DIM) -> unit-normalized (N, output_dim)."""
        z = pooled @ self.projection_weight.T + self.projection_bias
        norms = np.linalg.norm(z, axis=-1, keepdims=True)
        return z / np.clip(norms, 1e-6, None)

    @classmethod
    def fit(
        cls,
        pooled: np.ndarray,
        mhc_percentile: np.ndarray,
        hla_allele: np.ndarray,
        output_dim: int = 64,
        binder_percentile_threshold: float = 2.0,
    ) -> MhcHead:
        from sklearn.decomposition import PCA

        pca = PCA(n_components=output_dim, random_state=0).fit(pooled)
        weight = pca.components_  # (output_dim, EMBEDDING_DIM)
        bias = -pca.components_ @ pca.mean_
        head = cls(projection_weight=weight, projection_bias=bias)

        z = head.project(pooled)
        for allele in np.unique(hla_allele):
            mask = hla_allele == allele
            binder_mask = mask & (mhc_percentile <= binder_percentile_threshold)
            if binder_mask.sum() == 0:
                binder_mask = mask  # fall back to whole-allele mean if no strong binders present
            centroid = z[binder_mask].mean(axis=0)
            centroid = centroid / np.clip(np.linalg.norm(centroid), 1e-6, None)
            head.centroids[allele] = centroid

            cos_sim = z[mask] @ centroid
            target = np.clip(1 - mhc_percentile[mask] / 100, 0, 1)
            head.calibrators[allele] = IsotonicCalibrator.fit(cos_sim, target)
        return head

    def score(self, pooled: np.ndarray, hla_allele: str) -> dict:
        if hla_allele not in self.centroids:
            raise KeyError(
                f"no centroid calibrated for {hla_allele!r}; available: {list(self.centroids)}"
            )
        z = self.project(pooled[None, :])[0]
        cos_sim = float(z @ self.centroids[hla_allele])
        propensity = float(self.calibrators[hla_allele].predict(np.array([cos_sim]))[0])
        return {"cosine_similarity": cos_sim, "presentation_propensity": propensity}

    def save(self, path: Path) -> None:
        path.mkdir(parents=True, exist_ok=True)
        np.savez(
            path / "weights.npz",
            projection_weight=self.projection_weight,
            projection_bias=self.projection_bias,
            **{f"centroid__{a}": v for a, v in self.centroids.items()},
        )
        calibrators = {a: c.to_dict() for a, c in self.calibrators.items()}
        (path / "meta.json").write_text(json.dumps({"calibrators": calibrators}))

    @classmethod
    def load(cls, path: Path) -> MhcHead:
        arrays = np.load(path / "weights.npz")
        centroids = {
            k.removeprefix("centroid__"): arrays[k]
            for k in arrays.files
            if k.startswith("centroid__")
        }
        meta = json.loads((path / "meta.json").read_text())
        calibrators = {a: IsotonicCalibrator.from_dict(d) for a, d in meta["calibrators"].items()}
        return cls(
            projection_weight=arrays["projection_weight"],
            projection_bias=arrays["projection_bias"],
            centroids=centroids,
            calibrators=calibrators,
        )


# --------------------------------------------------------------------------
# Bundle + checkpoint manifest
# --------------------------------------------------------------------------


@dataclass
class ThreeHeadModel:
    cleavage: CleavageHead
    tap: TapHead
    mhc: MhcHead
    esm3_model_id: str = ESM3_MODEL_ID
    model_version: str = MODEL_VERSION
    dataset_version_hash: str = ""

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        self.cleavage.save(path / "cleavage")
        self.tap.save(path / "tap")
        self.mhc.save(path / "mhc")
        (path / "manifest.json").write_text(
            json.dumps(
                {
                    "esm3_model_id": self.esm3_model_id,
                    "model_version": self.model_version,
                    "dataset_version_hash": self.dataset_version_hash,
                },
                indent=2,
            )
        )

    @classmethod
    def load(cls, path: str | Path) -> ThreeHeadModel:
        path = Path(path)
        manifest = json.loads((path / "manifest.json").read_text())
        return cls(
            cleavage=CleavageHead.load(path / "cleavage"),
            tap=TapHead.load(path / "tap"),
            mhc=MhcHead.load(path / "mhc"),
            esm3_model_id=manifest["esm3_model_id"],
            model_version=manifest["model_version"],
            dataset_version_hash=manifest["dataset_version_hash"],
        )
