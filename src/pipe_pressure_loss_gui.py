# pipe_pressure_loss_gui.py
# Desktop GUI for Fluids Lab activity (Windows-friendly; PySide6 + Matplotlib)
#
# Student objective:
#   - Paste data copied from Excel
#   - Change only:
#       p_meas_err_delta_p   [%]
#       p_meas_err_V         [L]   (absolute)
#       p_meas_err_dt        [s]   (absolute)
#   - Re-run and inspect whether experimental uncertainty bounds cover theoretical bounds
#
# Packaging (Windows):
#   pip install pyside6 matplotlib pandas numpy pyinstaller
#   pyinstaller --noconfirm --onefile --windowed pipe_pressure_loss_gui.py
#
# Notes on pasted data:
#   Supported formats (tab-separated rows copied from Excel):
#   1) 4-column format (recommended):
#        D_mm    delta_p_mmH2O    V_accumulated_L    delta_t_s
#   2) 5-column format:
#        D_mm    start_V_L    end_V_L    delta_t_s    delta_p_mmH2O
#      (GUI computes V_accumulated = end - start; delta_p already in mmH2O)
#
# Theoretical f is computed from Colebrook-White with uncertainty propagated
# from V and dt only (matching your Example 7 concept), while experimental f
# uncertainty propagates delta_p (%), V (absolute L), and dt (absolute s).

from __future__ import annotations

import math
import sys
import traceback
from dataclasses import dataclass
from typing import Tuple, Optional

import numpy as np
import pandas as pd

from PySide6.QtCore import Qt, Signal, QObject
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QLabel,
    QPushButton,
    QPlainTextEdit,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QMessageBox,
    QGroupBox,
    QDoubleSpinBox,
    QSpinBox,
    QSlider,
    QSplitter,
    QFrame,
    QCheckBox,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


# ======================================================================================
# Core physics / uncertainty functions (adapted from your script, with modified model)
# ======================================================================================

def darcy_weisbach_f_lab_experimental(
    D_mm: float,
    L_m: float,
    delta_p_mmH2O: float,
    V_accumulated_L: float,
    delta_t_s: float,
    rho: float = 1000.0,
) -> float:
    """
    Experimental Darcy-Weisbach friction factor from measured pressure drop and flow.
    delta_p is mmH2O, V is L, D is mm, t is s.
    """
    if D_mm <= 0:
        raise ValueError("D must be > 0 mm.")
    if L_m <= 0:
        raise ValueError("L must be > 0 m.")
    if delta_p_mmH2O <= 0:
        raise ValueError("delta_p must be > 0 mmH2O.")
    if V_accumulated_L <= 0:
        raise ValueError("V_accumulated must be > 0 L.")
    if delta_t_s <= 0:
        raise ValueError("delta_t must be > 0 s.")
    if rho <= 0:
        raise ValueError("rho must be > 0 kg/m^3.")

    A = math.pi * ((D_mm / 1000.0) ** 2) / 4.0
    Q = (V_accumulated_L / 1000.0) / delta_t_s  # m3/s
    C = Q / A  # m/s

    # delta_p_mmH2O * 9.81 gives Pa (since mmH2O -> mH2O via /1000 and rho*g*h with rho=1000)
    return (2.0 * (delta_p_mmH2O * 9.81) * (D_mm / 1000.0)) / (rho * L_m * (C ** 2))


def _uniform_bounds_percent(center: float, percent: float) -> Tuple[float, float]:
    if percent < 0:
        raise ValueError("Percent uncertainty must be >= 0.")
    lo = center * (1.0 - percent / 100.0)
    hi = center * (1.0 + percent / 100.0)
    return lo, hi


def _uniform_bounds_absolute(center: float, abs_err: float) -> Tuple[float, float]:
    if abs_err < 0:
        raise ValueError("Absolute uncertainty must be >= 0.")
    lo = center - abs_err
    hi = center + abs_err
    return lo, hi


def _velocity_from_V_dt(D_mm: float, V_L, dt_s):
    """Velocity from accumulated volume and elapsed time. Supports scalar or numpy arrays."""
    A = math.pi * ((D_mm / 1000.0) ** 2) / 4.0
    Q = (np.asarray(V_L, dtype=float) / 1000.0) / np.asarray(dt_s, dtype=float)
    return Q / A


def _reynolds_from_velocity(D_mm: float, C, rho: float, mu: float):
    """
    Re = rho * V * D / mu
    D in mm (converted internally to m), V in m/s, rho in kg/m^3, mu in Pa.s
    Supports scalar/array C.
    """
    if rho <= 0:
        raise ValueError("rho must be > 0.")
    if mu <= 0:
        raise ValueError("mu must be > 0.")
    return (rho * np.asarray(C, dtype=float) * (D_mm / 1000.0)) / mu


def friction_factor_colebrook_white_scalar(
    Re: float,
    D_m: float,
    eps_m: float,
    max_iter: int = 50,
    tol: float = 1e-12,
) -> float:
    """
    Colebrook-White solved iteratively for Darcy friction factor.
    Handles turbulent/transitional-ish values; for laminar uses 64/Re.
    """
    if Re <= 0:
        raise ValueError("Re must be > 0.")
    if D_m <= 0:
        raise ValueError("D must be > 0.")
    if eps_m < 0:
        raise ValueError("eps must be >= 0.")

    if Re < 2300:
        return 64.0 / Re

    # Haaland initial guess (good starting point)
    rr = eps_m / D_m
    inv_sqrt_f = -1.8 * math.log10((rr / 3.7) ** 1.11 + 6.9 / Re)
    f = 1.0 / (inv_sqrt_f ** 2)

    for _ in range(max_iter):
        if f <= 0:
            f = 0.02
        rhs = -2.0 * math.log10(rr / 3.7 + 2.51 / (Re * math.sqrt(f)))
        f_new = 1.0 / (rhs ** 2)
        if abs(f_new - f) < tol:
            return f_new
        f = f_new

    return f


def friction_factor_colebrook_white_array(
    Re_arr,
    D_m: float,
    eps_m: float,
    max_iter: int = 25,
) -> np.ndarray:
    """
    Vectorized Colebrook-White using fixed-point updates.
    Laminar region uses 64/Re.
    """
    Re = np.asarray(Re_arr, dtype=float)
    if np.any(Re <= 0):
        raise ValueError("All Reynolds numbers must be > 0.")
    if D_m <= 0:
        raise ValueError("D_m must be > 0.")
    if eps_m < 0:
        raise ValueError("eps_m must be >= 0.")

    f = np.empty_like(Re, dtype=float)
    lam = Re < 2300.0
    turb = ~lam

    # Laminar exact
    f[lam] = 64.0 / Re[lam]

    if np.any(turb):
        Re_t = Re[turb]
        rr = eps_m / D_m

        # Haaland initial guess
        inv_sqrt_f = -1.8 * np.log10((rr / 3.7) ** 1.11 + 6.9 / Re_t)
        ft = 1.0 / (inv_sqrt_f ** 2)

        for _ in range(max_iter):
            ft = np.clip(ft, 1e-8, None)
            rhs = -2.0 * np.log10(rr / 3.7 + 2.51 / (Re_t * np.sqrt(ft)))
            ft = 1.0 / (rhs ** 2)

        f[turb] = ft

    return f


def mc_experimental_f_summary(
    D_mm: float,
    L_m: float,
    delta_p_mmH2O: float,
    V_accumulated_L: float,
    delta_t_s: float,
    p_meas_err_delta_p_pct: float,  # percent
    p_meas_err_V_abs_L: float,      # ABSOLUTE liters
    p_meas_err_dt_abs_s: float,     # ABSOLUTE seconds
    n_samples: int = 20000,
    seed: int = 1,
    rho: float = 1000.0,
) -> dict:
    """
    MC for experimental f with modified uncertainty model:
      - delta_p : percent
      - V       : absolute liters
      - dt      : absolute seconds
    """
    if n_samples < 100:
        raise ValueError("n_samples should be at least 100 for stable percentiles.")

    dp_lo, dp_hi = _uniform_bounds_percent(delta_p_mmH2O, p_meas_err_delta_p_pct)
    V_lo, V_hi = _uniform_bounds_absolute(V_accumulated_L, p_meas_err_V_abs_L)
    dt_lo, dt_hi = _uniform_bounds_absolute(delta_t_s, p_meas_err_dt_abs_s)

    if dp_lo <= 0:
        raise ValueError("delta_p lower bound <= 0. Reduce p_meas_err_delta_p (%).")
    if V_lo <= 0:
        raise ValueError("V_accumulated lower bound <= 0. Reduce p_meas_err_V (L).")
    if dt_lo <= 0:
        raise ValueError("delta_t lower bound <= 0. Reduce p_meas_err_dt (s).")

    rng = np.random.default_rng(seed)
    dp_s = rng.uniform(dp_lo, dp_hi, size=n_samples)
    V_s = rng.uniform(V_lo, V_hi, size=n_samples)
    dt_s = rng.uniform(dt_lo, dt_hi, size=n_samples)

    # Vectorized experimental f
    D_m = D_mm / 1000.0
    A = math.pi * (D_m ** 2) / 4.0
    Q = (V_s / 1000.0) / dt_s
    C = Q / A
    f_s = (2.0 * (dp_s * 9.81) * D_m) / (rho * L_m * (C ** 2))

    p05 = float(np.percentile(f_s, 5))
    p95 = float(np.percentile(f_s, 95))
    return {
        "f_nominal": float(darcy_weisbach_f_lab_experimental(
            D_mm=D_mm,
            L_m=L_m,
            delta_p_mmH2O=delta_p_mmH2O,
            V_accumulated_L=V_accumulated_L,
            delta_t_s=delta_t_s,
            rho=rho,
        )),
        "f_mean": float(np.mean(f_s)),
        "f_std": float(np.std(f_s, ddof=1)),
        "f_p05": p05,
        "f_p95": p95,
    }


def mc_theoretical_f_colebrook_from_V_dt_summary(
    D_mm: float,
    V_L: float,
    dt_s: float,
    p_meas_err_V_abs_L: float,
    p_meas_err_dt_abs_s: float,
    eps_m: float,
    rho: float,
    mu: float,
    n_samples: int = 20000,
    seed: int = 1,
) -> dict:
    """
    MC theoretical Darcy friction factor from Colebrook-White.
    Uncertainty is propagated from V and dt ONLY (matching your Example 7 logic),
    but now with ABSOLUTE uncertainties for V [L] and dt [s].
    """
    C0 = float(_velocity_from_V_dt(D_mm, V_L, dt_s))
    Re0 = float(_reynolds_from_velocity(D_mm, C0, rho=rho, mu=mu))
    f0 = float(friction_factor_colebrook_white_scalar(Re0, D_mm / 1000.0, eps_m))

    V_lo, V_hi = _uniform_bounds_absolute(V_L, p_meas_err_V_abs_L)
    dt_lo, dt_hi = _uniform_bounds_absolute(dt_s, p_meas_err_dt_abs_s)
    if V_lo <= 0:
        raise ValueError("Theoretical MC: V lower bound <= 0. Reduce p_meas_err_V (L).")
    if dt_lo <= 0:
        raise ValueError("Theoretical MC: dt lower bound <= 0. Reduce p_meas_err_dt (s).")

    if p_meas_err_V_abs_L == 0 and p_meas_err_dt_abs_s == 0:
        return {
            "Re_nominal": Re0,
            "f_th_nominal": f0,
            "f_th_mean": f0,
            "f_th_std": 0.0,
            "f_th_p05": f0,
            "f_th_p95": f0,
        }

    rng = np.random.default_rng(seed)
    V_s = rng.uniform(V_lo, V_hi, size=n_samples)
    dt_s_arr = rng.uniform(dt_lo, dt_hi, size=n_samples)

    C = _velocity_from_V_dt(D_mm, V_s, dt_s_arr)
    Re = _reynolds_from_velocity(D_mm, C, rho=rho, mu=mu)
    f_th = friction_factor_colebrook_white_array(Re, D_m=D_mm / 1000.0, eps_m=eps_m)

    return {
        "Re_nominal": Re0,
        "f_th_nominal": f0,
        "f_th_mean": float(np.mean(f_th)),
        "f_th_std": float(np.std(f_th, ddof=1)),
        "f_th_p05": float(np.percentile(f_th, 5)),
        "f_th_p95": float(np.percentile(f_th, 95)),
    }


def parse_pasted_table(text: str) -> pd.DataFrame:
    """
    Parse pasted Excel text (tab-separated) into canonical columns:
      D_mm, delta_p_mmH2O, V_accumulated_L, delta_t_s

    REQUIRED pasted format (exact order, tab-separated):
      D_mm    start_V_L    end_V_L    delta_t_s    delta_p_mmH2O

    Notes:
    - Blank lines are ignored.
    - A single header row is auto-skipped if non-numeric.
    - delta_p is read directly as mmH2O (no unit conversion)
    - V_accumulated_L = end_V_L - start_V_L
    """
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError("Paste box is empty.")

    rows = []
    for ln in lines:
        # Expect Excel copy-paste (tab-separated). Allow comma fallback only if user pasted CSV.
        parts = [p.strip() for p in (ln.split("\t") if "\t" in ln else ln.split(","))]
        rows.append(parts)

    # Require at least 5 columns now (strict format)
    rows = [r for r in rows if len(r) >= 5]
    if not rows:
        raise ValueError(
            "No valid rows found. Expected rows with 5 columns in this exact order:\n"
            "D_mm    start_V_L    end_V_L    delta_t_s    delta_p_mmH2O"
        )

    # Convert to DataFrame (ragged rows padded)
    max_cols = max(len(r) for r in rows)
    rows_pad = [r + [""] * (max_cols - len(r)) for r in rows]
    df_raw = pd.DataFrame(rows_pad)

    # Auto-skip first row if it looks like a header (non-numeric)
    def row_numeric_fraction(sr: pd.Series) -> float:
        vals = pd.to_numeric(sr, errors="coerce")
        return float(vals.notna().mean())

    if row_numeric_fraction(df_raw.iloc[0, :5]) < 0.8:
        df_raw = df_raw.iloc[1:].reset_index(drop=True)

    if df_raw.empty:
        raise ValueError("After removing header row, no data rows remain.")

    # Strictly parse first 5 columns as the required format
    first5 = df_raw.iloc[:, :5].apply(pd.to_numeric, errors="coerce")
    if not first5.notna().all(axis=1).all():
        raise ValueError(
            "Could not parse pasted data as numeric values.\n\n"
            "Required format (exact column order):\n"
            "D_mm    start_V_L    end_V_L    delta_t_s    delta_p_mmH2O"
        )

    D = first5.iloc[:, 0].astype(float)
    start_V = first5.iloc[:, 1].astype(float)
    end_V = first5.iloc[:, 2].astype(float)
    dt = first5.iloc[:, 3].astype(float)
    dp_mm = first5.iloc[:, 4].astype(float)

    df = pd.DataFrame({
        "D_mm": D,
        "delta_p_mmH2O": dp_mm,              # already mmH2O (no conversion)
        "V_accumulated_L": end_V - start_V,  # accumulated volume
        "delta_t_s": dt,
    })

    return _validate_input_df(df)


def _validate_input_df(df: pd.DataFrame) -> pd.DataFrame:
    required = ["D_mm", "delta_p_mmH2O", "V_accumulated_L", "delta_t_s"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    if df.empty:
        raise ValueError("No rows found.")
    for c in required:
        if (df[c] <= 0).any():
            bad_idx = df.index[df[c] <= 0].tolist()[:5]
            raise ValueError(f"Column {c} contains non-positive value(s). Example row index: {bad_idx}")

    return df.reset_index(drop=True)


# ======================================================================================
# Analysis workflow (Example-7 style)
# ======================================================================================

@dataclass
class RunConfig:
    L_m: float = 0.36
    eps_m: float = 0.0
    mu_Pas: float = 0.0009764
    rho_kgm3: float = 1000.0
    n_samples: int = 20000
    seed: int = 1


def run_activity_analysis(
    df_in: pd.DataFrame,
    p_meas_err_delta_p_pct: float,
    p_meas_err_V_abs_L: float,
    p_meas_err_dt_abs_s: float,
    cfg: RunConfig,
) -> tuple[pd.DataFrame, dict]:
    """
    Runs the Example-7 style workflow for all rows and computes coverage metrics.
    """
    rows = []
    for i, r in df_in.iterrows():
        Di = float(r["D_mm"])
        dpi = float(r["delta_p_mmH2O"])
        Vi = float(r["V_accumulated_L"])
        dti = float(r["delta_t_s"])

        exp_res = mc_experimental_f_summary(
            D_mm=Di,
            L_m=cfg.L_m,
            delta_p_mmH2O=dpi,
            V_accumulated_L=Vi,
            delta_t_s=dti,
            p_meas_err_delta_p_pct=p_meas_err_delta_p_pct,
            p_meas_err_V_abs_L=p_meas_err_V_abs_L,
            p_meas_err_dt_abs_s=p_meas_err_dt_abs_s,
            n_samples=cfg.n_samples,
            seed=cfg.seed,
            rho=cfg.rho_kgm3,
        )

        th_res = mc_theoretical_f_colebrook_from_V_dt_summary(
            D_mm=Di,
            V_L=Vi,
            dt_s=dti,
            p_meas_err_V_abs_L=p_meas_err_V_abs_L,
            p_meas_err_dt_abs_s=p_meas_err_dt_abs_s,
            eps_m=cfg.eps_m,
            rho=cfg.rho_kgm3,
            mu=cfg.mu_Pas,
            n_samples=cfg.n_samples,
            seed=cfg.seed,
        )

        # Coverage logic
        # "Theoretical interval inside experimental interval"
        theo_inside_exp = (
            (th_res["f_th_p05"] >= exp_res["f_p05"]) and
            (th_res["f_th_p95"] <= exp_res["f_p95"])
        )

        # Any interval overlap
        overlap = not (
            (exp_res["f_p95"] < th_res["f_th_p05"]) or
            (th_res["f_th_p95"] < exp_res["f_p05"])
        )

        rows.append({
            "meas_id": int(i),
            "D_mm": Di,
            "L_m": cfg.L_m,
            "Δp_mmH2O": dpi,
            "ΔV_L": Vi,
            "Δt_s": dti,

            "p_err_Δp_%": p_meas_err_delta_p_pct,     # percent
            "err_ΔV_abs_L": p_meas_err_V_abs_L,       # absolute liters
            "err_Δt_abs_s": p_meas_err_dt_abs_s,      # absolute seconds

            # Experimental
            "f_nominal": exp_res["f_nominal"],
            "f_mean": exp_res["f_mean"],
            "f_std": exp_res["f_std"],
            "f_p05": exp_res["f_p05"],
            "f_p95": exp_res["f_p95"],

            # Theoretical (Colebrook)
            "Re_nominal": th_res["Re_nominal"],
            "f_th_nominal": th_res["f_th_nominal"],
            "f_th_mean": th_res["f_th_mean"],
            "f_th_std": th_res["f_th_std"],
            "f_th_p05": th_res["f_th_p05"],
            "f_th_p95": th_res["f_th_p95"],

            # Coverage flags
            "theoretical_inside_experimental": int(theo_inside_exp),
            "any_overlap": int(overlap),
        })

    df_out = pd.DataFrame(rows)

    n = len(df_out)
    summary = {
        "n_rows": n,
        "n_inside": int(df_out["theoretical_inside_experimental"].sum()),
        "pct_inside": float(100.0 * df_out["theoretical_inside_experimental"].mean()) if n else 0.0,
        "n_overlap": int(df_out["any_overlap"].sum()),
        "pct_overlap": float(100.0 * df_out["any_overlap"].mean()) if n else 0.0,
    }
    return df_out, summary


# ======================================================================================
# GUI widgets
# ======================================================================================

class ValueSlider(QWidget):
    """
    Slider + numeric spinbox synced together.
    Supports integer or decimal values via scale factor.
    """
    valueChanged = Signal(float)

    def __init__(self, label: str, min_val: float, max_val: float, step: float, initial: float, suffix: str = ""):
        super().__init__()
        self.scale = int(round(1 / step)) if step < 1 else 1
        self.step = step

        layout = QGridLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.label = QLabel(label)
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(int(round(min_val * self.scale)))
        self.slider.setMaximum(int(round(max_val * self.scale)))
        self.slider.setValue(int(round(initial * self.scale)))

        if step >= 1:
            self.spin = QSpinBox()
            self.spin.setRange(int(min_val), int(max_val))
            self.spin.setSingleStep(int(step))
            self.spin.setValue(int(round(initial)))
            if suffix:
                self.spin.setSuffix(f" {suffix}")
        else:
            self.spin = QDoubleSpinBox()
            self.spin.setRange(min_val, max_val)
            self.spin.setSingleStep(step)
            self.spin.setDecimals(max(1, int(round(-math.log10(step)))))
            self.spin.setValue(initial)
            if suffix:
                self.spin.setSuffix(f" {suffix}")

        self.slider.valueChanged.connect(self._slider_to_spin)
        self.spin.valueChanged.connect(self._spin_to_slider)

        layout.addWidget(self.label, 0, 0, 1, 2)
        layout.addWidget(self.slider, 1, 0)
        layout.addWidget(self.spin, 1, 1)

    def _slider_to_spin(self, ivalue: int):
        val = ivalue / self.scale
        self.spin.blockSignals(True)
        self.spin.setValue(val)
        self.spin.blockSignals(False)
        self.valueChanged.emit(float(val))

    def _spin_to_slider(self, value):
        val = float(value)
        self.slider.blockSignals(True)
        self.slider.setValue(int(round(val * self.scale)))
        self.slider.blockSignals(False)
        self.valueChanged.emit(val)

    def value(self) -> float:
        return float(self.spin.value())


class MplCanvas(FigureCanvas):
    def __init__(self):
        self.fig = Figure(figsize=(7, 5), tight_layout=True)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Pipe Pressure Loss Uncertainty Activity (Desktop)")
        self.resize(1400, 900)

        self.df_input: Optional[pd.DataFrame] = None
        self.df_results: Optional[pd.DataFrame] = None

        central = QWidget()
        self.setCentralWidget(central)

        root_layout = QVBoxLayout(central)

        title = QLabel("Fluids Lab Activity for 'Pressure Loss in Pipes' : Calculation of Experimental vs Theoretical Friction Factor with Measurement Uncertainty")
        title_font = QFont()
        title_font.setPointSize(12)
        title_font.setBold(True)
        title.setFont(title_font)
        root_layout.addWidget(title)

        subtitle = QLabel(
            "   (1) Paste data from Excel, (2) adjust uncertainty assumptions, and (3) re-run. "
        )
        subtitle.setWordWrap(True)
        root_layout.addWidget(subtitle)

        splitter = QSplitter(Qt.Horizontal)
        root_layout.addWidget(splitter, 1)

        footer_row = QHBoxLayout()
        footer_row.addStretch(1)
        author_label = QLabel(
            "Developed by Mustafa Onur Onen | moonen1@sheffield.ac.uk\n"
            "MEE Fluids Lab - The University of Sheffield - 2026 • Built with OpenAI ChatGPT"
        )
        author_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        author_label.setStyleSheet("color: #555; padding-top: 4px; padding-right: 6px; font-size: 10pt;")
        author_label.setWordWrap(False)
        footer_row.addWidget(author_label)
        root_layout.addLayout(footer_row)


        # ---------------- LEFT PANEL: inputs ----------------
        left = QWidget()
        left_layout = QVBoxLayout(left)

        # Paste area
        paste_group = QGroupBox("1) Paste measurement table (from Excel)")
        paste_layout = QVBoxLayout(paste_group)

        self.paste_text = QPlainTextEdit()
        self.paste_text.setPlaceholderText(
            "Paste Excel rows in this exact format (tab-separated):\n"
            "D(mm)    start_V(L)    end_V(L)    Time(s)    ΔhL(mmH2O)\n\n"
            "Example:\n"
            "10    2.0    7.0    132    50.0"
        )
        paste_layout.addWidget(self.paste_text)

        btn_row = QHBoxLayout()
        self.btn_parse = QPushButton("Parse / Preview")
        self.btn_load_demo = QPushButton("Load Demo Data")
        self.btn_clear = QPushButton("Clear")
        btn_row.addWidget(self.btn_parse)
        btn_row.addWidget(self.btn_load_demo)
        btn_row.addWidget(self.btn_clear)
        paste_layout.addLayout(btn_row)

        left_layout.addWidget(paste_group)

        # Student controls
        control_group = QGroupBox("2) Student uncertainty inputs")
        control_layout = QVBoxLayout(control_group)

        self.s_dp = ValueSlider(
            label="Pressure loss uncertainty, ± (%)",
            min_val=0, max_val=50, step=0.1, initial=0.0, suffix="%"
        )
        self.s_V = ValueSlider(
            label="Volume uncertainty, ± (L) [ABSOLUTE]",
            min_val=0, max_val=1, step=0.01, initial=0.0, suffix="L"
        )
        self.s_dt = ValueSlider(
            label="Time uncertainty, ± (s) [ABSOLUTE]",
            min_val=0, max_val=20, step=0.1, initial=0.0, suffix="s"
        )
        control_layout.addWidget(self.s_dp)
        control_layout.addWidget(self.s_V)
        control_layout.addWidget(self.s_dt)

        self.btn_run = QPushButton("3) Run Analysis")
        self.btn_run.setMinimumHeight(36)
        control_layout.addWidget(self.btn_run)

        left_layout.addWidget(control_group)

        # Advanced controls
        adv_group = QGroupBox("Advanced settings (instructor)")
        adv_layout = QGridLayout(adv_group)

        self.sp_L = QDoubleSpinBox()
        self.sp_L.setRange(0.001, 100.0)
        self.sp_L.setDecimals(4)
        self.sp_L.setValue(0.36)
        self.sp_L.setSuffix(" m")

        self.sp_eps = QDoubleSpinBox()
        self.sp_eps.setRange(0.0, 0.01)
        self.sp_eps.setDecimals(6)
        self.sp_eps.setSingleStep(0.00001)
        self.sp_eps.setValue(0.0)
        self.sp_eps.setSuffix(" m")

        self.sp_mu = QDoubleSpinBox()
        self.sp_mu.setRange(1e-5, 1e-1)
        self.sp_mu.setDecimals(7)
        self.sp_mu.setSingleStep(0.0001)
        self.sp_mu.setValue(0.0009764)
        self.sp_mu.setSuffix(" Pa·s")

        self.sp_rho = QDoubleSpinBox()
        self.sp_rho.setRange(1.0, 2000.0)
        self.sp_rho.setDecimals(1)
        self.sp_rho.setValue(1000.0)
        self.sp_rho.setSuffix(" kg/m³")

        self.sp_nsamples = QSpinBox()
        self.sp_nsamples.setRange(1000, 200000)
        self.sp_nsamples.setSingleStep(1000)
        self.sp_nsamples.setValue(20000)

        self.sp_seed = QSpinBox()
        self.sp_seed.setRange(0, 999999)
        self.sp_seed.setValue(1)

        self.chk_show_means = QCheckBox("Show mean curves on plot")
        self.chk_show_means.setChecked(True)

        adv_layout.addWidget(QLabel("L"), 0, 0)
        adv_layout.addWidget(self.sp_L, 0, 1)
        adv_layout.addWidget(QLabel("eps"), 1, 0)
        adv_layout.addWidget(self.sp_eps, 1, 1)
        adv_layout.addWidget(QLabel("mu"), 2, 0)
        adv_layout.addWidget(self.sp_mu, 2, 1)
        adv_layout.addWidget(QLabel("rho"), 3, 0)
        adv_layout.addWidget(self.sp_rho, 3, 1)
        adv_layout.addWidget(QLabel("MC samples / row"), 4, 0)
        adv_layout.addWidget(self.sp_nsamples, 4, 1)
        adv_layout.addWidget(QLabel("Random seed"), 5, 0)
        adv_layout.addWidget(self.sp_seed, 5, 1)
        adv_layout.addWidget(self.chk_show_means, 6, 0, 1, 2)

        left_layout.addWidget(adv_group)

        # Status / summary
        summary_group = QGroupBox("Run summary")
        summary_layout = QVBoxLayout(summary_group)
        self.summary_text = QPlainTextEdit()
        self.summary_text.setReadOnly(True)
        self.summary_text.setMaximumHeight(180)
        summary_layout.addWidget(self.summary_text)

        btn_export_row = QHBoxLayout()
        self.btn_export_csv = QPushButton("Export Results CSV...")
        self.btn_export_csv.setEnabled(False)
        btn_export_row.addWidget(self.btn_export_csv)
        summary_layout.addLayout(btn_export_row)

        left_layout.addWidget(summary_group)
        left_layout.addStretch(1)

        splitter.addWidget(left)

        # ---------------- RIGHT PANEL: tables + plot ----------------
        right = QWidget()
        right_layout = QVBoxLayout(right)

        tabs = QTabWidget()

        # Input preview table
        self.table_input = QTableWidget()
        self.table_input.setAlternatingRowColors(True)
        tabs.addTab(self.table_input, "Parsed Input Preview")

        # Results table
        self.table_results = QTableWidget()
        self.table_results.setAlternatingRowColors(True)
        tabs.addTab(self.table_results, "Results Table")

        # Plot tab
        plot_container = QWidget()
        plot_layout = QVBoxLayout(plot_container)
        self.canvas = MplCanvas()
        plot_layout.addWidget(self.canvas)
        tabs.addTab(plot_container, "Uncertainty Plot")

        right_layout.addWidget(tabs)
        splitter.addWidget(right)

        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        # Connections
        self.btn_parse.clicked.connect(self.on_parse_clicked)
        self.btn_load_demo.clicked.connect(self.on_load_demo_clicked)
        self.btn_clear.clicked.connect(self.on_clear_clicked)
        self.btn_run.clicked.connect(self.on_run_clicked)
        self.btn_export_csv.clicked.connect(self.on_export_csv_clicked)

        self._write_summary("Ready.\nPaste data and click 'Parse / Preview'.")

    # ------------------------------------------------------------------
    # UI helpers
    # ------------------------------------------------------------------
    def _write_summary(self, text: str):
        self.summary_text.setPlainText(text)

    def _append_summary(self, text: str):
        cur = self.summary_text.toPlainText().rstrip()
        self.summary_text.setPlainText(cur + ("\n" if cur else "") + text)

    def _set_table_from_df(self, table: QTableWidget, df: pd.DataFrame, max_rows: Optional[int] = None):
        if df is None or df.empty:
            table.clear()
            table.setRowCount(0)
            table.setColumnCount(0)
            return

        show_df = df if max_rows is None else df.head(max_rows).copy()
        table.clear()
        table.setRowCount(len(show_df))
        table.setColumnCount(len(show_df.columns))
        table.setHorizontalHeaderLabels([str(c) for c in show_df.columns])

        for i in range(len(show_df)):
            for j, c in enumerate(show_df.columns):
                v = show_df.iloc[i, j]
                if isinstance(v, (float, np.floating)):
                    txt = f"{float(v):.6g}"
                else:
                    txt = str(v)
                item = QTableWidgetItem(txt)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                table.setItem(i, j, item)

        table.resizeColumnsToContents()
        table.resizeRowsToContents()

    def _get_config(self) -> RunConfig:
        return RunConfig(
            L_m=float(self.sp_L.value()),
            eps_m=float(self.sp_eps.value()),
            mu_Pas=float(self.sp_mu.value()),
            rho_kgm3=float(self.sp_rho.value()),
            n_samples=int(self.sp_nsamples.value()),
            seed=int(self.sp_seed.value()),
        )

    # ------------------------------------------------------------------
    # Events
    # ------------------------------------------------------------------
    def on_load_demo_clicked(self):
        demo = (
            "10.27\t1.0\t6.0\t71\t37.2\n"
            "10.27\t1.0\t6.0\t51\t73.7\n"
            "10.27\t1.0\t6.0\t40\t113.0\n"
            "10.27\t1.0\t6.0\t34\t154.5\n"
            "10.08\t1.0\t6.0\t74\t31.6\n"
            "10.08\t1.0\t6.0\t58\t58.0\n"
            "10.08\t1.0\t6.0\t42\t98.3\n"
            "10.08\t1.0\t6.0\t34\t146.5\n"
        )
        self.paste_text.setPlainText(demo)
        self._write_summary("Demo data loaded (required 5-column format). Click 'Parse / Preview'.")

    def on_clear_clicked(self):
        self.paste_text.clear()
        self.df_input = None
        self.df_results = None
        self._set_table_from_df(self.table_input, pd.DataFrame())
        self._set_table_from_df(self.table_results, pd.DataFrame())
        self.canvas.ax.clear()
        self.canvas.draw()
        self.btn_export_csv.setEnabled(False)
        self._write_summary("Cleared.")

    def on_parse_clicked(self):
        try:
            txt = self.paste_text.toPlainText()
            df = parse_pasted_table(txt)
            self.df_input = df
            self._set_table_from_df(self.table_input, df)
            self._write_summary(
                "Parsed successfully.\n"
                f"Rows: {len(df)}\n"
                "Accepted pasted format:\n"
                "  D(mm), Start_Vol(L), End_Vol(L), Time(s), Delta_hL(mmH2O)\n\n"
                "Now adjust uncertainty sliders and click 'Run Analysis'."
            )
        except Exception as e:
            QMessageBox.critical(self, "Parse error", str(e))
            self._append_summary(f"\nParse error:\n{e}")

    def on_run_clicked(self):
        try:
            if self.df_input is None or self.df_input.empty:
                # Try parsing automatically from current text
                self.df_input = parse_pasted_table(self.paste_text.toPlainText())
                self._set_table_from_df(self.table_input, self.df_input)

            cfg = self._get_config()
            p_dp = self.s_dp.value()        # percent
            p_V_abs = self.s_V.value()      # liters
            p_dt_abs = self.s_dt.value()    # seconds

            df_res, summary = run_activity_analysis(
                df_in=self.df_input,
                p_meas_err_delta_p_pct=p_dp,
                p_meas_err_V_abs_L=p_V_abs,
                p_meas_err_dt_abs_s=p_dt_abs,
                cfg=cfg,
            )

            self.df_results = df_res
            self._set_table_from_df(self.table_results, df_res)

            self._plot_results(df_res, show_means=self.chk_show_means.isChecked())
            self.btn_export_csv.setEnabled(True)

            self._write_summary(
                "Run completed.\n\n"
                f"Rows analyzed: {summary['n_rows']}\n"
                f"Theoretical interval fully inside experimental interval: "
                f"{summary['n_inside']} / {summary['n_rows']} "
                f"({summary['pct_inside']:.1f}%)\n"
                f"Any overlap between intervals: "
                f"{summary['n_overlap']} / {summary['n_rows']} "
                f"({summary['pct_overlap']:.1f}%)\n\n"
                "Current uncertainty settings:\n"
                f"  p_meas_err_delta_p = ±{p_dp:.3g} %\n"
                f"  p_meas_err_V       = ±{p_V_abs:.3g} L (absolute)\n"
                f"  p_meas_err_dt      = ±{p_dt_abs:.3g} s (absolute)\n\n"
                "Try changing only these three and rerun."
            )

        except Exception as e:
            tb = traceback.format_exc()
            QMessageBox.critical(self, "Run error", f"{e}\n\nDetails:\n{tb}")
            self._append_summary(f"\nRun error:\n{e}")

    def on_export_csv_clicked(self):
        try:
            if self.df_results is None or self.df_results.empty:
                QMessageBox.information(self, "No results", "No results to export yet.")
                return

            # Use a simple built-in path choice with Qt optional import
            from PySide6.QtWidgets import QFileDialog
            path, _ = QFileDialog.getSaveFileName(
                self, "Save Results CSV", "pipe_pressure_loss_results.csv", "CSV Files (*.csv)"
            )
            if not path:
                return
            self.df_results.to_csv(path, index=False)
            self._append_summary(f"\nResults exported:\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Export error", str(e))

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------
    def _plot_results(self, dfp: pd.DataFrame, show_means: bool = True):
        ax = self.canvas.ax
        ax.clear()

        dfp = dfp.sort_values("meas_id").reset_index(drop=True)
        x = dfp["meas_id"].to_numpy()

        # Experimental uncertainty band
        ax.fill_between(
            x,
            dfp["f_p05"].to_numpy(),
            dfp["f_p95"].to_numpy(),
            alpha=0.25,
            label="Experim. f. factor bounds",
        )

        # Theoretical uncertainty band
        ax.fill_between(
            x,
            dfp["f_th_p05"].to_numpy(),
            dfp["f_th_p95"].to_numpy(),
            alpha=0.25,
            label="Theoret. f. factor bounds",
        )

        if show_means:
            ax.plot(x, dfp["f_mean"].to_numpy(), marker="o", linewidth=1, label="Experim. f. factor (mean)")
            ax.plot(x, dfp["f_th_mean"].to_numpy(), marker="o", linewidth=1, label="Theoret. f. factor (mean)")

        # Highlight rows where theoretical interval is fully inside experimental
        inside_mask = dfp["theoretical_inside_experimental"].to_numpy().astype(bool)
        if inside_mask.any():
            ax.scatter(
                x[inside_mask],
                dfp.loc[inside_mask, "f_th_mean"].to_numpy(),
                marker="s",
                s=40,
                label="theoretical_inside_experimental",
            )

        ax.set_xlabel("meas_id")
        ax.set_ylabel("Darcy–Weisbach friction factor, f (or lambda) [-]")
        ax.set_title("Experimental vs theoretical uncertainty bounds")
        ax.legend()
        ax.grid(True, alpha=0.2)

        self.canvas.draw()


# ======================================================================================
# Main entry point
# ======================================================================================

def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()