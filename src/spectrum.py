"""Reproduce Figure 5 spectral luminosity curves from Joshi et al. (2014).

This module computes the spectral luminosity distribution for:
1. Pure Schwarzschild spacetime
2. JMN1 interior matched to Schwarzschild exterior

Physics model (kept unchanged from the working notebook/script):
- Thin-disc flux profile from metric functions A(r), B(r)
- Specific energy E(r), angular momentum L(r), angular velocity Omega(r)
- Observed redshift factor z(r)
- Spectral luminosity integrand integrated over radius for each h*nu/(kT*)

Reference:
P. S. Joshi, D. Malafarina, D. Narayan,
Class. Quantum Grav. 31 (2014) 015002.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import matplotlib
import numpy as np
import sympy as sp
from scipy.integrate import quad

# Use a non-interactive backend so scripts run in headless/local/CI/Colab contexts.
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ==============================
# Global physical/model settings
# ==============================
M = 1.0
MDOT = 1.0
ROUT = 1e4
L_FLOOR = 1e-20

# Schwarzschild-only model parameters
RIN_SCH_ONLY = 6.0

# JMN1 + Schwarzschild matched model parameters
M0 = 0.25
RB = 2 / M0
RIN_JMN = 1e-10
RIN_SCH_MATCHED = RB

# Dimensionless frequency grid definitions
HVKT_SCH = np.logspace(-6, 2, 120)
HVKT_JMN = np.logspace(-6, 2, 80)

# Symbolic radius variable
R_SYM = sp.symbols("r", positive=True)


# ==============================
# Symbolic metric toolkit
# ==============================
def build_metric_functions(a_expr: sp.Expr, b_expr: sp.Expr) -> dict[str, callable]:
    """Return numerical functions for circular geodesic disc quantities.

    Given metric functions ``A(r)`` and ``B(r)``, this builds the standard thin-disc
    quantities used in the script:
    - Metric determinant factor ``g = -r^2 A B``
    - Specific energy ``E(r)``
    - Specific angular momentum ``L(r)``
    - Angular velocity ``Omega(r)``
    - Radial derivatives ``dOmega/dr`` and ``dL/dr``

    The expressions are intentionally unchanged from the original implementation.
    """
    apr_expr = sp.diff(a_expr, R_SYM)
    e_expr = sp.sqrt((2 * a_expr**2) / (2 * a_expr - R_SYM * apr_expr))
    l_expr = sp.sqrt((R_SYM**3 * apr_expr) / (2 * a_expr - R_SYM * apr_expr))
    omega_expr = sp.sqrt(apr_expr / (2 * R_SYM))

    domega_expr = sp.diff(omega_expr, R_SYM)
    dl_expr = sp.diff(l_expr, R_SYM)
    gdet_expr = -R_SYM**2 * a_expr * b_expr

    # Lambdify exactly once per expression.
    return {
        "A": sp.lambdify(R_SYM, a_expr, "numpy"),
        "E": sp.lambdify(R_SYM, e_expr, "numpy"),
        "L": sp.lambdify(R_SYM, l_expr, "numpy"),
        "Omega": sp.lambdify(R_SYM, omega_expr, "numpy"),
        "dOmega": sp.lambdify(R_SYM, domega_expr, "numpy"),
        "dL": sp.lambdify(R_SYM, dl_expr, "numpy"),
        "gdet": sp.lambdify(R_SYM, gdet_expr, "numpy"),
    }


# ==============================
# Radiative quantities
# ==============================
def redshift_factor(radius: float, funcs: dict[str, callable]) -> float:
    """Return gravitational/Doppler redshift ``z(r)`` for disc emission.

    Formula kept unchanged:
    z = 1 / sqrt(-(-A + Omega^2 r^2)) - 1
    """
    return 1 / np.sqrt(-(-funcs["A"](radius) + funcs["Omega"](radius) ** 2 * radius**2)) - 1


def make_flux_function(
    funcs: dict[str, callable],
    rin: float,
    mdot: float,
    cache_size: int,
):
    """Build cached thin-disc flux ``F(r)`` with unchanged integrand/prefactor.

    Flux expression used:
    F(r) = - mdot / (4*pi*sqrt(-g)) * dOmega/dr / (E - Omega*L)^2
           * integral_{rin}^{r} (E - Omega*L) * dL/dr dr
    """

    @lru_cache(maxsize=cache_size)
    def flux(radius: float) -> float:
        integral_val, _ = quad(
            lambda x: (funcs["E"](x) - funcs["Omega"](x) * funcs["L"](x)) * funcs["dL"](x),
            rin,
            radius,
            limit=200,
        )

        prefactor = -(mdot / (4 * np.pi * np.sqrt(-funcs["gdet"](radius))))
        denom = (funcs["E"](radius) - funcs["Omega"](radius) * funcs["L"](radius)) ** 2

        if denom <= 0 or not np.isfinite(denom):
            return 0.0

        return prefactor * (funcs["dOmega"](radius) / denom) * integral_val

    return flux


def luminosity_integrand(
    radius: float,
    hvkt: float,
    funcs: dict[str, callable],
    flux_func,
) -> float:
    """Return radial integrand for observed spectral luminosity.

    This implements the same expression as the original script:
    dL_{nu,inf}/dr proportional to
      dL_inf(r) * ((1+z)^4 hvkt^4 / F) / (exp((1+z)hvkt/F^(1/4)) - 1) * 1/r

    where
    dL_inf(r) = 4*pi*r*sqrt(-g)*E(r)*F(r).
    """
    flux_val = flux_func(radius)
    if flux_val <= 0:
        return 0.0

    z_plus_one = 1 + redshift_factor(radius, funcs)
    expo = np.exp((z_plus_one * hvkt) / flux_val**0.25) - 1
    d_l_inf = 4 * np.pi * radius * np.sqrt(-funcs["gdet"](radius)) * funcs["E"](radius) * flux_val

    return (
        (15 / np.pi**4)
        * d_l_inf
        * ((z_plus_one**4 * hvkt**4) / flux_val)
        / expo
        * (1 / radius)
    )


# ==============================
# Spectrum builders
# ==============================
def compute_schwarzschild_spectrum() -> np.ndarray:
    """Compute ``(log10(hv/kT*), log10(nu L_nu,inf))`` for Schwarzschild spacetime."""
    a_sch_expr = 1 - 2 * M / R_SYM
    b_sch_expr = 1 / a_sch_expr
    sch_funcs = build_metric_functions(a_sch_expr, b_sch_expr)

    flux_sch = make_flux_function(
        funcs=sch_funcs,
        rin=RIN_SCH_ONLY,
        mdot=MDOT,
        cache_size=2000,
    )

    spectrum = []
    for hvkt in HVKT_SCH:
        val, _ = quad(
            lambda rr: luminosity_integrand(rr, hvkt, sch_funcs, flux_sch),
            RIN_SCH_ONLY,
            ROUT,
            limit=100,
        )
        spectrum.append((np.log10(hvkt), np.log10(val)))

    return np.array(spectrum)


def compute_jmn_plus_schwarzschild_spectrum() -> np.ndarray:
    """Compute matched JMN1 + Schwarzschild spectral luminosity curve.

    Interior (JMN1) is integrated on ``[RIN_JMN, RB]``.
    Exterior (Schwarzschild) is integrated on ``[RB, ROUT]``.
    Their contributions are summed per frequency bin.
    """
    # JMN1 interior metric
    a_jmn_expr = (1 - M0) * (R_SYM / RB) ** (M0 / (1 - M0))
    b_jmn_expr = 1 / (1 - M0)
    jmn_funcs = build_metric_functions(a_jmn_expr, b_jmn_expr)

    # Schwarzschild exterior metric
    a_sch_expr = 1 - 2 * M / R_SYM
    b_sch_expr = 1 / a_sch_expr
    sch_funcs = build_metric_functions(a_sch_expr, b_sch_expr)

    flux_jmn = make_flux_function(
        funcs=jmn_funcs,
        rin=RIN_JMN,
        mdot=MDOT,
        cache_size=50000,
    )
    flux_sch = make_flux_function(
        funcs=sch_funcs,
        rin=RIN_SCH_MATCHED,
        mdot=MDOT,
        cache_size=50000,
    )

    spectrum = []
    for hvkt in HVKT_JMN:
        l_jmn, _ = quad(
            lambda rr: luminosity_integrand(rr, hvkt, jmn_funcs, flux_jmn),
            RIN_JMN,
            RB,
            limit=100,
        )
        l_sch, _ = quad(
            lambda rr: luminosity_integrand(rr, hvkt, sch_funcs, flux_sch),
            RIN_SCH_MATCHED,
            ROUT,
            limit=100,
        )

        l_total = max(l_jmn + l_sch, L_FLOOR)
        spectrum.append((np.log10(hvkt), np.log10(l_total)))

    return np.array(spectrum)


# ==============================
# Plotting/output
# ==============================
def save_figure5(spec_sch: np.ndarray, spec_jmn: np.ndarray) -> Path:
    """Plot and save Figure 5-style comparison to ``figures/figure5.png``."""
    repo_root = Path(__file__).resolve().parent.parent
    output_dir = repo_root / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "figure5.png"

    plt.figure()
    plt.plot(spec_sch[:, 0], spec_sch[:, 1], "--", label="Schwarzschild")
    plt.plot(spec_jmn[:, 0], spec_jmn[:, 1], "--", label="JMN + Schwarzschild")
    plt.xlabel("log10(hν / kT*)")
    plt.ylabel("log10(ν Lν,∞)")
    plt.grid(True, linestyle=":")
    plt.xlim(-5, 1)
    plt.ylim(-10, 0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    return output_path


def main() -> None:
    """Compute both spectra and save Figure 5 output."""
    spec_sch = compute_schwarzschild_spectrum()
    spec_jmn = compute_jmn_plus_schwarzschild_spectrum()
    output_path = save_figure5(spec_sch, spec_jmn)
    print(f"Saved figure to: {output_path}")


if __name__ == "__main__":
    main()
