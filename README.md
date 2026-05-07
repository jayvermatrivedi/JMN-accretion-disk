# Accretion Disk Spectrum in JMN Naked Singularity Spacetime

## Project Overview
This repository provides a reproducible computational implementation of the spectral luminosity calculation presented in Figure 5 of:

P. S. Joshi, D. Malafarina, and R. Narayan, *Distinguishing black holes from naked singularities through their accretion disc properties*, Classical and Quantum Gravity **31** (2014) 015002.

The code computes and compares:
- The spectral luminosity distribution for a Schwarzschild spacetime.
- The spectral luminosity distribution for a matched JMN1 interior + Schwarzschild exterior spacetime.

The implementation is intended for graduate students and researchers working in relativistic astrophysics, compact objects, and accretion theory.

## Physics Background
Accretion disks around compact objects convert gravitational binding energy into radiation. The emitted spectrum depends on the underlying spacetime geometry because orbital energy, angular momentum, and redshift are metric-dependent.

This repository focuses on testing observational distinguishability between:
- A black hole spacetime (Schwarzschild), and
- A horizonless naked singularity model (JMN1 interior matched to Schwarzschild exterior).

The central observable is the dimensionless spectral luminosity profile plotted as
\[
\log_{10}\!\left(\nu L_{\nu,\infty}\right)
\quad\text{vs}\quad
\log_{10}\!\left(\frac{h\nu}{kT_*}\right).
\]

## JMN1 Metric (Interior)
The JMN1 interior metric used in the model is spherically symmetric with metric functions
\[
A(r) = (1-M_0)\left(\frac{r}{R_b}\right)^{\frac{M_0}{1-M_0}},
\qquad
B(r) = \frac{1}{1-M_0},
\]
with matching radius
\[
R_b = \frac{2}{M_0}.
\]

In the current setup, the parameter choice is
\[
M_0 = 0.25.
\]

The interior is matched to a Schwarzschild exterior at \(r=R_b\), enabling a piecewise spectrum construction from interior and exterior radial integrations.

## Schwarzschild Comparison (Exterior / Reference)
The Schwarzschild metric functions are
\[
A(r)=1-\frac{2M}{r},
\qquad
B(r)=\frac{1}{A(r)}.
\]

This spacetime is used in two ways:
- As the standalone black hole reference model.
- As the exterior region of the matched JMN1+Schwarzschild configuration.

The direct comparison between these spectra is the core scientific objective of the repository.

## Novikov-Thorne Thin Disk Formalism (Brief)
The implementation follows a standard relativistic thin-disk construction (Novikov-Thorne framework, specialized to static spherical metrics).

For circular equatorial geodesics, the code evaluates:
\[
E(r)=\sqrt{\frac{2A(r)^2}{2A(r)-rA'(r)}},
\]
\[
L(r)=\sqrt{\frac{r^3A'(r)}{2A(r)-rA'(r)}},
\]
\[
\Omega(r)=\sqrt{\frac{A'(r)}{2r}}.
\]

The flux profile is computed as
\[
F(r) = -\frac{\dot M}{4\pi\sqrt{-g}}\,\frac{\Omega_{,r}}{\left(E-\Omega L\right)^2}
\int_{r_{\rm in}}^{r}\left(E-\Omega L\right)L_{,r}\,dr,
\]
with
\[
g=-r^2A(r)B(r).
\]

The observed redshift factor is implemented as
\[
z(r)=\frac{1}{\sqrt{-\left(-A(r)+\Omega(r)^2r^2\right)}}-1.
\]

Using these ingredients, the code numerically integrates the spectral luminosity contribution over radius for each frequency bin.

## Installation
### Requirements
- Python 3.9+
- `numpy`
- `scipy`
- `sympy`
- `matplotlib`
- `jupyter`

Install with:

```bash
pip install -r requirements.txt
```

## Local Execution
From repository root:

```bash
python src/spectrum.py
```

On successful execution, the script generates:
- `figures/figure5.png`

The script is configured with a non-interactive Matplotlib backend, so it runs in terminal/headless environments and still saves the output figure.

## Google Colab Compatibility
This repository is runnable in Google Colab without code changes.

### Option 1: Upload repository archive
1. Upload the repository to Colab runtime.
2. Install dependencies:

```python
!pip install numpy scipy sympy matplotlib jupyter
```

3. Run:

```python
!python src/spectrum.py
```

4. Inspect generated file:

```python
from PIL import Image
Image.open('figures/figure5.png')
```

### Option 2: Clone from GitHub in Colab
```python
!git clone <your-repo-url>
%cd <repo-folder>
!pip install -r requirements.txt
!python src/spectrum.py
```

## Repository Structure
```text
.
├── figures/
│   └── figure5.png                # Generated spectrum figure
├── notebooks/                     # Optional exploratory notebooks
├── src/
│   └── spectrum.py                # Main computation and plotting script
├── requirements.txt               # Python dependencies
├── .gitignore                     # Scientific Python ignore rules
└── README.md
```

## Future Work
- JMN2 implementation for broader naked-singularity model comparison.
- Relativistic ray tracing from disk to observer (transfer functions, inclination effects).
- Extension to rotating compact objects and Kerr spacetime comparison.
- Photon sphere analysis and its imprint on spectral/observational diagnostics.

## Notes
- This repository does **not** include full-text paper PDFs.
- Scientific sources are cited in bibliography form below.

## Bibliography
1. P. S. Joshi, D. Malafarina, and R. Narayan, “Distinguishing black holes from naked singularities through their accretion disc properties,” *Classical and Quantum Gravity* **31** (2014) 015002.
2. I. D. Novikov and K. S. Thorne, “Astrophysics of black holes,” in *Black Holes (Les Astres Occlus)*, edited by C. DeWitt and B. S. DeWitt, Gordon and Breach, New York (1973), pp. 343-450.
