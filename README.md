# Accretion Disk Spectrum: JMN1 vs Schwarzschild

A compact, reproducible project to study the spectral luminosity distribution of accretion disks in:
- naked singularity spacetime (JMN1 interior matched to Schwarzschild exterior)
- black hole spacetime (Schwarzschild)

The implementation reproduces the Figure 5 style comparison from:

- P. S. Joshi, D. Malafarina, R. Narayan (2014)
- *Distinguishing black holes from naked singularities through their accretion disc properties*
- Class. Quantum Grav. **31** (2014) 015002

This repo computes and plots:
- Schwarzschild thin-disk spectrum
- JMN1 interior + Schwarzschild exterior matched spectrum

Output file:
- `figures/figure5.png`

---

## What this repo does

The main script [`src/spectrum.py`](src/spectrum.py) numerically evaluates the relativistic thin-disk flux and then integrates the observed spectral luminosity over radius.

The final plot is shown in dimensionless form:

$$
\log_{10}(\nu L_{\nu,\infty})
\quad \text{vs} \quad
\log_{10}\!\left(\frac{h\nu}{kT_*}\right)
$$

If math rendering is not enabled in your Markdown viewer, read it as:
- x-axis: `log10(h*nu/(k*T*))`
- y-axis: `log10(nu*L_nu,infinity)`

---

## Physics model (short version)

For a static spherical metric

$$
ds^2 = -A(r)dt^2 + B(r)dr^2 + r^2 d\Omega^2,
$$

the code uses standard circular-orbit quantities:

$$
E(r)=\sqrt{\frac{2A(r)^2}{2A(r)-rA'(r)}},
\quad
L(r)=\sqrt{\frac{r^3A'(r)}{2A(r)-rA'(r)}},
\quad
\Omega(r)=\sqrt{\frac{A'(r)}{2r}}.
$$

Flux is computed in the Novikov-Thorne thin-disk framework.

Metrics used:
- Schwarzschild: $A(r)=1-2M/r$, $B(r)=1/A(r)$
- JMN1 interior: $A(r)=(1-M_0)(r/R_b)^{M_0/(1-M_0)}$, $B(r)=1/(1-M_0)$, with $R_b=2/M_0$

Current default parameter in code:
- `M0 = 0.25`

---

## Quick start (terminal)

### 1) Clone the repo

```bash
git clone <your-repo-url>
cd <repo-folder>
```

### 2) Create and activate a virtual environment (recommended)

```bash
python -m venv .venv
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

### 3) Install dependencies

```bash
pip install -r requirements.txt
```

### 4) Run the spectrum script

```bash
python src/spectrum.py
```

### 5) Open the generated figure

Generated file:
- `figures/figure5.png`

---

## Run on Google Colab

### Option A: Clone directly from GitHub (recommended)

Create a new Colab notebook and run:

```python
!git clone <your-repo-url>
%cd <repo-folder>
!pip install -r requirements.txt
!python src/spectrum.py
```

Display the figure:

```python
from PIL import Image
Image.open('figures/figure5.png')
```

### Option B: Upload ZIP to Colab

1. Zip your local repo.
2. Upload it to Colab.
3. Extract and run:

```python
!unzip <repo-zip>.zip
%cd <repo-folder>
!pip install -r requirements.txt
!python src/spectrum.py
```

Then display:

```python
from PIL import Image
Image.open('figures/figure5.png')
```

---

## Project structure

```text
.
├── figures/
│   └── figure5.png
├── notebooks/
│   └── figure5_reproduction.ipynb
├── src/
│   └── spectrum.py
├── requirements.txt
└── README.md
```

---

## Reproducing with the notebook

If you prefer notebook workflow instead of running the script:

- Open `notebooks/figure5_reproduction.ipynb`
- Run all cells in order

---

## Notes

- The repository focuses on numerical reproduction of the spectrum plot.
- It does not include paper PDFs.
- For scientific context, please cite the original paper above.
