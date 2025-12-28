![Ledger](docs/ledger-logo.png)

**Epistemic-First Cosmological Inference**

`ledger` is the reference implementation for the **Parochial by Construction (PbC)** research framework. It provides the tooling to audit cosmological datasets for evidence of "Context Coupling"—the hypothesis that the history of the universe is not a pre-existing territory, but a dependency network fixed by the act of observation.

## 1. Theory & Intention

Standard cosmological inference (Temporal Evolution, or **TE**) assumes the past is fixed and independent of the observer. We measure the universe "as it is."

**Causal Evolution (CE)** proposes an alternative: Time is a bookkeeping metric for causal dependencies. In this view, the "Past" is a resource that is consumed to explain the "Present." Consequently, the specific constraints of the observer (the **Context**) impose a tax on the history that can be reconstructed.

**The Ledger Protocol** tests this by formalizing two entities:
1.  **The Record ($d_{obs}$):** The data we observe (e.g., CMB temperature/polarization).
2.  **The Context ($T_{ctx}$):** The cost of observation (e.g., scan strategy, exposure depth, masking).

If the universe behaves according to standard General Relativity, the Record and the Context should be statistically independent ($S_\gamma \approx 0$). If Causal Evolution is true, the Record will effectively "phase-lock" to the Context to minimize the computational cost of the history, producing a non-zero coupling signal.

### 1.1 Historical Motivation: The "Axis of Evil"

The PbC framework provides a theoretical basis for long-standing CMB anomalies.  The observed alignment of the quadrupole ($\ell=2$) and octopole ($\ell=3$) with the solar system's Ecliptic plane—often called the "Axis of Evil"—is interpreted here not as a fluke, but as a primary signal of **Phase-Locking (P1)** between the primordial record and the observer's scan context.

## 2. Research Goals

This repository targets the diagnostics defined in *Parochial by Construction: Dual Constructions for Cosmological Inference* [Draft, 2025].

* **P1 (Polarization Phase-Locking):** Tests if primordial CMB polarization phases ($\phi_{TE/EE}$) align with the anisotropy of the Planck scan strategy.
* **P5 (Target Selection Commutator):** Uses DESI `altmtl` realizations to test if target selection bookkeeping commutes with clustering inference.
* **P6 (Field Variance):** (Upcoming) Tests for context-dependent variance in JWST galaxy number counts.

## 3. Repository Structure

```text
ledger/
├── data/
│   ├── raw/                # Large FITS maps (Planck/DESI) - Ignored by git
│   ├── processed/          # Generated Context Templates (T_ctx) & Audit Results
│   └── mocks/              # Low-res synthetic data for CI/Testing
├── docs/
│   └── tex/                # LaTeX sources for Research Programs (P1, P2)
├── notebooks/              # Jupyter notebooks for prototyping and visualization
├── scripts/
│   ├── download_planck.sh  # Fetcher for public NPIPE/PR3 data
│   ├── generate_mock.py    # Generates synthetic data for testing
│   ├── 01_context_builder.py # Step 1: Raw Maps -> Context Vector (c)
│   ├── 02_planck_audit.py    # Step 2: Measure Phase-Locking (S_gamma)
│   └── run_injection.py    # Validation: Sensitivity testing
├── src/
│   └── pbc/                # Core Python Package
│       ├── context.py      # Logic for building T_ctx from exposure maps
│       ├── dual.py         # Causal Evolution posterior math (Woodbury identities)
│       ├── stats.py        # Estimators for P1-P6 diagnostics
│       └── utils/          # I/O and HEALPix wrappers
├── tests/
│   └── science/            # Physics validation (Null tests, Symmetry checks)
└── pyproject.toml          # Package definition and dependencies
```

## 4. Installation
**Prerequisites:** Python 3.9+

```bash
# 1. Clone the repository
git clone [https://github.com/your-org/ledger.git](https://github.com/your-org/ledger.git)
cd ledger

# 2. Create a virtual environment (Recommended)
python -m venv .venv
source .venv/bin/activate

# 3. Install in editable mode
# Includes 'numpy', 'healpy', 'tqdm'
pip install -e ".[dev]"
```

## 5. Quick Start (Simulation Mode)

If you cannot download the 2TB Planck archive, use the **Simulation Mode**. This generates synthetic "Standard Model" skies and "Dipole" scan strategies to verify the pipeline logic.

1. **Generate Mock Data:** Creates synthetic FITS files in `data/raw/`.

```bash
python scripts/generate_mock_data.py
```

2. **Build the Context Template:** Extracts the "Cost Vector" (c) from the mock scan strategy.

```bash
python scripts/01_context_builder.py --nside 512
```

3. **Audit the Ledger:** Cross-correlates the Mock Record with the Mock Context.

```bash
python scripts/02_planck_p1_audit.py --nside 512
# Expected: S_gamma ~ 0.0 (Null Result)
```

## 6. Real Data Workflow (Discovery Mode)

To run the actual scientific analysis on Planck 2018/PR4 data:

1. **Download Data:** Fetches the necessary Frequency and Component maps (requires ~2GB).

```bash
./scripts/download_planck.sh
```

2. Run High-Res Analysis:

```bash
# Step 1: Build the Context (NSIDE 2048)
python scripts/01_context_builder.py --nside 2048

# Step 2: Audit the Record
python scripts/02_planck_p1_audit.py --nside 2048
```

## 7. Validation & Safety

We adhere to a strict "Stop-Loss" protocol. Before claiming a result, the pipeline must pass Null (Bias) and Injection (Sensitivity) tests.

* **Null Test:** Verify *Sγ* is consistent with 0 on random skies.

```bash
pytest tests/science/test_null.py
```

* **Sensitivity Test:** Verify we can recover a hidden signal (λ=0.05).

```bash
python scripts/run_injection.py --lambda-inj 0.05 --nsims 100
```

## License

MIT

