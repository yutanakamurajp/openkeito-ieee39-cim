# IEEE 39-bus system with a WECC generic PV plant, from CIM to OpenKeito

This repository reproduces the IEEE 39-bus results of the paper

> Y. Nakamura, "OpenKeito: An Open Electromagnetic Transient Simulator for People and AI Agents That Starts
> From Power-Flow and CIM Data," *IEEE Access* (submitted).

It converts the CIM description of the IEEE 39-bus test system published with
[emthub / emtiop](https://github.com/temcdrm/emthub) into [OpenKeito](https://github.com/yutanakamurajp/OpenKeito)
models, runs them with OpenKeito, and makes the OpenKeito side of Figs. 14–18 and Tables 8–10 of the paper,
with every number quoted in the text (`figures/numbers.json`).

OpenKeito: https://github.com/yutanakamurajp/OpenKeito, DOI [10.5281/zenodo.24523022](https://doi.org/10.5281/zenodo.24523022).

## Contents

| Path | Content |
|---|---|
| `cim_to_openkeito.py` | The converter: reads `IEEE39.xml` and `IEEE39_ic.xml` with emtiop and writes one OpenKeito model per case and the models of the sequence-impedance scans |
| `cim/` | The input CIM files from emthub (unchanged, Apache-2.0, see `cim/NOTICE`) |
| `models/*.json` | The eight generated models (model format `"0.7"`) |
| `models/scan/` | Models for `openkeito scan` at 60 Hz and the impedances expected from the CIM (`expected_ohm.json`) |
| `run_cases.py` | Runs the eight cases and the scans with OpenKeito into `runs/` (Python standard library only) |
| `make_figures.py` | Figures, tables and numbers of the paper from `runs/` into `figures/` |
| `figures/` | The output of `make_figures.py` for OpenKeito v0.11.0, as used in the paper |
| `requirements.txt` | Python packages and the versions used |
| `CITATION.cff`, `LICENSE-MIT`, `LICENSE-APACHE` | Citation and licenses |

## Reproduce

1. **OpenKeito v0.11.0 or later.** The command `openkeito` is built from source with a
   [Rust toolchain](https://rustup.rs) (1.85 or later):

   ```sh
   git clone --branch v0.11.0 https://github.com/yutanakamurajp/OpenKeito.git
   cd OpenKeito
   cargo build --release -p openkeito-cli      # the command: target/release/openkeito
   python -m pip install .                     # the Python package (optional, not needed here)
   ```

   The Python package (`pip install .` in the OpenKeito folder) gives the same solver to Python, but the
   scripts here call the command `openkeito`.

2. **Python 3.10 or later** with the packages of `requirements.txt` (the results in `figures/` were made with
   Python 3.11 and these versions):

   ```sh
   python -m pip install -r requirements.txt
   ```

3. **(Optional) regenerate the models.** `models/` holds the converter's output; running it again gives the
   same files:

   ```sh
   python cim_to_openkeito.py cim/IEEE39.xml cim/IEEE39_ic.xml models
   ```

4. **Run the cases** (about 3 minutes on one desktop computer) and **make the figures**:

   ```sh
   python run_cases.py --openkeito /path/to/OpenKeito/target/release/openkeito
   python make_figures.py
   ```

   `run_cases.py` also takes the command from the environment variable `OPENKEITO` or the PATH. Every run must
   end with `status` `"completed"`; the warnings W183 (control blocks start from their own initial values) and
   W220 (the start is not exactly the phasor solution, because the PV references ramp up from 0.01 s) are
   expected.

5. **Compare.** `figures/manifest.csv` lists the sha256 of the CIM inputs and of every output. With OpenKeito
   v0.11.0 the figures (PDF), the CSV files, the tables and `numbers.json` are the same files as in this
   repository (checked with v0.11.0 on Windows; the PDFs also depend on the matplotlib version).

The computing time is not part of `figures/`. `python run_cases.py --timing` runs `nofault`, `slg16` and
`slg16_20s` three times each on one logical core with `--threads 1` and writes `runs/timing.json`; then
`make_figures.py` adds `openkeito_checks.compute_time` to `numbers.json`. For the paper (Intel Core Ultra 9
285K, Windows, one logical core, minimum of three runs) these were 8.3 s for `nofault` (6 s), 5.5 s for
`slg16` (3.5 s) and 33.4 s for `slg16_20s` (20 s), with a time step of 10 µs.

## Cases

| Case | Content |
|---|---|
| `nofault` | No fault, 6 s. The voltage PI of REEC_A starts at the bus voltage |
| `nofault_qreg0` | The same with the voltage PI starting at 0 (it is then held at its lower limit) |
| `slg16` | Single line to ground, bus 16 phase C, 0.001 Ω, 2.0–2.2 s (3.5 s) |
| `tpf02` | Three-phase fault to ground, bus 2, 3.0–3.1 s (4 s) |
| `slg16_20s` | 20 s with the fault of the emtiop example `IEEE39.prm` (bus 16 phase C, 12.0–12.2 s), voltage PI from 0, transformer magnetizing branches |
| `nofault_conditional_pi` | Sensitivity: PI blocks with OpenKeito's default conditional integration instead of `"anti_windup": "state"` |
| `slg16_20s_mag_transformer`, `slg16_20s_mag_shunt` | Sensitivity: magnetizing inside `three_phase_transformer`, or as shunts included in the phasor solution |

All cases: 10 µs, `2s-dirk`, steady-state initialization from the CIM initial conditions.

## Outputs and the paper

| File | Paper |
|---|---|
| `f14_oneline.pdf` | Fig. 14, one-line diagram drawn from the CIM diagram coordinates (overlapping buses moved as listed in `make_figures.py`) |
| `f15_steady_state.{pdf,csv}` | Fig. 15, phasor solution against the CIM initial conditions |
| `f16_pv_startup.{pdf,csv}` | Fig. 16, PV plant start-up (`nofault`) |
| `f17_fault_current.{pdf,csv}`, `f17_fault_current_cycles.csv` | Fig. 17, fault current of `slg16` and its AC and DC parts per cycle |
| `f18_response_20s.{pdf,csv}` | Fig. 18, PV reactive power and G35 speed over 20 s (`slg16_20s`) |
| `t08_mapping.tex` | Table 8, CIM classes and OpenKeito components |
| `t09_sequence_impedance.tex` | Table 9, sequence impedances of `openkeito scan` against the CIM |
| `t10_sensitivity.tex` | Table 10, sensitivity to the PI windup form and the magnetizing branches |
| `numbers.json` | The numbers of the text (with `openkeito_checks.definitions` for the derived indices) |
| `manifest.csv` | OpenKeito version and sha256 of the inputs and outputs |

Some figures and a table of the paper also show results of MATLAB/Simscape Electrical for comparison. Those
lines and numbers are not reproduced here: this repository holds neither the Simscape model nor its results,
only the OpenKeito side.

## How the CIM is represented

The module docstring of `cim_to_openkeito.py` lists every choice with its public source. In short:

* lines: `bergeron_line` with the positive- and zero-sequence parameters of the CIM;
* transformers: `three_phase_transformer`, the end-1 tap on the end-1 rated voltage (impedance seen from end 1 is
  (r + jx)·tap², the power-flow convention);
* loads: constant impedance at nominal voltage;
* synchronous machines: IEEE Std 1110 model 2.2 with the CIM open-circuit time constants; H = 3 s and the Xq''
  rule of emtiop's machine export (`AppendMachineDynamics`); the swing machine is a Thevenin source behind
  Ra + jXd';
* exciter and governor: the structure of the emtiop example `syncmach.pch` (ST1A-type exciter, KA from the CIM
  IEEET1; governor 1/R and Pmax from the CIM TGOV1);
* PV plant: WECC generic REGC_A, REEC_A and REPC_A (WECC Renewable Energy Modeling Task Force, 2014) with the
  parameters of the CIM dynamics as mapped in the emtiop example `ibr.pch`; a three-phase controlled current
  source, initializing source opened in all phases at 20 ms; limited PI controllers hold their integral
  (`"anti_windup": "state"`);
* magnetizing branches: linear R // L at the end-1 bus, excluded from the phasor solution
  (`"initialization": "exclude"`) so that the machines start from the CIM power flow.

Two treatments are not fixed by the data: the windup form of the limited PI controllers and how the
magnetizing branches enter the initial state. The sensitivity cases show their effect (Table 10).

## Inputs

`cim/IEEE39.xml` and `cim/IEEE39_ic.xml` are copied unchanged from the emthub repository,
https://github.com/temcdrm/emthub, folder `instances/` at commit `86d3aaf`:

| File | sha256 |
|---|---|
| `IEEE39.xml` | `45e472774bea6f8249c6637a9a32305f3c84563f7eae143b1fbdf549b8177ea6` |
| `IEEE39_ic.xml` | `cb7c26abb99b923a0ba3a72fde261a94cbe1312dcbba1b5041adcbb715d97a43` |

They can also be downloaded from
`https://raw.githubusercontent.com/temcdrm/emthub/86d3aaf/instances/IEEE39.xml` (and `IEEE39_ic.xml`).
The CIM is read with emtiop 1.0.1 (`pip install emtiop==1.0.1`), the PyPI package of emthub; the converter
relies on the table names of its `load_emt_dict`.

## Licenses

* The code and the models of this repository: MIT OR Apache-2.0 (`LICENSE-MIT`, `LICENSE-APACHE`), as OpenKeito.
* `cim/IEEE39.xml`, `cim/IEEE39_ic.xml`: emthub, Apache-2.0 (`cim/NOTICE`). The generated models contain data of
  these files; keep that notice when redistributing them.
* emtiop (Apache-2.0) is used as a library. The formulas named after its functions and example files were
  re-implemented here, not copied.

## Citation

See `CITATION.cff`. Please cite the paper above and OpenKeito (DOI 10.5281/zenodo.24523022).
