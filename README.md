# ang2-a14-analysis

Analysis code for *One side-chain atom switches the kinetic step governing angulin paralog
selectivity at tight junctions* (Kuzu, Niwa, Kimura & Kondoh).

The code here covers the computational analyses described in the paper's Methods:
contact occupancy in the molecular dynamics (MD) trajectories, solvent-accessible surface
area (SASA) of position 14, peak separation force from the steered MD (SMD) pulls, the
Brønsted / Φ-value analysis of the measured binding kinetics, and the exact-permutation
Spearman tests. It also contains the scripts that assemble the figure source-data workbook,
and the scripts that generated the SMD starting conformations.

The data these scripts read are deposited separately:

> Zenodo, **10.5281/zenodo.24103964** — https://doi.org/10.5281/zenodo.24103964

Structure predictions were produced with the AlphaFold Server
(https://alphafoldserver.com); no local prediction code was written or used.

---

## Layout

```
scripts/
  extract/        trajectory  ->  compressed per-replica arrays (.npz)
  analysis/       arrays + measured kinetics  ->  the numbers reported in the paper
  source_data/    arrays + kinetics  ->  the figure source-data workbook
  smd_setup/      MD trajectory  ->  SMD starting conformations and run packages
```

## What produces what

| Methods section | script | reads | writes / reports |
|---|---|---|---|
| Trajectory analysis — contact occupancy, formed-to-broken transitions | `extract/occ3.py` | production trajectories | `occ3/<tag>_r<n>.npz` |
| …statistics on the above | `analysis/occ_stats.py` | `occ3/` | per-residue occupancy table |
| Trajectory analysis — Y99 occupancy and contacting atom-pair counts | `extract/y99.py`, `analysis/y99_stats.py` | trajectories, `y99/` | Y99 occupancy 1.000 with no transitions; atom-pair counts |
| Solvent-accessible surface area | `extract/sasa_if.py` | trajectories (calls `gmx sasa`) | `sasa/bound_*.npz`, `sasa/free_*.npz` |
| …opening frequency and time-averaged exposure | `analysis/open_freq.py` | `sasa/` | `res14_open.npz`; 15.2 % vs 41.2 % at angulin-3, 28.0 % vs 35.2 % at angulin-1; wild-type time-averaged exposure 4.6 % and 20.3 % |
| …the separation statistic quoted in the text | `analysis/open_gap.py` | `res14_open.npz` | the gap of 2.2 SEM at angulin-3 |
| …further SASA statistics | `analysis/sasa_stats.py` | `sasa/` | |
| Steered MD — peak separation force | `extract/force_prof.py`, `analysis/force_stats.py` | the 216 pull force traces | `force_prof.npz`; per-system peak force |
| Brønsted analysis; Statistical analyses | `analysis/kin3.py` | the measured rate constants (tabulated in the script) | geometric means, ΔΔG‡, ΔΔG_eq, Φ with Monte-Carlo error propagation, Spearman ρ with exact permutation *p* |
| receptor-side negative control (source-data sheets 9 and 14) | `extract/recv_side.py`, `source_data/addsheet14.py` | trajectories, `recv/` | |
| figure source data | `source_data/mkxlsx.py`, then `addsheet12.py`, `addsheet14.py` | all of the above | `FIGURE_DATA_2026-10-24.xlsx` |
| Steered MD — starting conformations | `smd_setup/frames_var.sh` (with `count_ndx.py`, `write_chains.py`, `pmf_prep.py`) | production trajectories | 18 conformations per complex, with the 3.0–5.0 nm chain-separation check |
| Steered MD — run packaging | `smd_setup/mkpull_var.sh` | the conformations above | one pull job per conformation |

`kin3.py` carries the three independent SPR determinations per system as literal values, so
it runs on its own with no input files and reproduces every kinetic number in the paper.
It is the quickest way to check the analysis.

## Environment

```
python 3.10.20
numpy 2.2.6
scipy 1.15.2
MDAnalysis 2.9.0
openpyxl 3.1.5
GROMACS 2024.3   (built from source; used by sasa_if.py and by the SMD setup scripts)
```

`environment.yml` reproduces the Python side:

```bash
conda env create -f environment.yml
conda activate ang2
python scripts/analysis/kin3.py
```

## Paths

The scripts were run in place and carry the absolute paths of the machine they ran on
(`/home/yuki/...`). They are left exactly as run rather than rewritten, so that the code
published here is the code that produced the reported numbers. To run them elsewhere,
edit the path constants at the top of each file. The correspondence to the Zenodo
deposit is:

| path in the scripts | Zenodo deposit |
|---|---|
| `md_project_<tag>/`, `md_replicas/runs/<tag>_r<n>/`, `variant_replicas/runs/<tag>_r<n>/`, `free_replicas/runs/mono_<t>_r<n>/` | `05_trajectories/<system>/` |
| `pmf/smd_an/<tag>/prot.tpr`, `chainA_prot.ndx` | `02_md_inputs/<system>/` and `05_trajectories/<system>/<system>_protein.tpr` |
| `vpull/runs/<tag>/r<n>_<T>ns/pullf.xvg`, `pullx.xvg` | `04_source_data/smd_pull_traces/<system>/` |
| `confs_var/<tag>/r<n>_<T>ns/` | starting conformations; the corresponding run inputs are `03_smd_inputs/<system>/` |

System tags: `L` = angulin-1 complex, `i` = angulin-3 complex, `mono_` = free Ang2;
`wt, g, c, s, t, v` = wild-type Ang2, A14G, A14C, A14S, A14T, A14V. Replica numbers in the
tags are not always 1–3; the three replicas deposited for each system are the three on which
every reported number is based.

## Not included

The project ran a number of analyses that are not reported in the paper and have no
corresponding Methods section. Their code is not included here: the free-energy perturbation
calculations, the protein–protein association steered MD, and the exploratory analyses that
were examined and set aside. Sheets 10 and 11 of the deposited source-data workbook come
from two of those analyses and are therefore not regenerated by the scripts here; sheets
1–9 and 12–14 are. The docstring of `analysis/open_gap.py` refers to `grouping.py`, one of
the set-aside scripts, as the origin of the separation statistic it computes; `open_gap.py`
itself is self-contained and needs nothing from that file.

## A note on the console output

Several scripts print their tables and verdicts in Japanese, which is how they were written
and run. The numbers, labels and column headings are language-independent; only the prose
around them is not.

## License

MIT — see `LICENSE`.
