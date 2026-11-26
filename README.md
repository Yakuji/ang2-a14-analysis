# ang2-a14-analysis

Analysis code for *Steric bulk at one position of a tight junction modulator switches the kinetic
step governing angulin paralog selectivity* (Kuzu, Niwa, Kimura & Kondoh).

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

## Reading the output

Five scripts print their section headings, column names and verdicts in Japanese, which is
how they were written and run: `analysis/kin3.py`, `analysis/occ_stats.py`,
`analysis/y99_stats.py`, `analysis/sasa_stats.py` and `analysis/force_stats.py`. The code
itself — variable names, comments and the values printed — is in English or
language-independent. The scripts are published exactly as run rather than retranslated, so
this section gives the translation instead.

`analysis/open_freq.py` and `analysis/open_gap.py`, the extraction scripts and the
source-data scripts print in English and need nothing from this section.

### Glossary

| Japanese | English |
|---|---|
| 系 | system (one protein at one receptor) |
| 残基 / 名前 | residue / residue name |
| 対 | versus |
| 実測 | measured (as opposed to read off a figure) |
| 図読み | read off a published figure |
| 幾何平均 | geometric mean |
| 信頼区間 | confidence interval |
| 完全順列 / 順列下限 | exact permutation (test) / the floor of the exact permutation *p* |
| 単調・重なりなし | monotonic, with no overlap between groups |
| 重なりあり | the groups overlap |
| 全幅 | full range, highest ÷ lowest |
| 群 / 群間の比 | group (by γ-atom count) / ratio between adjacent groups |
| 判定 | verdict |
| 遷移状態で形成済み | already formed at the transition state (Φ ≈ 1) |
| 遷移状態では未形成 | not yet formed at the transition state (Φ ≈ 0) |
| 部分的 | partial (0.3 < Φ < 0.7) |
| KA がほぼ動かず Phi 不定 | K_A barely changes, so Φ is undefined |
| 変動する残基なし | no residue varies between systems |
| 誤差を伝播（モンテカルロ N 回） | error propagated by Monte Carlo, N draws |
| 会合障壁 / 複合体の不安定化 | association barrier / destabilization of the complex |
| 占有率 | contact occupancy |
| 両方 / 両受容体 | both / both receptors |
| ILDR2のみ / LSRのみ | angulin-3 only / angulin-1 only |
| 計 N 残基 | N residues in total |
| 補正 | correction (as in Holm correction) |
| 変異体 / 野生型 | mutant / wild type |
| 件 | count of items |
| 超 / 以上 | greater than / at least |
| 最良 | best (largest \|ρ\|) |
| 検定 N 残基 | N residues tested |
| 有意 | significant |
| 全系同一で除外 | excluded because the value is identical in every system |
| レプリカSD 対 幅 … 比 | replica-to-replica SD versus the spread between systems, and their ratio |
| ばらつきが覆う / 覆う | the scatter swamps the effect |
| 分離できている / 分離 | the systems separate cleanly |
| 順序は読めない | the ordering cannot be resolved |
| 断続性 | intermittency |
| 形成 / 破断の遷移回数 | number of formed-to-broken transitions |
| フレーム中 | out of N frames |
| 遊離状態 / 結合 | free state / bound |
| 側鎖 / 残基全体 | side chain / whole residue |
| 露出 / 露出率 | solvent exposure / fraction exposed |
| 埋没 / 埋没面積 | buried / buried area |
| 界面18残基 | the 18 interface residues |
| 重原子 / 重原子対 | heavy atom / heavy-atom pair |
| 最接近距離 | closest approach distance |
| 接触ペア数 | number of contacting pairs |
| 単独の余裕 | margin of that system on its own |
| 他5系で同じ計算をしたときの順位 | its rank when the same calculation is done for the other five systems |
| スイッチ | on–off switch of a contact |
| 飽和している | saturated |
| 平均 / 幅 / 比 | mean / spread / ratio |
| レプリカ間 | between replicas |
| なし | none |

### `analysis/kin3.py`

Section headings:

```
1. 実測 n=3 の幾何平均と 95% 信頼区間（対数正規として ln k 上で推定）
      Geometric means and 95% confidence intervals of the measured n = 3 determinations,
      estimated on ln k as log-normal
2. 図読み値との比較（これまでの全解析が使っていた値）
      Comparison with the figure-read values that the earlier analyses used
3. gamma 原子数との関係（実測値で再検定）
      Relation to the gamma-atom count, re-tested on the measured values
4. Phi 値 — 実測の誤差を伝播（モンテカルロ 200000 回）
      Phi values, with the measured error propagated by Monte Carlo, 200,000 draws
5. gamma 原子 1 個あたりの障壁上昇（ILDR2）と複合体の不安定化（LSR）
      Barrier increase per gamma atom at angulin-3, and destabilization of the complex
      at angulin-1
```

Example rows, with the Japanese parts translated:

```
対 ka  rho=-0.956 p=0.0222  単調・重なりなし  全幅 251倍
   vs ka, rho = -0.956, p = 0.0222, monotonic with no overlap, full range 251-fold
   群: 0個 2.87e+05 3.379e+05   1個 1.107e+04 7445   2個 2721 1344
      groups: 0 atoms ... | 1 atom ... | 2 atoms ...
   群間の比: 25.9倍 / 2.7倍
      ratios between adjacent groups: 25.9-fold then 2.7-fold

系    γ   ddG+ (kJ/mol)   ddG0 (kJ/mol)   Phi [68%区間]        判定
system                                         [68% interval]  verdict
A14C  1   +8.1 +- 1.7     +9.0 +- 2.3     0.89 [0.75, 1.09]    遷移状態で形成済み
                                                               already formed at the TS
A14G  0   -0.4 +- 1.0     -4.9 +- 1.3     0.08 [-0.14, 0.24]   遷移状態では未形成
                                                               not formed at the TS
```

`ddG+` is ΔΔG‡ and `ddG0` is ΔΔG_eq.

### `analysis/occ_stats.py`

```
1. 接触占有率が 0.2 を超える Ang-2 残基（野生型、3レプリカ平均）
      Ang2 residues whose contact occupancy exceeds 0.2 (wild type, mean of three replicas)
   columns:  残基 = residue number, 名前 = residue name,
             ILDR2 占有率 / LSR 占有率 = occupancy at angulin-3 / angulin-1,
             両受容体 = which receptors, with 両方 = both, ILDR2のみ / LSRのみ = only one
   計 23 残基  =  23 residues in total

2. 変異体と野生型の占有率の差（Holm 補正、受容体ごと 5変異体 x 23残基）
      Mutant-minus-wild-type occupancy differences (Holm correction, per receptor,
      5 mutants x 23 residues)
   3 SE 超 3 / 115、Holm p<0.05 0 件
      3 of 115 exceed 3 SE; none reaches Holm-corrected p < 0.05

3. 占有率 対 速度定数（6系、順列下限 2/720 = 0.0028）
      Occupancy versus the rate constants (six systems; the exact permutation floor is
      2/720 = 0.0028)
   最良 残基17 rho=+0.89 ... レプリカSD 0.189 対 幅 0.341 比 0.55 ばらつきが覆う
      best is residue 17 ...; replica SD 0.189 against a between-system spread of 0.341,
      ratio 0.55 — the scatter swamps the effect
   検定 22 残基、Holm有意 0 件、全系同一で除外 1 残基 [99]
      22 residues tested, none significant after Holm, 1 residue excluded because its value
      is identical in every system

4. 断続性: 形成/破断の遷移回数（野生型、1001フレーム中）
      Intermittency: number of formed-to-broken transitions (wild type, out of 1,001 frames)
```

### `analysis/y99_stats.py`

```
Tyr99 — 複合体MD、重原子 4.5 A 以内、レプリカ単位
   Tyr99 in the complex simulations, heavy atoms within 4.5 A, replica as the unit
（参考: 接触占有率は12系×3レプリカすべてで 1.000、スイッチ 0 回。飽和している）
   for reference: contact occupancy is 1.000 in all 12 systems x 3 replicas, with zero
   on-off switches — it is saturated
Tyr99 と受容体の接触ペア数 [重原子対]
   number of Tyr99-receptor contacting pairs, in heavy-atom pairs
Tyr99 が接触する受容体残基の数 [残基]
   number of receptor residues Tyr99 contacts
Tyr99 の受容体への最接近距離 [A]
   closest approach of Tyr99 to the receptor, in angstrom
レプリカSD / 変異体間の幅 = 0.56  <<< 0.5 超、順序は読めない
   replica SD divided by the spread between systems is 0.56 — above 0.5, so the ordering
   cannot be resolved
A14G 単独の余裕 -0.12 SE   他5系で同じ計算をしたときの順位 6 / 6
   A14G on its own has a margin of -0.12 SE, and ranks 6 of 6 when the same calculation is
   done for the other five systems
```

The first two lines are the statement quoted in the paper: Y99 is in contact in every frame
of every system, with no on–off transitions.

### `analysis/sasa_stats.py`

```
1. 遊離状態での位置14の溶媒露出（nm^2、3レプリカ）
      Solvent exposure of position 14 in the free state (nm^2, three replicas)
   columns:  側鎖 SASA = side-chain SASA, 残基全体 SASA = whole-residue SASA,
             側鎖の露出率 = fraction of the side chain exposed
   gamma原子数 対 側鎖SASA : rho=+0.717  完全順列 p=0.1778
      gamma-atom count versus side-chain SASA, exact permutation p
   レプリカSD平均 0.0044 対 変異体間の幅 0.0742 比 0.06 分離できている
      mean replica SD against the between-system spread, ratio 0.06 — they separate cleanly
   側鎖が露出しているか: 最小 wt 0.0189 nm^2、最大 A14T 0.0742 nm^2
      is the side chain exposed: smallest is wild type, largest is A14T

2. 遊離状態の SASA 対 速度定数（6系、順列下限 2/720 = 0.0028）
      Free-state SASA versus the rate constants
   row labels:  位置14 側鎖 = position 14 side chain, 位置14 残基全体 = whole residue 14,
                界面18残基 = the 18 interface residues
   trailing word:  分離 = the systems separate, 覆う = the scatter swamps the effect

3. 結合で埋まる面積（遊離 - 複合体、nm^2）
      Area buried on binding (free minus complex, nm^2)
   埋没面積 対 速度定数、および野生型との差（Holm 補正、受容体ごと 5変異体 x 3量）
      Buried area versus the rate constants, and the difference from wild type
      (Holm correction, per receptor, 5 mutants x 3 quantities)
   row labels:  側鎖 = side chain, 残基 = whole residue, 界面 = interface
```

### `analysis/force_stats.py`

```
系ごとのレプリカ平均（3レプリカの平均 ± レプリカ間SD）
   Per-system replica averages (mean of three replicas +- SD between replicas)
対 ka : / 対 kd :
   versus ka / versus kd
各変異体 対 野生型（SE倍、レプリカ間ばらつきから）
   Each mutant against wild type, in multiples of SE derived from the replica-to-replica
   scatter
3 SE 以上:
   those reaching at least 3 SE:
なし
   none
```

Column names are already in English: `Fmax` is the peak separation force, the quantity
reported in Figures 7F and 7G.

## License

MIT — see `LICENSE`.
