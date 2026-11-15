# Design study: alternatives evaluated against the finalized design

This note documents the alternatives that were implemented and benchmarked against the
finalized design, and what the numbers show — including the honest negative result that
the estimator is already near-optimal, and the trade-off that led to the feedforward
variant **not** being adopted. Everything here is reproducible from the repository
(`run_headtohead.m`, `run_rigor.m`, `run_band_sweep.m`).

The **finalized design is the default** in `turret_params.m`; every other variant is one
frozen configuration reachable by an override, evaluated with no per-case retuning.

| selector | values (default **bold**) | meaning |
|---|---|---|
| `P.est`  | **`CA`** · `adaptQ` · `IMM` | estimator / predictor |
| `P.ctrl` | **`lag`** · `lagff` · `pi` | lag compensator on the Tp-predicted reference; `lagff` = + feedforward that aims at the *current* estimate; `pi` = PI compensator |
| `P.acq`  | **`dualmode`** · `linear` | dual-mode acquisition; `linear` = LTI tracker throughout |

```matlab
simulate_case('C', struct('acq','linear'))     % baseline column
simulate_case('C', struct('ctrl','lagff'))     % feedforward variant (not adopted)
simulate_case('C', struct('est','IMM'))        % IMM estimator
```

---

## 1. Estimator and prediction

- **IMM {CV, CA, CT}** — implemented: constant-velocity, constant-acceleration and a
  5-state coordinated-turn EKF (turn rate as a filtered state), blended by Markov-mixed
  model probabilities.
- **Adaptive-Q (NIS-triggered)** — implemented: a single CA filter that raises the jerk
  PSD when the normalized innovation squared exceeds a chi-square gate.
- **PCRLB + NEES/NIS** (`run_rigor.m`, Part A) — implemented. On a maneuver-plus-loss
  trajectory (Monte-Carlo M = 300): filter RMS = **0.945 x sqrt(PCRLB)** (near-optimal),
  mean NIS = 1.87 (2 dof), mean NEES = 3.65 (n_x = 6, deliberately conservative).

**Result.** Neither IMM nor adaptive-Q beats the single CA filter, and the PCRLB explains
why: the CA filter is already within ~5% of the information floor, so there is no
estimation headroom to recover. On the stress suite the IMM is clearly worse (RMS miss
2.06 m vs 1.22 m, in-turn-loss prediction error 6.40 m vs 4.26 m).

## 2. Controller and acquisition

- **Finalized: lag-matched compensator** `Gc(s) = Kp/(1+s/p)` on the Tp-predicted
  reference, with `Kp = 1/(Tp*Kdc)` and `p = 2/Tp` derived so the closed-loop group delay
  equals Tp (predictor lead and loop lag cancel). Margins: PM 65.2 deg, GM 38.0 dB.
- **PI compensator** — evaluated for comparison (`design_analysis.m`): type-2, so zero
  ramp lag; combined with the predictor it *leads* the target by Tp (111 m at 800 m in a
  parabolic-angle test), and its gain-margin frequency sits on the plant resonance.
- **Dual-mode acquisition** — implemented: a max-rate slew while the pointing error is
  large, then a bumpless hand-off to the linear tracker inside a band. `run_band_sweep.m`
  checks the 2.5 deg band: 1.0 deg chatters (468 slew re-entries), 1.5 deg re-enters 16
  times, 2.0 deg and above are clean.
- **Feedforward variant (`lagff`)** — implemented, *not adopted.* The model-inverse
  `Vff = thetadot/Kdc + (Rm*Jeq/(N*Kt))*thetaddot` is derived from the plant, not tuned.
  Controller-only (`run_rigor.m`, Part B; true kinematics fed in; peak LOS rate 0.08 rad/s):

  | LOS freq [rad/s] | 0.2 | 0.3 | 0.5 | 0.7 | 1.0 | 1.4 | 2.0 |
  |---|---|---|---|---|---|---|---|
  | lag only, RMS pointing [deg] | 0.0014 | 0.0064 | 0.0176 | 0.0338 | 0.0673 | 0.1282 | 0.2407 |
  | with feedforward [deg] | 0.0003 | 0.0003 | 0.0005 | 0.0008 | 0.0011 | 0.0019 | 0.0030 |

  **Why it is not the finalized design.** The feedforward cancels the loop delay, so it
  must be driven by the *current* line of sight and the Tp-ahead prediction is no longer
  in the control path. The brief's architecture requires the turret to track the predicted
  elevation reference. An earlier version of this project had this variant as the default
  while the paper described the predicted reference driving the loop; the mismatch was
  found by comparing the logged reference with `atan2` of the predicted state, and the
  lag design was made the default.

## 3. Head-to-head (`results/headtohead.csv`)

One frozen design per column, no retuning. Base = lag compensator + linear acquisition;
Final = + dual-mode; FF = Final + feedforward variant; aQ / IMM = estimator swapped in Final.

| metric (mean) | suite | Base | **Final** | FF variant | aQ | IMM |
|---|---|---|---|---|---|---|
| Success rate [%]       | Orig-40   | 100 | **100** | 100 | 100 | 100 |
| Interception time [s]  | Orig-40   | 2.87 | **2.67** | 2.63 | 2.67 | 2.67 |
| RMS pointing [deg]     | Orig-40   | 0.098 | **0.074** | 0.073 | 0.075 | 0.089 |
| RMS miss [m]           | Orig-40   | 0.98 | **0.75** | 0.59 | 0.79 | 1.28 |
| Success rate [%]       | Stress-24 | 100 | **95.8** | 100 | 95.8 | 100 |
| Interception time [s]  | Stress-24 | 3.74 | **3.91** | 3.55 | 3.98 | 4.40 |
| RMS pointing [deg]     | Stress-24 | 0.110 | **0.103** | 0.088 | 0.094 | 0.120 |
| RMS miss [m]           | Stress-24 | 1.37 | **1.22** | 0.89 | 1.18 | 2.06 |

The dual-mode front-end improves interception time and miss distance on the 40 original
cases (-7% and -23%). The finalized design is weaker than the feedforward variant in two
respects, both reported in the paper: it intercepts 23/24 stress cases, and in one of the
40 original cases (seed 116) the beam drifts 25-36 m off the target in the last two
seconds after four hits (so that case has no defined acquisition time and is excluded
from the post-acquisition statistics). Both are consistent with prediction error over the
0.5 s horizon, since the variant that does not use the prediction does not show them.

## 4. How to reproduce

```matlab
run_headtohead     % -> results/headtohead.{csv,tex,mat}, figures/fig_headtohead*
run_rigor          % -> figures/fig_pcrlb.*, fig_consistency.*, fig_feedforward.*, results/rigor.txt
run_band_sweep     % -> results/band_sweep.txt
run_scenarios ; run_final_evaluation
```

`headtohead_figure.m` redraws the ablation figures from `results/headtohead.mat` without
re-simulating.
