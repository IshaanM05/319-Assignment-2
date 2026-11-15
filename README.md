# Predictive Tracking and Turret Elevation Control for a Maneuvering Aerial Target

ME319 Control Systems (Autumn 2026, S3), IIT Bombay — Assignment 2.
Author: Ishaan Mondal. Repository: https://github.com/IshaanM05/319-Assignment-2

One fixed estimator/predictor/controller design (never retuned) for the elevation axis
of the Guardian point-defense turret: a constant-acceleration (CA) Kalman filter in
Cartesian coordinates, a fixed 0.5 s kinematic prediction that gives the predicted
elevation reference, and a conventional LTI lag compensator `Gc(s) = Kp/(1+s/p)` that
tracks it. `Kp` and `p` are derived so the closed-loop group delay equals the prediction
horizon (predictor lead and loop lag cancel, so the beam lies on the target's current
line of sight when a pulse is issued). A bounded dual-mode acquisition front-end handles
the initial slew. The estimator is *justified*, not assumed: a posterior Cramér–Rao
bound analysis shows it is within ~5% of the information-theoretic floor, so IMM and
adaptive-process-noise alternatives (both implemented) give no gain and are not adopted.
A model-inverse feedforward variant is also implemented and benchmarked, but is **not**
the finalized design because it aims at the current estimate and so takes the prescribed
prediction out of the control path. Full derivation, Bode/margin analysis, scenario
results, a 40-case evaluation sweep, a 24-case stress suite and an ablation are in
[`paper/paper.pdf`](paper/paper.pdf) (source: [`paper/paper.tex`](paper/paper.tex)).

**Headline result:** 100% interception across the three standardized scenarios and across
40 randomized bounded maneuvers with randomized measurement-loss intervals — one fixed
design, no retuning. On the 40-case suite: mean post-acquisition RMS miss 0.75 m, mean
interception time 2.67 s. On the harder 24-case stress suite the design intercepts 23/24
(95.8%). **Known limitations** (reported in the paper): in one of the 40 cases (seed 116)
the beam drifts 25-36 m off the target in the last two seconds after four hits, so that
case is excluded from the post-acquisition statistics (acquisition time undefined); and
the feedforward variant, which does not use the prediction, reaches 100% on both suites
with 0.59 m / 0.89 m RMS miss.

## Repository layout

```
turret_plant.m              Elevation-drive ODEs, states [theta_t, theta_dot_t, i]
turret_controller.m         CA-KF + 0.5 s prediction + LTI lag controller on the predicted reference
                            + dual-mode acquisition + engagement-pulse logic
target_motion_generator.m   Scenarios A/B/C + randomized + stress generators (20 Hz radar)
turret_params.m             Single source of truth; the defaults ARE the finalized design
simulate_case.m             Closed-loop RK4 simulation driver + metrics
design_analysis.m           Plant derivation check, Bode/margins, anti-windup demo
run_scenarios.m             Reproduces all Scenario A/B/C(+C2) figures and tables
run_final_evaluation.m      Reproduces the 40-case evaluation table/figure
run_headtohead.m            Ablation: baseline -> +dual-mode (final), feedforward variant, estimator swaps
run_band_sweep.m            Sensitivity of the dual-mode hand-off band
headtohead_figure.m         Draws the ablation figures from the saved head-to-head results
run_rigor.m                 PCRLB optimality + NEES/NIS consistency + controller-only feedforward study
paper_figure_style.m        Forces print-ready (white, black-text) MATLAB figures
IMPROVEMENTS.md             Design study: research menu, what was built, findings
figures/                    All generated .png/.pdf figures used in the paper
results/                    All generated .csv/.txt/.mat numeric results
paper/paper.tex, paper.pdf  The IEEE-format paper (6 pages + 1-page AI appendix)
```

## Design selectors and ablation

The **defaults in `turret_params.m` are the finalized design**. Every element is a
selectable flag, so each ablation column is one frozen design (no per-case retuning):

| selector | values (default **bold**) | meaning |
|---|---|---|
| `P.est`  | **`CA`** · `adaptQ` · `IMM` | estimator: CA, NIS-adaptive Q, or CV/CA/coordinated-turn IMM |
| `P.ctrl` | **`lag`** · `lagff` · `pi` | lag compensator on the Tp-predicted reference; `lagff` = + feedforward, aims at the *current* estimate (comparison only); `pi` = PI compensator |
| `P.acq`  | **`dualmode`** · `linear` | bounded dual-mode acquisition; `linear` = LTI tracker throughout |

```matlab
L = simulate_case('C');                                       % the finalized design
L = simulate_case('C', struct('acq','linear'));               % baseline: linear acquisition
L = simulate_case('C', struct('ctrl','lagff'));               % feedforward variant (not adopted)
L = simulate_case(struct('stress',3));                        % a high-rate-turn stress case
run_headtohead      % -> results/headtohead.{csv,tex,mat}, figures/fig_headtohead*
run_rigor           % -> results/rigor.txt, figures/fig_pcrlb.*, fig_consistency.*, fig_feedforward.*
run_band_sweep      % -> results/band_sweep.txt
```

**Why each element.** Starting from the linear baseline (lag compensator, linear
acquisition: 100% on the 40 cases, 0.98 m mean RMS miss), the **dual-mode** front-end cuts
the interception time from 2.87 s to 2.67 s and the RMS miss to **0.75 m** (-23%), with
100% interception preserved on the 40 cases (stress suite: 23/24). Swapping the CA estimator for
**adaptive-Q or IMM does not help** — the PCRLB shows CA is already within ~5% of the
information bound — so the single-model CA filter is kept. The **feedforward variant**
(`lagff`) is more accurate (0.59 m / 0.89 m, 100% on both suites) but aims at the current
estimate, so the Tp prediction no longer drives the turret; it is therefore reported as a
comparison and not adopted.
See [`IMPROVEMENTS.md`](IMPROVEMENTS.md) for the design study behind these choices.

## Requirements

- MATLAB (tested on R2026a) with the **Control System Toolbox** (`tf`, `margin`, `bode`,
  `feedback`, `stepinfo`, `bandwidth`, `lsim` are used by `design_analysis.m`; the
  simulation/estimation code itself — `simulate_case.m`, `turret_controller.m`,
  `turret_plant.m`, `target_motion_generator.m` — uses no toolbox and also runs under
  plain MATLAB or GNU Octave).
- (Optional, to rebuild the paper) a LaTeX distribution with the `IEEEtran` class and
  TikZ (e.g. MiKTeX or TeX Live): `pdflatex paper/paper.tex` (run twice).

## Reproducing every figure and table

From the repository root in MATLAB:

```matlab
design_analysis        % -> figures/fig_bode.*, fig_closedloop.*, fig_antiwindup.*
                       %    results/design_numbers.txt  (all Bode/margin/step numbers
                       %    quoted in Secs. III-IV)
run_scenarios          % -> figures/fig_scenario_A.*, _B.*, _C.*, _C2.*, fig_cov.*
                       %    results/scenario_metrics.csv/.tex, results/scenarios.mat
                       %    (Table III, Fig. 5 of the paper)
run_final_evaluation   % -> figures/fig_final_eval*.*
                       %    results/final_eval.csv, final_eval_summary.txt
                       %    (Fig. 6 of the paper; ~50 s for 40 cases)
run_headtohead         % -> results/headtohead.{csv,tex,mat}
                       %    figures/fig_headtohead.* and fig_headtohead_full.*
                       %    (Table IV of the paper; supplementary ablation figures;
                       %     ~5 min for 320 closed-loop simulations)
run_rigor              % -> figures/fig_pcrlb.*, fig_consistency.*, fig_feedforward.*
                       %    results/rigor.txt (Fig. 4 of the paper)
run_band_sweep         % -> results/band_sweep.txt (hand-off band sensitivity, Sec. IV-D)
```

Each script is self-contained (clears persistent state at the top) and writes directly
into `figures/` and `results/` at the repository root, so they can be run in any order.
All figures are exported with a white background by `paper_figure_style.m`, independent
of the local MATLAB theme.

To rebuild the paper PDF:

```bash
cd paper
pdflatex paper.tex
pdflatex paper.tex   # second pass resolves cross-references
```

## Design summary

| | |
|---|---|
| Plant | `Gp(s) = 3.125 / [s(s^2 + 10.47s + 3129.7)]`, type-1, `ωn = 55.94` rad/s, `ζ = 0.094` |
| Estimator | 6-state CA Kalman filter, Cartesian, Joseph-form covariance, `qj = 200` m²/s⁵ (within ~5% of the PCRLB) |
| Prediction | `x_pred = F(Tp) x`, `Tp = 0.5` s (fixed), `θref = atan2(y_pred, x_pred)` |
| Feedback | `Gc(s) = Kp/(1+s/p)`, `Kp = 2003` V/rad, `p = 4` rad/s — derived so the closed-loop group delay `= Tp` |
| Acquisition | dual-mode: rate-limited slew, then a bumpless hand-off to the LTI tracker (band 2.5°) |
| Margins | `PM = 65.2°`, `GM = 38.0 dB @ 29.4` rad/s, crossover `1.82` rad/s |
| Pulse logic | fire when `d⊥ ≤ 1.5` m and cross-beam `σ⊥ ≤ 1.0` m, `1.0` s recharge, current (not predicted) estimate |
| Result | 100% on 3 scenarios and 40/40 randomized cases (23/24 stress), mean RMS miss 0.75 m, mean interception time 2.67 s, no retuning |

See `paper/paper.pdf` Secs. III–VI for the full derivation, the PCRLB justification, and
the ablation (including why the PI compensator, the IMM/adaptive-Q estimators and the
feedforward variant were evaluated and not adopted).

## Function interfaces (exact signatures required by the brief)

```matlab
[r_meas, theta_meas, meas_valid] = target_motion_generator(t, scenario)
[V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret)
[theta_dot_turret, theta_ddot_turret, i_dot] = turret_plant(t, V_m, theta_turret, theta_dot_turret, i)
```

`turret_controller.m` additionally returns an optional 3rd output `dbg` (internal
estimator/predictor state, used only for plotting by `simulate_case.m`) and
`target_motion_generator.m` an optional 4th output `truth` (ground truth, used only for
scoring by `simulate_case.m`); neither is required by or breaks the mandated interface.
`scenario` accepts `'A'`/`'B'`/`'C'` (standardized scenarios), a struct
`struct('seed', n)` for a randomized bounded-maneuver + randomized-loss case (used by
`run_final_evaluation.m`), `struct('stress', k)` for a sustained high-rate coordinated-turn
case with a loss inside the turn (used by `run_headtohead.m`), or
`struct('base','B','loss',[t0 dur; ...])` for a fixed motion with custom loss timing
(case C2). Any other case identifier is mapped deterministically to one randomized case of
the same envelope-respecting family, so an external evaluation script that supplies its
own case names still runs the closed loop.

Persistent state in `turret_controller.m` (Kalman filter state/covariance, compensator
state, last-pulse time) is cleared automatically whenever `t` goes backwards (i.e. a new
run starts), and can also be cleared manually with `clear turret_controller`.

## AI usage

Per the assignment brief, AI usage was permitted and used throughout (derivation, MATLAB
implementation, simulation, figure generation, manuscript drafting). Full disclosure,
including what was accepted, corrected and rejected, is in the 1-page AI Appendix at the
end of `paper/paper.pdf`.
