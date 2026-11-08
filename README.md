# Predictive Tracking and Elevation Pointing Control for a Maneuvering Object

ME319 Control Systems (Autumn 2026, S3), IIT Bombay — bonus assignment.
Author: Ishaan Mondal.

One fixed estimator/predictor/controller design (never retuned) for the elevation
axis of a single-axis pointing mount: a constant-acceleration (CA) Kalman filter in
Cartesian coordinates, a fixed 0.5 s kinematic prediction, and a single lag
compensator `Gc(s) = Kp/(1+s/p)` whose pole is derived so the closed loop's own lag
equals the prediction horizon. Full derivation, Bode/margin analysis, scenario
results and a 40-case randomized robustness sweep are in [`paper/paper.pdf`](paper/paper.pdf)
(source: [`paper/paper.tex`](paper/paper.tex)).

**Headline result:** 100% interception success (40/40) across 40 randomized bounded
maneuvers with randomized sensor-loss intervals, one fixed design, no retuning.

## Repository layout

```
turret_plant.m              Plant ODEs (state [theta_turret, theta_dot_turret, i])
turret_controller.m         CA-KF + 0.5 s prediction + lag compensator + pulse logic
target_motion_generator.m   Scenarios A/B/C + randomized bounded-maneuver generator
turret_params.m             Single source of truth for every design constant
simulate_case.m             Closed-loop RK4 simulation driver + metrics
design_analysis.m           Plant derivation check, Bode/margins, anti-windup demo
run_scenarios.m             Reproduces all Scenario A/B/C(+C2) figures and tables
run_final_evaluation.m      Reproduces the 40-case robustness-sweep table/figure
figures/                    All generated .png/.pdf figures used in the paper
results/                    All generated .csv/.txt/.mat numeric results
paper/paper.tex, paper.pdf  The IEEE-format paper (incl. 1-page AI appendix)
```

## Requirements

- MATLAB (tested on R2026a) with the **Control System Toolbox** (`tf`, `margin`,
  `bode`, `feedback`, `stepinfo`, `bandwidth`, `lsim` are used by `design_analysis.m`;
  the simulation/estimation code itself, `simulate_case.m` / `turret_controller.m` /
  `turret_plant.m` / `target_motion_generator.m`, uses no toolbox and will also run
  under plain MATLAB or GNU Octave).
- (Optional, to rebuild the paper) a LaTeX distribution with the `IEEEtran` class
  (e.g. MiKTeX or TeX Live) — `pdflatex paper/paper.tex` (run twice).

## Reproducing every figure and table

From the repository root in MATLAB:

```matlab
design_analysis        % -> figures/fig_bode.*, fig_closedloop.*, fig_antiwindup.*
                        %    results/design_numbers.txt  (all Bode/margin/step numbers
                        %    quoted in the paper, Sec. IV)
run_scenarios           % -> figures/fig_scenario_A.*, _B.*, _C.*, _C2.*, fig_cov.*
                         %    results/scenario_metrics.csv/.tex, results/scenarios.mat
                         %    (Table II, Fig. 5-6 of the paper)
run_final_evaluation    % -> figures/fig_final_eval.*
                         %    results/final_eval.csv, final_eval_summary.txt
                         %    (Table III, Fig. 7 of the paper; ~20 s for 40 cases)
```

Each script is self-contained (clears persistent state at the top) and writes
directly into `figures/` and `results/` at the repository root, so they can be run
in any order. Total run time for all three is well under two minutes.

To rebuild the paper PDF:

```bash
cd paper
pdflatex paper.tex
pdflatex paper.tex   # second pass resolves cross-references
```

## Design summary

| | |
|---|---|
| Plant | `Gp(s) = 3.125 / [s(s^2 + 10.47s + 3129.7)]`, type-1, `ωn=55.94` rad/s, `ζ=0.094` |
| Estimator | 6-state CA Kalman filter, Cartesian, Joseph-form covariance, `qj=200` m²/s⁵ |
| Prediction | `x_pred = F(tp) x`, `tp = 0.5` s (fixed), `θref = atan2(y_pred, x_pred)`, clamped |
| Controller | `Gc(s) = Kp/(1+s/p)`, `Kp=2003` V/rad, `p=4` rad/s — derived so closed-loop lag `= tp` |
| Margins | `PM = 65.2°`, `GM = 38.0 dB @ 29.4 rad/s`, crossover `1.82` rad/s |
| Pulse logic | fire when `d_perp ≤ 1.5` m and cross-beam `σ ≤ 1.0` m, `1.0` s refractory, current (not predicted) estimate |

See `paper/paper.pdf` Sec. III–IV for the full derivation and justification
(including why the brief's suggested PI compensator was evaluated and rejected).

## Function interfaces (exact signatures required by the brief)

```matlab
[r_meas, theta_meas, meas_valid] = target_motion_generator(t, scenario)
[V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret)
[theta_dot_turret, theta_ddot_turret, i_dot] = turret_plant(t, V_m, theta_turret, theta_dot_turret, i)
```

`turret_controller.m` additionally returns an optional 3rd output `dbg` (internal
estimator/predictor state, used only for plotting by `simulate_case.m`) and
`target_motion_generator.m` an optional 4th output `truth` (ground truth, used only
for scoring by `simulate_case.m`); neither is required by or breaks the mandated
interface. `scenario` accepts `'A'`/`'B'`/`'C'` (standardized scenarios) or a struct
`struct('seed', n)` for a randomized bounded-maneuver + randomized-loss case (used by
`run_final_evaluation.m`), or `struct('base','B','loss',[t0 dur; ...])` for a fixed
motion with custom loss timing (used for the illustrative case C2).

Persistent state in `turret_controller.m` (Kalman filter state/covariance,
compensator state, last-pulse time) is cleared automatically whenever `t` goes
backwards (i.e. a new run starts), and can also be cleared manually with
`clear turret_controller`.

## AI usage

Per the assignment brief, AI usage was permitted and used throughout (derivation,
MATLAB implementation, simulation, figure generation, manuscript drafting). Full
disclosure, including what was accepted, modified and rejected, is in the
1-page AI Appendix at the end of `paper/paper.pdf`.
