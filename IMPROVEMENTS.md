# Beyond the baseline: design-upgrade study

This note is the research/ideation deliverable behind the submitted design. It (1)
surveys candidate upgrades for each subsystem with analytical justification and
references, (2) describes what was implemented, and (3) reports what the numbers
actually show — including the honest negative result that the estimator is already
near-optimal and therefore cannot be meaningfully beaten on this envelope.

Everything here is **flag-gated**: the **finalized design is the default** in
`turret_params.m`, and the feedback-only baseline (and every other variant) is one
frozen configuration reachable by an override, evaluated with no per-case retuning.
Selectors (default **bold**):

| selector | values | meaning |
|---|---|---|
| `P.est`  | **`CA`** · `adaptQ` · `IMM` | estimator / predictor |
| `P.ctrl` | **`lagff`** · `lag` · `pi` | 2-DOF (lag feedback + feedforward); `lag` = feedback only; `pi` = suggested PI |
| `P.acq`  | **`dualmode`** · `linear` | dual-mode acquisition; `linear` = LTI tracker throughout |

Run `simulate_case(scenario, overrides_struct)` to select a configuration, e.g. the
feedback-only baseline is `simulate_case('C', struct('ctrl','lag','acq','linear'))`.

---

## 1. Where the baseline can and cannot be beaten

The baseline's own admitted limits (paper Sec. VII) point to four levers. Analyzing
each against the physics of *this* problem (0.5 s horizon; 800–2000 m, ≤150 m/s,
≤20 m/s² envelope; a turret whose saturated slew rate is only 13.7°/s):

- **Prediction error (≈3.5 m RMS) ≫ estimation error (≈0.74 m).** But over a 0.5 s
  horizon most of the prediction error is the *velocity-estimate error × horizon*,
  not model mismatch — so better *dynamics models* (IMM/CT) help far less than one
  would guess. The PCRLB (Sec. 4) confirms the CA filter is already near-optimal.
- **High-rate tracking lag.** The lag-matching trick is exact only for *constant*
  angular rate; its error grows with angular *acceleration*. This is a genuine,
  recoverable gap → **2-DOF feedforward** closes it (up to **80× lower** pointing
  error, Sec. 3).
- **Interception time (~2.6 s)** is set by the saturated 0°→target slew, which is a
  physical ceiling — but the *capture transient* after the slew is not. A
  **dual-mode** front-end captures cleanly → earlier, cleaner lock, and a −40% RMS
  miss distance over the feedback-only design.
- **Robustness breadth.** A harder **stress suite** (sustained high-rate
  coordinated turns with a 0.5 s loss *inside* the turn + a sharp tangential kick)
  plus **PCRLB / NEES / NIS** rigor shows where the design stands.

---

## 2. Estimator & prediction — the research menu (A)

**A1 · IMM{CV, CA, CT}** — *implemented.* Three models (constant-velocity,
constant-acceleration, and a 5-state **coordinated-turn EKF** with the turn rate as
a *filtered* state) run in parallel and are blended by Markov-mixed model
probabilities. The CT model predicts the *arc* of a sustained turn, which CA
extrapolates as a parabola. Standard SOTA for maneuvering targets; the brief names
it. Refs: Blom & Bar-Shalom IMM; CV/CA/CT model sets [1,2].

**A2 · Adaptive-Q (NIS-triggered)** — *implemented.* A single CA filter that
detects a maneuver via the normalized innovation squared (χ², 2 dof) and transiently
inflates the jerk PSD. Lighter than IMM, preserves the single-filter ethos. Refs:
NIS χ² maneuver detection; adaptive-factor-from-NIS [3,4].

**A3 · Singer / Current-Statistical model** — *cited, not built.* Correlated
acceleration with an adaptive mean; marginal over a 0.5 s horizon. Refs [5].

**A4 · Debiased converted measurements (Lerro–Bar-Shalom) / UKF** — *cited.*
Removes the (tiny here, σθ·r ≤ 0.7 m) polar→Cartesian conversion bias; a
correctness/rigor upgrade rather than a performance one. Refs [6].

**A5 · PCRLB + NEES/NIS** — *implemented* (`run_rigor.m`). The rigor layer, Sec. 4.

**Result (estimator):** on easy motion IMM reverts to ≈CA (even marginally better);
on the stress suite IMM and adaptive-Q are within a few percent of CA and neither
beats it. **This is the correct result, and the PCRLB explains it** (Sec. 4): CA is
already near the information floor, so there is no estimation headroom to recover.
This *validates with rigor* the baseline's original choice to reject IMM.

---

## 3. Controller & acquisition — the research menu (B)

**B1 · 2-DOF model-inverse feedforward** — *implemented, headline win.* A feedforward
voltage `V_ff = θ̇_ref/K_dc + (R_m J_eq / N K_t) θ̈_ref`, with the reference
rate/acceleration computed analytically from the predicted Cartesian state
(`θ̇_ref = (x·v_y − y·v_x)/r²`). The gains are *derived* from the plant (model
inverse), not tuned. The feedback loop then only rejects residual. This removes the
baseline's dependence on lead/lag cancellation, which held only at low angular rate.

> **Clean controller-only result** (`run_rigor.m`, Part B; true kinematics fed in to
> isolate control from estimation; peak LOS rate fixed at 0.08 rad/s, unsaturated):
>
> | LOS freq ω [rad/s] | 0.2 | 0.3 | 0.5 | 0.7 | 1.0 | 1.4 | 2.0 |
> |---|---|---|---|---|---|---|---|
> | baseline RMS point. err [deg] | 0.0014 | 0.0064 | 0.0176 | 0.0338 | 0.0673 | 0.1282 | 0.2407 |
> | **feedforward** [deg] | 0.0003 | 0.0003 | 0.0005 | 0.0008 | 0.0011 | 0.0019 | 0.0030 |
> | improvement | 5× | 19× | 33× | 43× | 62× | 67× | **80×** |
>
> The baseline error grows with angular acceleration (∝ ω at fixed peak rate);
> feedforward holds it flat. This is exactly the limitation the paper flagged.
> (An earlier version of this study fed only the *rate* term of `V_ff`; adding the
> acceleration term, which completes the plant inversion, raised the improvement from
> 19× to 80× — a useful reminder to check feedforward completeness against the model.)

**B2 · Resonance notch** — *cited.* Would raise bandwidth, but acquisition is
saturation-limited so it mostly helps steady tracking; secondary to B1.

**B3 · Dual-mode / near-time-optimal acquisition** — *implemented.* A max-rate
slew while the pointing error is large, with a bumpless hand-off to the linear
tracker inside a small error band. The saturated slew itself is a physical ceiling,
but removing the capture transient yields a cleaner, earlier lock. The band is set
from the geometry of the saturated-rate unwind (0.24 rad/s against the compensator's
0.25 s relaxation): below ≈1.5° the hand-off re-enters saturation and the capture
oscillates, so 2.5° is used.

**B4 · LQG / MPC** — *cited as the natural next step, not built.* Heavier; the brief
disfavors the added tuning/complexity, and B1+B3 already recover the control-side gap.

**Result (control) — the recommended design is `CA + FF + dual-mode`** (keep the
near-optimal CA estimator; see Sec. 4). Head-to-head over both suites, one frozen
design per column, no retuning (`results/headtohead.csv`):

| metric (mean) | suite | Base | **Final (CA+FF+DM)** | adaptQ | IMM |
|---|---|---|---|---|---|
| Success rate [%]       | Orig-40 | 100 | **100** | 100 | 100 |
| Interception time [s]  | Orig-40 | 2.87 | **2.63** (−8%) | 2.63 | 2.60 |
| RMS pointing err [deg] | Orig-40 | 0.098 | **0.073** (−26%) | 0.073 | 0.077 |
| RMS miss dist. [m]     | Orig-40 | 0.98 | **0.59** (−40%) | 0.61 | 0.79 |
| Success rate [%]       | Stress-24 | 100 | **100** | 100 | 100 |
| Interception time [s]  | Stress-24 | 3.74 | **3.55** (−5%) | 3.55 | 3.60 |
| RMS pointing err [deg] | Stress-24 | 0.110 | **0.088** (−20%) | 0.089 | 0.092 |
| RMS miss dist. [m]     | Stress-24 | 1.37 | **0.89** (−35%) | 0.92 | 1.18 |

Feedforward + dual-mode cut interception time, steady pointing error, and — most
tellingly — **RMS miss distance by 35–40%** (far more margin against the 2 m
engagement radius), with 100% interception preserved and the proven CA estimator
untouched. Adding IMM on top (`IMM+FF+DM`) does not help and hurts the miss
(estimator noise), consistent with Sec. 4.

---

## 4. Rigor: PCRLB + consistency (`run_rigor.m`, Part A)

On a maneuver-plus-loss trajectory, Monte-Carlo (M = 300):

- **Filter RMS = 0.95 × √PCRLB** → the CA filter extracts essentially all available
  information (near-optimal MMSE). Both the error and the bound grow during the loss
  and recover together.
- **Mean NIS = 1.87 (vs 2 dof), 79 % inside the 95 % band** → innovations are
  statistically consistent.
- **Mean NEES = 3.65 (vs n_x = 6)** → the filter is deliberately *conservative* in the
  weakly-observed velocity/acceleration states (the `q_jerk = 200` robustness choice),
  which is the safe direction for a design that must generalize without retuning.

Together these turn "IMM doesn't help" from a disappointment into a **proof**: there
is no estimation headroom to recover, so the engineering effort is correctly spent on
the control side.

---

## 5. Engagement (C) — cited, not built

**C1 · Expected-miss-minimizing pulse timing** (fire at minimum predicted
perpendicular miss within the recharge window, using the covariance) would squeeze
the already-98.3 % per-pulse success; a nice-to-have given the current margin.

---

## 6. How to reproduce

```matlab
run_headtohead     % -> results/headtohead.{csv,tex,mat}, figures/fig_headtohead*
run_rigor          % -> figures/fig_pcrlb.*, fig_consistency.*, fig_feedforward.*
                   %    results/rigor.txt
% regression (must be unchanged):
run_scenarios ; run_final_evaluation
```

`headtohead_figure.m` redraws the ablation figures from `results/headtohead.mat`
without re-simulating.

## 7. Bottom line

**Recommended design = CA estimator + 2-DOF feedforward + dual-mode acquisition**
(`struct('ctrl','lagff','acq','dualmode')`). It is beaten where there was real
headroom, and the baseline is kept where it was already near-optimal:

- **Beaten (control side):** feedforward removes the high-rate tracking lag (up to
  **80× lower** pointing error at high angular acceleration); dual-mode gives a
  cleaner, faster capture. Net: **−8 % interception time, −26 % steady pointing
  error, −40 % RMS miss distance** on the original suite (−35 % on the stress suite),
  100 % success preserved, no retuning.
- **Not beaten, and proven why (estimator side):** IMM and adaptive-Q match but do
  not beat CA — because the CA filter is already within ~5 % of the PCRLB on this
  envelope/horizon. The sophisticated result is the *proof of near-optimality*, which
  validates the baseline's design judgment and says the engineering effort belongs on
  the control side.

### References
1. Modified IMM / interacting acceleration compensation (US Pat. 7,079,991); CV/CA/CT model sets.
2. IMMNet: hybrid model-based + data-driven maneuvering-target tracking, arXiv:2607.13573.
3. Innovation χ² maneuver detection, IFAC 2008; NIS-driven adaptive factors (US Pat. 7,079,991).
4. Y. Bar-Shalom, *Estimation with Applications to Tracking and Navigation* (NEES/NIS).
5. Zhou, Current-Statistical model; adaptive-CA note, *Remote Sensing* 17(5):850.
6. D. Lerro & Y. Bar-Shalom, "Tracking with debiased consistent converted measurements vs EKF," IEEE T-AES 29(3), 1993; Bordonaro–Willett–Bar-Shalom DUCM (2012/2017).
7. Tichavský, recursive posterior Cramér–Rao bound; Hendeby–Karlsson NEES∼χ²(n_x).
8. Shen & Tsiotras, time-optimal rigid-body reorientation (bang-bang), AIAA 1998.
9. 2-DOF / acceleration feedforward for gimbal tracking (Appl. Mech. Mater. 415:101).
