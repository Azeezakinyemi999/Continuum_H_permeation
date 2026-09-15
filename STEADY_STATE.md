# Steady State, Transients, and the Accumulation Model

**Why the existing model needs no time axis, what it already computes, and the two conditions under which that would stop being true.**

Verified against `ACTIVE_STUDY = 'Guo_etal_2025_316L'` at the default operating point (873 K, 1 bar upstream, 1 mm wall, 48 nm Cr₂O₃).

---

## 1. The three model classes

Hydrogen does not damage a wall by being present in the gas. It has to get in, build up somewhere, and then degrade something. Three model classes describe those three links, and they are routinely conflated:

**gas → flux → inventory → property loss**

| | Question it answers | Output | Needs time? |
|---|---|---|---|
| **Permeation** | How fast does H cross? | $J$, $\Phi$, PRF, regime | No |
| **Accumulation** | How much H is where? | $C_L$, $C_T$, inventory, $C_{max}$ | Only if unsteady |
| **Embrittlement index** | How much strength is lost? | EI, RRA, $K_{IH}$, $\sigma_c(\theta)$ | No |

The permeation model sets the boundary condition, the accumulation model sets the driving variable, and the embrittlement index sets the failure criterion. A coating that cuts permeation by 10³ has *not* cut embrittlement risk by 10³ — inventory may already sit above the critical concentration.

---

## 2. What the model has today

Everything in [`calculations/`](calculations/) is **steady state**. That was verified, not assumed:

- No time-integration machinery anywhere. A search for `solve_ivp`, `odeint`, `scipy.integrate`, `RK45`, Euler stepping, `time_step`, `erfc`, `t_lag` and `breakthrough` returns **zero hits in executable code**. The word "transient" appears twice, both times in docstrings as a caveat about what the model cannot do.
- All 26 solver call sites are `brentq` / `root_scalar` — **root finders, not integrators**. They solve algebraic flux-balance equations ("find the interface pressure at which flux in = flux out"), which is the definition of a steady state.

### Beware: several outputs look temporal but are not

| Output | What it actually counts |
|---|---|
| `convergence_history` | root-solver iterations |
| `iteration`, `iterations` | nonlinear-solve index |
| `level4_iterations` | $D_{eff}$ self-consistency loop |
| `max_iterations = 15` | solver iteration cap |
| `converged` | did the algebraic solve converge |

Plotting `convergence_history` produces a settling curve, but its x-axis is solver iterations, not seconds. It carries no physical time information. Likewise `profiles['x']` is a **spatial** axis (0 → 1 mm), not time.

### Steady state does not mean zero spatial structure

This distinction matters. The model is a genuine 1-D **spatial** model:

| Quantity | Status |
|---|---|
| $C_L$ lattice concentration — `C_up`, `C_down`, `profiles['C']` | implemented |
| $C_T$ trapped concentration, per trap population | implemented |
| $\theta(x)$ trap occupancy profile | implemented |
| $C_L(x)$ spatial profile across the wall | implemented |
| $J$, $\Phi$, PRF, regime | implemented |
| $D_{eff}$ from grain boundaries + trapping | implemented |
| $\partial C/\partial t$ — inventory history, filling transient | **absent** |
| $\nabla\sigma_h$ — stress-driven drift | **absent** |
| $t_{lag}$, breakthrough time | **absent** |

At the default operating point the partition is:

$$
C_L = 12.33\ \mathrm{mol/m^3}, \qquad
\sum_i C_{T,i} = 4.26\ \mathrm{mol/m^3}, \qquad
C_{total} = 16.60\ \mathrm{mol/m^3}
$$

so **25.7 % of the hydrogen in the wall is trapped.** Cross-check: the trapped fraction computed this way agrees with `mobile_fraction = 0.7383` obtained independently from the $D_{eff}$ route, to within 0.5 %.

> One correction worth recording. Although $D_{eff}$ is evaluated at every grid point, the Oriani form as implemented, $D_{eff} = D_L/\left(1+\sum_i (N_{T,i}/N_L)e^{E_{b,i}/RT}\right)$, has **no $C_L$ dependence** — only $N_T$, $E_b$, $N_L$ and $T$. So `profiles['D']` comes out flat (2.069e-11 m²/s across the wall) while `profiles['theta']` varies. That is correct physics in the dilute-trap limit, but the comment at [`interface_solver.py:497`](calculations/interface_solver.py#L497) claiming "$D_{eff}$ depends on concentration" does not match the implementation.

---

## 3. The argument: at steady state, time is not needed

The claim to test is whether the wall actually *is* at steady state for the duration that matters. Four characteristic times govern this.

**Diffusion lag** — how long a concentration gradient takes to establish across a layer:

$$
t_{\text{lag}} = \frac{L^{2}}{6D}
$$

**Inventory fill time** — how long the steady flux takes to deposit the stored hydrogen:

$$
t_{\text{fill}} = \frac{C \cdot L}{J}
$$

Evaluated at the default operating point, with $J = 3.197\times10^{-6}$ mol m⁻² s⁻¹:

| Process | Expression | Value |
|---|---|---|
| Oxide diffusion lag | $L_{ox}^2/6D_{ox}$ | 28 s |
| Metal diffusion lag | $L_m^2/6D_{eff}$ | 16 min |
| Fill trap inventory | $\sum C_T \cdot L_m / J$ | 22 min |
| Fill lattice inventory | $\bar{C_L} \cdot L_m / J$ | 32 min |
| **Time to reach steady state** | max of the above | **≈ 32 min** |

Compared against operating duration:

| Duration | Transient as a fraction |
|---|---|
| 1 day | 2.2 % |
| 1 year | 0.0061 % |
| 40-year design life | **0.00015 %** |

**The system is at steady state for essentially all of its life.** A steady-state model is therefore the *correct* model here, not an approximation requiring apology.

### This also disposes of the accumulation question

At steady state, $J$ and $C$ are constant in time, so the integrals collapse to multiplications:

$$
\text{inventory} = \left(C_L + \sum_i C_{T,i}\right) \times \text{volume}
$$

$$
\text{cumulative throughput} = J \times A \times t
$$

Every term on the right already exists. Total trapped inventory is $4.26\ \mathrm{mol/m^3} \times 1\ \mathrm{mm} = 4.3\ \mathrm{mmol/m^2}$ — obtainable today, with no new code and no time integration.

---

## 4. The two conditions that would break the argument

Both are checkable with what exists.

**Condition 1 — a duty cycle shorter than ~30 min.** Start–stop cycling, thermal transients or pressure pulses on that timescale mean the wall never reaches steady state, and a steady-state flux would then over- or under-predict. For continuous high-temperature operation this does not apply.

**Condition 2 — lower temperature, or a much thicker wall.** Both lag terms scale as $L^2/D$, and $D$ is Arrhenius. Dropping from 873 K to 573 K reduces $D_{eff}$ by roughly two orders of magnitude, pushing the metal lag from 16 min toward a day or more. At the bottom of the configured `T_range` (623 K) the transient stops being negligible.

**Recommendation:** run the timescale calculation at the coldest and thickest condition of interest rather than assuming the 873 K result generalises. The deliverable is not the assumption but the number attached to it — *"the transient is 0.0002 % of design life at 873 K"* is a far stronger statement in a paper than *"steady state is assumed."*

---

## 5. When a time axis would still be worth adding

Only for **validation**, not for the physics.

| Goal | What to build |
|---|---|
| Flux, permeability, PRF, regime, inventory | **Nothing** — already steady state, already correct |
| Validate the trapping block against TDS | **Nothing** — needs `reduction_factor`, which exists |
| Predict $t_{lag}$, breakthrough, transient curves | Transient 1-D: add $\partial C/\partial t$ to the existing spatial grid |
| Total H inventory over service life | Transient + trap state variable |
| Crack-tip $C_{max}$ → embrittlement index | Full accumulation model with stress coupling |

The transient 1-D option is a convenience for comparing against measured permeation curves — it would let the model produce $D_{app}$ the way an experimentalist does, from a lag time:

$$
\frac{\partial C}{\partial t} = \frac{\partial}{\partial x}\left(D_{eff}\frac{\partial C}{\partial x}\right)
$$

It is a modest extension, because the spatial grid, the boundary conditions and $D_{eff}$ all already exist; only a time loop is added. It would **not change a single steady-state prediction.**

The full accumulation model is a different scope entirely. The stress-driven drift term

$$
\frac{\partial C_L}{\partial t} + \frac{\partial C_T}{\partial t}
= \nabla\cdot\left(D\nabla C_L - \frac{D\,C_L V_H}{RT}\nabla\sigma_h\right)
$$

requires a coupled finite-element stress solve. That is a different codebase from an analytical resistance network, and it only earns its cost if the work proceeds all the way to an embrittlement criterion.

---

## 6. Bottom line

For the questions this model is built to answer — steady permeation flux, effective permeability, barrier effectiveness, rate-limiting regime, and hydrogen inventory — **time is not needed, and adding it would not change any answer.**

The honest way to state this in writing is not "steady state is assumed" but:

> The wall reaches steady state in approximately 32 minutes at 873 K, which is $1.5\times10^{-6}$ of the 40-year design life. A steady-state formulation is therefore exact to within the accuracy of every other input.
