# Steady State, Transients, and Why Trapping Is Not in the Model

**Why the model needs no time axis, why hydrogen trapping therefore drops out of it, and the conditions under which that would stop being true.**

Verified against the `Guo_etal_2025` study (316L / Cr₂O₃) at its default operating point: 873 K, 1 bar upstream, 1 mm metal, 1 µm Cr₂O₃.

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

**This repository is a permeation model.** Accumulation (how much hydrogen the wall holds, and where) and embrittlement are future work.

---

## 2. What the model has today

Everything in [`calculations/`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/tree/main/calculations/) is **steady state**:

- There is no time integration anywhere: no `solve_ivp`, `odeint`, time stepping, `t_lag` or breakthrough time.
- Every solver call is a root finder (`brentq` / `root_scalar`) on an algebraic flux balance ("find the interface pressure, or the surface coverage, at which flux in = flux out"), which is the definition of a steady state.

The model is still a 1-D **spatial** model: the lattice concentration falls across each layer, and the interface pressure and surface coverage are solved for. What it does not have:

| Quantity | Status |
|---|---|
| $J$, $\Phi$, PRF, regime | implemented |
| $C_L$ at each face, interface pressure, surface coverage $\theta$ | implemented |
| $D_{metal}$ with grain-boundary fast paths | implemented |
| $C_T$, trapped hydrogen; total inventory | **not modelled** — see below |
| $\partial C/\partial t$, $t_{lag}$, breakthrough, filling transient | **absent** |
| $\nabla\sigma_h$, stress-driven drift | **absent** |

### Why hydrogen trapping is not in the model

Hydrogen in an alloy sits either in ordinary interstitial lattice sites, where it is mobile, or in traps (dislocations, grain boundaries, vacancies, carbide interfaces), where it is held by a binding energy $E_b$ and does not move. Only lattice hydrogen carries flux:

$$
J = -D_L\,\frac{\partial C_L}{\partial x},
\qquad
\frac{\partial (C_L + C_T)}{\partial t} = -\frac{\partial J}{\partial x}.
$$

Traps add storage. If they are in local equilibrium with the lattice (Oriani), $C_T$ is a function of $C_L$ and the balance becomes a diffusion equation for $C_L$ with an **effective** diffusivity

$$
\frac{\partial C_L}{\partial t} = D_{eff}\,\frac{\partial^2 C_L}{\partial x^2},
\qquad
D_{eff} = \frac{D_L}{1 + \mathrm{d}C_T/\mathrm{d}C_L}.
$$

$D_{eff}$ governs how fast the wall *fills*: the time lag $t_{lag} = L^2/6D_{eff}$ and the build-up of inventory. At steady state, $\partial C/\partial t = 0$, the traps are full to their equilibrium occupancy and exchange no net hydrogen, and the flux is

$$
J = \frac{D_L\,(C_{L,up} - C_{L,down})}{L}.
$$

The same result holds with saturable traps, because $J = -D_{eff}\,\partial C_{total}/\partial x$ integrated across the wall gives $J L = \int D_{eff}\,\mathrm{d}C_{total} = \int D_L\,\mathrm{d}C_L = D_L\,\Delta C_L$. It also holds for slow (McNabb–Foster) trapping, since once the profile stops changing the net trapping rate is zero. **The steady-state permeation flux does not depend on trap density, binding energy or occupancy.**

Until this was corrected, the code computed $J = D_{eff}\,\Delta C_L/L$, pairing the trap-reduced $D_{eff}$ with the *lattice* concentrations. That counts trapping twice. At the 316L default point it lowered the Level 5 flux by about 26 % (4.30 → 3.18 µmol m⁻² s⁻¹), and at the top of the sensitivity range for vacancy density it lowered the bare-metal flux by up to 35×. It also gave the trap inputs a spurious place among the sensitivity rankings.

Trapping was therefore removed, not switched off. It is still the right physics for **inventory** and **time lag**, which matter for hydrogen embrittlement and tritium accounting. Those are accumulation questions (§1) and outside this project's scope. If they are taken up, trapping belongs in a separate step that takes the steady $C_L(x)$ profile from this model and computes $C_T(x)$, the inventory and $t_{lag}$. It does not belong in the flux.

**Two caveats.**

- Grain boundaries are both traps and fast paths. The fast-path part does change the steady flux, and stays in the model as the Level 4 grain-boundary enhancement.
- The diffusivity *data* can carry trapping. A diffusivity measured from a permeation time lag in a trap-rich specimen is an apparent $D_{eff}$, and a solubility computed as $\Phi/D_{app}$ is an apparent total solubility. Their product, the permeability $\Phi$, is still right, so Levels 1–5 are unaffected. The Level 6 surface coupling, however, uses $D$ and $K_s$ separately and should be given lattice values.

---

## 3. The argument: at steady state, time is not needed

The claim to test is whether the wall actually *is* at steady state for the duration that matters. Three characteristic times govern this.

**Diffusion lag** — how long a concentration gradient takes to establish across a layer:

$$
t_{\text{lag}} = \frac{L^{2}}{6D}
$$

**Lattice fill time** — how long the steady flux takes to deposit the mobile hydrogen in the metal:

$$
t_{\text{fill}} = \frac{\bar{C}_L \cdot L}{J}
$$

Evaluated at the default operating point, with $J = 4.30\times10^{-6}$ mol m⁻² s⁻¹:

| Process | Expression | Value |
|---|---|---|
| Oxide diffusion lag | $L_{ox}^2/6D_{ox}$ | 3.3 h |
| Metal diffusion lag | $L_m^2/6D_{metal}$ | 12 min |
| Fill lattice inventory | $\bar{C}_L \cdot L_m / J$ | 35 min |
| **Time to reach steady state** | max of the above | **≈ 3.3 h** |

These are lower bounds. Traps slow the filling: the metal lag uses $D_{eff} \le D_L$, and the traps must be filled too. With the trap densities previously configured, that lengthened the metal terms by about a third, which is still far below the oxide lag.

Compared against operating duration:

| Duration | Transient as a fraction |
|---|---|
| 1 day | 14 % |
| 1 year | 0.04 % |
| 40-year design life | **0.001 %** |

**The system is at steady state for essentially all of its life.** A steady-state model is therefore the *correct* model here, not an approximation requiring apology.

Cumulative throughput is then simply $J \times A \times t$.

---

## 4. The conditions that would break the argument

**Condition 1 — a duty cycle shorter than a few hours.** Start–stop cycling, thermal transients or pressure pulses on that timescale mean the wall never reaches steady state, and a steady-state flux would then over- or under-predict. For continuous high-temperature operation this does not apply.

**Condition 2 — lower temperature, a thicker wall or coating, or heavy trapping.** Every lag term scales as $L^2/D$, and $D$ is Arrhenius. Dropping from 873 K to 573 K reduces the diffusivities by roughly two orders of magnitude, pushing the lag from hours toward weeks. Strong, dense traps lengthen the metal lag by the factor $1 + \mathrm{d}C_T/\mathrm{d}C_L$.

**Recommendation:** run the timescale calculation at the coldest and thickest condition of interest rather than assuming the 873 K result generalises. The deliverable is not the assumption but the number attached to it — *"the transient is 0.001 % of design life at 873 K"* is a far stronger statement in a paper than *"steady state is assumed."*

---

## 5. When a time axis would be worth adding

| Goal | What to build |
|---|---|
| Flux, permeability, PRF, regime | **Nothing** — already steady state, already correct |
| Hydrogen inventory (lattice + trapped) | Post-processing of the steady $C_L(x)$ with a trap model — no time axis, but trap data |
| Predict $t_{lag}$, breakthrough, transient curves | Transient 1-D: $\partial C_L/\partial t = \partial_x(D_{eff}\,\partial_x C_L)$ on the spatial grid |
| Inventory build-up over service life | Transient + trap state variable |
| Crack-tip $C_{max}$ → embrittlement index | Full accumulation model with stress coupling |

The transient 1-D option is a convenience for comparing against measured permeation curves. It would let the model produce $D_{app}$ the way an experimentalist does, from a lag time, and it is where a trap model would first be needed. It would **not change a single steady-state prediction.**

The full accumulation model is a different scope entirely. The stress-driven drift term

$$
\frac{\partial C_L}{\partial t} + \frac{\partial C_T}{\partial t}
= \nabla\cdot\left(D\nabla C_L - \frac{D\,C_L V_H}{RT}\nabla\sigma_h\right)
$$

requires a coupled finite-element stress solve. That is a different codebase from an analytical resistance network, and it only earns its cost if the work proceeds all the way to an embrittlement criterion.

---

## 6. Bottom line

For the questions this model is built to answer — steady permeation flux, effective permeability, barrier effectiveness and rate-limiting regime — **time is not needed, trapping does not enter, and adding either would not change any answer.**

The honest way to state this in writing is not "steady state is assumed" but:

> The wall reaches steady state in about 3 hours at 873 K, which is $10^{-5}$ of the 40-year design life. A steady-state formulation is therefore exact to within the accuracy of every other input, and at steady state hydrogen trapping does not affect the permeation flux.
