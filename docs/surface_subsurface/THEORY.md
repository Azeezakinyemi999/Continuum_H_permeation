---
title: "Surface, subsurface and bulk hydrogen"
subtitle: "How the surface kinetics, the MLMD and the continuum model connect, and what a kinetic absorption step would add"
date: "October 2026"
---

# 0. What this note is for

This note sets out the theory before anything is built. It does not change the
existing code bases. It answers four questions:

1. Can the surface kinetics (RMG and Cantera, with no bulk) be computed separately
   from the subsurface and the bulk (MLMD or literature), and combined later in the
   continuum model?
2. What exactly does the continuum model (`MHI_permeation`) assume about the
   surface, the subsurface and recombination?
3. What would a kinetic surface-to-subsurface step (absorption) add, and how could it
   be built as a separate module that calls the existing code?
4. How can each assumption of each level, in both continuum models, be tested?

**Short answers.**

1. **Yes.** The equilibria separate exactly. The steady state separates whenever
   every step exchanges hydrogen much faster than the permeation flux. The rate
   constants can always be computed separately; only the final steady-state solve
   couples them.
2. In the continuum model, $\theta$ is the coverage of the **gas-facing surface**.
   Recombination is the **surface** step $2\mathrm{H^*} \to \mathrm{H_2}$. The
   subsurface has no coverage: it is a concentration, fixed by local equilibrium
   with the surface. This is the "coverage-modified Sieverts" law.
3. A kinetic absorption step gives back today's model when absorption is fast. It
   adds an absorption-limited regime and fits the code's existing one-unknown solve.
4. Each assumption gets a limit check, a relaxation test and an experimental
   signature that could disprove it (section 9). One of them, trapping in the Level 4
   steady-state flux, conflicts with standard theory unless $D$ and $K_s$ are
   defined in a particular way (section 9.2).

# 1. Notation and counting

| Symbol | Meaning | Units |
|---|---|---|
| $P$, $P^\circ$ | $\mathrm{H_2}$ pressure; standard pressure ($10^5$ Pa) | Pa |
| $\Gamma$ | surface site density | mol m$^{-2}$ |
| $\theta_H$, $\theta_v$, $\theta_X$ | fraction of surface sites holding H, empty, holding another species | -- |
| $c_s$, $c$ | H concentration just below the surface; in the bulk | mol H m$^{-3}$ |
| $N_L$, $d$ | subsurface interstitial sites; subsurface layer spacing | mol m$^{-3}$; m |
| $k_d$ | dissociative adsorption coefficient | mol $\mathrm{H_2}$ m$^{-2}$ s$^{-1}$ Pa$^{-1}$ |
| $k_r^s$ | surface recombination coefficient | m$^2$ mol$^{-1}$ s$^{-1}$ |
| $K$; $K_{eq} = K/P^\circ$ | surface equilibrium constant (Cantera form; continuum form) | --; Pa$^{-1}$ |
| $k_a$ | absorption coefficient, surface to subsurface | mol H m$^{-2}$ s$^{-1}$ |
| $k_b$ | resurfacing coefficient, subsurface to surface | m s$^{-1}$ |
| $c^* = k_a/k_b$ | surface-subsurface partition | mol H m$^{-3}$ |
| $K_s$ | Sieverts constant | mol H m$^{-3}$ Pa$^{-1/2}$ |
| $D$, $L$ | diffusivity; wall thickness | m$^2$ s$^{-1}$; m |

**Counting rule.** Rates on the gas side are counted in mol $\mathrm{H_2}$. Anything
inside the solid is counted in mol H. One $\mathrm{H_2}$ delivers two H.

- The surface package reports fluxes in mol $\mathrm{H_2}$.
- MLMD reports the permeability in mol H.
- The continuum model does not declare a convention; its equation catalogue calls
  the surface flux "the net rate of atoms".
- Table 5 (`Hydrogen-Surface-kinetics/docs/literature/`) counts the dissolved
  concentration in mol $\mathrm{H_2}$; its $k_r$ is in mol $\mathrm{H_2}^{-1}$ m$^4$ s$^{-1}$.

# 2. The four surface steps

The surface exchanges hydrogen with the gas (two steps) and with the metal beneath
it (two steps). Diffusion then carries it through the wall.

$$
\begin{aligned}
&(1)\ \ \mathrm{H_2(g)} + 2\,* \to 2\,\mathrm{H^*} &\qquad r_1 &= k_d\,P\,\theta_v^2 &&[\mathrm{mol\,H_2\,m^{-2}\,s^{-1}}]\\
&(2)\ \ 2\,\mathrm{H^*} \to \mathrm{H_2(g)} + 2\,* & r_2 &= k_r^s\,(\Gamma\,\theta_H)^2 &&[\mathrm{mol\,H_2\,m^{-2}\,s^{-1}}]\\
&(3)\ \ \mathrm{H^*} + V_{sub} \to \mathrm{H_{sub}} + * & r_3 &= k_a\,\theta_H &&[\mathrm{mol\,H\,m^{-2}\,s^{-1}}]\\
&(4)\ \ \mathrm{H_{sub}} + * \to \mathrm{H^*} + V_{sub} & r_4 &= k_b\,c_s\,\theta_v &&[\mathrm{mol\,H\,m^{-2}\,s^{-1}}]\\
&(5)\ \ \text{diffusion into the wall} & J_D &= \tfrac{D}{L}\,(c_s - c_L) &&[\mathrm{mol\,H\,m^{-2}\,s^{-1}}]
\end{aligned}
$$

- $\theta_v = 1 - \theta_H - \sum_X \theta_X$ counts every occupied site, impurities included.
- Steps 3 and 4 assume a dilute subsurface: the factor $(1 - c_s/N_L)$ is dropped.
- $\theta_v$ appears in step 4 because an H atom coming up needs an empty surface site.

The absorption coefficients come from the MLMD Hop A rates $\nu_{fwd}$ and $\nu_{rev}$ (s$^{-1}$):

$$
k_a = \Gamma\,\nu_{fwd}, \qquad k_b = d\,\nu_{rev}.
$$

Here $\Gamma\theta_H$ is the H per area on the surface, and $c_s d$ is the H per area
in the first subsurface layer.

**Balances.** Hydrogen accumulates on the surface or in the subsurface only if
the steps around it are unbalanced:

$$
\Gamma\,\frac{d\theta_H}{dt} = 2(r_1 - r_2) - (r_3 - r_4),
\qquad
d\,\frac{dc_s}{dt} \approx (r_3 - r_4) - J_D .
$$

At steady state one flux passes through every step:

$$
\boxed{\;2(r_1 - r_2) \;=\; r_3 - r_4 \;=\; J_D \;\equiv\; J\;}
$$

| Step | RMG-Cantera | MLMD | Continuum, Model 1 (Levels 1-5) | Continuum, Model 2 (Level 6 family) |
|---|---|---|---|---|
| (1)(2) dissociation, recombination | **computed** | dissociation barrier (Part 1) | equilibrium | **kinetic** |
| (3)(4) absorption, resurfacing | absent | **Hop A** (Part 2) | equilibrium | equilibrium |
| (5) diffusion | absent | **$D$** (Part 3) | Fick and Sieverts | Fick and Sieverts |
| exit surface | absent | -- | equilibrium | equilibrium |

# 3. Equilibrium: the four steps separate exactly

At equilibrium $J = 0$. The steps are in series, so zero net flux means zero net
rate in **every** step. Each pair balances on its own:

$$
r_1 = r_2 \;\Rightarrow\; \frac{\theta_H}{\theta_v} = \sqrt{\frac{K\,P}{P^\circ}},
\qquad \frac{K}{P^\circ} = \frac{k_d}{k_r^s\,\Gamma^2}
\qquad\text{(this is what Cantera computes)}
$$

$$
r_3 = r_4 \;\Rightarrow\; c_s = c^*\,\frac{\theta_H}{\theta_v},
\qquad c^* = \frac{k_a}{k_b} = \frac{\Gamma}{d}\,\frac{\nu_{fwd}}{\nu_{rev}}
\qquad\text{(this is what MLMD Hop A gives)}
$$

Together they give Sieverts' law:

$$
c_s = c^*\sqrt{K/P^\circ}\;\sqrt{P} = K_s\sqrt{P},
\qquad
\boxed{\;K_s = c^*\sqrt{K/P^\circ}\;}
$$

So **Cantera's surface-only equilibrium is exactly the surface state of the full
four-step equilibrium.** Leaving the subsurface out of Cantera loses nothing at
equilibrium.

The same relation in enthalpies, per mol H:

$$
\boxed{\;\Delta H_s = \tfrac{1}{2}\,\Delta H_{ads} + \Delta H_{ss}\;}
$$

where $\Delta H_{ads}$ is the adsorption enthalpy per mol $\mathrm{H_2}$ (negative), and
$\Delta H_{ss}$ is the surface-to-subsurface enthalpy (the Hop A reaction enthalpy,
positive on Ni).

**Consequence.** $K_s$, $K$ and $c^*$ are not independent: any two fix the third.
If $K$ comes from Cantera and $K_s$ from the bulk, then $c^* = K_s/\sqrt{K/P^\circ}$
follows, and it must agree with MLMD's Hop A. This is the consistency check for
mixing sources.

# 4. Steady state with permeation: when the separation holds

With $J \neq 0$ the pairs no longer balance individually. Each step carries $J$ by
moving a little away from equilibrium. Close to equilibrium, a step behaves like a
resistor whose conductance is its **exchange rate**, the one-way rate at equilibrium:

$$
R_{0,surf} = r_1 = r_2\big|_{eq},
\qquad
R_{0,sub} = r_3 = r_4\big|_{eq} = k_a\,\theta_{H,eq},
\qquad
J_{D,max} = \frac{D\,K_s\sqrt{P}}{L}.
$$

For a step with net rate $J$, the chemical-potential drop is about
$\Delta\mu/RT \approx J/R_0$. Two dimensionless numbers decide everything:

$$
W_s = \frac{J_{D,max}}{2\,R_{0,surf}}, \qquad W_a = \frac{J_{D,max}}{R_{0,sub}}.
$$

| Condition | Meaning | What to do |
|---|---|---|
| $W_s \ll 1$ and $W_a \ll 1$ | both surface pairs stay at equilibrium | compute every step separately; Sieverts applies at the surface |
| $W_s \sim 1$ | dissociation and recombination limit | solve the surface pair together with the flux (Level 6 does this) |
| $W_a \sim 1$ | absorption limits | solve the subsurface pair together with the flux (no model does this yet; section 7) |

Impurities lower $R_{0,surf}$ through $\theta_v$. Oxide surfaces can lower both
exchange rates. Those are the cases where the ratios approach 1.

# 5. One equation for all three descriptions

Assume absorption is fast, so that $c_s = c^*\theta_H/\theta_v$. Let impurities
adsorb in their own equilibrium with the gas, $\theta_X = b\,\theta_v$ with
$b = \sum_X K_X p_X$. Because all site fractions add up to one,

$$
\theta_v = \frac{1}{1 + b + c_s/c^*},
\qquad
\theta_H = \frac{c_s/c^*}{1 + b + c_s/c^*}.
$$

Putting these into $r_1 - r_2$ gives the **master equation** (in mol $\mathrm{H_2}$):

$$
\boxed{\;J_{\mathrm{H_2}} = \frac{k_d\,P - k_r^b\,c_s^2}{\left(1 + b + c_s/c^*\right)^2},
\qquad k_r^b = \frac{k_r^s\,\Gamma^2}{c^{*2}}\;}
\qquad [k_r^b] = \mathrm{m^4\,mol^{-1}\,s^{-1}}
$$

| Limit | Equation | Used by |
|---|---|---|
| clean surface, low coverage: $b = 0$, $c_s \ll c^*$ | $J = k_d P - k_r^b c_s^2$ | Table 5 and the permeation literature |
| written with coverages instead of $c_s$ | $J = k_d P\theta_v^2 - k_r^s\Gamma^2\theta_H^2$ | RMG-Cantera |
| clean surface, any coverage: $b = 0$ | $J = (k_d P - k_r^b c_s^2)/(1 + c_s/c^*)^2$ | the continuum's Level 6 |

The denominator describes the surface filling up. An experiment fitted with
$J = k_d P - k_r c^2$ has assumed it equals 1.

**Bulk equilibrium.** Setting $J = 0$ removes the denominator:

$$
c_s = \sqrt{k_d/k_r^b}\,\sqrt{P}
\;\Rightarrow\;
K_s = \sqrt{k_d/k_r^b}
\;\Leftrightarrow\;
k_r^b = k_d/K_s^2 .
$$

Impurities ($b$) and crowding ($c_s/c^*$) change how fast equilibrium is reached,
never where it lies.

**Reading Table 5.**

1. **The hidden solubility.** Each row contains a Sieverts constant:
   $K_s = \sqrt{k_d/k_r}$ and $\Delta H_s = (E_d - E_r)/2$.
2. **$k_d$** is the same coefficient as Cantera's clean-surface $k_d$. The
   "oxidized" and "activated" rows are effective values, $k_d/(1+b)^2$.
3. **$k_r$** is not comparable with $k_r^s$. Converting needs $c^*$:
   $k_r^s = k_r^b c^{*2}/\Gamma^2$. With $c^* = K_s/\sqrt{K/P^\circ}$, this reduces to
   $k_r^s = k_d P^\circ/(K\Gamma^2)$, which is only surface detailed balance.
   **Permeation experiments therefore determine $k_d$ and $K_s$, but not $K$ or $c^*$
   separately**, unless the surface saturates and the denominator becomes visible.
4. **Units.** Table 5 counts $c$ in mol $\mathrm{H_2}$ m$^{-3}$. With $c$ in mol H,
   $k_r^b = k_r^{\mathrm{Table\,5}}/4$ and $K_s(\mathrm{H}) = 2\,K_s(\mathrm{H_2})$.

# 6. What the continuum model does today

## 6.1 Two models

| | Model 1 (Levels 1-5) | Model 2 (Level 6 family) |
|---|---|---|
| $\mathrm{H_2}$ dissociates | at the buried oxide/metal interface, instantly | at the gas-facing surface, at a finite rate |
| the oxide carries | $\mathrm{H_2}$ molecules (Henry, $C = K_{ox}P$) | H atoms (Sieverts, $C = K_{ox}\sqrt{P}$) |
| steps (1)(2) | equilibrium | kinetic |
| steps (3)(4) | equilibrium | equilibrium |
| surface data that can enter the flux | none | $k_{diss}$, $K_{eq}$ for the oxide; $k_{diss,metal}$, $K_{eq,metal}$ for bare metal and pinholes |

Model 1 is the master equation with $k_d \to \infty$: the flux is pure Sieverts, and
neither $K$, $b$ nor $c^*$ appears. Surface data can only test whether Model 1's
equilibrium assumption holds ($W_s \ll 1$), at $P_{up}$ on bare metal or at $P_{int}$
under an oxide.

## 6.2 The coverage is the surface coverage

$$
\mu_{H,ads} = \mu^\circ_{H,ads} + RT\ln\frac{\theta}{1-\theta},
\qquad
J_{surf} = k_{diss}\left[P(1-\theta)^2 - \frac{\theta^2}{K_{eq}}\right]
$$

$\theta$ is the fraction of adsorption sites on the gas-facing surface that hold H.
Only H adsorbs, and there is one type of site. The subsurface has no coverage; it is
the concentration

$$
C = \frac{K}{\sqrt{K_{eq}}}\,\frac{\theta}{1-\theta} = K\,g(\theta),
\qquad K = K_s\ \text{(metal) or}\ K_{ox}\ \text{(oxide)},
$$

computed inside `g_theta` (`calculations/surface_kinetics.py`). This is section 3's
$c_s = c^*\theta_H/\theta_v$ with $c^* = K/\sqrt{K_{eq}}$: steps 3 and 4 are taken to be
at equilibrium.

## 6.3 Recombination is surface recombination

$k_{recomb}\theta^2$ with $k_{recomb} = k_{diss}/K_{eq}$ is $2\mathrm{H^*} \to \mathrm{H_2(g)}$.
It needs two neighbouring filled sites. No subsurface term appears. Because $C$ is tied to
$\theta$, the same term can be written in the bulk form,

$$
k_{recomb}\,\theta^2 = \frac{k_{diss}}{K^2}\,C^2\,(1-\theta)^2
\;\longrightarrow\; k_r^b\,C^2 \quad (\theta \ll 1),
$$

which is why it matches Table 5 at low coverage.

## 6.4 The hidden partition, and the 316L constants

$c^* = K/\sqrt{K_{eq}}$ is never an input; it follows from two constants taken from
separate sources. In the 316L study:

- $K_{eq,metal,ref} = 1.52\times10^{-3}$ Pa$^{-1}$, which is numerically $K_{s,ref}^2 = 0.039^2$;
- $H_{eq,metal} = 19{,}297 \approx 2H_s = 2\times9{,}648$ J mol$^{-1}$.

So $c^* = 1.000$ mol m$^{-3}$ at every temperature, and $\Delta H_{ss} = H_s - \tfrac12 H_{eq} \approx 0$.
The metal's $k_{recomb} = 1.729\times10^{-8}$ equals the bulk-form $k_d/K_s^2 = 1.728\times10^{-8}$,
and $E_{recomb} = 35.7$ kJ mol$^{-1} = E_d - 2H_s$.

The 316L metal "coverage" is therefore the subsurface concentration in disguise:
$\theta/(1-\theta) = C/(1\ \mathrm{mol\,m^{-3}})$. The other two studies label these
constants PLACEHOLDER.

## 6.5 Putting Cantera results into Model 2 (bare metal)

| Continuum parameter | From Cantera |
|---|---|
| `k_diss_metal` | $k_d$ (times 2 if the continuum counts atoms) |
| `K_eq_metal` | $K/P^\circ$ |
| `E_diss_metal` | intrinsic $E_a$ of $k_d$ |
| `H_eq_metal` | $\Delta H^\circ$ (**negative**, exothermic) |
| impurity case | `k_diss_metal` $\times f$ **and** `K_eq_metal` $\times f$, where $f = 1/(1+b)^2 = \theta_v^2/(1-\theta_H)^2$ |

- $1 + b = (1-\theta_H)/\theta_v$ comes from the Cantera steady state of that case. If the
  mechanism has two-site species, take $\theta_v$ from the rate, $\sqrt{k_{d,eff}/k_d}$.
- With this mapping the continuum's rate law equals Cantera's term by term, and
  $g(\theta)$ becomes $\theta_H/(\theta_v\sqrt{K_{eq}})$.
- **Do not use the observed $\Delta E_a$ and $A/A_{ref}$.** They describe $k_d\theta_v^2$,
  which already contains the site blocking, so the continuum's $(1-\theta)^2$ would
  count it twice. Changing `k_diss` alone also changes $k_{recomb}$, which impurities
  do not affect.

## 6.6 Other points in the current code

- **The exit surface is ideal.** $J = \beta(g(\theta) - \sqrt{P_{down}})$: recombination on the
  downstream side is taken to be instant.
- **The resistance split mixes units.** `smooth_surface_resistance` returns
  $1/(k_{diss}P(1-\theta)^2)$ in m$^2$ s mol$^{-1}$, while $1/\alpha$ and $1/\beta$ are in
  Pa$^{1/2}$ m$^2$ s mol$^{-1}$, and the code adds them. In the docs' L1+L6 example the
  code reports 2.8% surface ("metal-limited"); on a consistent $\sqrt{P}$ basis the
  surface carries 89%, which reproduces the observed 9.16$\times$ drop in flux. The
  fluxes themselves are correct.
- **The Cr$_2$O$_3$ gas-facing surface needs its own constants.** The Ni, Fe and Cu
  mechanisms cannot supply them.

# 7. Adding a kinetic absorption step

## 7.1 Bare metal

The unknowns are $\theta$ and $c_s$. Counting in mol H, so $k_{diss} = 2k_d$:

$$
\begin{aligned}
J_{surf} &= k_{diss}\left[P_{up}(1-\theta)^2 - \theta^2/K_{eq}\right]\\
J_{abs}  &= k_a\,\theta - k_b\,c_s\,(1-\theta), \qquad k_a/k_b = c^* = K_s/\sqrt{K_{eq}}\ \ \text{(detailed balance)}\\
J_D      &= \tfrac{D}{L}\,(c_s - c_{down})
\end{aligned}
$$

Setting $J_{abs} = J_D$ gives $c_s$ in closed form, just as $g(\theta)$ is today:

$$
c_s(\theta) = \frac{k_a\,\theta + (D/L)\,c_{down}}{k_b\,(1-\theta) + D/L},
\qquad
\text{then solve}\quad J_{surf}(\theta) = \tfrac{D}{L}\left(c_s(\theta) - c_{down}\right),\ \ \theta\in(0,1).
$$

**Fast limit.** For $k_a, k_b \to \infty$ at fixed $c^*$, $c_s \to c^*\theta/(1-\theta) = K_s\,g(\theta)$,
which is today's code exactly.

## 7.2 Oxide surface (Model 2)

Write $u = c_{s,ox}/K_{ox}$, and $\gamma = \alpha\beta/(\alpha+\beta)$ for the oxide and metal in
series, with $\alpha = D_{ox}K_{ox}/L_{ox}$ and $\beta = D_mK_s/L_m$:

$$
k_a/k_b = K_{ox}/\sqrt{K_{eq}},
\qquad
u(\theta) = \frac{k_a\,\theta + \gamma\sqrt{P_{down}}}{k_b\,K_{ox}(1-\theta) + \gamma},
\qquad
J_{surf}(\theta) = \gamma\left(u(\theta) - \sqrt{P_{down}}\right).
$$

The fast limit gives $u \to g(\theta)$, today's L2b+L6.

## 7.3 Impurities in the extension

Step 4 needs a truly empty site: $\theta_v = (1-\theta)/(1+b)$. So, in addition to section
6.5, $k_b \to k_b/(1+b)$, with $k_a$ unchanged. Detailed balance is preserved:
$k_a(1+b)/k_b = c^*(1+b) = K/\sqrt{K_{eq}/(1+b)^2}$.

## 7.4 Resistances in consistent units

Every resistance is a drop in $\sqrt{P}$-equivalent divided by $J$, so they add:

$$
R_{surf} = \frac{\sqrt{P_{up}} - g(\theta)}{J},\quad
R_{abs} = \frac{g(\theta) - c_s/K_s}{J},\quad
R_{metal} = \frac{c_s/K_s - \sqrt{P_{down}}}{J},\quad
\sum R = \frac{\sqrt{P_{up}} - \sqrt{P_{down}}}{J}.
$$

For the oxide: $R_{abs} = (g - u)/J$, $R_{ox} = (u - \sqrt{P_{int}})/J$, and
$R_{metal} = (\sqrt{P_{int}} - \sqrt{P_{down}})/J$.

## 7.5 How to recognise the limiting step

| Limiting step | Pressure exponent $n$ | Flux vs thickness $L$ | Coverage $\theta$ |
|---|---|---|---|
| diffusion | 0.5 | $\propto 1/L$ | $\approx \theta_{eq}$ |
| dissociation | 1 | independent | $\ll \theta_{eq}$ (starved) |
| absorption | 0.5, falling at high coverage | independent | $\approx \theta_{eq}$ (fills up) |

**Absorption-limited permeation has the same pressure exponent as diffusion-limited
permeation.** Only the thickness dependence separates them, so an exponent of 0.5 on its
own does not prove the coverage-modified Sieverts law.

**Error of the coverage-modified Sieverts law.** In the checks run so far, the flux error
is about $\tfrac12\,J/(k_a\theta)$: about 5% when $J/(k_a\theta) = 0.08$, and negligible
below 0.01.

## 7.6 Inputs and decisions

| Input | Bare metal | Oxide |
|---|---|---|
| $k_d$, $K$ | Cantera (`arrhenius_intrinsic.csv`) | MLMD Parts 1-2 on Cr$_2$O$_3$, or the "oxidized" rows of Table 5 |
| $k_a = \Gamma\nu_{fwd}$, $k_b = d\,\nu_{rev}$ | MLMD `part2_tst_rates.csv` (Hop A, 400-1200 K) | MLMD Hop A on Cr$_2$O$_3$ (not run yet) |
| $D$, $K_s$ or $K_{ox}$ | MLMD Part 3, or the literature | MLMD Part 3 (`Cr_oxide_supercell` has run) |

**Decisions to make.**

1. **Where $K_s$ comes from.** Either $K_s = (k_a/k_b)\sqrt{K/P^\circ}$, fully atomistic, or
   $K_s$ from the literature with $k_b = k_a\sqrt{K/P^\circ}/K_s$. For a paper, show both.
2. **How to combine Hop A sites.** Use the fastest site, a Boltzmann-weighted average,
   or report the range.
3. **Counting.** Choose atoms or molecules once, and convert every input at the boundary.

# 8. Conditions for the transfer to be valid

1. **Steady state**, or a surface that is fast compared with changes in the bulk.
2. **Absorption at local equilibrium** ($W_a \ll 1$), for the coverage-modified Sieverts law.
   Otherwise use section 7.
3. **A dilute subsurface and bulk** ($c \ll N_L$), so that $c^*$ is a constant and Sieverts applies.
4. **The same ideal Langmuir surface in both models:** one type of site, no H-H
   interaction, and dissociation needing two empty sites.
5. **Low coverage**, for Table 5's linear form: $c_s \ll c^*$ and $b \ll 1$. Otherwise keep the denominator.
6. **Impurities in equilibrium with the gas**, independent of the hydrogen flux, so that $b$ is a constant.
7. **The same surface state.** Cantera's $k_d$ is for the RMG facet, a clean single crystal;
   Table 5 was measured on real alloys.
8. **The same conventions:** H or $\mathrm{H_2}$ counting, $P^\circ = 1$ bar, and Cantera's
   $\Gamma$ wherever $k_r^s$ appears.
9. **A closed thermodynamic cycle:** $\Delta H_s - \tfrac12\Delta H_{ads} \approx \Delta H_{ss}$ (Hop A).

# 9. Validating the assumptions

The continuum model combines the surface, the oxide and the metal through a set of
assumptions. None of them can be proven in general. For each one, three kinds of
evidence can be shown:

1. **Limit checks.** A level reduces to the level below it when its new physics is
   switched off. These already exist, in the "Verified limit checks" section of each
   level page in `docs/theory/`.
2. **Relaxation tests.** The same wall is solved again with the assumption removed,
   and $|\Delta J/J|$ is reported across the operating window. The assumption is
   justified wherever the error stays below a stated tolerance, for example 5%.
3. **Experimental signatures.** Each assumption predicts a trend that can be measured.
   A measured trend that contradicts it disproves the assumption.

## 9.1 Assumption register

**Shared by all levels.**

| Assumption | Relaxation test | Holds when | Disproved by |
|---|---|---|---|
| steady state | time lag $t_{lag} = L^2/(6D_{eff})$ | $t_{lag} \ll$ exposure time | the flux still rising when it is measured |
| dilute hydrogen, Sieverts' law | -- | $C_L/N_L \ll 1$ | $n \neq 0.5$ on a diffusion-limited wall |
| absorption at local equilibrium | four-step solver (section 7) | $W_a \ll 1$ | $n = 0.5$, but $J$ independent of $L$ |
| ideal exit surface | the surface pair added on the exit side | exit exchange rate $\gg J$ | $J$ changes with the state of the exit face |

**Model 1 (Levels 1-5).**

| Level | Assumption | Relaxation test | Disproved by |
|---|---|---|---|
| L1 | dissociation and recombination at equilibrium | L1 against L1+L6 | $n \to 1$, or $J$ independent of $L$, on thin walls |
| L2a, L2b | the oxide carries $\mathrm{H_2}$ molecules (Henry) | the Model 2 oxide | $n \approx 0.5$ on an oxide-limited wall |
| L2b | oxide/metal interface at Sieverts equilibrium | an interface kinetic step, with its rate swept | $J$ independent of oxide thickness although the oxide limits |
| L3 | defect paths in parallel and independent (no spreading under a pinhole) | a 2D unit cell: one pinhole above a metal slab | coating effectiveness vs coverage disagrees with the area fractions |
| L4 | grain boundaries as a parallel path; Oriani trapping | 1D diffusion with traps, run to steady state | see section 9.2 |
| L5 | each defect path sits over its own patch of metal | the 2D unit cell, with $D_{eff}$ in the metal | -- |

**Model 2 (the Level 6 family).**

| Assumption | Relaxation test | Disproved by |
|---|---|---|
| ideal Langmuir surface: one site type, $(1-\theta)^2$, only H adsorbs | the Cantera microkinetics, mapped with $f$ (section 6.5) | the pressure dependence of $J$ disagrees with the Cantera coverages |
| coverage-modified Sieverts law | four-step solver (section 7) | a high-coverage, low-flux, thickness-independent regime |
| the oxide carries H atoms (Sieverts) | the Model 1 oxide | $n \approx 1$ with $J \propto 1/L_{ox}$ |
| the pinhole metal surface is kinetic | `use_sieverts_pinhole` on and off | -- |

The `use_sieverts_pinhole` flag is read only by `Surface_proposal.ipynb`, not by the
model wrappers, so the last test needs the flag passed through first.

## 9.2 Level 4: trapping and the steady-state flux

At steady state, traps in local equilibrium with the lattice do not change the flux.
With $\partial C/\partial t = 0$, the diffusion-trapping equation reduces to

$$
D_L\,\frac{\partial^2 C_L}{\partial x^2} = 0
\quad\Rightarrow\quad
J = \frac{D_L\left(C_{L,up} - C_{L,down}\right)}{L},
$$

which contains no trap parameter. Traps change how long the wall takes to reach this
state (the time lag), not the state itself.

Level 4 instead computes (`calculate_defective_metal_flux` in
`calculations/permeation_calc.py`):

$$
J = \frac{D_{eff}\,K_s\left(\sqrt{P_{up}} - \sqrt{P_{down}}\right)}{L},
\qquad
D_{eff} = \frac{D_L}{1 + \sum_i N_{T,i}K_i/N_L}.
$$

The two agree only if $D_L$ (`D_ref`) is the trap-free lattice diffusivity **and** $K_s$
is the total solubility, lattice plus trapped. Then $D_{eff}K_s = D_L K_{s,L}$ in the
dilute-trap limit.

If $D$ and $K_s$ both come from a time-lag permeation experiment, both are apparent
values, and their product is already the trap-free permeability. The code then applies
trapping a second time. In the example on the Level 4 docs page the trapping term is
0.354, so the flux comes out lower by a factor of 1.354 (26%).

**Tests.**

- **Numerical:** run 1D diffusion with Oriani traps to steady state, with and without
  traps. The steady flux should be identical; only the time lag should differ.
- **Experimental:** annealed and cold-worked samples of the same alloy should give the
  same steady flux but different time lags.
- **First check:** what the source of each study reports for `D_ref` and `K_s_ref`: a
  lattice or an apparent diffusivity, and a lattice or a total solubility.

## 9.3 Decisive experiments

| Measurement | What it decides |
|---|---|
| $J$ vs $P$, bare metal | L1 vs L1+L6: $n = 0.5$ (diffusion) or towards 1 (surface) |
| $J$ vs $P$ **and** $J$ vs oxide thickness, oxidized metal | Model 1 vs Model 2. On its own, $n = 1$ fits both a molecular oxide (Model 1, $J \propto 1/L_{ox}$) and a slow surface (Model 2, $J$ independent of $L_{ox}$); the thickness dependence separates them |
| $J$ vs metal thickness | diffusion-limited ($J \propto 1/L$) vs surface- or absorption-limited ($J$ independent of $L$) |
| apparent activation energy | $E_D + \Delta H_s$ (bulk) vs $E_{diss}$ (surface) |
| impurity on and off | whether the surface takes part at all |
| steady $J$ and time lag, different microstructures | Level 4 trapping (section 9.2) |

**Use data that was not used for fitting.** The study constants were calibrated on the
sources behind the `Guo_etal_2025`, `Fuerst_et_al_2024` and `incoloy802_cr2o3` folders,
so reproducing those fluxes is circular. A different pressure range, wall thickness or
oxide state, or another dataset, is needed.

## 9.4 The numerical experiments as a notebook

Each relaxation test is one block of `assumption_tests.ipynb`. It solves the wall with
and without the assumption over a grid of $T$, $P$, $L$, $L_{ox}$ and $b$, and reports
a validity window.

| Test | Relaxed model | Status |
|---|---|---|
| surface equilibrium (L1) | L1+L6 | exists |
| Henry vs Sieverts oxide (L2, L6) | Model 1 vs Model 2 oxide | exists |
| absorption equilibrium (all levels) | four-step solver (section 7) | new |
| ideal exit surface (all levels) | exit surface pair | new |
| interface equilibrium (L2b) | interface kinetic step | new |
| independent parallel paths (L3, L5) | 2D pinhole unit cell | new |
| trapping in the steady flux (L4) | 1D diffusion-trapping transient | new |
| ideal Langmuir surface (L6) | Cantera coverages with the $f$ mapping | new |

The output has one row per assumption: the largest $|\Delta J/J|$ over the grid, and the
part of the grid where it exceeds the tolerance.

# 10. A separate module that calls the existing code

Proposed layout for the module, with no edits to `calculations/`. This note itself sits
in `docs/surface_subsurface/`; the code would go in a folder of its own, for example
`surface_subsurface/` at the repository root.

```text
surface_subsurface/
  four_step.py            bare-metal and oxide solvers with the absorption pair;
                          resistances in sqrt(P) units
  relaxed.py              the other relaxed models: exit surface pair, interface step,
                          1D trap transient, 2D pinhole cell
  inputs.py               read the Cantera and MLMD outputs; unit and counting conversions
  checks.py               fast-limit tests, cycle closure, W_s and W_a
  assumption_tests.ipynb  the relaxation tests of section 9.4
  analysis.ipynb          figures and tables
```

- **Existing functions it calls, read-only.** From `calculations.surface_kinetics`:
  `surface_flux`, `g_theta`, `solve_steady_state_flux_L1L6` and
  `solve_steady_state_flux_direct`, as the fast-limit references. From
  `calculations.config.model_config`: `get_surface_kinetics_from_config`, for the
  existing study constants.
- **A study must be named before importing `calculations`.** This folder has no
  `study_config.py`, so set `MHI_STUDY=Application/study_config.py` (or another study
  folder) before the import, or it raises.
- **Branch.** The folder works on `main` or on a branch. A branch keeps `main` clean
  until the results show whether anything should move into `calculations/`.

# Appendix A. Checks already run

These were run with throwaway scripts against the existing code and data. The module in
section 10 should reproduce them as its tests.

**A.1 The continuum's Level 6 equals the master equation.** The relative difference
is $2\times10^{-16}$.

**A.2 The impurity mapping is exact.** For $\theta_H = 0.01$, 0.3 and 0.8, the full
competitive-Langmuir rate equals the continuum rate with $f$, and $c_s$ agrees to all
printed digits.

**A.3 The bare-metal absorption step reduces to today's L1+L6** (docs example: $\theta = 0.5736$, $J = 4.722\times10^{-7}$):

| $k_b$ (m s$^{-1}$) | $\theta$ | $J/J_{today}$ |
|---|---|---|
| $10^{6}$ | 0.5736 | 1.000 |
| $10^{-2}$ | 0.5736 | 1.000 |
| $10^{-5}$ | 0.583 | 0.957 |
| $10^{-6}$ | 0.646 | 0.682 |
| $10^{-7}$ | 0.816 | 0.164 |

**A.4 The oxide absorption step reduces to today's L2b+L6:** $J/J_{today} = 1.0000$,
0.9998, 0.828 and 0.046 for $k_b = 10^{6}$, $10^{-6}$, $10^{-9}$ and $10^{-11}$ m s$^{-1}$.

**A.5 Absorption-limited signature** (fast surface, fast diffusion, slow entry):

| $P$ (Pa) | $\theta$ | $n$ | $J(0.1\,\mathrm{mm})/J(1\,\mathrm{mm})$ |
|---|---|---|---|
| $10^2$ | 0.000 | 0.500 | 1.000 |
| $10^4$ | 0.003 | 0.498 | 1.000 |
| $10^6$ | 0.031 | 0.485 | 1.000 |
| $10^8$ | 0.240 | 0.380 | 1.000 |

**A.6 Ni separation numbers** (1 bar, 1 mm, MLMD vibrational route; Hop A range across
sites; mol H m$^{-2}$ s$^{-1}$):

| $T$ (K) | $J_{D,max}$ | $2R_{0,surf}$ | $R_{0,sub}$ | largest $W$ |
|---|---|---|---|---|
| 800 | $6.3\times10^{-6}$ | 310 | 0.56 to 280 | $1\times10^{-5}$ |
| 1000 | $3.0\times10^{-5}$ | 490 | 4.4 to 800 | $7\times10^{-6}$ |
| 1200 | $8.4\times10^{-5}$ | 520 | 16 to 1460 | $5\times10^{-6}$ |

**A.7 Ni thermodynamic cycle.** Cantera gives $\tfrac12\Delta H_{ads} = -46.7$ kJ mol$^{-1}$.
MLMD gives $\Delta H_s = 12.5$ (octahedral) and 32.8 (tetrahedral) kJ mol$^{-1}$. The implied
$\Delta H_{ss}$ is 59.2 and 79.4 kJ mol$^{-1}$, against MLMD Hop A values of 63.8 and 84.0.
The cycle closes to within about 5 kJ mol$^{-1}$.

**A.8 Table 5, H$_2$-Ni row.** $\Delta H_s = (30 - 2.5)/2 = 13.75$ kJ mol$^{-1}$, and
$K_s(1000\ \mathrm{K}) = 0.089$ mol $\mathrm{H_2}$ m$^{-3}$ Pa$^{-1/2}$ (0.178 in mol H).
