# Validating the Trapping Block Against TDS

**A zero-free-parameter test of `calculate_effective_diffusivity_trapping`, plus a typeset reference for the permeation / accumulation / embrittlement equation set.**

Code: [`calculations/defective_metal.py`](calculations/defective_metal.py) · Preview this file with ⇧⌘V in VS Code to render the math.

---

## 1. The mechanism

Hydrogen dissolved in an alloy does not sit in one kind of site. Most of it occupies ordinary interstitial positions in the lattice, where it hops freely and carries flux. Some fraction is caught at **traps** — dislocation cores, carbide and oxide interfaces, grain boundaries, vacancies — which are energetically deeper by a binding energy $E_b$ of roughly 20–60 kJ/mol. Hydrogen sitting in a trap is immobile on the timescale of a hop, so it contributes to *inventory* but not to *transport*.

Two consequences follow, and they pull in opposite directions:

- Traps **raise** the total hydrogen content of the material (they add storage capacity).
- Traps **lower** the apparent diffusivity (they remove carriers from the mobile population).

The second effect is what the trapping block computes, and it is what makes an alloy's measured diffusivity depend on its cold work and heat treatment rather than on composition alone.

---

## 2. The problem this test solves

The trapping block takes **two inputs per trap population** and returns **one number**:

$$
\underbrace{N_{T,i},\; E_{b,i}}_{\text{inputs}}
\;\longrightarrow\;
\frac{D_{\text{eff}}}{D_{\text{lattice}}}
= \frac{1}{1+\displaystyle\sum_i \frac{N_{T,i}}{N_L}\exp\!\left(\frac{E_{b,i}}{RT}\right)}
$$

Most papers obtain $N_T$ and $E_b$ **by fitting them to the permeation curve**. Once that is done, the statement "the model reproduces the permeation data" carries no information — two adjustable knobs were tuned to hit one number. The comparison is circular.

**Thermal desorption spectroscopy (TDS) breaks the circularity.** It measures $N_T$ and $E_b$ from an entirely different physical signal: the rate at which hydrogen leaves a charged specimen during a temperature ramp. That signal has nothing to do with steady-state flux. So TDS supplies the inputs for free, and the ratio $D_{\text{eff}}/D_{\text{lattice}}$ becomes a genuine prediction with no free parameters left to adjust.

---

## 3. The four steps

### Step 1 — TDS gives you `trap_list`

Charge a specimen to saturation, then heat it at a constant rate $\beta$ while recording the desorption flux against temperature. The spectrum shows distinct peaks; each peak is one trap population.

**Trap density, from peak area.** Integrate the desorption rate across the peak:

$$
C_{T,i} = \int \frac{\mathrm{d}C}{\mathrm{d}t}\,\mathrm{d}t
\qquad\Longrightarrow\qquad
N_{T,i} \approx C_{T,i}\,N_A \quad (\text{for } \theta \to 1)
$$

This is the `density` field.

**Binding energy, from peak temperature.** Repeat the ramp at 3–4 heating rates, record each peak temperature $T_p$, and apply a **Kissinger** analysis:

$$
\ln\!\left(\frac{\beta}{T_p^{2}}\right) = -\frac{E_a}{R\,T_p} + \text{const}
$$

Plotting $\ln(\beta/T_p^2)$ against $1/T_p$ gives a straight line of slope $-E_a/R$.

> ⚠️ **Do not stop here.** Kissinger returns the **detrapping** activation energy, not the binding energy:
>
> $$
> E_a = E_b + E_m \qquad\Longrightarrow\qquad E_b = E_a - E_m
> $$
>
> where $E_m$ is the lattice migration energy (the Arrhenius slope of $D_{\text{lattice}}$). The code expects $E_b$. Feeding $E_a$ in directly overestimates $E_b$ by 10–20 kJ/mol, and since it enters through $\exp(E_b/RT)$ the error is exponential. This is the single most common mistake in this comparison.

### Step 2 — the model predicts the ratio

Pass `trap_list` in at the operating temperature and read off `reduction_factor`. No parameter has been touched.

```python
from calculations.defective_metal import calculate_effective_diffusivity_trapping

traps = [{'name': 'dislocations', 'binding_energy': 30e3, 'density': 1e26}]
r = calculate_effective_diffusivity_trapping(
        D_lattice=1e-11, temperature=573, trap_list=traps,
        lattice_concentration=1e-3, lattice_density=1e29)
r['reduction_factor']   # -> the prediction
```

### Step 3 — permeation measures the same ratio, independently

You need the numerator and the denominator separately, both from lag times:

$$
D_{\text{app}} = \frac{L^{2}}{6\,t_{\text{lag}}}
$$

$$
D_{\text{app}} = \frac{L^{2}}{6\,t_{\text{lag}}(\text{low } T)}
\qquad\qquad
D_{\text{lattice}} = \frac{L^{2}}{6\,t_{\text{lag}}(\text{high } T)}
$$

At high temperature the exponential collapses, $\exp(E_b/RT)\to 1$, and trapping switches itself off — so the *same specimen* hands you the denominator. That is the clean route: one specimen, one apparatus, systematic errors cancel between numerator and denominator. The alternatives are a fully annealed specimen or a literature single-crystal value, both of which reintroduce specimen-to-specimen scatter.

### Step 4 — compare two dimensionless numbers

TDS-predicted against permeation-measured. For a model of this class, agreement within a factor of ~2 is a genuine pass.

---

## 4. Why the temperature sweep is the real test

A match at one temperature could be luck. The **shape** cannot be.

Run at $E_b = 30$ kJ/mol and $N_T = 10^{26}\,\mathrm{m^{-3}}$ (so $N_T/N_L = 10^{-3}$, i.e. cold-worked):

| $T$ [K] | $K = \exp(E_b/RT)$ | $(N_T/N_L)\,K$ | $D_{\text{eff}}/D_{\text{lattice}}$ |
|--------:|-------------------:|---------------:|------------------------------------:|
| 373 | $1.590\times10^{4}$ | 15.898 | 0.0592 |
| 473 | $2.056\times10^{3}$ | 2.056 | 0.3272 |
| 573 | $5.431\times10^{2}$ | 0.543 | 0.6480 |
| 673 | $2.131\times10^{2}$ | 0.213 | 0.8244 |
| 773 | $1.065\times10^{2}$ | 0.106 | 0.9038 |
| 873 | $6.238\times10^{1}$ | 0.062 | 0.9413 |

A **17× suppression at 373 K, decaying to essentially nothing by 873 K.** The two parameters control different features of that curve:

- $N_T/N_L$ sets the **depth** of the suppression.
- $E_b$ sets the **temperature at which it turns on** — the position of the knee.

So an Arrhenius plot of the *measured* $D_{\text{app}}$ shows a kink:

$$
\text{above the knee:}\quad \frac{\mathrm{d}\ln D_{\text{app}}}{\mathrm{d}(1/T)} = -\frac{E_m}{R}
\qquad\qquad
\text{below the knee:}\quad \frac{\mathrm{d}\ln D_{\text{app}}}{\mathrm{d}(1/T)} = -\frac{E_m + E_b}{R}
$$

```
  ln D_app

    │
    │ ●●●                      slope = − E_m / R          (above the knee:
    │    ●●●                                               trap-free lattice)
    │       ●
    │        ⌄  ← knee, position set by E_b
    │         ●●
    │           ●●●            slope = − (E_m + E_b) / R   (below the knee:
    │              ●●●                                      trapping active)
    └─────────────────────────  1/T
```

TDS fixes both the knee location and the slope change **before** you look at the permeation data. That is falsifiable in a way a single flux number never is.

---

## 5. Two caveats before trusting a pass

**Oriani assumes low occupancy.** The code warns at $\theta > 0.9$ ([`defective_metal.py:1134`](calculations/defective_metal.py#L1134)). Heed it. Near saturation the local-equilibrium assumption fails and you need McNabb–Foster trapping kinetics instead. Note the tension: TDS charging is normally done *at* saturation, which is precisely the regime where the Oriani formula is least valid — so check $\theta$ at your **operating** concentration, not at the charging concentration.

**TDS and permeation must see the same specimen state.** Traps are microstructural, not compositional. Different heat treatment, different cold work, different oxide $\Rightarrow$ different $N_T$. Cut both specimens from the same plate with the same processing history, or the comparison means nothing.

### If the two disagree by more than ~3×

Usual culprits, in order of likelihood:

1. $E_m$ was not subtracted from $E_a$ (see the warning in Step 1).
2. Peak deconvolution assigned area to the wrong trap population.
3. A trap population TDS could not see, because it desorbs above the top temperature of the ramp.

---

# Appendix A — Symbols

| Symbol | Meaning | Units |
|---|---|---|
| $N_{T,i}$ | trap density, population $i$ | m⁻³ |
| $N_L$ | lattice interstitial site density (FCC $\approx 10^{29}$) | m⁻³ |
| $E_{b,i}$ | trap binding energy | J/mol |
| $E_a$ | detrapping activation energy (Kissinger output) | J/mol |
| $E_m$ | lattice migration energy | J/mol |
| $\theta_i$ | trap occupancy — Oriani valid only for $\theta \lesssim 0.9$ | – |
| $K_i$ | trap equilibrium constant, $\exp(E_{b,i}/RT)$ | – |
| $C_L,\;C_{T}$ | lattice (mobile) and trapped H concentration | mol/m³ |
| $D_{\text{lattice}}$ | intrinsic lattice diffusivity | m²/s |
| $D_{\text{eff}},\;D_{\text{app}}$ | effective / apparent diffusivity | m²/s |
| $K_S$ | Sieverts solubility constant | mol m⁻³ Pa⁻¹ᐟ² |
| $\Phi$ | permeability | mol m⁻¹ s⁻¹ Pa⁻¹ᐟ² |
| $V_H$ | partial molar volume of H | m³/mol |
| $\sigma_h$ | hydrostatic stress | Pa |
| $\beta$ | TDS heating rate | K/s |
| $R$ | gas constant, 8.314 | J mol⁻¹ K⁻¹ |

---

# Appendix B — Equation reference

## B.1 Trapping (Oriani local equilibrium)

Ratio of trapped to mobile hydrogen at each trap population:

$$
\frac{C_{T,i}}{C_L} = \frac{N_{T,i}}{N_L}\,K_i,
\qquad K_i = \exp\!\left(\frac{E_{b,i}}{RT}\right)
$$

Mobile fraction, from the mass balance $C_{\text{total}} = C_L + \sum_i C_{T,i}$:

$$
f_{\text{mobile}} = \frac{C_L}{C_{\text{total}}}
= \frac{1}{1+\displaystyle\sum_i \frac{N_{T,i}}{N_L}K_i}
$$

Only mobile hydrogen carries flux, so $D_{\text{eff}} = D_{\text{lattice}}\,f_{\text{mobile}}$:

$$
\boxed{\;D_{\text{eff}} = \frac{D_{\text{lattice}}}{1+\displaystyle\sum_i \frac{N_{T,i}}{N_L}\exp\!\left(\frac{E_{b,i}}{RT}\right)}\;}
$$

$$
\frac{D_{\text{eff}}}{D_{\text{lattice}}} \equiv \texttt{reduction\_factor}
$$

Critical trap density, at which $D_{\text{eff}} = D_{\text{lattice}}/2$:

$$
N_T^{*} = \frac{N_L}{K}
$$

Parallel grain-boundary path (Level 4, combined with trapping):

$$
D_{\text{eff}} = \frac{(1-f_{gb})D_{\text{bulk}} + f_{gb}\,\alpha\,D_{\text{bulk}}}{1+\sum_i N_{T,i}K_i/N_L}
$$

## B.2 Permeation

$$
\text{Fick:}\quad J = -D\,\frac{\partial C}{\partial x}
\qquad\qquad
\text{Sieverts:}\quad C_0 = K_S\sqrt{p_{\mathrm{H_2}}}
$$

$$
\text{Steady flux:}\quad J = \frac{\Phi\left(\sqrt{p_{\text{up}}}-\sqrt{p_{\text{down}}}\right)}{L}
$$

$$
\text{Series stack:}\quad R_{\text{tot}} = \sum_i \frac{L_i}{\Phi_i},
\qquad J = \frac{\Delta P_{\text{eff}}}{R_{\text{tot}}}
$$

$$
\text{Barrier value:}\quad \mathrm{PRF} = \frac{J_{\text{bare}}}{J_{\text{coated}}}
\qquad\qquad
\text{Lag time:}\quad t_{\text{lag}} = \frac{L^{2}}{6D}
$$

$$
\text{Regime test:}\quad J \propto p^{\,n},
\qquad
n =
\begin{cases}
0.5 & \text{diffusion-limited} \\
\to 1 & \text{surface-limited}
\end{cases}
$$

## B.3 Accumulation (Sofronis–McMeeking)

$$
\underbrace{\frac{\partial C_L}{\partial t} + \frac{\partial C_T}{\partial t}}_{\text{lattice + trapped storage}}
= \nabla\cdot\left(
\underbrace{D\nabla C_L}_{\text{Fickian}}
- \underbrace{\frac{D\,C_L V_H}{RT}\nabla\sigma_h}_{\text{stress-driven drift}}
\right)
$$

Hydrogen migrates *up* the hydrostatic stress gradient, concentrating at notch roots and crack tips.

## B.4 Embrittlement index

Empirical susceptibility ratios from paired tests in hydrogen and in air:

$$
\mathrm{EI} = \frac{\varepsilon_{\text{air}}-\varepsilon_{\mathrm{H}}}{\varepsilon_{\text{air}}}
\qquad\qquad
\mathrm{RRA} = \frac{\mathrm{RA}_{\text{air}}-\mathrm{RA}_{\mathrm{H}}}{\mathrm{RA}_{\text{air}}}
$$

$$
\text{Empirical fit:}\quad \mathrm{EI} = a + b\ln C_{\mathrm{H}}
\qquad (C_{\mathrm{H}} > C_{\text{th}})
$$

Mechanism-based decohesion law (Serebrinsky–Carter–Ortiz), cohesive strength falling with H coverage $\theta$:

$$
\sigma_c(\theta) = \sigma_{c0}\left(1 - 1.0467\,\theta + 0.1687\,\theta^{2}\right)
$$

---

# Appendix C — Where each quantity lives

| Quantity | Role | Source |
|---|---|---|
| $D_{\text{lattice}}$, $K_S$, $\Phi_{\text{oxide}}$ | **input** (intrinsic) | [`data/material_data.py`](data/material_data.py), [`data/oxide_properties.py`](data/oxide_properties.py) |
| $N_T$, $E_b$ | **input** (microstructural) | TDS — or fitted, which forfeits the test |
| $D_{\text{eff}}$, `reduction_factor` | **output** (effective) | [`calculations/defective_metal.py`](calculations/defective_metal.py) |
| $\Phi_{\text{eff}}$, PRF, regime | **output** (system) | [`calculations/classify_regime.py`](calculations/classify_regime.py) |
| parameter ranges | SA only | [`data/sensitivity_parameters.csv`](data/sensitivity_parameters.csv) |

The distinction in the first column is the important one: **intrinsic properties go in, effective properties come out.** The model does not derive $D_{\text{lattice}}$ or $K_S$ — it consumes them and predicts what a real, trapped, coated, defective wall will appear to have.
