# Foundations: chemical potential and the driving force

Every quantity the model tracks — pressure, surface coverage, dissolved
concentration — is a different way of writing the same thing: the chemical
potential of hydrogen at some point in the wall. Flux flows down that potential,
and each level of the hierarchy is a different assumption about which step in the
chain is slow enough to matter.

This chapter establishes that chain, because it is what licenses the $\sqrt{P}$
driving force the code actually uses.

## 1. The gas phase sets the boundary condition

For molecular hydrogen in the gas,

$$\mu_{H_2} = \mu^{\circ}_{H_2} + RT \ln(P/P_0)$$

with $\mu^{\circ}_{H_2}$ the standard chemical potential at reference pressure
$P_0$, and $R = 8.314$ J/mol/K. Raising the pressure raises the thermodynamic
availability of H₂ to react — it is the supply side of the whole problem.

## 2. Surface coverage and the Langmuir isotherm

Hydrogen adsorbs as atoms, so an adsorbed atom is in equilibrium with *half* a
gas molecule:

$$\mu_{H,\text{ads}} = \tfrac{1}{2}\mu_{H_2}$$

Statistical thermodynamics gives the adsorbed atom's potential in terms of the
fractional coverage $\theta$:

$$\mu_{H,\text{ads}} = \mu^{\circ}_{H,\text{ads}} + RT \ln\!\left[\frac{\theta}{1-\theta}\right]$$

The $\theta/(1-\theta)$ activity is what makes a surface saturate. As
$\theta \to 0$ there are empty sites everywhere and the potential is low; as
$\theta \to 1$ the last few sites become very expensive and the potential
diverges.

Equating the two expressions and collecting the standard-state terms into an
equilibrium constant $K_{eq}$ yields the Langmuir isotherm in the form the code
uses:

$$\frac{\theta}{1-\theta} = \sqrt{K_{eq} P}
\qquad\Longleftrightarrow\qquad
\theta = \frac{\sqrt{K_{eq}P}}{1 + \sqrt{K_{eq}P}}$$

The square root is not cosmetic. It is the signature of a diatomic molecule
dissociating into two atoms, and it propagates all the way to the flux law.

### Virtual pressure

Inverting the isotherm gives the pressure that *would* be in equilibrium with a
given coverage:

$$P_{\text{virtual}}(\theta) = \frac{1}{K_{eq}}\left[\frac{\theta}{1-\theta}\right]^2$$

This is the single most useful construct in the surface-kinetics level. It
converts a coverage — which has no pressure units — into an equivalent pressure
that can be compared directly against the gas:

- $P > P_{\text{virtual}}$: net adsorption; the surface is being filled
- $P < P_{\text{virtual}}$: net desorption
- $P = P_{\text{virtual}}$: surface equilibrium, and surface kinetics have
  dropped out of the problem

In the code this appears as `g_theta`, which returns $\sqrt{P_{\text{virtual}}}$
rather than $P_{\text{virtual}}$ because every downstream flux law wants the
square root anyway:

$$g(\theta) = \frac{\theta}{(1-\theta)\sqrt{K_{eq}}} = \sqrt{P_{\text{virtual}}(\theta)}$$

## 3. Dissolved concentration and the Sieverts connection

At the surface–subsurface boundary, dissolved hydrogen equilibrates with adsorbed
hydrogen, $\mu_{H,\text{dissolved}} = \mu_{H,\text{ads}}$. Writing the dissolved
potential in the ideal-dilute form $\mu^{\circ} + RT\ln(C/C_0)$ and solving for
$C$ collects the standard states into the solubility constant:

$$C_{ox,\text{up}} = \frac{K_{ox}}{\sqrt{K_{eq}}}\left[\frac{\theta}{1-\theta}\right]$$

which is exactly $K_{ox}\sqrt{P_{\text{virtual}}}$:

$$K_{ox}\sqrt{P_{\text{eff}}}
= K_{ox}\sqrt{\frac{1}{K_{eq}}\left[\frac{\theta}{1-\theta}\right]^2}
= \frac{K_{ox}}{\sqrt{K_{eq}}}\left[\frac{\theta}{1-\theta}\right]$$

So the surface-kinetics levels do not replace Sieverts' law — they replace the
*pressure* that Sieverts' law is evaluated at. The interface concentration is
still $K\sqrt{P}$; it is just $P_{\text{virtual}}$ rather than $P_{\text{gas}}$.
That substitution is the whole content of Level 6, and it is why the two
equilibrium descriptions coexist without contradiction — see
{doc}`equilibrium-models`.

| Condition | $P_{\text{eff}}$ vs $P_{\text{gas}}$ | Meaning |
|---|---|---|
| equilibrium | $P_{\text{eff}} = P_{\text{gas}}$ | surface kinetics not limiting |
| surface-limited, upstream | $P_{\text{eff}} < P_{\text{gas}}$ | adsorption cannot keep up; surface starved |
| surface-limited, downstream | $P_{\text{eff}} > P_{\text{gas}}$ | desorption cannot keep up; hydrogen backs up |

## 4. The surface rate law as a potential difference

The net dissociation flux is the difference of forward and reverse rates,

$$J_{\text{surface}} = k_{\text{diss}} P (1-\theta)^2 - k_{\text{recomb}} \theta^2$$

which is `surface_flux` in the code. The $(1-\theta)^2$ and $\theta^2$ factors
are combinatorial: dissociation needs two adjacent empty sites, recombination
needs two adjacent filled ones.

Detailed balance ties the two rate constants to the isotherm. Setting
$J_{\text{surface}} = 0$ recovers the Langmuir isotherm exactly, which requires

$$K_{eq} = \frac{k_{\text{diss}}}{k_{\text{recomb}}}$$

and the code enforces this rather than storing $k_{\text{recomb}}$
independently — `surface_flux` computes `k_recomb = k_diss / K_eq` internally, so
the isotherm and the rate law cannot drift apart.

Factoring the rate law reveals its thermodynamic content:

$$J_{\text{surface}} = k_{\text{recomb}}(1-\theta)^2 K_{eq}\left[P - P_{\text{virtual}}(\theta)\right]$$

The flux is proportional to the *departure from surface equilibrium*, measured as
a pressure difference — the Onsager form $J = L\,\Delta\mu$ with
$L_{\text{surf}}$ absorbing the kinetic prefactor. When the surface is fast,
$\theta$ adjusts until $P_{\text{virtual}} \to P$ and this term vanishes from the
resistance budget.

## 5. The equilibrium constant is the thermodynamic anchor

$K_{eq}$ carries the temperature dependence of the surface:

$$K_{eq}(T) = K_{eq,\text{ref}}\exp\left[\frac{-\Delta H_{\text{ads}}}{R}\left(\frac{1}{T} - \frac{1}{T_{\text{ref}}}\right)\right]$$

This is the reference-temperature Arrhenius form used throughout the project;
`arrhenius` implements it, and every temperature-dependent property goes through
it rather than through a bare $\exp(-E/RT)$ with a pre-exponential.

```{important}
$T_{\text{ref}}$ is a property of the **material and the study**, not a global
constant. The values in the configuration are:

| Material | $T_{\text{ref}}$ [K] |
|---|---|
| 316L metal surface (`Guo_etal_2025_316L`) | 673 |
| Hastelloy N metal surface (`fuerst_etal_2024`) | 965 |
| Incoloy 802 metal surface (`incoloy802_cr2o3`) | 965 |
| Cr₂O₃ oxide surface (all three studies) | 1623 |

Reference-temperature Arrhenius is used precisely because a measured value at a
stated temperature is more trustworthy than an extrapolated pre-exponential, so
quoting the wrong $T_{\text{ref}}$ silently rescales the property.
```

## 6. The chain, end to end

```text
μ_H₂(gas, upstream)          ← overall driving force →      μ_H₂(gas, downstream)
        │                                                             ▲
        ▼  finite-rate dissociation (Level 6)                         │
μ_H₂(virtual, θ_up)          ← gradient through oxide →     μ_H₂(virtual, θ_down)
        │                                                             ▲
        ▼  local equilibrium (always assumed)                         │
μ_H(dissolved, C_up)         ← Fickian diffusion →          μ_H(dissolved, C_down)
```

Each variable in the code is one line of that diagram:

| Code quantity | Chemical potential it represents |
|---|---|
| `P_upstream`, `P_downstream` | $\mu_{H_2} = \mu^{\circ} + RT\ln(P/P_0)$ |
| `theta` | $\mu_{H,\text{ads}} = \mu^{\circ} + RT\ln[\theta/(1-\theta)]$ |
| `C_up`, `C_down` | $\mu_{H,\text{diss}} = \mu^{\circ} + RT\ln(C/C_0)$ |
| `k_diss`, `k_recomb` | kinetic coefficients in $J = L\,\Delta\mu$ |
| `K_eq` | $\exp(-\Delta G^{\circ}/RT)$ — sets equilibrium coverage |

## 7. Why the code can use a √P driving force

The total potential drop across the wall is

$$\Delta\mu_{\text{total}} = \frac{RT}{2}\ln\!\left(\frac{P_{\text{up}}}{P_{\text{down}}}\right)
= RT \ln\!\left(\frac{\sqrt{P_{\text{up}}}}{\sqrt{P_{\text{down}}}}\right)$$

and it decomposes across the series steps:

$$\Delta\mu_{\text{total}} = \Delta\mu_{\text{surface,up}} + \Delta\mu_{\text{oxide}} + \Delta\mu_{\text{metal}} + \Delta\mu_{\text{surface,down}}$$

with the surface terms present only when dissociation is finite-rate. For a
Sieverts-law solid the concentration is $K\sqrt{P}$, so a gradient in
$\sqrt{P}$ *is* a gradient in concentration, and hence proportional to the
gradient in $\mu$. This is why the metal flux laws are written with
$\sqrt{P_{\text{up}}} - \sqrt{P_{\text{down}}}$ rather than $P_{\text{up}} -
P_{\text{down}}$.

The oxide is the interesting case, because **which law applies to it depends on
the level**. Below Level 6 there is no dissociation step, so hydrogen is assumed
to dissolve and diffuse through the oxide as intact molecules — Henry's law,
linear in $P$. Once Level 6 models dissociative adsorption explicitly, the
species entering the oxide is atomic, and its transport becomes Sieverts-like in
$\sqrt{P}$. The species assumption changes with the level, and so does the
pressure scaling. {doc}`equilibrium-models` sets out both forms, where each
applies, and an important unresolved consequence for the shared solubility
constant.
