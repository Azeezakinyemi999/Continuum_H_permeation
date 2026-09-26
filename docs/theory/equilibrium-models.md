# Equilibrium models: Sieverts, Henry and Langmuir

The model uses three different equilibrium descriptions at once. That looks
contradictory until you notice that each applies at a **different location** or
to a **different chemical species**, and that is the whole content of this
chapter.

| Description | Where it applies | Species | Concentration |
|---|---|---|---|
| Langmuir kinetics | gas–oxide surface | H₂ → 2 H(ads) | $\theta/(1-\theta) = \sqrt{K_{eq}P}$ |
| Henry's law | inside the oxide, **below Level 6** | H₂ (intact) | $C = K_{ox} P$ |
| Sieverts' law | inside the oxide **with Level 6**, and always in the metal | H (atomic) | $C = K\sqrt{P}$ |

## Why Langmuir and Sieverts do not conflict

The question that prompts this chapter is whether it is legitimate to assume
Sieverts' *equilibrium* — a thermodynamic condition — while simultaneously
imposing Langmuir surface coverage as a *non-equilibrium, rate-limited*
condition. They are not in conflict, because they sit at different interfaces.

**At the gas-facing surface**, dissociative adsorption is explicitly *not* in
equilibrium. Coverage $\theta$ is not given by the isotherm evaluated at
$P_{\text{up}}$; it is whatever value makes the fluxes balance,

$$J_{\text{surface}} = J_{\text{oxide}} = J_{\text{metal}}$$

at steady state. This matters because dissociative adsorption on oxides carries a
high activation barrier, so it is a plausible rate-limiting step. Modelling it as
instantaneous would discard the physics Level 6 exists to capture.

**At the buried oxide–metal interface**, local equilibrium *is* assumed:

$$C_{\text{int}} = K\sqrt{P_{\text{int}}}$$

The justification is a rate argument, not a thermodynamic one. Atomic hydrogen
exchange across a buried solid–solid interface is fast compared with bulk
diffusion through either layer, so that interface never limits the flux and there
is nothing to gain from modelling its kinetics. Only interfaces that can be slow
get a kinetic treatment; the gas-facing surface can be slow, the buried one
cannot.

So the resistance chain reads:

```text
gas (P_up)
  │  Langmuir kinetics — NOT in equilibrium, can be rate-limiting
  ▼
surface coverage θ
  │  diffusion through the oxide
  ▼
oxide–metal interface (P_int)
  │  Sieverts equilibrium — assumed fast, never limiting
  ▼
diffusion through the metal
  │
  ▼
gas (P_down)
```

Each arrow is a resistance in series. The model's job is to find which one
dominates, and `regime` reports the answer.

## The oxide changes species with the level

Below Level 6 there is no dissociation step in the model at all. Hydrogen is
therefore assumed to dissolve and diffuse through the oxide **as intact H₂
molecules**, which is Henry's law: concentration linear in pressure, and so flux
linear in the pressure *difference*.

$$C = K_{ox}P \qquad\Longrightarrow\qquad
J_{ox} = \frac{D_{ox}K_{ox}}{L_{ox}}\left(P_{\text{up}} - P_{\text{int}}\right)$$

This is `molecular_diffusion_flux`, used by Levels 2a, 2b and 3.

Once Level 6 is switched on, the gas-facing surface dissociates H₂ explicitly, so
what enters the oxide is **atomic hydrogen**. Transport then takes the Sieverts
form, linear in $\sqrt{P}$:

$$C = K_{ox}\sqrt{P} \qquad\Longrightarrow\qquad
J_{ox} = \alpha\left(\sqrt{P_{\text{eff}}} - \sqrt{P_{\text{int}}}\right),
\qquad \alpha \equiv \frac{D_{ox}K_{ox}}{L_{ox}}$$

This is `oxide_flux`, used by every Level 6 combination, with
$\sqrt{P_{\text{eff}}} = g(\theta)$ from {doc}`foundations`.

The species assumption is a deliberate consequence of what the level models. It
is not an inconsistency in itself: a model with no dissociation step has no
atomic hydrogen to transport, and a model with one has no intact molecules left.

## Where each law is enforced in the code

| Function | Layer | Law |
|---|---|---|
| `sieverts_concentration` | metal surface | $C = K_s\sqrt{P}$ |
| `calculate_simple_metal_flux` | metal bulk | Sieverts + Fick |
| `molecular_diffusion_flux` | oxide, no L6 | Henry + Fick |
| `oxide_flux` | oxide, with L6 | Sieverts, via $g(\theta)$ |
| `g_theta` | gas–oxide surface | Langmuir isotherm, inverted |
| `surface_flux` | gas–oxide surface | Langmuir–Hinshelwood rate law |
| `solve_interface_pressure` | oxide–metal interface | Sieverts equilibrium + flux balance |

## The pinhole special case: a third surface to worry about

A pinhole is a hole clean through the oxide, so it exposes **bare metal directly
to the gas**. That creates a second gas-facing surface, with its own dissociation
kinetics, which need not be as slow as dissociation on the oxide. Level 6
therefore carries a separate pair of metal-surface constants,
`k_diss_metal` and `K_eq_metal`, distinct from the oxide's `k_diss` and `K_eq`.

Whether to model that surface kinetically at all is a configuration choice:

```python
OXIDE_DEFECTS['use_sieverts_pinhole']   # False in all three studies
```

- `False` — the pinhole's metal surface gets finite-rate kinetics from
  `METALS[...]['surface_kinetics']`.
- `True` — metal dissociation is assumed fast, so the pinhole sits at the
  Sieverts limit and `k_diss_metal` is not needed.

```{note}
**Neither model wrapper reads this flag.** `level5_model_wrapper` and
`level5L6_model_wrapper` do not take it from the parameter dictionary, and it is
absent from `DEFAULT_PARAMS_LEVEL5L6` for that reason. The only consumer is
`Surface_proposal.ipynb`, which reads it from
`build_simulation_config()['oxide_defects']`.

So changing it has **no effect on the sensitivity analysis or on any wrapper
call** — only on that notebook. If you need the Sieverts-limit pinhole in a
wrapper-driven study, it has to be plumbed through the parameter dictionary
first.
```
