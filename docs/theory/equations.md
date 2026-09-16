# Equation catalogue

The theory chapters derive each level's equations in the context that motivates
them. This page is the complement: a single consolidated list of every equation
in the model, each tagged with the file and function that implements it.

## The compiled catalogue

[**Model_Equations.pdf**](../_static/Model_Equations.pdf) — the full typeset
catalogue, organised as:

| Section | Contents |
|---|---|
| A | Bulk transport, Levels 1–5 |
| B | Surface kinetics, Level 6 |

Every equation carries its `file → function` location, so it can be traced
directly to the implementation.

The LaTeX source is
[`latex/Model_Equations.tex`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/latex/Model_Equations.tex)
and remains the single source of truth for the typeset form; the PDF here is a
build of it. Regenerate the PDF from that source rather than editing it.

```{note}
The catalogue and the theory chapters are generated from different processes —
the chapters are re-verified against the running code by
`docs/_tools/worked_values.py`, while the PDF is compiled from LaTeX by hand.
Where the two disagree about an implementation detail, the chapters are the more
current, because their numbers fail a check when the code moves and the PDF's do
not.
```

## Quick index of the central results

| Level | Result |
|---|---|
| 1 | $J = \dfrac{D K_s}{L}\left(\sqrt{P_{\text{up}}} - \sqrt{P_{\text{down}}}\right)$ |
| 2a | $J = \dfrac{D_{ox}K_{ox}}{L_{ox}}\left(P_{\text{up}} - P_{\text{down}}\right)$ |
| 2b | $\alpha(P_{\text{up}} - P_{\text{int}}) = \beta(\sqrt{P_{\text{int}}} - \sqrt{P_{\text{down}}})$ |
| 3 | $J_{\text{total}} = \sum_i j^{(i)} f^{(i)}$ |
| 4 | $D_{\text{eff}} = \left[(1-f_{gb})D_L + f_{gb}\alpha D_L\right]\big/\left(1 + \sum_i N_{T,i}K_i/N_L\right)$ |
| 5 | Level 3 paths, each solved with Level 4's $D_{\text{eff}}$ |
| 6 | $J_{\text{surf}} = k_{\text{diss}}P(1-\theta)^2 - k_{\text{recomb}}\theta^2$, with $\sqrt{P_{\text{int}}}(\theta)$ closed-form |

Symbols and units are tabulated in {doc}`../PARAMETERS`; the physical reasoning
behind each is in the corresponding theory chapter.
