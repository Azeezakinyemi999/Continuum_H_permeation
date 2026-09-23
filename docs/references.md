# References

## Foundational theory

| Source | Used for |
|---|---|
| Sieverts (1929) | Square-root pressure dependence of dissolved hydrogen |
| Richardson & Antill (1955) | Permeation through metals |
| Langmuir (1918) | Surface adsorption isotherm |
| Oriani (1970) | Local-equilibrium trapping model — {doc}`theory/level4-microstructure` |
| McNabb & Foster (1963) | Kinetic trapping; the alternative when Oriani's low-occupancy assumption fails |
| Strehlow & Savage (1974) | Parallel-path transport through defective coatings — {doc}`theory/level3-defective-oxide` |
| Hart (1957), Underwood (1970) | Grain-boundary diffusion and stereology |
| Palumbo & Aust (1990) | Grain-boundary structure dependence |

## Material property sources

| Source | Supplies |
|---|---|
| Guo et al. (2025) | 316L diffusivity, solubility, permeability — the active study |
| Fuerst et al. (2024) | Hastelloy N permeation data |
| Schmidt et al. (1985) | Incoloy 802 (X40 NiCrAlTi 31/19) diffusivity and solubility |
| Nemanic et al. (2023) | Cr₂O₃ transport properties |
| Stover (1986) | Cr₂O₃ activation energies |
| Grant et al. (1988) | Surface dissociation rate constants |
| Lu et al. (2022) | Trap binding energies from thermal desorption spectroscopy |
| Zhu et al. (2021) | Grain size and dislocation density from EBSD |
| Young et al. (1997) | M₆C carbide trap density |

## Sensitivity analysis methodology

| Source | Used for |
|---|---|
| Pianosi & Wagener (2018) | PAWN index, and the dummy-parameter significance test |
| Borgonovo (2007) | The δ moment-independent importance measure |
| Puy, Lo Piano & Saltelli (2020) | Documented failure modes of PAWN — see {doc}`caveats` |

## Where the numbers actually live

Material values are not stored in this documentation. They live in the study
configuration modules, each entry carrying a `reference` field naming its source:

- [`calculations/config/studies/`](https://github.com/Azeezakinyemi999/MHI_permeation/tree/main/calculations/config/studies/)

The full bibliography in BibTeX form, 32 entries, is
[`latex/references.bib`](https://github.com/Azeezakinyemi999/MHI_permeation/blob/main/latex/references.bib).
