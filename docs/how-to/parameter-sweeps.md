# How to read the parameter sweeps

Four notebooks in each study folder sweep one parameter at a time through every level
of both models, L1 to L5 and L1L6 to L5L6, and fit the result. This page says what
each one fits, how to read its figures and tables, and in particular the difference
between the two pressure exponents the tables report: the **regime $n$** and the
**net $n$**.

The notebooks are `effect-of-temperature.ipynb`, `effect-of-pressure.ipynb`,
`effect-of-thickness.ipynb` and `effect-of-pinholes.ipynb`, in `Guo_etal_2025/`,
`Fuerst_et_al_2024/` and `incoloy802_cr2o3/`. Every value they use comes from the
study's own `study_config.py`, so the same notebook gives each study's curves.

## What the four notebooks sweep

| Notebook | Swept | Held at | Fit | Exponent in the tables |
|---|---|---|---|---|
| `effect-of-temperature` | $T$ over `T_range` | `P_upstream` | $J = A\,e^{-E_a/RT}$ | $E_a$ (kJ/mol) |
| `effect-of-pressure` | $P_{up}$ over `P_range` | `T_operating` | $J = A\,P^{\,n}$ | $n$ |
| `effect-of-thickness` | $L_{ox}$, $L_m$, or both scaled | `T_operating`, `P_upstream` | $J = A\,L^{\,m}$ | $m$ |
| `effect-of-pinholes` | pinhole area fraction $f$ over `defect_fraction_range` | `T_operating`, `P_upstream` | $J = A\,f^{\,n}$ | $n$ |

Each sweep is done at the study's operating point and again at the low, middle and
high end of the other ranges: three pressures for the temperature sweep, three
temperatures for the pressure sweep, both sets for the thickness and pinhole sweeps.

The apparent permeability $\Phi_{app} = J\,L/\Delta\sqrt{P}$ ($J\,L/\Delta P$ for L2a,
whose oxide follows Henry's law) is computed from the same flux and fitted the same
way. It is a material property only where thickness-proportional Sieverts transport
controls; see {doc}`../caveats`.

## Regimes, and why each fit covers one regime

Every simulated point is labelled with its **regime**: the path that carries most of
the flux (intact oxide, pinholes, or the bare wall) and the step that limits it
(surface dissociation, oxide diffusion, metal diffusion, lattice or trap-controlled).
The labels come from `calculations.wall_labels`. A point where no single step limits
is labelled "no single limiting step (not fitted)": it is drawn, never fitted.

A power law or an Arrhenius law describes one controlling step. Across a change of
regime the slope changes, so one line through every point describes none of them.
The fits are therefore made **per regime**, over the window where that regime holds,
and only when the window has at least three points. These are the lines the figures
draw, coloured like the points they fit.

## The pressure exponent: regime n and net n

Write the **local exponent** of the flux at a pressure as

$$n(P) = \frac{d\ln J}{d\ln P}.$$

For a clean power law it is constant: $n \approx 1$ where a Henry-law oxide (Model 1)
or the surface rate (Model 2) controls, $n = 0.5$ where Sieverts transport controls.
Across the full `P_range` it changes wherever the regime changes. The pressure tables
give two summaries of it.

### The regime n: the slope inside one regime

For each regime, $\ln J = \ln A + n\,\ln P$ is fitted by least squares over that
regime's pressure window. Its slope $n$ names the controlling step, which is why it is
the slope to compare with the reference lines of the figure, $n = 1$ and $n = 0.5$.
Each regime fit is drawn, and it is a row of the table with its prefactor $A$.

### The net n: the trend over the whole range

The table's row **net (whole range)** gives, for each temperature,

$$n_{net} = \frac{\ln\!\big(J(P_{max})/J(P_{min})\big)}{\ln\!\big(P_{max}/P_{min}\big)}
          = \frac{1}{\ln(P_{max}/P_{min})}\int_{P_{min}}^{P_{max}} n(P)\,d\ln P ,$$

the net change of $\ln J$ over the net change of $\ln P$ between the two ends of the
range. The second form is the same number written as an average: $n_{net}$ is the
**local exponent averaged over the range, every decade of pressure counting equally**.
It answers "over the whole range studied, by how much does the flux rise per decade of
pressure?"

It is not a fitted line. It uses only the two ends, assumes nothing about the shape in
between, and has no prefactor: the table leaves $A$ empty on that row, and the figure
does not draw it.

### Why not one least-squares line over the whole range

A least-squares line over every point has a slope too, but it is a different average.
On a grid uniform in $\ln P$, its slope weights the local exponent by a parabola,
largest in the middle of the range and zero at the two ends, so it under-counts the
regimes at either end; and the line itself need not pass near any point. Where one
regime holds over the whole range all three numbers agree. Where the regime changes,
the net $n$ is the one with a plain meaning, which is why the tables report it instead.

### Reading the two together

- **The net $n$ equals one regime $n$**: that regime holds over (nearly) the whole range.
- **The net $n$ lies between the regime values**: the range is split between regimes,
  and $n_{net}$ is close to their exponents weighted by the decades of pressure each one
  spans. The short transition stretches between windows count too, which is why the
  match is close rather than exact.
- **The net $n$ changes with temperature**: the regimes move. Compare the rows of the
  three-temperature table (`<level>_pressure_3T_fits.csv`).

Quote the regime $n$ to name a mechanism, and the net $n$ to state the overall pressure
dependence over the range, for instance when comparing with a measurement that reports
one exponent across a wide pressure span. Neither is a material constant: both belong to
the operating point and the range.

The same two summaries are given for $\Phi_{app}$. Its exponent is $n - 0.5$ (or $n - 1$
for L2a), so a flat $\Phi_{app}$, slope 0, is the Sieverts regime.

## The temperature sweep: regime fits and the single fit

The figures draw one Arrhenius fit per regime, on $\ln J$ against $1/T$; the slope of each
gives that regime's $E_a$. The table also has a **single fit**: one least-squares
Arrhenius line over every temperature. It is reported for comparison and not drawn. Where
the regime changes with temperature, which it does for most defective and surface levels,
it averages over those regimes; quote the regime values there. $\Phi_{app}$ is taken at a
fixed pressure, so its activation energy is that of $J$.

## The thickness sweep: one figure per swept layer

Three sweeps, each its own figure, flux on the left and $\Phi_{app}$ on the right:

1. **oxide thickness** over `L_oxide_range`, the metal fixed at `L_metal`;
2. **metal thickness** over `L_metal_range`, the oxide fixed at `L_oxide`;
3. **the whole wall scaled**: both layers multiplied by one factor, the ratio fixed.

A level gets only the sweeps of the layers it has. A one-layer level, L1, L4 and L1L6
(metal) or L2a and L2aL6 (oxide), gets just its own layer, because scaling its whole wall
would repeat that sweep exactly.

The fit is $J \propto L^{\,m}$ per regime. $m = -1$ when the swept thickness controls
the flux, $m \approx 0$ when something else does: the other layer, or the surface. The
whole-wall sweep separates bulk from surface: every diffusion step scales as $1/L$, so a
bulk-controlled wall gives exactly $m = -1$ and a flat $\Phi_{app}$, while surface kinetics,
which do not depend on thickness, pull $m$ toward 0 and make $\Phi_{app}$ grow with $L$.

## The pinhole sweep: PRF and the half-flux fraction

With cracks and oxide grain boundaries set to zero, the wall is two paths side by side,
intact oxide and bare metal under the pinholes, so $J = (1-f)\,J_i + f\,J_p$. One number
sets the curve, the permeation reduction factor $\mathrm{PRF} = J_p/J_i$, and pinholes
carry half the flux at

$$f_{1/2} = \frac{1}{1 + \mathrm{PRF}} .$$

A PRF below 1 means the pinhole path passes less than the intact oxide, as when
dissociation on the bare metal is slow; the flux then falls with $f$. The fits are
$J \propto f^{\,n}$ per regime; the paths table, `<level>_pinholes_paths.csv`, gives every
point's flux through each path, its share, its limiting step, the PRF and $f_{1/2}$.

## Reading a figure

- **Panels**: the flux above or left, $\Phi_{app}$ below or right, on log axes.
- **The legend is outside the panels, on the right.** *Simulation*, level with the panels:
  one marker per regime. Colour is the limiting step; a hollow marker means a defect path
  carries most of the flux; a diamond means the metal is trap-controlled. *Lines*, below it
  at the bottom: the fit line, the reference slopes, and in a three-condition figure the line
  style of each condition.
- **Reference slopes** (pressure figures): dashed lines through the middle point for
  $n = 1$ (long dashes) and $n = 0.5$ (short dashes), each labelled on the panel. Compare
  a regime fit's slope with them.
- **Three-condition figures**: one curve per condition, named by its line style and by
  the label at its right-hand end.

## Reading a table

Every figure has a fits table beside it, `<figure>_fits.csv`, also shown in the notebook.

| Column | Meaning |
|---|---|
| `quantity` | `J`, the flux, or `Phi_app`, the apparent permeability |
| `condition` | where the sweep was made: the operating point, or one of the three conditions |
| `sweep` | thickness only: `oxide`, `metal` or `total` (the whole wall) |
| `fit` | the regime fitted, by its dominant path and limiting step; or `single fit` (temperature), or `net (whole range)` (pressure) |
| window, `points` | the range the row covers and how many simulated points it contains |
| exponent | $E_a$, $n$ or $m$ for the flux; for $\Phi_{app}$, its own slope |
| `A` | the prefactor of the fit, SI units; empty on a `net (whole range)` row |

The tables give the fitted values only; they carry no goodness-of-fit column.

## Files

Each notebook writes into the folder named after it, `<notebook>-results/`:

| Notebook | Figures | Tables |
|---|---|---|
| temperature | `<level>_arrhenius`, `<level>_arrhenius_3P` | `<figure>_fits.csv` |
| pressure | `<level>_pressure`, `<level>_pressure_3T` | `<figure>_fits.csv` |
| thickness | `<level>_thickness_<sweep>`, `<level>_thickness_<sweep>_3P`, `_3T`, with `<sweep>` = `oxide`, `metal`, `total` | `<level>_thickness_fits.csv`, `_3P_fits.csv`, `_3T_fits.csv` |
| pinholes | `<level>_pinholes`, `<level>_pinholes_3P`, `_3T` | `<figure>_fits.csv`, `<figure>_paths.csv` |

Figures are saved as PNG and PDF. The study report, `python -m calculations.report
<study>`, has one section per notebook built from these files.

## Running the sweeps

Open a notebook in its study folder and run it top to bottom; the four are independent
of each other and of the other notebooks. To run all of them headless and then build the
study report:

```bash
python -m calculations.report Guo_etal_2025 --run-notebooks effect-of-temperature,effect-of-pressure,effect-of-thickness,effect-of-pinholes
```

On a laptop a sweep notebook takes about 5 to 13 minutes, the thickness sweep the longest.
