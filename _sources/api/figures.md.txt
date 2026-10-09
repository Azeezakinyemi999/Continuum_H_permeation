# `figures` — where a notebook's figures go

Each notebook writes into one folder beside it, named after it, created relative to
the working directory — the same rule the study lookup and `RESULTS_DIR` follow:

```python
from calculations.figures import results_dir, save_figure

FIGS = results_dir("Proposal")            # -> Proposal-results/
save_figure(fig, "Level1_Perfect_Metal", FIGS)
```

The sensitivity plotters in {doc}`sensitivity` all return their figure, so they wrap
directly — including the ones that return `None` for a plot they declined to draw,
which `save_figure` treats as a no-op rather than an error.

```{eval-rst}
.. automodule:: calculations.figures
   :members:
   :member-order: bysource
```
