"""
Analysis functions for metabolite_set. Port of analysis.R.

The numeric core is where tidypolars earns its place:
  - transform_metabolites   -> mutate() with if_else / min-of-nonzero / log
  - summarize_correlations  -> pivot_longer() + group_by() + summarise()

Two deliberate departures from the R original:
  1. RNG. R saved and restored `.Random.seed` to avoid clobbering the global
     stream. Here we use a *local* numpy Generator, so there is no global state
     to protect and the save/restore dance disappears entirely.
  2. Random values will NOT match R's runif() stream byte-for-byte (different
     RNG). The design intent is preserved (reproducible under a seed; only zero
     rows imputed; non-zeros untouched), and every deterministic result
     -- fixed_, the logs, and any zero-free correlation -- matches R exactly.
"""

from __future__ import annotations

import warnings
from typing import Optional

import numpy as np
import polars as pl
import tidypolars as tp

from .metabolite_set import metabolite_set


# --- derived column names for one metabolite (single source of truth) --------
def derived_names(col: str) -> dict[str, str]:
    return {
        "fixed": f"fixed_{col}",
        "random": f"random_{col}",
        "log_fixed": f"log_fixed_{col}",
        "log_random": f"log_random_{col}",
    }


# =============================================================================
# transform_metabolites
# =============================================================================
def transform_metabolites(x: metabolite_set, seed: Optional[int] = 123) -> metabolite_set:
    """
    For each metabolite column, append four derived columns:
      fixed_       zeros -> min(nonzero)/3
      random_      zeros -> that value + SD-scaled uniform jitter
      log_fixed_   natural log of fixed_
      log_random_  natural log of random_

    All-zero columns are skipped (no derived columns), matching the R behaviour.
    Returns a new metabolite_set; metabolite_cols is unchanged.
    """
    rng = np.random.default_rng(seed)  # local stream; global RNG untouched
    data = x.data

    # eligible = has at least one non-zero value
    eligible = []
    for col in x.metabolite_cols:
        vals = data.pull(col)
        if (vals > 0).sum() > 0:
            eligible.append(col)

    if not eligible:
        return metabolite_set(data, metabolites=x.metabolite_cols)

    # --- fixed_ : one mutate, all eligible columns at once -------------------
    fixed_exprs = {}
    for col in eligible:
        nm = derived_names(col)
        floor = tp.col(col).filter(tp.col(col) > 0).min() / 3
        fixed_exprs[nm["fixed"]] = tp.if_else(tp.col(col) == 0, floor, tp.col(col))
    data = data.mutate(**fixed_exprs)

    # --- log_fixed_ : references the freshly-made fixed_ columns --------------
    logf_exprs = {
        derived_names(c)["log_fixed"]: tp.col(derived_names(c)["fixed"]).log()
        for c in eligible
    }
    data = data.mutate(**logf_exprs)

    # --- random_ : jitter drawn from the local numpy Generator ---------------
    # Built as concrete Series (only zero rows imputed) then attached in one go.
    random_series = {}
    for col in eligible:
        nm = derived_names(col)
        vals = data.pull(col).to_numpy().astype(float)
        floor = vals[vals > 0].min() / 3.0
        sd = np.std(vals[~np.isnan(vals)], ddof=1)  # R sd(): sample SD, na.rm
        zero_mask = vals == 0
        out = vals.copy()
        jitter = rng.uniform(0.0, 0.2, size=int(zero_mask.sum()))
        out[zero_mask] = floor + sd * jitter
        random_series[nm["random"]] = pl.Series(nm["random"], out)
    data = data.mutate(**{name: tp.lit(s) for name, s in random_series.items()})

    # --- log_random_ ---------------------------------------------------------
    logr_exprs = {
        derived_names(c)["log_random"]: tp.col(derived_names(c)["random"]).log()
        for c in eligible
    }
    data = data.mutate(**logr_exprs)

    return metabolite_set(data, metabolites=x.metabolite_cols)


# =============================================================================
# calculate_metabolite_correlations
# =============================================================================
_PAIRS = [
    ("original_vs_fixed", "original", "fixed"),
    ("original_vs_random", "original", "random"),
    ("original_vs_log_fixed", "original", "log_fixed"),
    ("original_vs_log_random", "original", "log_random"),
    ("fixed_vs_random", "fixed", "random"),
    ("fixed_vs_log_fixed", "fixed", "log_fixed"),
    ("fixed_vs_log_random", "fixed", "log_random"),
    ("random_vs_log_fixed", "random", "log_fixed"),
    ("random_vs_log_random", "random", "log_random"),
    ("log_fixed_vs_log_random", "log_fixed", "log_random"),
]

CORR_TYPES = [name for name, _, _ in _PAIRS]


def calculate_metabolite_correlations(x: metabolite_set) -> Optional[tp.tibble]:
    """
    Pairwise Pearson correlations among the five variants of each metabolite.
    One row per metabolite, one column per variant pair. Returns None (with a
    warning) if the set has not been transformed yet.
    """
    data = x.data
    names = set(data.names)
    pdf = data.as_polars()

    rows = []
    for col in x.metabolite_cols:
        nm = derived_names(col)
        if not set(nm.values()).issubset(names):
            warnings.warn(
                f"Skipping '{col}': derived columns missing. "
                "Run transform_metabolites() first.",
                stacklevel=2,
            )
            return None

        variant_col = {
            "original": col,
            "fixed": nm["fixed"],
            "random": nm["random"],
            "log_fixed": nm["log_fixed"],
            "log_random": nm["log_random"],
        }

        def cc(a: str, b: str) -> float:
            ca, cb = variant_col[a], variant_col[b]
            d = pdf.select([ca, cb]).drop_nulls()  # use = "complete.obs"
            return d.select(pl.corr(pl.col(ca), pl.col(cb))).item()

        row = {"metabolite": col}
        for name, a, b in _PAIRS:
            row[name] = cc(a, b)
        rows.append(row)

    return tp.from_polars(pl.DataFrame(rows))


# =============================================================================
# summarize_correlations
# =============================================================================
def summarize_correlations(correlation_results: tp.tibble) -> tp.tibble:
    """
    mean / median / min / max / sd per correlation type, computed by reshaping
    long and grouping -- the idiomatic tidyverse move.
    """
    long = correlation_results.pivot_longer(
        cols=CORR_TYPES,
        names_to="correlation_type",
        values_to="value",
    )
    summ = long.summarize(
        mean=tp.col("value").mean(),
        median=tp.col("value").median(),
        min=tp.col("value").min(),
        max=tp.col("value").max(),
        sd=tp.col("value").std(),
        _by="correlation_type",
    )
    # deterministic ordering to match the pair definition order
    order = {name: i for i, name in enumerate(CORR_TYPES)}
    ordered = summ.as_polars().sort(
        by=pl.col("correlation_type").replace_strict(order, return_dtype=pl.Int64)
    )
    return tp.from_polars(ordered)


# =============================================================================
# run_metabolite_analysis  (full workflow)
# =============================================================================
def run_metabolite_analysis(x: metabolite_set, seed: Optional[int] = 123) -> dict:
    """Transform, correlate, summarise. Returns a dict of the three results."""
    transformed = transform_metabolites(x, seed=seed)
    correlations = calculate_metabolite_correlations(transformed)
    return {
        "transformed": transformed,
        "correlations": correlations,
        "summary": summarize_correlations(correlations),
    }
