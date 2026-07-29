"""
metabolite_set : a small class holding a wide metabolomics table plus the
explicitly-resolved set of metabolite (feature) columns.

Port of met_cor.R / metabolite_set.R. Identification happens ONCE, at
construction, and is validated there. Every downstream function reads
`ms.metabolite_cols` instead of re-guessing.

The R original used S3 (an object + generics dispatching on class). Python's
natural equivalent is a plain class with methods, so the S3 `UseMethod`
plumbing simply disappears.
"""

from __future__ import annotations

import re
from typing import Iterable, Optional, Sequence

import polars as pl
import tidypolars as tp


def _as_tibble(data) -> tp.tibble:
    """Accept a tidypolars tibble, polars/pandas DataFrame, or dict."""
    if isinstance(data, tp.tibble):
        return data
    if isinstance(data, pl.DataFrame):
        return tp.from_polars(data)
    # pandas / dict / anything polars can construct from
    try:
        import pandas as pd  # noqa
        if isinstance(data, pd.DataFrame):
            return tp.from_pandas(data)
    except ImportError:
        pass
    return tp.from_polars(pl.DataFrame(data))


def resolve_metabolite_cols(
    data: tp.tibble,
    metabolites: Optional[Sequence[str]] = None,
    pattern: Optional[str] = None,
    id_cols: Optional[Sequence[str]] = None,
) -> list[str]:
    """
    Resolve which columns are metabolites via exactly one of three strategies:

    - `metabolites`: an explicit list of column names.
    - `pattern`:     a regex matched against the column names (e.g. r"^HMDB").
    - `id_cols`:     names of metadata columns; every *other* numeric column is
                     treated as a metabolite. Usually the most robust choice.
    """
    supplied = [metabolites is not None, pattern is not None, id_cols is not None]
    if sum(supplied) == 0:
        raise ValueError("Specify exactly one of `metabolites`, `pattern`, or `id_cols`.")
    if sum(supplied) > 1:
        raise ValueError("Specify only one of `metabolites`, `pattern`, or `id_cols`.")

    names = list(data.names)
    schema = dict(zip(data.names, data.dtypes))  # tibble blocks .schema

    if metabolites is not None:
        missing = [m for m in metabolites if m not in names]
        if missing:
            raise ValueError("Metabolite columns not found in data: " + ", ".join(missing))
        cols = list(metabolites)
    elif pattern is not None:
        cols = [n for n in names if re.search(pattern, n)]
    else:  # id_cols
        missing = [c for c in id_cols if c not in names]
        if missing:
            raise ValueError("`id_cols` not found in data: " + ", ".join(missing))
        candidate = [n for n in names if n not in set(id_cols)]
        cols = [c for c in candidate if schema[c].is_numeric()]

    if len(cols) == 0:
        raise ValueError("No metabolite columns matched the given specification.")
    return cols


class metabolite_set:
    """
    A wide metabolomics table + the resolved metabolite column names.

    Construct with exactly one of `metabolites`, `pattern`, or `id_cols`.
    `.data` is a tidypolars tibble; `.metabolite_cols` is a list[str].
    """

    def __init__(
        self,
        data,
        metabolites: Optional[Sequence[str]] = None,
        pattern: Optional[str] = None,
        id_cols: Optional[Sequence[str]] = None,
    ):
        tib = _as_tibble(data)
        cols = resolve_metabolite_cols(tib, metabolites, pattern, id_cols)
        self.data: tp.tibble = tib
        self.metabolite_cols: list[str] = cols
        self._validate()

    def _validate(self) -> "metabolite_set":
        names = set(self.data.names)
        missing = [c for c in self.metabolite_cols if c not in names]
        if missing:
            raise ValueError("metabolite_cols refer to absent columns: " + ", ".join(missing))
        schema = dict(zip(self.data.names, self.data.dtypes))
        non_numeric = [c for c in self.metabolite_cols if not schema[c].is_numeric()]
        if non_numeric:
            raise ValueError(
                "Metabolite columns must be numeric; these are not: " + ", ".join(non_numeric)
            )
        return self

    # --- printing: mirrors format.metabolite_set / print.metabolite_set ------
    def __repr__(self) -> str:
        eg = self.metabolite_cols[:3]
        more = len(self.metabolite_cols) - len(eg)
        eg_str = ", ".join(eg)
        if more > 0:
            eg_str += f", ... ({more} more)"
        return (
            "<metabolite_set>\n"
            f"  samples:     {self.data.nrow}\n"
            f"  metabolites: {len(self.metabolite_cols)}\n"
            f"  columns:     {eg_str}"
        )
