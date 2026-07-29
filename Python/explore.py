"""
Reproduces the R testthat fixture (helper-fixtures.R) and checks the port
against the hand-computed numbers, then prints the workflow output so you can
poke at it interactively.

    glucose: min nonzero = 3  -> zero imputed to 3/3 = 1
    lactate: min nonzero = 2  -> zero imputed to 2/3
    citrate: no zeros         -> fixed_ == original, so cor(original, fixed) = 1
"""
import math
import polars as pl
from corrstability import (
    metabolite_set,
    transform_metabolites,
    calculate_metabolite_correlations,
    summarize_correlations,
    run_metabolite_analysis,
)


def make_metab_df():
    return pl.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3", "s4"],
            "group": ["ctrl", "case", "ctrl", "case"],
            "glucose": [0.0, 3, 6, 9],
            "lactate": [2.0, 0, 4, 8],
            "citrate": [1.0, 2, 3, 4],
        }
    )


def make_metab_set():
    return metabolite_set(make_metab_df(), id_cols=["sample_id", "group"])


def approx(a, b, tol=1e-9):
    return abs(a - b) < tol


# ---- construction -----------------------------------------------------------
ms = make_metab_set()
assert set(ms.metabolite_cols) == {"glucose", "lactate", "citrate"}
print(ms)
print()

# id_cols drops non-numeric non-id columns
df2 = make_metab_df().with_columns(pl.Series("note", ["a", "b", "c", "d"]))
ms2 = metabolite_set(df2, id_cols=["sample_id", "group"])
assert "note" not in ms2.metabolite_cols

# pattern strategy
dfp = make_metab_df().rename({"glucose": "HMDB0001", "lactate": "HMDB0002"})
msp = metabolite_set(dfp, pattern=r"^HMDB")
assert set(msp.metabolite_cols) == {"HMDB0001", "HMDB0002"}

# error paths
for bad in [
    lambda: metabolite_set(make_metab_df()),                                   # none
    lambda: metabolite_set(make_metab_df(), pattern="^gl", id_cols="sample_id"),# two
    lambda: metabolite_set(make_metab_df(), metabolites=["glucose", "nope"]),  # missing
    lambda: metabolite_set(make_metab_df(), pattern="ZZZ_no_match"),           # no match
]:
    try:
        bad(); raise AssertionError("expected an error")
    except ValueError:
        pass

# ---- transform: deterministic fixed_ ---------------------------------------
ms_t = transform_metabolites(make_metab_set(), seed=123)
assert ms_t.metabolite_cols == ms.metabolite_cols
fg = ms_t.data.pull("fixed_glucose").to_list()
fl = ms_t.data.pull("fixed_lactate").to_list()
fc = ms_t.data.pull("fixed_citrate").to_list()
assert fg == [1, 3, 6, 9],           f"fixed_glucose={fg}"
assert all(approx(a, b) for a, b in zip(fl, [2, 2/3, 4, 8])), f"fixed_lactate={fl}"
assert fc == ms_t.data.pull("citrate").to_list()
# log_fixed_ is the log of fixed_
lfg = ms_t.data.pull("log_fixed_glucose").to_list()
assert all(approx(a, math.log(b)) for a, b in zip(lfg, fg))

# random_ only touches the zero rows; non-zeros untouched
rg = ms_t.data.pull("random_glucose").to_list()
assert rg[1:] == [3, 6, 9], f"random_glucose nonzeros changed: {rg}"

# reproducible under the same seed; different across seeds
a = transform_metabolites(make_metab_set(), seed=1)
b = transform_metabolites(make_metab_set(), seed=1)
c = transform_metabolites(make_metab_set(), seed=2)
assert a.data.pull("random_glucose").to_list() == b.data.pull("random_glucose").to_list()
assert a.data.pull("fixed_glucose").to_list() == c.data.pull("fixed_glucose").to_list()
assert a.data.pull("random_glucose").to_list() != c.data.pull("random_glucose").to_list()

# all-zero column is skipped
dfz = make_metab_df().with_columns(pl.Series("empty", [0.0, 0, 0, 0]))
msz = metabolite_set(dfz, id_cols=["sample_id", "group"])
mszt = transform_metabolites(msz, seed=123)
assert "empty" in mszt.metabolite_cols
assert "fixed_empty" not in mszt.data.names

# ---- correlations -----------------------------------------------------------
corrs = calculate_metabolite_correlations(ms_t)
assert corrs.nrow == 3
assert set(corrs.pull("metabolite").to_list()) == {"glucose", "lactate", "citrate"}
citrate_row = corrs.as_polars().filter(pl.col("metabolite") == "citrate")
assert approx(citrate_row["original_vs_fixed"].item(), 1.0), "zero-free -> corr 1"

# untransformed set warns and returns None
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    assert calculate_metabolite_correlations(make_metab_set()) is None

# ---- summary ----------------------------------------------------------------
summ = summarize_correlations(corrs)
assert summ.nrow == 10
assert set(summ.names) == {"correlation_type", "mean", "median", "min", "max", "sd"}
mins = summ.pull("min").to_list()
maxs = summ.pull("max").to_list()
assert all(mn <= mx for mn, mx in zip(mins, maxs))

# ---- full workflow == steps by hand ----------------------------------------
res = run_metabolite_analysis(make_metab_set(), seed=123)
assert set(res.keys()) == {"transformed", "correlations", "summary"}
manual_c = calculate_metabolite_correlations(transform_metabolites(make_metab_set(), seed=123))
assert res["correlations"].as_polars().equals(manual_c.as_polars())

print("ALL CHECKS PASSED\n")
print("=== summary ===")
print(res["summary"].as_polars())
print("\n=== correlations ===")
print(res["correlations"].as_polars())
