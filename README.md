# corrStability

Quantifying how much the choice of **zero-imputation strategy** perturbs
metabolite correlations in untargeted metabolomics data.

Metabolomics tables are full of zeros (below detection limit). Before analysis
those zeros must be filled in — but *how* you fill them is a judgement call, and
different choices can quietly change your downstream results. `corrStability`
makes that sensitivity measurable: it generates several imputed variants of each
metabolite and reports how strongly they agree. Correlations near 1 mean the
imputation choice barely matters; lower or more variable correlations flag
metabolites whose signal depends on how their zeros were handled.

## The idea, as pseudocode

```
GOAL: measure how sensitive each metabolite is to the zero-filling method.

for each metabolite column m:

    nonzeros ← values of m that are > 0
    floor    ← min(nonzeros) / 3            # a detection-limit surrogate

    # --- two strategies for replacing the zeros ---
    fixed[m]  ← m, with every zero set to  floor
    random[m] ← m, with every zero set to  floor + jitter,
                where  jitter = sd(m) * Uniform(0, 0.2)

    # --- log-transformed counterparts of each ---
    log_fixed[m]  ← log(fixed[m])
    log_random[m] ← log(random[m])

    # --- how closely do the five variants track each other? ---
    variants ← { original, fixed, random, log_fixed, log_random }
    for each pair (a, b) in variants:
        corr[m, a·b] ← correlation( a[m], b[m] )

# --- roll up across all metabolites ---
for each variant pair (a · b):
    report  mean, median, min, max, sd  of  corr[·, a·b]
```

The result is one correlation row per metabolite (every variant pair) plus a
summary table describing, for each pair, how stable that relationship is across
the whole dataset.

## Installation

```r
# install.packages("devtools")
devtools::install_github("prashantbhandari87/corrStability")
```

## Usage

The one thing you must tell it is **which columns are metabolites** — there is
no reliance on a naming convention. Pick whichever of the three strategies fits
your data:

```r
library(corrStability)

# A) name your metadata columns; the numeric remainder are metabolites
ms <- metabolite_set(data, id_cols = c("sample_id", "group"))

# B) the metabolites share a pattern
ms <- metabolite_set(data, pattern = "^HMDB")

# C) list them explicitly
ms <- metabolite_set(data, metabolites = c("glucose", "lactate", "citrate"))

ms
#> <metabolite_set>
#>   samples:     48
#>   metabolites: 212
#>   columns:     glucose, lactate, citrate, ... (209 more)
```

Always print the object to confirm it grabbed the metabolite count you expect —
if that number is wrong, your column specification is wrong, and this is the
place to catch it.

Then run the full workflow:

```r
res <- run_metabolite_analysis(ms, seed = 123)

res$summary        # mean / median / min / max / sd per correlation type
res$correlations   # one row per metabolite, one column per variant pair
res$transformed    # the metabolite_set with fixed_/random_/log_ columns added

saveRDS(res, "samples_correlations.rds")
```

Or run the steps individually, e.g. to inspect the imputed columns before
correlating:

```r
ms_t  <- transform_metabolites(ms, seed = 123)   # returns a metabolite_set
corrs <- calculate_metabolite_correlations(ms_t) # data.frame
summ  <- summarize_correlations(corrs)           # data.frame

head(ms_t$data)   # peek at the fixed_/random_/log_fixed_/log_random_ columns
```

The `seed` makes the random imputation reproducible without touching your
global RNG. Pass `seed = NULL` if you want fresh randomness on each run.

