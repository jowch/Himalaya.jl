# HimalayaDB.jl

Read-only, programmatic access to HimalayaUI's curated SAXS annotations.

```julia
using HimalayaDB
db = connect()                       # HIMALAYA_DB_PATH, or ~/.himalaya/himalaya.db
exps = exposures(db; sample=1)
peaks = curated_peaks(db, exps[1].id)   # auto ∪ adds − excludes; excluded peaks
                                         # are returned tagged (excluded ∈ {0,1}),
                                         # not removed — the truly-effective set
                                         # is filter(p -> p.excluded == 0, peaks)
cands = index_candidates(db, exps[1].id)  # every candidate index; populated whenever indices exist

using DataFrames
dataframe(peaks)                     # -> DataFrame

idx = reconstruct_index(db, cands[1].id)  # -> Himalaya.Index{Pn3m}
q, I, σ = load_trace(db, exps[1].id)      # opt-in .dat loading

# The exposure's durable indexed assignment (the "confirmed index" HimalayaUI uses).
# Empty when the assignment state is form_factor/null. Note: on upgraded DBs this can
# be auto-seeded, so a non-empty result is not by itself proof of a human decision.
confirmed = confirmed_indices(db, exps[1].id)
```
