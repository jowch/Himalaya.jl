# HimalayaDB.jl

Read-only, programmatic access to HimalayaUI's curated SAXS annotations.

```julia
using HimalayaDB
db = connect()                       # HIMALAYA_DB_PATH, or ~/.himalaya/himalaya.db
exps = exposures(db; sample=1)
peaks = curated_peaks(db, exps[1].id)   # auto ∪ adds − excludes
idxs  = confirmed_indices(db, exps[1].id)

using DataFrames
dataframe(peaks)                     # -> DataFrame

using SparseArrays                   # opt-in typed reconstruction
idx = reconstruct_index(db, idxs[1].id)   # -> Himalaya.Index{Pn3m}
q, I, σ = load_trace(db, exps[1].id)      # opt-in .dat loading
```
