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

## Design

HimalayaDB reads the HimalayaUI database file directly, with no server running. These rules keep it safe to use next to the app:

- **It only reads.** It never writes, migrates or chmods the file. `connect` opens with `PRAGMA query_only = ON`, so any write or DDL fails loudly. SQLite.jl has no read-only open flag, so the pragma is the guard.
- **It does not own the schema.** HimalayaUI creates and migrates every table. HimalayaDB only `SELECT`s, and on a DB that predates a migration it warns at `connect` rather than healing anything.
- **It duplicates as little as possible.** The only copied logic is the effective-peaks query (auto peaks minus excludes, plus adds) and the confirmed-index query. `test/test_contract.jl` runs both packages' readers on the same fixture DB and fails if they disagree, so schema drift in HimalayaUI shows up here.
- **"Confirmed" means the `assignments` model.** `confirmed_indices` returns the exposure's indices when its `assignments.state` is `indexed` (the default when no row exists). On upgraded DBs that can be auto-seeded, so it is not proof of a human decision.
- **Reading while the app writes.** The DB runs in WAL mode. A reader opened while HimalayaUI is writing needs write access to the DB's directory, because SQLite creates the `-wal` and `-shm` sidecars there. Reading a copied `.db` with the app stopped needs no such access.
- **Traces are opt-in.** Annotations come from the DB alone. `load_trace` reads the `.dat` file from disk and errors clearly if it isn't reachable.
