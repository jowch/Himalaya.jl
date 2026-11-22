# HimalayaDB.jl

Read-only, programmatic access to HimalayaUI's curated SAXS annotations.

```julia
using HimalayaDB
db = connect()                       # HIMALAYA_DB_PATH, or ~/.himalaya/himalaya.db

exp = experiments(db)[1]             # one per beamtime
sts = setups(db; experiment=exp.id)  # one per setup folder: dirs, patterns, energy_kev,
                                     # flight_path_m, beam centre, pixel_size_um, q_units
smp = samples(db; experiment=exp.id)[1]
exps = exposures(db; sample=smp.id)  # each row has setup_id; join it to `sts` for geometry

call = phase_call(db, smp.id)        # the sample's phase call, as the Call list shows it:
                                     # (; state, phase, call_phases, members)
sample_labels(db, smp.id)            # confirmed and suggested labels (Buffer, Calibrant, …)
subtraction_background(db, smp.id)   # the background sample's id, or nothing
bridged_samples(db, smp.id)          # samples a user marked as the same specimen

peaks = curated_peaks(db, exps[1].id)   # auto ∪ adds − excludes; excluded peaks
                                         # are returned tagged (excluded ∈ {0,1}),
                                         # not removed — the truly-effective set
                                         # is filter(p -> p.excluded == 0, peaks)
cands = index_candidates(db, exps[1].id)  # every candidate index
confirmed = confirmed_indices(db, exps[1].id)  # this exposure's members of the call

using DataFrames
dataframe(peaks)                     # -> DataFrame

idx = reconstruct_index(db, cands[1].id)  # -> Himalaya.Index{Pn3m}
q, I, σ = load_trace(db, exps[1].id)      # opt-in .dat loading, from the exposure's setup
```

## Design

HimalayaDB reads the HimalayaUI database file directly, with no server running. These rules keep it safe to use next to the app:

- **It only reads.** It never writes, migrates or chmods the file. `connect` opens with `PRAGMA query_only = ON`, so any write or DDL fails loudly. SQLite.jl has no read-only open flag, so the pragma is the guard.
- **It does not own the schema.** HimalayaUI creates and migrates every table. HimalayaDB only `SELECT`s, and on a DB that predates a migration it warns at `connect` rather than healing anything.
- **It duplicates as little as possible.** The only copied logic is the effective-peaks query (auto peaks minus excludes, plus adds) and the confirmed-index query. `test/test_contract.jl` runs both packages' readers on the same fixture DB and fails if they disagree, so schema drift in HimalayaUI shows up here.
- **"Confirmed" means the `assignments` model.** `confirmed_indices` returns the exposure's indices when its `assignments.state` is `indexed` (the default when no row exists). On upgraded DBs that can be auto-seeded, so it is not proof of a human decision.
- **Reading while the app writes.** The DB runs in WAL mode. A reader opened while HimalayaUI is writing needs write access to the DB's directory, because SQLite creates the `-wal` and `-shm` sidecars there. Reading a copied `.db` with the app stopped needs no such access.
- **Traces are opt-in.** Annotations come from the DB alone. `load_trace` reads the `.dat` file from disk and errors clearly if it isn't reachable.

HimalayaUI calls an exposure a "frame". The database and this package keep
"exposure".

## Upgrading from 0.1

0.2 follows the one-workspace-per-beamtime schema. Directories, file patterns
and geometry moved from the experiment to its setups, so `experiments` no
longer returns `data_dir`, `analysis_dir`, `experiment_type`, `energy_kev` or
`flight_path_m`. Read them from `setups` and join on `setup_id`: two setups of
one beamtime can have different flight paths. 0.1 fails on any database
HimalayaUI has opened since the setups migration.

## Upgrading to 0.4

Samples bridged as the same sample share one call, so `phase_call` of a bridged
sample now includes the indices confirmed on the other members' exposures. A
member's `exposure_id` says which sample's frame holds the evidence.

## Upgrading to 0.5

HimalayaUI dropped `samples.merged_into_id`, a column only the retired merge
route ever wrote. `samples` and the annotation readers no longer filter on it.
0.4 fails on any database HimalayaUI has opened since then.

## Keeping it in step with HimalayaUI

The package re-implements a few HimalayaUI reads. `test/test_contract.jl` runs
both on one fixture database, and CI runs this suite, so a schema change in
HimalayaUI that breaks a reader fails the build.
