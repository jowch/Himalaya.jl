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

HimalayaUI calls an exposure a "frame". The database and this package keep
"exposure".

## Upgrading from 0.1

0.2 follows the one-workspace-per-beamtime schema. Directories, file patterns
and geometry moved from the experiment to its setups, so `experiments` no
longer returns `data_dir`, `analysis_dir`, `experiment_type`, `energy_kev` or
`flight_path_m`. Read them from `setups` and join on `setup_id`: two setups of
one beamtime can have different flight paths. 0.1 fails on any database
HimalayaUI has opened since the setups migration.

## Keeping it in step with HimalayaUI

The package re-implements a few HimalayaUI reads. `test/test_contract.jl` runs
both on one fixture database, and CI runs this suite, so a schema change in
HimalayaUI that breaks a reader fails the build.
