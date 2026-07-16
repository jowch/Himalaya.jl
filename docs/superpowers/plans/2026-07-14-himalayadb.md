# HimalayaDB.jl Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone, read-only Julia package that reads HimalayaUI's curated SAXS annotations (effective peak set, index candidates, confirmed phase) directly from the SQLite file, offline, returning Tables.jl rows for downstream analysis and plotting.

**Architecture:** New package at `packages/HimalayaDB/`, sibling to `HimalayaUI`. It opens the central SQLite DB read-only (never migrates, never writes), and re-implements only the small load-bearing query logic (the effective-peaks UNION + the confirmed-group join). A contract test cross-checks its results against HimalayaUI on a shared fixture DB, so the duplicated SQL can't silently drift. Core `Himalaya` is a dependency only for opt-in `Index{P}` reconstruction; `DataFrames` is a weakdep extension.

**Tech Stack:** Julia, SQLite.jl 1.8.0, DBInterface, Tables, SparseArrays (stdlib), core `Himalaya`; `DataFrames` (weakdep); `HimalayaUI` + `Test` (test-only).

## Global Constraints

- **Read-only, always.** No `CREATE`/`INSERT`/`UPDATE`/`DELETE`/migrations/`chmod` in `src/`. Enforced at runtime by `PRAGMA query_only = ON`. (Direct INSERTs are allowed in `test/` fixture code only.)
- **No HTTP / Oxygen dependency.** Never `import Oxygen`.
- **Does not own the schema.** `SELECT` only; never `CREATE TABLE`. The schema is created by HimalayaUI.
- **DB path resolution:** `ENV["HIMALAYA_DB_PATH"]`, else `~/.himalaya/himalaya.db`. Never create directories.
- **Julia ≥ 1.9** (package extensions / weakdeps).
- **Return type:** Tables.jl rows (`Tables.rowtable(...)` → `Vector{<:NamedTuple}`), SQL `NULL` → `missing`, matching HimalayaUI's convention.
- **v1 read surface only:** experiments/samples/exposures + `curated_peaks`, `index_candidates`, `confirmed_indices`, `reconstruct_index`, `load_trace`. No series/comparisons/tags/user_actions.
- **Test framework:** stdlib `Test`.

---

## File Structure

```
packages/HimalayaDB/
  Project.toml            # deps + [weakdeps]/[extensions]; [sources] for Himalaya path
  src/HimalayaDB.jl       # module: includes, exports, `dataframe` stub
  src/connect.jl          # connect(), default_himalaya_db_path()
  src/queries.jl          # experiments/samples/exposures, curated_peaks, index_candidates, confirmed_indices
  src/reconstruct.jl      # resolve_phase(), reconstruct_index()
  src/trace.jl            # load_dat(), load_trace()
  ext/DataFramesExt.jl    # dataframe(rows) = DataFrame(rows)
  test/
    Project.toml          # test env: HimalayaUI, Himalaya, Test, DataFrames (+ [sources])
    runtests.jl           # includes the test files
    fixture.jl            # build_fixture(path, analysis_dir) -> ids (uses HimalayaUI.open_db)
    test_connect.jl
    test_queries.jl
    test_reconstruct.jl
    test_trace.jl
    test_dataframes.jl
    test_contract.jl      # cross-check curated_peaks/index_candidates vs HimalayaUI
```

---

### Task 1: Package skeleton, `connect`, and the shared test fixture

**Files:**
- Create: `packages/HimalayaDB/Project.toml`
- Create: `packages/HimalayaDB/src/HimalayaDB.jl`
- Create: `packages/HimalayaDB/src/connect.jl`
- Create: `packages/HimalayaDB/test/Project.toml`
- Create: `packages/HimalayaDB/test/runtests.jl`
- Create: `packages/HimalayaDB/test/fixture.jl`
- Test: `packages/HimalayaDB/test/test_connect.jl`

**Interfaces:**
- Produces: `HimalayaDB.connect(path=default_himalaya_db_path()) -> SQLite.DB` (read-only); `HimalayaDB.default_himalaya_db_path() -> String`; test helper `build_fixture(path, analysis_dir) -> NamedTuple` with fields `experiment_id, sample_id, exposure_id, exposure2_id, index_id, group_id, auto_peak_ids::Vector{Int}` (`exposure2_id` is an uncurated exposure: active auto group, no custom group).

- [ ] **Step 1: Mint a UUID and read the two path-dep UUIDs**

Run:
```bash
julia -e 'using UUIDs; println(uuid4())'
grep -m1 '^uuid' /home/jonathanchen/projects/Himalaya.jl/Project.toml
grep -m1 '^uuid' /home/jonathanchen/projects/Himalaya.jl/packages/HimalayaUI/Project.toml
```
Record the generated UUID (for `HimalayaDB`) and the `Himalaya` / `HimalayaUI` UUIDs. Also check whether `packages/HimalayaUI/Project.toml` uses a `[sources]` block or a dev'd Manifest to reference `Himalaya`; mirror whichever it uses. The blocks below assume `[sources]` (Julia ≥1.11); if the repo predates it, `Pkg.develop(path=...)` the local packages instead.

- [ ] **Step 2: Write `Project.toml`**

`packages/HimalayaDB/Project.toml` (replace `<HIMALAYADB_UUID>` and `<HIMALAYA_UUID>` with the values from Step 1):
```toml
name = "HimalayaDB"
uuid = "<HIMALAYADB_UUID>"
authors = ["Jonathan Chen"]
version = "0.1.0"

[deps]
SQLite = "0aa819cd-b072-5ff4-a722-6bc24af294d9"
DBInterface = "a10d1c49-ce27-4219-8d33-6db1a4562965"
Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Himalaya = "<HIMALAYA_UUID>"

[weakdeps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"

[extensions]
DataFramesExt = "DataFrames"

[sources]
Himalaya = {path = "../.."}

[compat]
SQLite = "1.8"
DBInterface = "2"
Tables = "1"
DataFrames = "1"
julia = "1.9"

[extras]

[targets]
```

- [ ] **Step 3: Write the module file**

`packages/HimalayaDB/src/HimalayaDB.jl`:
```julia
module HimalayaDB

using SQLite, DBInterface, Tables, SparseArrays
import Himalaya

export connect,
    experiments, samples, exposures,
    curated_peaks, index_candidates, confirmed_indices,
    reconstruct_index, load_trace, dataframe

"""
    dataframe(rows) -> DataFrame

Convert Tables.jl rows to a DataFrame. Method provided by the DataFrames
weakdep extension — call `using DataFrames` to activate it.
"""
function dataframe end

include("connect.jl")
include("queries.jl")
include("reconstruct.jl")
include("trace.jl")

end # module
```
Note: `queries.jl`, `reconstruct.jl`, `trace.jl` don't exist yet; create empty stubs so the module loads:
```bash
: > packages/HimalayaDB/src/queries.jl
: > packages/HimalayaDB/src/reconstruct.jl
: > packages/HimalayaDB/src/trace.jl
```

- [ ] **Step 4: Write `connect.jl`**

`packages/HimalayaDB/src/connect.jl`:
```julia
"""
    default_himalaya_db_path() -> String

Resolve the DB path from `HIMALAYA_DB_PATH`, else `~/.himalaya/himalaya.db`.
Never creates directories (unlike HimalayaUI's `default_db_path`).
"""
default_himalaya_db_path() =
    get(ENV, "HIMALAYA_DB_PATH", joinpath(homedir(), ".himalaya", "himalaya.db"))

"""
    connect(path = default_himalaya_db_path()) -> SQLite.DB

Open the HimalayaUI database read-only. Errors if the file is missing.
Never runs migrations and never chmods the file.
"""
function connect(path::AbstractString = default_himalaya_db_path())
    isfile(path) || throw(ArgumentError(
        "HimalayaDB.connect: no database at $path (set HIMALAYA_DB_PATH?)"))
    db = SQLite.DB(path)
    # ponytail: SQLite.jl 1.8 opens RW at the OS layer (no readonly flag on DB());
    # query_only=ON makes every write/DDL fail loudly, which is all a reader needs.
    # Upgrade path if a truly read-only filesystem is ever required: drop to
    # SQLite.C.sqlite3_open_v2 with SQLITE_OPEN_READONLY.
    DBInterface.execute(db, "PRAGMA query_only = ON")
    return db
end
```

- [ ] **Step 5: Write the shared test fixture**

`packages/HimalayaDB/test/fixture.jl`:
```julia
using SQLite, DBInterface, Tables
import HimalayaUI

"""
    build_fixture(path, analysis_dir) -> NamedTuple

Create a schema-correct DB at `path` (via HimalayaUI.open_db) populated with a
deterministic sample: 3 auto peaks (q = 0.10, 0.1414, 0.1732), one `exclude`
curation on the middle peak, one `add` curation at q = 0.20, one candidate
`Pn3m` index supported by two auto peaks, and a confirmed custom index group.
Also creates a second, UNCURATED exposure (`exposure2_id`): an `active` auto
group with a member but no custom group — used to prove `confirmed_indices`
filters on `kind='custom'`, not `active=1`.
Direct INSERTs are fine here — a reader only cares that the view tables are populated.
"""
function build_fixture(path::AbstractString, analysis_dir::AbstractString)
    db = HimalayaUI.open_db(path)  # creates schema + migrations
    try
        experiment_id = HimalayaUI.create_experiment!(db; name="fixture",
            path=analysis_dir, data_dir=analysis_dir, analysis_dir=analysis_dir)
        sample_id = HimalayaUI.create_sample!(db; experiment_id=experiment_id,
            name="s1", display_name="Sample 1")
        exposure_id = HimalayaUI.create_exposure!(db; sample_id=sample_id,
            filename="s1", kind="file", status="accepted")

        auto_peak_ids = Int[]
        for (q, sharp, fidx) in [(0.10, 1.0, 10), (0.1414, 0.8, 20), (0.1732, 0.6, 30)]
            r = DBInterface.execute(db,
                "INSERT INTO auto_peaks (exposure_id, q, intensity, prominence, sharpness, findpeaks_index) VALUES (?,?,?,?,?,?)",
                [exposure_id, q, 100.0, 5.0, sharp, fidx])
            push!(auto_peak_ids, Int(DBInterface.lastrowid(r)))
        end

        DBInterface.execute(db,
            "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)",
            [exposure_id, 0.1414])
        DBInterface.execute(db,
            "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)",
            [exposure_id, 0.20])

        ri = DBInterface.execute(db,
            "INSERT INTO indices (exposure_id, phase, basis, score, kind, status) VALUES (?, 'Pn3m', ?, ?, 'auto', 'candidate')",
            [exposure_id, 0.10, 0.9])
        index_id = Int(DBInterface.lastrowid(ri))
        # supporting peaks: ratio positions 1 and 2 -> first two auto peaks
        for (pos, pid) in [(1, auto_peak_ids[1]), (2, auto_peak_ids[3])]
            DBInterface.execute(db,
                "INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position) VALUES (?, ?, 'auto', ?)",
                [index_id, pid, pos])
        end

        rg = DBInterface.execute(db,
            "INSERT INTO index_groups (exposure_id, kind, active) VALUES (?, 'custom', 1)",
            [exposure_id])
        group_id = Int(DBInterface.lastrowid(rg))
        DBInterface.execute(db,
            "INSERT INTO index_group_members (group_id, index_id) VALUES (?, ?)",
            [group_id, index_id])

        # Second, UNCURATED exposure: an active auto group with a member, but NO
        # custom group. confirmed_indices must return EMPTY here (the curator
        # committed nothing) — this pins the filter to kind='custom', not active=1.
        exposure2_id = HimalayaUI.create_exposure!(db; sample_id=sample_id,
            filename="s2", kind="file", status="accepted")
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, intensity, prominence, sharpness, findpeaks_index) VALUES (?,?,?,?,?,?)",
            [exposure2_id, 0.11, 100.0, 5.0, 1.0, 12])
        ri2 = DBInterface.execute(db,
            "INSERT INTO indices (exposure_id, phase, basis, score, kind, status) VALUES (?, 'Im3m', ?, ?, 'auto', 'candidate')",
            [exposure2_id, 0.11, 0.5])
        index2_id = Int(DBInterface.lastrowid(ri2))
        rg2 = DBInterface.execute(db,
            "INSERT INTO index_groups (exposure_id, kind, active) VALUES (?, 'auto', 1)",
            [exposure2_id])
        group2_id = Int(DBInterface.lastrowid(rg2))
        DBInterface.execute(db,
            "INSERT INTO index_group_members (group_id, index_id) VALUES (?, ?)",
            [group2_id, index2_id])

        return (; experiment_id, sample_id, exposure_id, exposure2_id,
                  index_id, group_id, auto_peak_ids)
    finally
        close(db)
    end
end
```

- [ ] **Step 6: Write `test/Project.toml` and `test/runtests.jl`**

`packages/HimalayaDB/test/Project.toml` (replace UUIDs from Step 1):
```toml
[deps]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
SQLite = "0aa819cd-b072-5ff4-a722-6bc24af294d9"
DBInterface = "a10d1c49-ce27-4219-8d33-6db1a4562965"
Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Himalaya = "<HIMALAYA_UUID>"
HimalayaUI = "<HIMALAYAUI_UUID>"

[sources]
Himalaya = {path = "../../.."}
HimalayaUI = {path = "../../HimalayaUI"}
```
`packages/HimalayaDB/test/runtests.jl`:
```julia
using Test
using HimalayaDB
include("fixture.jl")

@testset "HimalayaDB" begin
    include("test_connect.jl")
end
```
(Later tasks add their `include(...)` lines here.)

- [ ] **Step 7: Write the failing connect test**

`packages/HimalayaDB/test/test_connect.jl`:
```julia
@testset "connect" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    build_fixture(dbpath, dir)

    db = HimalayaDB.connect(dbpath)
    @test db isa SQLite.DB
    # read works
    @test !isempty(Tables.rowtable(DBInterface.execute(db, "SELECT id FROM experiments")))
    # writes are blocked by query_only=ON
    @test_throws Exception DBInterface.execute(db, "CREATE TABLE t (x INTEGER)")
    close(db)

    # missing file errors, does not create it
    missing_path = joinpath(dir, "nope.db")
    @test_throws ArgumentError HimalayaDB.connect(missing_path)
    @test !isfile(missing_path)
end
```

- [ ] **Step 8: Instantiate and run — verify pass**

Run:
```bash
julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.instantiate()'
julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20
```
Expected: the `connect` testset passes. (If `Pkg.test` can't see the test `[sources]` paths, adjust the relative paths in `test/Project.toml` to match the actual worktree depth.)

- [ ] **Step 9: Commit**

```bash
git add packages/HimalayaDB
git commit -m "feat(HimalayaDB): package skeleton + read-only connect + test fixture"
```

---

### Task 2: Navigation readers (`experiments`, `samples`, `exposures`)

**Files:**
- Modify: `packages/HimalayaDB/src/queries.jl`
- Modify: `packages/HimalayaDB/test/runtests.jl`
- Test: `packages/HimalayaDB/test/test_queries.jl`

**Interfaces:**
- Consumes: `connect`, `build_fixture`.
- Produces: `experiments(db) -> rows`; `samples(db; experiment=nothing) -> rows`; `exposures(db; sample=nothing) -> rows`. Rows are `Vector{<:NamedTuple}`.

- [ ] **Step 1: Write the failing test**

`packages/HimalayaDB/test/test_queries.jl`:
```julia
@testset "navigation" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    exps = experiments(db)
    @test length(exps) == 1
    @test exps[1].id == ids.experiment_id
    @test exps[1].name == "fixture"

    smps = samples(db; experiment=ids.experiment_id)
    @test length(smps) == 1
    @test smps[1].id == ids.sample_id

    exps_all = samples(db)  # no filter
    @test length(exps_all) == 1

    exs = exposures(db; sample=ids.sample_id)
    @test length(exs) == 2                 # exposure_id + the uncurated exposure2_id
    @test exs[1].id == ids.exposure_id
    @test exs[1].filename == "s1"

    close(db)
end
```
Add `include("test_queries.jl")` to `runtests.jl`.

- [ ] **Step 2: Run — verify it fails**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: FAIL — `UndefVarError: experiments not defined`.

- [ ] **Step 3: Implement**

Append to `packages/HimalayaDB/src/queries.jl`:
```julia
"""
    experiments(db) -> Vector{<:NamedTuple}

All experiments, ordered by id.
"""
experiments(db::SQLite.DB) = Tables.rowtable(DBInterface.execute(db,
    "SELECT id, name, path, data_dir, analysis_dir, experiment_type, energy_kev, flight_path_m, created_at FROM experiments ORDER BY id"))

"""
    samples(db; experiment=nothing) -> Vector{<:NamedTuple}

Samples, optionally filtered to one experiment id.
"""
function samples(db::SQLite.DB; experiment::Union{Integer,Nothing}=nothing)
    if experiment === nothing
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, experiment_id, name, notes FROM samples ORDER BY id"))
    else
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, experiment_id, name, notes FROM samples WHERE experiment_id = ? ORDER BY id",
            [Int(experiment)]))
    end
end

"""
    exposures(db; sample=nothing) -> Vector{<:NamedTuple}

Exposures, optionally filtered to one sample id.
"""
function exposures(db::SQLite.DB; sample::Union{Integer,Nothing}=nothing)
    if sample === nothing
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, sample_id, filename, kind, selected, status, image_path FROM exposures ORDER BY id"))
    else
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, sample_id, filename, kind, selected, status, image_path FROM exposures WHERE sample_id = ? ORDER BY id",
            [Int(sample)]))
    end
end
```

- [ ] **Step 4: Run — verify pass**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaDB/src/queries.jl packages/HimalayaDB/test
git commit -m "feat(HimalayaDB): experiments/samples/exposures navigation readers"
```

---

### Task 3: `curated_peaks` (the effective-peaks UNION)

**Files:**
- Modify: `packages/HimalayaDB/src/queries.jl`
- Test: `packages/HimalayaDB/test/test_queries.jl`

**Interfaces:**
- Produces: `curated_peaks(db, exposure_id) -> rows` with columns `id, exposure_id, q, intensity, prominence, sharpness, source, excluded`, sorted by `q`. `source ∈ {"auto","manual"}`; `excluded ∈ {0,1}` (only auto peaks matched to an `exclude` curation within tolerance are `1`).

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaDB/test/test_queries.jl`:
```julia
@testset "curated_peaks" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    peaks = curated_peaks(db, ids.exposure_id)
    # 3 auto + 1 manual add = 4 rows
    @test length(peaks) == 4
    autos = filter(p -> p.source == "auto", peaks)
    manuals = filter(p -> p.source == "manual", peaks)
    @test length(autos) == 3
    @test length(manuals) == 1
    @test manuals[1].q == 0.20
    # the middle auto peak (0.1414) is excluded
    excluded = filter(p -> p.excluded == 1, peaks)
    @test length(excluded) == 1
    @test excluded[1].q == 0.1414
    # sorted by q
    @test issorted([p.q for p in peaks])

    close(db)
end
```

- [ ] **Step 2: Run — verify it fails**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: FAIL — `UndefVarError: curated_peaks not defined`.

- [ ] **Step 3: Implement (port of `HimalayaUI.get_peaks_for_exposure`, pipeline.jl:600)**

Append to `packages/HimalayaDB/src/queries.jl`:
```julia
"""
    curated_peaks(db, exposure_id) -> Vector{<:NamedTuple}

The effective (curated) peak set for an exposure: `auto_peaks ∪ adds − excludes`.
Each row is tagged `source ∈ {"auto","manual"}` and `excluded ∈ {0,1}`.
Excludes are matched to auto peaks by q within `MAX(1e-6, ABS(q)*0.001)`.

Mirrors `HimalayaUI.get_peaks_for_exposure`; the contract test guards drift.
"""
curated_peaks(db::SQLite.DB, exposure_id::Integer) = Tables.rowtable(DBInterface.execute(db, """
    SELECT a.id, a.exposure_id, a.q, a.intensity, a.prominence, a.sharpness,
           'auto' AS source,
           CASE WHEN c.q IS NOT NULL THEN 1 ELSE 0 END AS excluded
    FROM auto_peaks a
    LEFT JOIN peak_curations c
        ON c.exposure_id = a.exposure_id
       AND c.kind = 'exclude'
       AND ABS(c.q - a.q) <= MAX(1e-6, ABS(a.q) * 0.001)
    WHERE a.exposure_id = ?
    UNION ALL
    SELECT id, exposure_id, q, NULL AS intensity, NULL AS prominence, NULL AS sharpness,
           'manual' AS source, 0 AS excluded
    FROM peak_curations
    WHERE exposure_id = ? AND kind = 'add'
    ORDER BY q
""", [Int(exposure_id), Int(exposure_id)]))
```

- [ ] **Step 4: Run — verify pass**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaDB/src/queries.jl packages/HimalayaDB/test/test_queries.jl
git commit -m "feat(HimalayaDB): curated_peaks effective-peaks UNION"
```

---

### Task 4: `index_candidates` and `confirmed_indices`

> See "Post-implementation corrections" below — `confirmed_indices` shipped against the `assignments`/`assignment_members` model, not the `custom` `index_groups` model described here.

**Files:**
- Modify: `packages/HimalayaDB/src/queries.jl`
- Test: `packages/HimalayaDB/test/test_queries.jl`

**Interfaces:**
- Produces: `index_candidates(db, exposure_id) -> rows` (all `indices`, columns `id, exposure_id, phase, basis, score, r_squared, lattice_d, status, kind, inputs_hash`, sorted `score DESC`). `confirmed_indices(db, exposure_id) -> rows` (same columns; only indices that are members of the exposure's `custom` group, sorted `score DESC`).

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaDB/test/test_queries.jl`:
```julia
@testset "indices" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    cands = index_candidates(db, ids.exposure_id)
    @test length(cands) == 1
    @test cands[1].id == ids.index_id
    @test cands[1].phase == "Pn3m"
    @test cands[1].kind == "auto"

    conf = confirmed_indices(db, ids.exposure_id)
    @test length(conf) == 1
    @test conf[1].id == ids.index_id
    @test conf[1].phase == "Pn3m"

    # uncurated exposure2: has a candidate index in an ACTIVE auto group, but no
    # custom group -> confirmed_indices must be empty (filter is kind='custom',
    # not active=1). index_candidates still sees the index.
    @test length(index_candidates(db, ids.exposure2_id)) == 1
    @test isempty(confirmed_indices(db, ids.exposure2_id))

    close(db)
end
```

- [ ] **Step 2: Run — verify it fails**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: FAIL — `UndefVarError: index_candidates not defined`.

- [ ] **Step 3: Implement**

Append to `packages/HimalayaDB/src/queries.jl`:
```julia
const _INDEX_COLS = "id, exposure_id, phase, basis, score, r_squared, lattice_d, status, kind, inputs_hash"

"""
    index_candidates(db, exposure_id) -> Vector{<:NamedTuple}

All candidate index-choices for an exposure (auto + speculative), sorted by score.
Mirrors `HimalayaUI.get_indices_for_exposure`.
"""
index_candidates(db::SQLite.DB, exposure_id::Integer) = Tables.rowtable(DBInterface.execute(db,
    "SELECT $_INDEX_COLS FROM indices WHERE exposure_id = ? ORDER BY score DESC",
    [Int(exposure_id)]))

"""
    confirmed_indices(db, exposure_id) -> Vector{<:NamedTuple}

The human-confirmed indices: members of the exposure's `kind='custom'` index
group — what the curator committed to. Sorted by score. Empty when the curator
has never touched the exposure (no custom group exists yet).

Verified against HimalayaUI's write/read paths: human curation always lands in
the on-demand `kind='custom'` group (routes_analysis.jl `ensure_custom_group!`),
and `kind='custom' ⟹ active=1`, so filtering on `active=1` alone would wrongly
include the pre-curation auto group. No single HimalayaUI getter returns this
shape, so the members join is composed here.
"""
confirmed_indices(db::SQLite.DB, exposure_id::Integer) = Tables.rowtable(DBInterface.execute(db, """
    SELECT i.id, i.exposure_id, i.phase, i.basis, i.score, i.r_squared,
           i.lattice_d, i.status, i.kind, i.inputs_hash
    FROM indices i
    JOIN index_group_members m ON m.index_id = i.id
    JOIN index_groups g        ON g.id = m.group_id
    WHERE g.exposure_id = ? AND g.kind = 'custom'
    ORDER BY i.score DESC
""", [Int(exposure_id)]))
```
Note: every column is qualified with `i.` because `index_groups` shares the
column names `exposure_id` and `kind` with `indices` — an unqualified SELECT
list would throw `ambiguous column name`.

- [ ] **Step 4: Run — verify pass**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaDB/src/queries.jl packages/HimalayaDB/test/test_queries.jl
git commit -m "feat(HimalayaDB): index_candidates + confirmed_indices readers"
```

---

### Task 5: `reconstruct_index` (opt-in core `Himalaya.Index{P}`)

**Files:**
- Modify: `packages/HimalayaDB/src/reconstruct.jl`
- Modify: `packages/HimalayaDB/test/runtests.jl`
- Test: `packages/HimalayaDB/test/test_reconstruct.jl`

**Interfaces:**
- Produces: `reconstruct_index(db, index_id) -> Himalaya.Index{P}` where `P <: Himalaya.Phase`. Throws `ArgumentError` on unknown phase or missing index.

- [ ] **Step 1: Write the failing test**

`packages/HimalayaDB/test/test_reconstruct.jl`:
```julia
@testset "reconstruct_index" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    idx = reconstruct_index(db, ids.index_id)
    @test idx isa Himalaya.Index
    @test Himalaya.phase(idx) === Himalaya.Pn3m
    @test Himalaya.basis(idx) == 0.10
    # two supporting peaks at ratio positions 1 and 2
    @test Himalaya.numpeaks(idx) == 2

    @test_throws ArgumentError reconstruct_index(db, 999999)
    close(db)
end
```
Add `include("test_reconstruct.jl")` to `runtests.jl`.

- [ ] **Step 2: Run — verify it fails**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: FAIL — `UndefVarError: reconstruct_index not defined`.

- [ ] **Step 3: Implement**

`packages/HimalayaDB/src/reconstruct.jl`:
```julia
# Local copy of HimalayaUI's phase-name→type resolver (speculative.jl:19).
# Tiny and stable; reimplemented rather than depending on a HimalayaUI internal.
function resolve_phase(name::AbstractString)
    bare = last(split(String(name), '.'))
    P = try
        getfield(Himalaya, Symbol(bare))
    catch
        return nothing
    end
    (P isa Type && P <: Himalaya.Phase) ? P : nothing
end

"""
    reconstruct_index(db, index_id) -> Himalaya.Index{P}

Rebuild a core `Himalaya.Index{P}` from an `indices` row and its `index_peaks`
supporting peaks. Peaks/sharpness are sparse vectors indexed by ratio position.
Throws `ArgumentError` if the index is missing or its phase name is unknown.
"""
function reconstruct_index(db::SQLite.DB, index_id::Integer)
    irows = Tables.rowtable(DBInterface.execute(db,
        "SELECT phase, basis FROM indices WHERE id = ?", [Int(index_id)]))
    isempty(irows) && throw(ArgumentError("reconstruct_index: no index $index_id"))
    irow = irows[1]

    P = resolve_phase(irow.phase)
    P === nothing && throw(ArgumentError(
        "reconstruct_index: unknown phase $(irow.phase) for index $index_id"))

    prows = Tables.rowtable(DBInterface.execute(db, """
        SELECT ip.ratio_position,
               COALESCE(ap.q, pc.q) AS q_observed,
               ap.sharpness AS sharpness
        FROM index_peaks ip
        LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
        LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
        WHERE ip.index_id = ?
        ORDER BY ip.ratio_position
    """, [Int(index_id)]))

    n = isempty(prows) ? 0 : maximum(Int(r.ratio_position) for r in prows)
    peaks = spzeros(Float64, n)
    sharp = spzeros(Float64, n)
    for r in prows
        pos = Int(r.ratio_position)
        peaks[pos] = Float64(r.q_observed)
        sharp[pos] = r.sharpness === missing ? 0.0 : Float64(r.sharpness)
    end

    return Himalaya.Index{P}(Float64(irow.basis), peaks, sharp)
end
```

- [ ] **Step 4: Run — verify pass**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaDB/src/reconstruct.jl packages/HimalayaDB/test/test_reconstruct.jl packages/HimalayaDB/test/runtests.jl
git commit -m "feat(HimalayaDB): reconstruct_index -> Himalaya.Index{P}"
```

---

### Task 6: `load_trace` (opt-in `.dat` loading)

**Files:**
- Modify: `packages/HimalayaDB/src/trace.jl`
- Modify: `packages/HimalayaDB/test/runtests.jl`
- Test: `packages/HimalayaDB/test/test_trace.jl`

**Interfaces:**
- Produces: `load_dat(path) -> (q, I, σ)` (three `Vector{Float64}`); `load_trace(db, exposure_id; pattern=nothing) -> (q, I, σ)` (default pattern read from the experiment's `integration_pattern` column, falling back to `"{name}.dat"` when NULL). Throws with an actionable message if the exposure or `.dat` file is missing.

- [ ] **Step 1: Write the failing test**

`packages/HimalayaDB/test/test_trace.jl`:
```julia
@testset "load_trace" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    # fixture filename is "s1"; default pattern -> "s1.dat" in analysis_dir (=dir)
    open(joinpath(dir, "s1.dat"), "w") do io
        println(io, "# q I sigma")
        println(io, "0.10 100.0 1.0")
        println(io, "0.20 50.0 0.5")
    end

    db = HimalayaDB.connect(dbpath)
    q, I, σ = load_trace(db, ids.exposure_id)
    @test q == [0.10, 0.20]
    @test I == [100.0, 50.0]
    @test σ == [1.0, 0.5]

    # explicit pattern= override
    open(joinpath(dir, "custom_s1.dat"), "w") do io
        println(io, "0.30 10.0 0.1")
    end
    qo, Io, σo = load_trace(db, ids.exposure_id; pattern="custom_{name}.dat")
    @test qo == [0.30]

    # missing exposure errors
    @test_throws ArgumentError load_trace(db, 999999)
    close(db)
end
```
Add `include("test_trace.jl")` to `runtests.jl`.

Note: the fixture's experiment has `integration_pattern = NULL`, so `load_trace(db, exposure_id)` with no keyword exercises the NULL→`"{name}.dat"` fallback and resolves `s1.dat`; the second call exercises the explicit-override branch.

- [ ] **Step 2: Run — verify it fails**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: FAIL — `UndefVarError: load_trace not defined`.

- [ ] **Step 3: Implement**

`packages/HimalayaDB/src/trace.jl`:
```julia
"""
    load_dat(path) -> (q, I, σ)

Parse a whitespace-delimited `.dat` trace file: columns q, I, σ. Blank lines and
`#` comments are skipped. Errors if a data row has fewer than 3 columns.
"""
function load_dat(path::AbstractString)
    q = Float64[]; I = Float64[]; σ = Float64[]
    for ln in eachline(path)
        s = strip(ln)
        (isempty(s) || startswith(s, '#')) && continue
        cols = split(s)
        length(cols) >= 3 || error("$path: expected ≥3 columns, got $(length(cols))")
        push!(q, parse(Float64, cols[1]))
        push!(I, parse(Float64, cols[2]))
        push!(σ, parse(Float64, cols[3]))
    end
    return (q, I, σ)
end

"""
    load_trace(db, exposure_id; pattern=nothing) -> (q, I, σ)

Resolve an exposure's on-disk `.dat` path and parse it. The path is
`joinpath(experiments.analysis_dir, replace(pat, "{name}" => exposures.filename))`,
where `pat` is the `pattern` keyword if given, else the experiment's
`integration_pattern` column, else `"{name}.dat"` when that column is NULL.

Pass `pattern=` to override (e.g. `"{name}_tot.dat"`). The file-not-found error
names the path it tried, which reveals a pattern mismatch.
"""
function load_trace(db::SQLite.DB, exposure_id::Integer;
                    pattern::Union{AbstractString,Nothing}=nothing)
    rows = Tables.rowtable(DBInterface.execute(db, """
        SELECT e.filename, x.analysis_dir, x.integration_pattern
        FROM exposures e
        JOIN samples s     ON s.id = e.sample_id
        JOIN experiments x ON x.id = s.experiment_id
        WHERE e.id = ?
    """, [Int(exposure_id)]))
    isempty(rows) && throw(ArgumentError("load_trace: no exposure $exposure_id"))
    row = rows[1]
    row.filename === missing && throw(ArgumentError(
        "load_trace: exposure $exposure_id has no filename"))
    pat = pattern !== nothing ? String(pattern) :
          (row.integration_pattern === missing ? "{name}.dat" : String(row.integration_pattern))
    path = joinpath(String(row.analysis_dir),
                    replace(pat, "{name}" => String(row.filename)))
    isfile(path) || throw(ArgumentError(
        "load_trace: trace file not found at $path (wrong `pattern=`, or data dir not present?)"))
    return load_dat(path)
end
```

- [ ] **Step 4: Run — verify pass**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaDB/src/trace.jl packages/HimalayaDB/test/test_trace.jl packages/HimalayaDB/test/runtests.jl
git commit -m "feat(HimalayaDB): load_trace opt-in .dat loading"
```

---

### Task 7: DataFrames weakdep extension

**Files:**
- Create: `packages/HimalayaDB/ext/DataFramesExt.jl`
- Modify: `packages/HimalayaDB/test/runtests.jl`
- Test: `packages/HimalayaDB/test/test_dataframes.jl`

**Interfaces:**
- Produces: `HimalayaDB.dataframe(rows) -> DataFrame` (method added when `DataFrames` is loaded).

- [ ] **Step 1: Write the failing test**

`packages/HimalayaDB/test/test_dataframes.jl`:
```julia
using DataFrames

@testset "dataframe extension" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    df = dataframe(curated_peaks(db, ids.exposure_id))
    @test df isa DataFrame
    @test nrow(df) == 4
    @test "source" in names(df)
    close(db)
end
```
Add `include("test_dataframes.jl")` to `runtests.jl`.

- [ ] **Step 2: Run — verify it fails**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: FAIL — `MethodError: no method matching dataframe(...)` (the stub has no methods until the extension loads).

- [ ] **Step 3: Implement the extension**

`packages/HimalayaDB/ext/DataFramesExt.jl`:
```julia
module DataFramesExt

using HimalayaDB, DataFrames, Tables

HimalayaDB.dataframe(rows) = DataFrame(Tables.columntable(rows))

end # module
```

- [ ] **Step 4: Run — verify pass**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: PASS. (The extension auto-loads because `DataFrames` is in the test env; `using DataFrames` in the test triggers it.)

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaDB/ext/DataFramesExt.jl packages/HimalayaDB/test/test_dataframes.jl packages/HimalayaDB/test/runtests.jl
git commit -m "feat(HimalayaDB): DataFrames weakdep extension"
```

---

### Task 8: Contract test — cross-check vs HimalayaUI (the drift guard)

**Files:**
- Modify: `packages/HimalayaDB/test/runtests.jl`
- Test: `packages/HimalayaDB/test/test_contract.jl`

**Interfaces:**
- Consumes: `build_fixture`, `connect`, `curated_peaks`, `index_candidates`; `HimalayaUI.open_db`, `HimalayaUI.get_peaks_for_exposure`, `HimalayaUI.get_indices_for_exposure`.

- [ ] **Step 1: Write the contract test**

`packages/HimalayaDB/test/test_contract.jl`:
```julia
import HimalayaUI

@testset "contract vs HimalayaUI" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)

    # HimalayaUI's own readers (read-write open) on the same DB
    uidb = HimalayaUI.open_db(dbpath)
    ui_peaks = HimalayaUI.get_peaks_for_exposure(uidb, ids.exposure_id)
    ui_indices = HimalayaUI.get_indices_for_exposure(uidb, ids.exposure_id)
    close(uidb)

    db = HimalayaDB.connect(dbpath)
    db_peaks = curated_peaks(db, ids.exposure_id)
    db_indices = index_candidates(db, ids.exposure_id)
    close(db)

    # Same rows, same order, same effective-peaks logic.
    @test [(p.q, p.source, p.excluded) for p in db_peaks] ==
          [(p.q, p.source, p.excluded) for p in ui_peaks]
    @test [(i.id, i.phase, i.score) for i in db_indices] ==
          [(i.id, i.phase, i.score) for i in ui_indices]
end
```
Add `include("test_contract.jl")` to `runtests.jl`.

- [ ] **Step 2: Run — verify it passes**

Run: `julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: PASS. (This confirms the duplicated SQL matches HimalayaUI. If it ever fails after a HimalayaUI schema/logic change, that's the drift guard doing its job — reconcile `curated_peaks`/`index_candidates` with the changed HimalayaUI query.)

- [ ] **Step 3: Add a README**

Create `packages/HimalayaDB/README.md` with a short usage example:
```markdown
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
```

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaDB/test/test_contract.jl packages/HimalayaDB/test/runtests.jl packages/HimalayaDB/README.md
git commit -m "test(HimalayaDB): contract cross-check vs HimalayaUI + README"
```

---

## Self-Review

**Spec coverage** (against `docs/superpowers/specs/2026-07-14-himalayadb-design.md`):
- §Connection (read-only, no migrations/chmod, WAL risk) → Task 1 (`connect`, `query_only=ON`, `isfile` guard, ponytail ceiling note). ✓
- §Read API (`experiments`/`samples`/`exposures`, `curated_peaks`, `index_candidates`, `confirmed_indices`) → Tasks 2–4. ✓ (`confirmed_indices` composes the members join per the extraction correction — noted in-task.)
- §Opt-in `reconstruct_index` → Task 5. ✓
- §Opt-in `load_trace` → Task 6. ✓
- §DataFrames weakdep → Task 7. ✓
- §Contract/drift guard → Task 8 (scoped to `curated_peaks` + `index_candidates`, the two with direct HimalayaUI equivalents; `confirmed_indices` unit-tested in Task 4, since HimalayaUI has no equivalent getter). ✓
- §Return types (Tables.jl rows) → all query tasks. ✓
- §Non-goals (no writes/HTTP/schema-ownership/series/tags) → Global Constraints + scope. ✓

**Placeholder scan:** No "TBD/TODO". The two `<..._UUID>` tokens are real values resolved in Task 1 Step 1 (not placeholders — the plan says exactly how to obtain them). The `[sources]` relative paths carry an explicit "adjust to actual worktree depth" instruction.

**Type consistency:** `connect`/`curated_peaks`/`index_candidates`/`confirmed_indices`/`reconstruct_index`/`load_trace`/`load_dat`/`dataframe`/`resolve_phase`/`build_fixture` names and signatures are consistent across tasks and match the exports in Task 1. `curated_peaks` output columns (`id, exposure_id, q, intensity, prominence, sharpness, source, excluded`) match both the Task 3 SQL and the Task 8 contract assertion.

---

## Post-implementation corrections

The plan above was written from a data-model map extracted against the *main clone*, which sat on an older branch (`speculative-peak-durability`) with a different schema than the worktree (`main`) the package was actually built and tested against. Two corrections were applied during execution; the **shipped `packages/HimalayaDB/` code is the source of truth**, not the code blocks above:

1. **`samples` reader** — dropped `display_name` (the `migrate_samples_name_collapse!` migration renames `display_name`→`name` on every fresh DB, so no `display_name` column exists). Also: `create_exposure!` now requires `experiment_id`; `create_sample!` takes only `name` (no `display_name`). The fixture was adjusted accordingly.

2. **`confirmed_indices` (Critical, caught by the final whole-branch review)** — the plan's `index_groups`/`index_group_members` `kind='custom'` path was **retired by HimalayaUI's D-10 plotting redesign**. Confirmation now lives in `assignments(exposure_id, state)` + `assignment_members(exposure_id, index_id)`. The shipped query reads `assignment_members JOIN indices` gated by `COALESCE(assignments.state,'indexed')='indexed'`, `ORDER BY i.score DESC NULLS LAST, i.id ASC` — mirroring HimalayaUI's own confirmed-index read. It returns the exposure's durable *indexed assignment* (which may be auto-seeded on migrated DBs), not strictly a human-only confirmation. The fixture confirms through the real `assignments`/`assignment_members` tables (and uses a `state='form_factor'` exposure to exercise the state gate); the `group_id` fixture field was removed.

3. **`reconstruct_index` sizing** — `SparseVector` length uses `length(Himalaya.phaseratios(P))` (the codebase convention) rather than `maximum(ratio_position)`, so position-indexed access and `Base.show` are correct.

**Known follow-ups (non-blocking):** `confirmed_indices` has no contract test (HimalayaUI exposes no single equivalent getter — its confirmed read lives inside `compute_member_snapshot`); assorted low-risk coverage gaps and defensive-polish Minors are recorded in the SDD progress ledger.
