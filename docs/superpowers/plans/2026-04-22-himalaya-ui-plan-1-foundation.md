# HimalayaUI Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the `HimalayaUI` sub-package in the monorepo with a working SQLite-backed analysis pipeline that can ingest a manifest CSV, run peak-finding + indexing on `.dat` files, persist results, and be driven from the command line.

**Architecture:** `HimalayaUI` lives at `packages/HimalayaUI/` with its own `Project.toml` depending on the root `Himalaya` package via a local path dev dependency. The package is split into focused files: `db.jl` owns all SQLite interaction, `manifest.jl` owns CSV parsing, `datfile.jl` owns `.dat` parsing, `pipeline.jl` orchestrates analysis, and `cli.jl` exposes `himalaya init/analyze/show`. No web server in this plan — that is Plan 2.

**Tech Stack:** Julia 1.9+, SQLite.jl, CSV.jl, DelimitedFiles (stdlib), ArgParse.jl, Himalaya.jl (local)

**This is Plan 1 of 6. Subsequent plans:** Plan 2 — REST API (Oxygen.jl); Plan 3 — Frontend shell (Vite + layout + sample list); Plan 4 — Trace viewer; Plan 5 — Phase panel + Miller plot; Plan 6 — Properties panel.

---

## File Map

| File | Responsibility |
|---|---|
| `packages/HimalayaUI/Project.toml` | Package manifest, deps |
| `packages/HimalayaUI/src/HimalayaUI.jl` | Module entry, includes |
| `packages/HimalayaUI/src/db.jl` | SQLite schema, CRUD helpers |
| `packages/HimalayaUI/src/datfile.jl` | Parse three-column `.dat` files |
| `packages/HimalayaUI/src/manifest.jl` | Parse Google Sheets CSV manifest |
| `packages/HimalayaUI/src/pipeline.jl` | Orchestrate analysis, auto-group, persist |
| `packages/HimalayaUI/src/cli.jl` | `himalaya init / analyze / show` |
| `packages/HimalayaUI/test/runtests.jl` | Test orchestrator |
| `packages/HimalayaUI/test/test_db.jl` | Schema + CRUD tests |
| `packages/HimalayaUI/test/test_datfile.jl` | `.dat` parser tests |
| `packages/HimalayaUI/test/test_manifest.jl` | Manifest parser tests |
| `packages/HimalayaUI/test/test_pipeline.jl` | Pipeline + auto-group tests |

---

### Task 1: Create the HimalayaUI sub-package

**Files:**
- Create: `packages/HimalayaUI/Project.toml`
- Create: `packages/HimalayaUI/src/HimalayaUI.jl`
- Create: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Create directory structure**

```bash
mkdir -p packages/HimalayaUI/src packages/HimalayaUI/test
```

- [ ] **Step 2: Write `Project.toml`**

```toml
name = "HimalayaUI"
uuid = "a1b2c3d4-e5f6-7890-abcd-ef1234567890"
version = "0.1.0"

[deps]
ArgParse = "c7e460c6-2fb9-53a9-8c5b-16f535851c63"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DBInterface = "a10d1c49-ce27-4219-8d33-6db1a4562965"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
Himalaya = "00000000-0000-0000-0000-000000000000"
JSON3 = "0f8b85d8-7e73-4b5c-badb-23f01f67f2e7"
SQLite = "0aa819cd-b072-5ff4-a722-6bc24af294d9"
Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[compat]
julia = "1.9"
```

> Note: the Himalaya UUID above is a placeholder — the actual UUID will be filled in by `Pkg.develop` in Step 4. Leave it as-is for now.

- [ ] **Step 3: Write the module entry point**

`packages/HimalayaUI/src/HimalayaUI.jl`:
```julia
module HimalayaUI

include("db.jl")
include("datfile.jl")
include("manifest.jl")
include("pipeline.jl")
include("cli.jl")

end
```

- [ ] **Step 4: Wire up the local Himalaya dependency**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.develop(path="../.."); Pkg.add(["SQLite", "CSV", "ArgParse", "JSON3", "Tables", "DBInterface"]); Pkg.instantiate()'
cd ../..
```

- [ ] **Step 5: Write the test orchestrator**

`packages/HimalayaUI/test/runtests.jl`:
```julia
using Test

@testset "HimalayaUI" begin
    include("test_db.jl")
    include("test_datfile.jl")
    include("test_manifest.jl")
    include("test_pipeline.jl")
end
```

- [ ] **Step 6: Verify the package loads**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI; println("ok")'
```

Expected output: `ok`

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/
git commit -m "feat: scaffold HimalayaUI sub-package in monorepo"
```

---

### Task 2: SQLite schema and CRUD helpers

**Files:**
- Create: `packages/HimalayaUI/src/db.jl`
- Create: `packages/HimalayaUI/test/test_db.jl`

- [ ] **Step 1: Write the failing test**

`packages/HimalayaUI/test/test_db.jl`:
```julia
using Test, SQLite, DBInterface
using HimalayaUI: create_schema!, create_experiment!, create_sample!,
                  create_exposure!, get_experiment, get_samples, get_exposures

@testset "db schema" begin
    db = SQLite.DB()  # in-memory
    create_schema!(db)

    tables = Set(r[1] for r in DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table'"))

    for t in ["users", "experiments", "samples", "sample_tags",
              "exposures", "exposure_sources", "exposure_tags",
              "peaks", "indices", "index_peaks",
              "index_groups", "index_group_members", "user_actions"]
        @test t in tables
    end
end

@testset "db CRUD" begin
    db = SQLite.DB()
    create_schema!(db)

    exp_id = create_experiment!(db;
        name        = "TestRun",
        path        = "/data/exp1",
        data_dir    = "/data/exp1/data",
        analysis_dir = "/data/exp1/analysis/automatic_analysis")
    @test exp_id == 1

    exp = get_experiment(db, exp_id)
    @test exp.name == "TestRun"
    @test exp.path == "/data/exp1"

    s_id = create_sample!(db; experiment_id = exp_id, label = "D1", name = "UX1")
    @test s_id == 1

    samples = get_samples(db, exp_id)
    @test length(samples) == 1
    @test first(samples).label == "D1"

    e_id = create_exposure!(db; sample_id = s_id, filename = "JC001", kind = "file")
    @test e_id == 1

    exposures = get_exposures(db, s_id)
    @test length(exposures) == 1
    @test first(exposures).filename == "JC001"
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project=packages/HimalayaUI -e '
using Pkg; Pkg.test("HimalayaUI"; test_args=["test_db"])'
```

Expected: error — `create_schema!` not defined.

- [ ] **Step 3: Implement `db.jl`**

`packages/HimalayaUI/src/db.jl`:
```julia
using SQLite, DBInterface, Tables

const SCHEMA = """
CREATE TABLE IF NOT EXISTS users (
    id       INTEGER PRIMARY KEY,
    username TEXT UNIQUE NOT NULL
);

CREATE TABLE IF NOT EXISTS experiments (
    id             INTEGER PRIMARY KEY,
    name           TEXT,
    path           TEXT NOT NULL,
    data_dir       TEXT NOT NULL,
    analysis_dir   TEXT NOT NULL,
    manifest_path  TEXT,
    created_at     DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS samples (
    id            INTEGER PRIMARY KEY,
    experiment_id INTEGER REFERENCES experiments(id),
    label         TEXT,
    name          TEXT,
    notes         TEXT
);

CREATE TABLE IF NOT EXISTS sample_tags (
    id        INTEGER PRIMARY KEY,
    sample_id INTEGER REFERENCES samples(id),
    key       TEXT NOT NULL,
    value     TEXT NOT NULL,
    source    TEXT DEFAULT 'manual'
);

CREATE TABLE IF NOT EXISTS exposures (
    id        INTEGER PRIMARY KEY,
    sample_id INTEGER REFERENCES samples(id),
    filename  TEXT,
    kind      TEXT DEFAULT 'file',
    selected  BOOLEAN DEFAULT FALSE
);

CREATE TABLE IF NOT EXISTS exposure_sources (
    averaged_exposure_id INTEGER REFERENCES exposures(id),
    source_exposure_id   INTEGER REFERENCES exposures(id),
    role                 TEXT DEFAULT 'signal',
    PRIMARY KEY (averaged_exposure_id, source_exposure_id)
);

CREATE TABLE IF NOT EXISTS exposure_tags (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    key         TEXT NOT NULL,
    value       TEXT NOT NULL,
    source      TEXT DEFAULT 'manual'
);

CREATE TABLE IF NOT EXISTS peaks (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    q           REAL NOT NULL,
    intensity   REAL,
    prominence  REAL,
    sharpness   REAL,
    source      TEXT DEFAULT 'auto'
);

CREATE TABLE IF NOT EXISTS indices (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    phase       TEXT NOT NULL,
    basis       REAL NOT NULL,
    score       REAL,
    r_squared   REAL,
    lattice_d   REAL,
    status      TEXT DEFAULT 'candidate'
);

CREATE TABLE IF NOT EXISTS index_peaks (
    index_id       INTEGER REFERENCES indices(id),
    peak_id        INTEGER REFERENCES peaks(id),
    ratio_position INTEGER,
    residual       REAL,
    PRIMARY KEY (index_id, peak_id)
);

CREATE TABLE IF NOT EXISTS index_groups (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    kind        TEXT NOT NULL DEFAULT 'auto',
    active      BOOLEAN DEFAULT FALSE,
    created_by  INTEGER REFERENCES users(id),
    created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS index_group_members (
    group_id  INTEGER REFERENCES index_groups(id),
    index_id  INTEGER REFERENCES indices(id),
    PRIMARY KEY (group_id, index_id)
);

CREATE TABLE IF NOT EXISTS user_actions (
    id          INTEGER PRIMARY KEY,
    user_id     INTEGER REFERENCES users(id),
    timestamp   DATETIME DEFAULT CURRENT_TIMESTAMP,
    action      TEXT,
    entity_type TEXT,
    entity_id   INTEGER,
    note        TEXT
);
"""

function create_schema!(db::SQLite.DB)
    for stmt in split(SCHEMA, ";")
        s = strip(stmt)
        isempty(s) && continue
        DBInterface.execute(db, s)
    end
end

function create_experiment!(db::SQLite.DB;
        name::Union{String,Nothing} = nothing,
        path::String,
        data_dir::String,
        analysis_dir::String,
        manifest_path::Union{String,Nothing} = nothing)
    DBInterface.execute(db,
        "INSERT INTO experiments (name, path, data_dir, analysis_dir, manifest_path)
         VALUES (?, ?, ?, ?, ?)",
        [name, path, data_dir, analysis_dir, manifest_path])
    Int(DBInterface.lastrowid(db))
end

function create_sample!(db::SQLite.DB;
        experiment_id::Int,
        label::Union{String,Nothing} = nothing,
        name::Union{String,Nothing}  = nothing,
        notes::Union{String,Nothing} = nothing)
    DBInterface.execute(db,
        "INSERT INTO samples (experiment_id, label, name, notes) VALUES (?, ?, ?, ?)",
        [experiment_id, label, name, notes])
    Int(DBInterface.lastrowid(db))
end

function create_exposure!(db::SQLite.DB;
        sample_id::Int,
        filename::Union{String,Nothing} = nothing,
        kind::String                    = "file",
        selected::Bool                  = false)
    DBInterface.execute(db,
        "INSERT INTO exposures (sample_id, filename, kind, selected) VALUES (?, ?, ?, ?)",
        [sample_id, filename, kind, Int(selected)])
    Int(DBInterface.lastrowid(db))
end

function get_experiment(db::SQLite.DB, id::Int)
    rows = collect(DBInterface.execute(db,
        "SELECT * FROM experiments WHERE id = ?", [id]))
    isempty(rows) && error("experiment $id not found")
    first(rows)
end

function get_samples(db::SQLite.DB, experiment_id::Int)
    collect(DBInterface.execute(db,
        "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id", [experiment_id]))
end

function get_exposures(db::SQLite.DB, sample_id::Int)
    collect(DBInterface.execute(db,
        "SELECT * FROM exposures WHERE sample_id = ? ORDER BY id", [sample_id]))
end

function open_db(experiment_path::String)::SQLite.DB
    db_path = joinpath(experiment_path, "himalaya.db")
    db = SQLite.DB(db_path)
    create_schema!(db)
    db
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: `db schema` and `db CRUD` pass; other testsets error (files not included yet — that is fine since `runtests.jl` wraps all in one `@testset`). To run only `test_db.jl` in isolation:

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_db.jl")'
```

Expected: all `test_db.jl` tests pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_db.jl
git commit -m "feat(HimalayaUI): SQLite schema and CRUD helpers"
```

---

### Task 3: `.dat` file parser

**Files:**
- Create: `packages/HimalayaUI/src/datfile.jl`
- Create: `packages/HimalayaUI/test/test_datfile.jl`

- [ ] **Step 1: Write the failing test**

`packages/HimalayaUI/test/test_datfile.jl`:
```julia
using Test
using HimalayaUI: load_dat

const EXAMPLE_DAT = joinpath(@__DIR__, "..", "..", "..", "..", "test", "data", "example_tot.dat")

@testset "load_dat" begin
    q, I, σ = load_dat(EXAMPLE_DAT)

    @test length(q) == 922
    @test length(I) == 922
    @test length(σ) == 922

    @test q[1] ≈ 6.258469e-03 atol=1e-8
    @test I[1] ≈ 4.527260e+04 atol=1.0
    @test σ[1] ≈ 2.127736e+02 atol=1e-2

    @test all(q .> 0)
    @test all(I .> 0)
    @test all(σ .> 0)
    @test issorted(q)
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_datfile.jl")'
```

Expected: `load_dat` not defined.

- [ ] **Step 3: Implement `datfile.jl`**

`packages/HimalayaUI/src/datfile.jl`:
```julia
using DelimitedFiles

"""
    load_dat(path) -> (q, I, σ)

Parse a three-column whitespace-separated SAXS integration file.
Returns vectors of scattering vector q, intensity I, and uncertainty σ.
No header rows expected.
"""
function load_dat(path::AbstractString)
    data = readdlm(path, Float64)
    size(data, 2) >= 3 || error("$path: expected ≥3 columns, got $(size(data,2))")
    data[:, 1], data[:, 2], data[:, 3]
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_datfile.jl")'
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/datfile.jl packages/HimalayaUI/test/test_datfile.jl
git commit -m "feat(HimalayaUI): .dat file parser"
```

---

### Task 4: Manifest CSV parser

**Files:**
- Create: `packages/HimalayaUI/src/manifest.jl`
- Create: `packages/HimalayaUI/test/test_manifest.jl`

- [ ] **Step 1: Write the failing test**

`packages/HimalayaUI/test/test_manifest.jl`:
```julia
using Test
using HimalayaUI: expand_filename_range, parse_manifest, ManifestSample

@testset "expand_filename_range" begin
    # prefix-only range: JC001-004
    @test expand_filename_range("JC001-004") == ["JC001", "JC002", "JC003", "JC004"]
    # full-prefix range: JC013-JC016
    @test expand_filename_range("JC013-JC016") == ["JC013", "JC014", "JC015", "JC016"]
    # single filename (no dash)
    @test expand_filename_range("JC001") == ["JC001"]
end

const MANIFEST_CSV = """
#\tSample\tName\tType\tTime(s)\t\t#\t\tFilename(s)\tNotes (Sample)\tNotes (Exposure)
\tFlight Path: DNA, 0.7 m, Capillaries\t\t\t\t\t\t\t\t\t
1\tD1\tUX1\tControl\t20\t\t\t\tJC001-004\tclear\t
2\tD2\tUX2\tControl\t20\t\t\t\tJC005-008\tclear\t
3\tD3\tUL1\tControl\t20\t\t\t\tJC009-JC012\tclear\t
4\tD4\tUL2\tSample\t20\t\t\t\tJC013-JC016\tcondensed\tsq
"""

@testset "parse_manifest" begin
    samples = parse_manifest(IOBuffer(MANIFEST_CSV))

    @test length(samples) == 4

    s1 = samples[1]
    @test s1.label == "D1"
    @test s1.name  == "UX1"
    @test s1.notes_sample == "clear"
    @test s1.filenames == ["JC001", "JC002", "JC003", "JC004"]

    s3 = samples[3]
    @test s3.filenames == ["JC009", "JC010", "JC011", "JC012"]

    s4 = samples[4]
    @test s4.label == "D4"
    @test s4.notes_sample   == "condensed"
    @test s4.notes_exposure == "sq"
    @test s4.filenames == ["JC013", "JC014", "JC015", "JC016"]
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_manifest.jl")'
```

Expected: `expand_filename_range` not defined.

- [ ] **Step 3: Implement `manifest.jl`**

`packages/HimalayaUI/src/manifest.jl`:
```julia
using CSV, Tables

struct ManifestSample
    label          ::String
    name           ::String
    notes_sample   ::String
    notes_exposure ::String
    filenames      ::Vector{String}
end

"""
    expand_filename_range(s) -> Vector{String}

Expand a filename range like "JC001-004" or "JC013-JC016" to individual
filenames. Returns a single-element vector for plain filenames.
"""
function expand_filename_range(s::AbstractString)::Vector{String}
    m = match(r"^([A-Za-z]+)(\d+)-(?:[A-Za-z]*)(\d+)$", s)
    m === nothing && return [s]
    prefix, start_s, stop_s = m[1], m[2], m[3]
    width  = length(start_s)
    start  = parse(Int, start_s)
    stop   = parse(Int, stop_s)
    [string(prefix, lpad(i, width, '0')) for i in start:stop]
end

"""
    parse_manifest(io_or_path) -> Vector{ManifestSample}

Parse a tab-separated manifest CSV exported from Google Sheets.
Skips section header rows (rows where the first column is non-numeric or empty).
"""
function parse_manifest(source)::Vector{ManifestSample}
    lines = readlines(source)
    # Skip the column header row (first line)
    samples = ManifestSample[]
    for line in lines[2:end]
        cols = split(line, '\t')
        length(cols) < 9 && continue
        num_str = strip(cols[1])
        # Skip section header rows (non-numeric first column)
        tryparse(Int, num_str) === nothing && continue

        label          = strip(get(cols, 2,  ""))
        name           = strip(get(cols, 3,  ""))
        notes_sample   = strip(get(cols, 10, ""))
        notes_exposure = strip(get(cols, 11, ""))
        filename_str   = strip(get(cols, 9,  ""))

        isempty(filename_str) && continue
        filenames = expand_filename_range(filename_str)

        push!(samples, ManifestSample(label, name, notes_sample, notes_exposure, filenames))
    end
    samples
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_manifest.jl")'
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/manifest.jl packages/HimalayaUI/test/test_manifest.jl
git commit -m "feat(HimalayaUI): manifest CSV parser with filename range expansion"
```

---

### Task 5: Auto-group algorithm

**Files:**
- Create: `packages/HimalayaUI/src/pipeline.jl` (partial — auto_group only)
- Create: `packages/HimalayaUI/test/test_pipeline.jl` (partial)

- [ ] **Step 1: Write the failing test**

`packages/HimalayaUI/test/test_pipeline.jl`:
```julia
using Test
using Himalaya: indexpeaks, Index, peaks, score
using HimalayaUI: auto_group

@testset "auto_group" begin
    # Three peaks that could match Im3m or Pn3m.
    # Use real indexpeaks to get a realistic set of candidates.
    qs    = [0.1000, 0.1414, 0.2000]   # ratios ≈ 1 : √2 : 2 → Im3m-ish
    proms = [1.0, 0.8, 0.6]

    candidates = indexpeaks(qs, proms)

    group = auto_group(candidates)

    # Group must be non-empty if any candidates exist
    if !isempty(candidates)
        @test !isempty(group)
        # No two indices in the group share a peak
        peak_sets = [Set(peaks(idx)) for idx in group]
        for i in eachindex(peak_sets), j in eachindex(peak_sets)
            i == j && continue
            @test isempty(intersect(peak_sets[i], peak_sets[j]))
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_pipeline.jl")'
```

Expected: `auto_group` not defined.

- [ ] **Step 3: Implement `auto_group` in `pipeline.jl`**

`packages/HimalayaUI/src/pipeline.jl`:
```julia
using Himalaya
using SQLite, DBInterface

"""
    auto_group(indices) -> Vector{Index}

Greedily select a non-overlapping set of indices by descending score.
An index is added to the group only if none of its peaks are already
claimed by a previously selected index.
"""
function auto_group(indices::Vector{<:Himalaya.Index})::Vector{<:Himalaya.Index}
    isempty(indices) && return indices
    sorted    = sort(indices; by = score, rev = true)
    claimed   = Set{Float64}()
    group     = eltype(indices)[]

    for idx in sorted
        idx_peaks = Set(peaks(idx))
        isempty(intersect(idx_peaks, claimed)) || continue
        push!(group, idx)
        union!(claimed, idx_peaks)
    end
    group
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_pipeline.jl")'
```

Expected: `auto_group` testset passes.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/pipeline.jl packages/HimalayaUI/test/test_pipeline.jl
git commit -m "feat(HimalayaUI): greedy auto_group for non-overlapping index selection"
```

---

### Task 6: Persist analysis results to SQLite

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl` (add `persist_analysis!`)
- Modify: `packages/HimalayaUI/test/test_pipeline.jl` (add persistence tests)

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_pipeline.jl`:
```julia
using HimalayaUI: create_schema!, create_experiment!, create_sample!,
                  create_exposure!, persist_analysis!, get_peaks_for_exposure,
                  get_indices_for_exposure, get_groups_for_exposure

@testset "persist_analysis!" begin
    db = SQLite.DB()
    create_schema!(db)
    exp_id  = create_experiment!(db; path="/tmp", data_dir="/tmp/data",
                                     analysis_dir="/tmp/analysis")
    s_id    = create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    e_id    = create_exposure!(db; sample_id=s_id, filename="example_tot.dat")

    dat_path = joinpath(@__DIR__, "..", "..", "..", "..", "test", "data", "example_tot.dat")
    q, I, σ  = HimalayaUI.load_dat(dat_path)
    peaks_result  = Himalaya.findpeaks(q, I, σ)
    candidates    = Himalaya.indexpeaks(peaks_result.q, peaks_result.prominence)
    group_indices = HimalayaUI.auto_group(candidates)

    persist_analysis!(db, e_id, q, I, peaks_result, candidates, group_indices)

    stored_peaks   = get_peaks_for_exposure(db, e_id)
    stored_indices = get_indices_for_exposure(db, e_id)
    stored_groups  = get_groups_for_exposure(db, e_id)

    @test length(stored_peaks) == length(peaks_result.q)
    @test length(stored_indices) == length(candidates)
    @test length(stored_groups) == 1
    @test stored_groups[1].kind   == "auto"
    @test stored_groups[1].active == 1
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_pipeline.jl")'
```

Expected: `persist_analysis!` not defined.

- [ ] **Step 3: Implement persistence helpers in `pipeline.jl`**

Append to `packages/HimalayaUI/src/pipeline.jl`:
```julia
function persist_analysis!(db::SQLite.DB, exposure_id::Int,
                            q_full::Vector{Float64},
                            I_full::Vector{Float64},
                            peaks_result::NamedTuple,
                            candidates::Vector{<:Himalaya.Index},
                            group_indices::Vector{<:Himalaya.Index})

    # Remove old auto peaks and indices for this exposure (preserve manual peaks)
    DBInterface.execute(db,
        "DELETE FROM index_group_members WHERE group_id IN
         (SELECT id FROM index_groups WHERE exposure_id = ? AND kind = 'auto')",
        [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM index_groups WHERE exposure_id = ? AND kind = 'auto'",
        [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM index_peaks WHERE index_id IN
         (SELECT id FROM indices WHERE exposure_id = ?)", [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM indices WHERE exposure_id = ?", [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM peaks WHERE exposure_id = ? AND source = 'auto'", [exposure_id])

    # Persist auto peaks
    q_to_peak_id = Dict{Float64, Int}()
    for (i, qval) in enumerate(peaks_result.q)
        # Find intensity at this q in the full trace
        full_idx = findfirst(==(peaks_result.indices[i]), eachindex(q_full))
        intensity = full_idx !== nothing ? I_full[full_idx] : NaN
        DBInterface.execute(db,
            "INSERT INTO peaks (exposure_id, q, intensity, prominence, sharpness, source)
             VALUES (?, ?, ?, ?, ?, 'auto')",
            [exposure_id, qval, intensity,
             peaks_result.prominence[i], peaks_result.sharpness[i]])
        q_to_peak_id[qval] = Int(DBInterface.lastrowid(db))
    end

    # Persist candidate indices
    candidate_to_db_id = Dict{Int, Int}()
    for (ci, idx) in enumerate(candidates)
        P          = Himalaya.phase(idx)
        fit_result = Himalaya.fit(idx)
        s          = Himalaya.score(idx)
        DBInterface.execute(db,
            "INSERT INTO indices (exposure_id, phase, basis, score, r_squared, lattice_d, status)
             VALUES (?, ?, ?, ?, ?, ?, 'candidate')",
            [exposure_id, string(P), Himalaya.basis(idx),
             s, fit_result.R², fit_result.d])
        db_id = Int(DBInterface.lastrowid(db))
        candidate_to_db_id[ci] = db_id

        # Persist index_peaks
        peak_vals       = Himalaya.peaks(idx)
        ratio_positions = findnz_positions(idx)
        for (rpos, qval) in zip(ratio_positions, peak_vals)
            peak_id = get(q_to_peak_id, qval, nothing)
            peak_id === nothing && continue
            ideal  = Himalaya.phaseratios(P; normalize=true)[rpos] * Himalaya.basis(idx)
            resid  = abs(qval - ideal)
            DBInterface.execute(db,
                "INSERT OR IGNORE INTO index_peaks (index_id, peak_id, ratio_position, residual)
                 VALUES (?, ?, ?, ?)",
                [db_id, peak_id, rpos, resid])
        end
    end

    # Persist auto group
    DBInterface.execute(db,
        "INSERT INTO index_groups (exposure_id, kind, active) VALUES (?, 'auto', 1)",
        [exposure_id])
    group_db_id = Int(DBInterface.lastrowid(db))

    # Find DB ids for group members by matching candidates list
    group_set = Set(group_indices)
    for (ci, idx) in enumerate(candidates)
        idx in group_set || continue
        db_id = candidate_to_db_id[ci]
        DBInterface.execute(db,
            "INSERT INTO index_group_members (group_id, index_id) VALUES (?, ?)",
            [group_db_id, db_id])
    end
end

# Helper: get the ratio positions of assigned peaks from an Index sparse vector
function findnz_positions(idx::Himalaya.Index)
    peak_vec = idx.peaks
    nz_pos, _ = SparseArrays.findnz(peak_vec)
    nz_pos
end

function get_peaks_for_exposure(db::SQLite.DB, exposure_id::Int)
    collect(DBInterface.execute(db,
        "SELECT * FROM peaks WHERE exposure_id = ? ORDER BY q", [exposure_id]))
end

function get_indices_for_exposure(db::SQLite.DB, exposure_id::Int)
    collect(DBInterface.execute(db,
        "SELECT * FROM indices WHERE exposure_id = ? ORDER BY score DESC", [exposure_id]))
end

function get_groups_for_exposure(db::SQLite.DB, exposure_id::Int)
    collect(DBInterface.execute(db,
        "SELECT * FROM index_groups WHERE exposure_id = ? ORDER BY id", [exposure_id]))
end
```

Also add `using SparseArrays` to the top of `pipeline.jl`.

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_pipeline.jl")'
```

Expected: all testsets pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/pipeline.jl packages/HimalayaUI/test/test_pipeline.jl
git commit -m "feat(HimalayaUI): persist peaks, indices, and auto group to SQLite"
```

---

### Task 7: Top-level `analyze_exposure!` and `init_experiment!`

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl` (add top-level functions)
- Modify: `packages/HimalayaUI/test/test_pipeline.jl` (add integration test)

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_pipeline.jl`:
```julia
using HimalayaUI: init_experiment!, analyze_exposure!, open_db

@testset "init_experiment!" begin
    tmp = mktempdir()
    data_dir     = joinpath(tmp, "data")
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(data_dir)
    mkpath(analysis_dir)

    db = open_db(tmp)
    exp_id = init_experiment!(db;
        name         = "TestExp",
        path         = tmp,
        data_dir     = data_dir,
        analysis_dir = analysis_dir)

    @test exp_id == 1
    exp = HimalayaUI.get_experiment(db, exp_id)
    @test exp.name == "TestExp"
end

@testset "analyze_exposure! integration" begin
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)

    # Copy a real .dat file into analysis_dir
    src = joinpath(@__DIR__, "..", "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db     = open_db(tmp)
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    e_id   = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)

    @test length(get_peaks_for_exposure(db, e_id))   > 0
    @test length(get_indices_for_exposure(db, e_id)) > 0
    @test length(get_groups_for_exposure(db, e_id))  == 1
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_pipeline.jl")'
```

Expected: `init_experiment!` / `analyze_exposure!` not defined.

- [ ] **Step 3: Implement in `pipeline.jl`**

Append to `packages/HimalayaUI/src/pipeline.jl`:
```julia
"""
    init_experiment!(db; kwargs...) -> experiment_id

Create the experiment row. Thin wrapper over `create_experiment!`.
"""
function init_experiment!(db::SQLite.DB; kwargs...)
    create_experiment!(db; kwargs...)
end

"""
    analyze_exposure!(db, exposure_id, analysis_dir)

Load the .dat file for `exposure_id`, run findpeaks + indexpeaks,
auto-group results, and persist everything to the DB.
The .dat filename is taken from `exposures.filename` with `.dat` appended.
"""
function analyze_exposure!(db::SQLite.DB, exposure_id::Int, analysis_dir::String)
    row = only(DBInterface.execute(db,
        "SELECT filename FROM exposures WHERE id = ?", [exposure_id]))
    dat_path = joinpath(analysis_dir, row[:filename] * ".dat")
    isfile(dat_path) || error("dat file not found: $dat_path")

    q, I, σ      = load_dat(dat_path)
    peaks_result = Himalaya.findpeaks(q, I, σ)
    candidates   = Himalaya.indexpeaks(peaks_result.q, peaks_result.prominence)
    group        = auto_group(candidates)

    persist_analysis!(db, exposure_id, q, I, peaks_result, candidates, group)
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_pipeline.jl")'
```

Expected: all testsets pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/pipeline.jl packages/HimalayaUI/test/test_pipeline.jl
git commit -m "feat(HimalayaUI): analyze_exposure! and init_experiment! top-level API"
```

---

### Task 8: CLI — `himalaya init` and `himalaya analyze`

**Files:**
- Create: `packages/HimalayaUI/src/cli.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl` (export `main`)

- [ ] **Step 1: Implement `cli.jl`**

`packages/HimalayaUI/src/cli.jl`:
```julia
using ArgParse

function parse_cli(args)
    s = ArgParseSettings(prog = "himalaya")

    @add_arg_table! s begin
        "command"
            help     = "init | analyze | show"
            required = true
    end

    # Parse just the command first, then dispatch
    parsed = parse_args(args, s; as_symbols = true)
    parsed[:command]
end

function cli_init(args)
    s = ArgParseSettings(prog = "himalaya init")
    @add_arg_table! s begin
        "experiment_path"
            help     = "path to experiment directory"
            required = true
        "--manifest", "-m"
            help     = "path to manifest CSV"
            default  = nothing
        "--beamline"
            help     = "beamline profile name (default: default)"
            default  = "default"
        "--name"
            help     = "experiment name"
            default  = nothing
    end
    p = parse_args(args, s; as_symbols = true)

    exp_path     = p[:experiment_path]
    manifest_path = p[:manifest]

    # Resolve default beamline paths
    data_dir     = joinpath(exp_path, "data")
    analysis_dir = joinpath(exp_path, "analysis", "automatic_analysis")

    db     = open_db(exp_path)
    exp_id = init_experiment!(db;
        name          = something(p[:name], basename(exp_path)),
        path          = exp_path,
        data_dir      = data_dir,
        analysis_dir  = analysis_dir,
        manifest_path = manifest_path)

    if manifest_path !== nothing && isfile(manifest_path)
        samples = parse_manifest(manifest_path)
        for ms in samples
            s_id = create_sample!(db;
                experiment_id  = exp_id,
                label          = ms.label,
                name           = ms.name,
                notes          = ms.notes_sample)
            for filename in ms.filenames
                create_exposure!(db; sample_id = s_id, filename = filename)
            end
            if !isempty(ms.notes_exposure)
                # store exposure note as a tag on the first exposure
                e_id = Int(only(DBInterface.execute(db,
                    "SELECT id FROM exposures WHERE sample_id = ? LIMIT 1",
                    [s_id]))[1])
                DBInterface.execute(db,
                    "INSERT INTO exposure_tags (exposure_id, key, value, source)
                     VALUES (?, 'note', ?, 'manifest')",
                    [e_id, ms.notes_exposure])
            end
        end
        println("Imported $(length(samples)) samples from manifest.")
    end

    println("Initialized experiment #$exp_id at $exp_path")
end

function cli_analyze(args)
    s = ArgParseSettings(prog = "himalaya analyze")
    @add_arg_table! s begin
        "experiment_path"
            required = true
        "--sample", "-s"
            help    = "analyze only this sample label (e.g. D1)"
            default = nothing
    end
    p = parse_args(args, s; as_symbols = true)

    db  = open_db(p[:experiment_path])
    exp = get_experiment(db, 1)

    sample_filter = p[:sample]
    samples       = get_samples(db, 1)
    sample_filter !== nothing && filter!(s -> s[:label] == sample_filter, samples)

    for sample in samples
        exposures = get_exposures(db, sample[:id])
        for exp_row in exposures
            e_id = exp_row[:id]
            print("  Analyzing $(sample[:label]) / $(exp_row[:filename]) ... ")
            try
                analyze_exposure!(db, e_id, exp[:analysis_dir])
                println("done")
            catch e
                println("SKIP ($(e.msg))")
            end
        end
    end
end

function cli_show(args)
    s = ArgParseSettings(prog = "himalaya show")
    @add_arg_table! s begin
        "experiment_path"
            required = true
        "--sample", "-s"
            help     = "sample label"
            required = true
    end
    p = parse_args(args, s; as_symbols = true)

    db      = open_db(p[:experiment_path])
    samples = get_samples(db, 1)
    sample  = findfirst(r -> r[:label] == p[:sample], samples)
    sample === nothing && error("sample $(p[:sample]) not found")
    sample_row = samples[sample]

    exposures = get_exposures(db, sample_row[:id])
    for exp_row in exposures
        e_id    = exp_row[:id]
        pks     = get_peaks_for_exposure(db, e_id)
        idxs    = get_indices_for_exposure(db, e_id)
        groups  = get_groups_for_exposure(db, e_id)

        println("\nExposure: $(exp_row[:filename])")
        println("  Peaks ($(length(pks))):")
        for pk in pks
            @printf "    q=%.4f  prom=%.3f  sharp=%.3f  [%s]\n" pk[:q] pk[:prominence] pk[:sharpness] pk[:source]
        end
        println("  Indices ($(length(idxs))):")
        for idx in idxs
            @printf "    %-6s  basis=%.4f  score=%.3f  R²=%.4f  d=%.2f\n" idx[:phase] idx[:basis] idx[:score] idx[:r_squared] idx[:lattice_d]
        end
        active = findfirst(g -> g[:active] == 1, groups)
        if active !== nothing
            println("  Active group: $(groups[active][:kind])")
        end
    end
end

function main(args = ARGS)
    isempty(args) && (println("Usage: himalaya <command> [args]"); return)
    cmd  = popfirst!(args)

    if cmd == "init"
        cli_init(args)
    elseif cmd == "analyze"
        cli_analyze(args)
    elseif cmd == "show"
        cli_show(args)
    else
        println("Unknown command: $cmd. Available: init, analyze, show")
    end
end
```

Add `using Printf` to the top of `cli.jl`.

- [ ] **Step 2: Update module entry to export `main`**

`packages/HimalayaUI/src/HimalayaUI.jl`:
```julia
module HimalayaUI

include("db.jl")
include("datfile.jl")
include("manifest.jl")
include("pipeline.jl")
include("cli.jl")

export main

end
```

- [ ] **Step 3: Smoke test the CLI end-to-end**

```bash
# Create a temp experiment directory
mkdir -p /tmp/himalaya_test/data
mkdir -p /tmp/himalaya_test/analysis/automatic_analysis
cp test/data/example_tot.dat /tmp/himalaya_test/analysis/automatic_analysis/example_tot.dat

# Init
julia --project=packages/HimalayaUI -e '
using HimalayaUI
main(["init", "/tmp/himalaya_test", "--name", "TestRun"])'
```

Expected output:
```
Initialized experiment #1 at /tmp/himalaya_test
```

```bash
# Manually add a sample + exposure via Julia REPL (no manifest in this smoke test)
julia --project=packages/HimalayaUI -e '
using HimalayaUI, SQLite
db = open_db("/tmp/himalaya_test")
s_id = create_sample!(db; experiment_id=1, label="D1", name="UX1")
create_exposure!(db; sample_id=s_id, filename="example_tot")'

# Analyze
julia --project=packages/HimalayaUI -e '
using HimalayaUI
main(["analyze", "/tmp/himalaya_test"])'
```

Expected:
```
  Analyzing D1 / example_tot ... done
```

```bash
# Show
julia --project=packages/HimalayaUI -e '
using HimalayaUI
main(["show", "/tmp/himalaya_test", "--sample", "D1"])'
```

Expected: printed table of peaks and indices for D1.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/src/cli.jl packages/HimalayaUI/src/HimalayaUI.jl
git commit -m "feat(HimalayaUI): CLI — himalaya init / analyze / show"
```

---

### Task 9: Full test suite pass

- [ ] **Step 1: Run the complete HimalayaUI test suite**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: all testsets pass with zero failures.

- [ ] **Step 2: Run the root Himalaya test suite to confirm nothing regressed**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: all pass.

- [ ] **Step 3: Commit if any fixes were needed**

```bash
git add -p
git commit -m "fix(HimalayaUI): full test suite fixes"
```

---

## Self-Review

**Spec coverage check:**

| Spec section | Covered by task |
|---|---|
| Monorepo / HimalayaUI sub-package | Task 1 |
| SQLite schema (all tables) | Task 2 |
| `.dat` file parser | Task 3 |
| Manifest CSV parser, filename range expansion | Task 4 |
| Auto-group algorithm (non-overlapping, greedy) | Task 5 |
| Persist peaks, indices, index_peaks, index_groups | Task 6 |
| `analyze_exposure!`, `init_experiment!` | Task 7 |
| CLI: `himalaya init`, `analyze`, `show` | Task 8 |
| Beamline config (default paths) | Covered in `cli_init` defaults |
| `.dat` format (whitespace-separated, no header) | Task 3 |
| `manifest_path` stored in `experiments` table | Task 2 schema + Task 8 |
| Manual peaks preserved across re-analysis | Task 6 `persist_analysis!` |
| Tags from manifest stored with `source='manifest'` | Task 8 `cli_init` |

**Not in this plan (correct — deferred to later plans):**
- REST API (Plan 2)
- Frontend (Plans 3–6)
- `himalaya serve` and `himalaya export` CLI commands (Plan 2)
- User identification (Plan 2)
- Index group modification via API (Plan 2)
