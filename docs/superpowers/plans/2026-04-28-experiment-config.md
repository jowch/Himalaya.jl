# Experiment Config System Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace all hardcoded file I/O conventions in HimalayaUI with a config-driven system driven by a per-experiment `experiment.toml` file, enabling support for multiple beamlines and experiment types without code changes.

**Architecture:** A new `config.jl` module owns an `ExperimentConfig` struct, TOML loading, filesystem-based file resolution, and a config-driven manifest parser. The config TOML blob is stored in the `experiments` DB row so the system is self-contained; `experiment.toml` on disk is the source of truth and can be re-ingested at any time. Experiment directories are read-only at runtime — Himalaya only reads from them, never writes.

**Tech Stack:** Julia 1.9+, TOML stdlib, SQLite.jl, ArgParse.jl, existing HimalayaUI conventions (`@testset` / `@test`, Oxygen.jl routes)

---

## File Map

| Action | Path | Responsibility |
|--------|------|----------------|
| Create | `packages/HimalayaUI/src/config.jl` | `ExperimentConfig` struct, `load_config`, `config_from_db`, `resolve_files`, `parse_manifest` |
| Create | `packages/HimalayaUI/configs/simple.toml` | Built-in template encoding current hardcoded behavior |
| Create | `packages/HimalayaUI/test/test_config.jl` | Unit tests for config.jl |
| Modify | `packages/HimalayaUI/Project.toml` | Add `TOML` stdlib dep |
| Modify | `packages/HimalayaUI/src/db.jl` | Add `config`, `experiment_type`, `energy_kev`, `flight_path_m` columns; update `create_experiment` |
| Modify | `packages/HimalayaUI/src/manifest.jl` | Delegate to `config.jl`; keep `ManifestSample` struct |
| Modify | `packages/HimalayaUI/src/pipeline.jl` | Use `config_from_db` for path resolution in `analyze_exposure!` |
| Modify | `packages/HimalayaUI/src/cli.jl` | Add `config` subcommand; update `init` to read `experiment.toml`; add `reingest`; remove `--beamline` |
| Modify | `packages/HimalayaUI/src/server.jl` | Add `POST /api/experiments/:id/reingest` route |
| Modify | `packages/HimalayaUI/test/runtests.jl` | Include `test_config.jl` |

---

## Task 1: Add TOML dep + `ExperimentConfig` struct + `load_config`

**Files:**
- Create: `packages/HimalayaUI/src/config.jl`
- Create: `packages/HimalayaUI/test/test_config.jl`
- Modify: `packages/HimalayaUI/Project.toml`
- Modify: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Add `TOML` to Project.toml**

In `packages/HimalayaUI/Project.toml`, add to `[deps]`:
```toml
TOML = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
```
And to `[compat]`:
```toml
TOML = "1"
```

- [ ] **Step 2: Write the failing test**

Create `packages/HimalayaUI/test/test_config.jl`:
```julia
using Test, HimalayaUI
import HimalayaUI: load_config, ExperimentConfig

@testset "load_config" begin
    mktempdir() do dir
        toml = joinpath(dir, "experiment.toml")
        write(toml, """
        [experiment]
        name        = "TestRun/ExpA"
        description = "test"
        manifest    = "manifest.csv"

        [beamline]
        energy_kev    = 12.0
        flight_path_m = 2.5

        [manifest]
        delimiter      = "\\t"
        skip_rows      = 1
        header_row     = 0
        sample_id      = 1
        label          = 2
        name           = 3
        filenames      = 9
        notes_sample   = 10
        notes_exposure = 11

        [layout]
        data_dir      = "data"
        analysis_dir  = "analysis/automatic_analysis"
        exposure_type = "simple"

        [files]
        integration = "{name}.dat"
        image       = "{name}.tiff"
        """)
        cfg = load_config(toml)
        @test cfg.name == "TestRun/ExpA"
        @test cfg.energy_kev == 12.0
        @test cfg.flight_path_m == 2.5
        @test cfg.manifest_file == "manifest.csv"
        @test cfg.delimiter == "\t"
        @test cfg.skip_rows == 1
        @test cfg.header_row == 0
        @test cfg.col_sample_id == 1
        @test cfg.col_label == 2
        @test cfg.col_name == 3
        @test cfg.col_filenames == 9
        @test cfg.col_notes_sample == 10
        @test cfg.col_notes_exposure == 11
        @test cfg.data_dir == "data"
        @test cfg.analysis_dir == "analysis/automatic_analysis"
        @test cfg.exposure_type == "simple"
        @test cfg.integration_pattern == "{name}.dat"
        @test cfg.image_pattern == "{name}.tiff"
    end
end

@testset "load_config validates path traversal" begin
    mktempdir() do dir
        toml = joinpath(dir, "experiment.toml")
        write(toml, """
        [experiment]
        name = "X" 
        description = ""
        manifest = "manifest.csv"
        [beamline]
        energy_kev = 0.0
        flight_path_m = 0.0
        [manifest]
        delimiter = "\\t"
        skip_rows = 0
        header_row = 0
        sample_id = 1
        label = 2
        name = 3
        filenames = 9
        notes_sample = 10
        notes_exposure = 11
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "../{name}.dat"
        image = "{name}.tiff"
        """)
        @test_throws ErrorException load_config(toml)
    end
end
```

- [ ] **Step 3: Run test to confirm it fails**

```bash
cd /path/to/Himalaya.jl
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: `UndefVarError: ExperimentConfig not defined` or similar.

- [ ] **Step 4: Create `packages/HimalayaUI/src/config.jl`**

```julia
module ConfigModule

using TOML

struct ExperimentConfig
    # [experiment]
    name               ::String
    description        ::String
    manifest_file      ::String        # relative path to manifest CSV
    # [beamline]
    energy_kev         ::Union{Float64,Nothing}
    flight_path_m      ::Union{Float64,Nothing}
    # [manifest]
    delimiter          ::String
    skip_rows          ::Int
    header_row         ::Int
    col_sample_id      ::Union{Int,String}
    col_label          ::Union{Int,String}
    col_name           ::Union{Int,String}
    col_filenames      ::Union{Int,String}
    col_notes_sample   ::Union{Int,String}
    col_notes_exposure ::Union{Int,String}
    # [layout]
    data_dir           ::String
    analysis_dir       ::String
    exposure_type      ::String
    # [files]
    integration_pattern::String
    image_pattern      ::String
end

function _validate_pattern(pattern::String, field::String)
    startswith(pattern, "/") && error("$field must not be an absolute path: $pattern")
    contains(pattern, "..") && error("$field must not traverse upward (..): $pattern")
    contains(pattern, "{name}") || error("$field must contain {name}: $pattern")
end

function load_config(path::AbstractString)::ExperimentConfig
    isfile(path) || error("experiment.toml not found: $path")
    d = TOML.parsefile(path)

    exp    = get(d, "experiment", Dict())
    bl     = get(d, "beamline",   Dict())
    mf     = get(d, "manifest",   Dict())
    layout = get(d, "layout",     Dict())
    files  = get(d, "files",      Dict())

    integration = get(files, "integration", "{name}.dat")
    image       = get(files, "image",       "{name}.tiff")
    _validate_pattern(integration, "files.integration")
    _validate_pattern(image,       "files.image")

    ExperimentConfig(
        get(exp, "name",        ""),
        get(exp, "description", ""),
        get(exp, "manifest",    "manifest.csv"),
        get(bl,  "energy_kev",    nothing),
        get(bl,  "flight_path_m", nothing),
        get(mf,  "delimiter",      "\t"),
        get(mf,  "skip_rows",      1),
        get(mf,  "header_row",     0),
        get(mf,  "sample_id",      1),
        get(mf,  "label",          2),
        get(mf,  "name",           3),
        get(mf,  "filenames",      9),
        get(mf,  "notes_sample",   10),
        get(mf,  "notes_exposure", 11),
        get(layout, "data_dir",      "data"),
        get(layout, "analysis_dir",  "analysis/automatic_analysis"),
        get(layout, "exposure_type", "simple"),
        integration,
        image,
    )
end

end  # module
```

- [ ] **Step 5: Export from HimalayaUI module**

In `packages/HimalayaUI/src/HimalayaUI.jl`, add:
```julia
include("config.jl")
using .ConfigModule: ExperimentConfig, load_config
export ExperimentConfig, load_config
```
(Add the `include` and `using`/`export` lines alongside the existing includes.)

- [ ] **Step 6: Add test file to runtests.jl**

In `packages/HimalayaUI/test/runtests.jl`, add:
```julia
include("test_config.jl")
```

- [ ] **Step 7: Run tests to confirm they pass**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: `test_config.jl` tests pass, all existing tests still pass.

- [ ] **Step 8: Commit**

```bash
git add packages/HimalayaUI/src/config.jl \
        packages/HimalayaUI/test/test_config.jl \
        packages/HimalayaUI/test/runtests.jl \
        packages/HimalayaUI/Project.toml \
        packages/HimalayaUI/src/HimalayaUI.jl
git commit -m "feat(config): ExperimentConfig struct + load_config"
```

---

## Task 2: Built-in `simple.toml` template

**Files:**
- Create: `packages/HimalayaUI/configs/simple.toml`

- [ ] **Step 1: Create the configs directory and simple.toml**

Create `packages/HimalayaUI/configs/simple.toml`:
```toml
[experiment]
name        = ""
description = ""
manifest    = "manifest.csv"

[beamline]
energy_kev    = 0.0
flight_path_m = 0.0

[manifest]
delimiter      = "\t"
skip_rows      = 1
header_row     = 0
sample_id      = 1
label          = 2
name           = 3
filenames      = 9
notes_sample   = 10
notes_exposure = 11

[layout]
data_dir      = "data"
analysis_dir  = "analysis/automatic_analysis"
exposure_type = "simple"

[files]
integration = "{name}.dat"
image       = "{name}.tiff"
```

- [ ] **Step 2: Add `configs_dir()` helper to config.jl**

In `packages/HimalayaUI/src/config.jl`, add inside `module ConfigModule`:
```julia
configs_dir() = joinpath(@__DIR__, "..", "configs")

function list_config_types()::Vector{String}
    dir = configs_dir()
    isdir(dir) || return String[]
    [splitext(f)[1] for f in readdir(dir) if endswith(f, ".toml")]
end

function load_builtin_config(type_name::AbstractString)::ExperimentConfig
    path = joinpath(configs_dir(), type_name * ".toml")
    isfile(path) || error("Unknown config type '$type_name'. Available: $(join(list_config_types(), ", "))")
    load_config(path)
end
```

Export from `HimalayaUI.jl`: `list_config_types, load_builtin_config`

- [ ] **Step 3: Write test for builtin loading**

In `packages/HimalayaUI/test/test_config.jl`, add:
```julia
@testset "load_builtin_config simple" begin
    cfg = HimalayaUI.load_builtin_config("simple")
    @test cfg.delimiter == "\t"
    @test cfg.col_sample_id == 1
    @test cfg.data_dir == "data"
    @test cfg.analysis_dir == "analysis/automatic_analysis"
    @test cfg.integration_pattern == "{name}.dat"
    @test cfg.image_pattern == "{name}.tiff"
end

@testset "list_config_types" begin
    types = HimalayaUI.list_config_types()
    @test "simple" in types
end
```

- [ ] **Step 4: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/configs/simple.toml \
        packages/HimalayaUI/src/config.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/test/test_config.jl
git commit -m "feat(config): simple.toml built-in template + list/load helpers"
```

---

## Task 3: `resolve_files` — prefix-based filesystem scan

**Files:**
- Modify: `packages/HimalayaUI/src/config.jl`
- Modify: `packages/HimalayaUI/test/test_config.jl`

- [ ] **Step 1: Write the failing tests**

Add to `packages/HimalayaUI/test/test_config.jl`:
```julia
@testset "resolve_files" begin
    mktempdir() do dir
        # Create fake .dat files
        for name in ["JC001", "JC002", "JC003", "AB001"]
            write(joinpath(dir, name * ".dat"), "")
        end

        cfg = HimalayaUI.load_builtin_config("simple")

        # Single stem prefix: JC001 → finds JC001.dat only (exact prefix)
        @test HimalayaUI.resolve_files(cfg, dir, "JC001", cfg.integration_pattern) == ["JC001"]

        # Range prefix: JC002 → finds JC002.dat
        @test HimalayaUI.resolve_files(cfg, dir, "JC002", cfg.integration_pattern) == ["JC002"]

        # Broad prefix: JC → finds all JC*.dat sorted
        @test HimalayaUI.resolve_files(cfg, dir, "JC", cfg.integration_pattern) == ["JC001", "JC002", "JC003"]

        # Non-existent prefix → empty
        @test HimalayaUI.resolve_files(cfg, dir, "ZZ", cfg.integration_pattern) == String[]

        # Non-existent directory → empty
        @test HimalayaUI.resolve_files(cfg, "/no/such/dir", "JC", cfg.integration_pattern) == String[]
    end
end

@testset "resolve_files with subdir pattern" begin
    mktempdir() do dir
        subdir = joinpath(dir, "integrated")
        mkpath(subdir)
        write(joinpath(subdir, "SA001.dat"), "")
        write(joinpath(subdir, "SA002.dat"), "")

        cfg = HimalayaUI.load_builtin_config("simple")
        results = HimalayaUI.resolve_files(cfg, dir, "SA", "integrated/{name}.dat")
        @test results == ["SA001", "SA002"]
    end
end
```

- [ ] **Step 2: Run test to confirm failure**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: `UndefVarError: resolve_files not defined`.

- [ ] **Step 3: Implement `resolve_files` in config.jl**

Add inside `module ConfigModule` in `config.jl`:
```julia
function resolve_files(
    ::ExperimentConfig,
    base_dir::AbstractString,
    prefix::AbstractString,
    file_pattern::String,
)::Vector{String}
    parts = split(file_pattern, "{name}"; limit=2)
    length(parts) == 2 || error("file pattern must contain exactly one {name}: $file_pattern")
    before, after = String.(parts)

    scan_subdir    = dirname(before)
    file_prefix_local = basename(before) * prefix
    scan_dir = isempty(scan_subdir) ? base_dir : joinpath(base_dir, scan_subdir)

    isdir(scan_dir) || return String[]

    matches = filter(readdir(scan_dir)) do f
        startswith(f, file_prefix_local) && endswith(f, after)
    end
    sort!(matches)
    # Strip suffix to return bare stems
    [m[1:end-length(after)] for m in matches]
end
```

Export `resolve_files` from `HimalayaUI.jl`.

- [ ] **Step 4: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/config.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/test/test_config.jl
git commit -m "feat(config): resolve_files prefix-based filesystem scan"
```

---

## Task 4: Config-driven `parse_manifest`

**Files:**
- Modify: `packages/HimalayaUI/src/config.jl`
- Modify: `packages/HimalayaUI/src/manifest.jl`
- Modify: `packages/HimalayaUI/test/test_config.jl`

The current `manifest.jl` exports `ManifestSample` and `parse_manifest(path)`. The new version adds an overload `parse_manifest(config, path)` and the old signature becomes a thin wrapper using `load_builtin_config("simple")` for backward compatibility.

- [ ] **Step 1: Understand the existing `ManifestSample` struct**

Read `packages/HimalayaUI/src/manifest.jl`. Note that `ManifestSample` has fields:
`label::String`, `name::String`, `notes_sample::String`, `notes_exposure::String`, `filenames::Vector{String}`

Do not change this struct — `config.jl` is a translation layer that produces it.

- [ ] **Step 2: Write failing tests**

Add to `packages/HimalayaUI/test/test_config.jl`:
```julia
@testset "parse_manifest positional columns" begin
    mktempdir() do dir
        csv = joinpath(dir, "manifest.csv")
        # Tab-separated, 1 skip row, no header row
        # Columns: 1=id, 2=label, 3=name, 4-8=ignored, 9=filenames, 10=notes_sample, 11=notes_exposure
        write(csv, join([
            "header row to skip",
            "1\tD1\tUX1\t\t\t\t\t\tJC001-003\tnote_s\tnote_e",
            "2\tD2\tUX2\t\t\t\t\t\tJC004\t\t",
            "not a number\tskip\tme\t\t\t\t\t\tJC999\t\t",
        ], "\n"))

        cfg = HimalayaUI.load_builtin_config("simple")
        samples = HimalayaUI.parse_manifest(cfg, csv)

        @test length(samples) == 2
        @test samples[1].label == "D1"
        @test samples[1].name  == "UX1"
        @test samples[1].filenames == ["JC001", "JC002", "JC003"]
        @test samples[1].notes_sample   == "note_s"
        @test samples[1].notes_exposure == "note_e"
        @test samples[2].filenames == ["JC004"]
        @test samples[2].notes_sample   == ""
        @test samples[2].notes_exposure == ""
    end
end

@testset "parse_manifest named columns" begin
    mktempdir() do dir
        csv = joinpath(dir, "manifest.csv")
        # header_row=1, skip_rows=0, comma-separated, named columns
        write(csv, join([
            "#,Sample,Name,Type,Time,x,y,z,Filename(s),Notes (Sample),Notes (Exposure)",
            "1,D1,UX1,,,,,,AB001-002,s_note,e_note",
            "2,D2,UX2,,,,,,AB003,,",
        ], "\n"))

        cfg_toml = """
        [experiment]
        name = "X"
        description = ""
        manifest = "manifest.csv"
        [beamline]
        energy_kev = 0.0
        flight_path_m = 0.0
        [manifest]
        delimiter = ","
        skip_rows = 0
        header_row = 1
        sample_id = "#"
        label = "Sample"
        name = "Name"
        filenames = "Filename(s)"
        notes_sample = "Notes (Sample)"
        notes_exposure = "Notes (Exposure)"
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "{name}.dat"
        image = "{name}.tiff"
        """
        toml_path = joinpath(dir, "experiment.toml")
        write(toml_path, cfg_toml)
        cfg = HimalayaUI.load_config(toml_path)

        samples = HimalayaUI.parse_manifest(cfg, csv)
        @test length(samples) == 2
        @test samples[1].label    == "D1"
        @test samples[1].filenames == ["AB001", "AB002"]
        @test samples[1].notes_sample   == "s_note"
        @test samples[1].notes_exposure == "e_note"
        @test samples[2].filenames == ["AB003"]
    end
end
```

- [ ] **Step 3: Run to confirm failure**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: `UndefVarError: parse_manifest not defined` (config overload).

- [ ] **Step 4: Implement config-driven `parse_manifest` in config.jl**

Add inside `module ConfigModule` in `config.jl`:
```julia
using ..HimalayaUICore: ManifestSample   # adjust to actual module path if different

function _resolve_col(col::Union{Int,String}, header_map::Dict{String,Int}, field::String)::Int
    if col isa String
        idx = get(header_map, col, nothing)
        idx === nothing && @warn "Column '$col' not found in header; field $field will be empty"
        return idx === nothing ? 0 : idx
    end
    return col
end

function _expand_range(entry::AbstractString)::Vector{String}
    # "JC001-004" → ["JC001","JC002","JC003","JC004"]
    # "JC013-JC016" → ["JC013","JC014","JC015","JC016"]
    m = match(r"^([A-Za-z]+)(\d+)-(?:[A-Za-z]*)(\d+)$", entry)
    if m !== nothing
        prefix, s_start, s_end = m.captures
        n_start, n_end = parse(Int, s_start), parse(Int, s_end)
        width = length(s_start)
        return [prefix * lpad(string(i), width, '0') for i in n_start:n_end]
    end
    return [entry]
end

function parse_manifest(cfg::ExperimentConfig, path::AbstractString)::Vector{ManifestSample}
    lines = readlines(path)
    isempty(lines) && return ManifestSample[]

    # Skip leading rows
    lines = lines[cfg.skip_rows+1:end]

    # Build header map if header_row specified
    header_map = Dict{String,Int}()
    data_start = 1
    if cfg.header_row > 0
        # header_row is 1-based relative to the *original* file;
        # after skipping skip_rows, it's at position (header_row - skip_rows)
        hdr_idx = cfg.header_row - cfg.skip_rows
        if hdr_idx >= 1 && hdr_idx <= length(lines)
            cols = split(lines[hdr_idx], cfg.delimiter)
            for (i, c) in enumerate(cols)
                header_map[strip(c)] = i
            end
            data_start = hdr_idx + 1
        end
    end

    samples = ManifestSample[]
    for line in lines[data_start:end]
        cols = split(line, cfg.delimiter)

        # Resolve column indices
        function getcol(col)
            idx = _resolve_col(col, header_map, string(col))
            idx == 0 && return ""
            idx <= length(cols) ? strip(cols[idx]) : ""
        end

        id_str = getcol(cfg.col_sample_id)
        tryparse(Int, id_str) === nothing && continue

        filenames_str = getcol(cfg.col_filenames)
        isempty(filenames_str) && continue

        filenames = _expand_range(filenames_str)

        push!(samples, ManifestSample(
            getcol(cfg.col_label),
            getcol(cfg.col_name),
            getcol(cfg.col_notes_sample),
            getcol(cfg.col_notes_exposure),
            filenames,
        ))
    end
    samples
end
```

> **Note on `ManifestSample`:** Both `config.jl` and `manifest.jl` are included into the same top-level `HimalayaUI` module via `HimalayaUI.jl`. As long as `manifest.jl` is included before `config.jl`, `ManifestSample` is already in scope — no import needed. Verify the include order in `HimalayaUI.jl` and ensure `manifest.jl` comes first.

Export `parse_manifest` from `HimalayaUI.jl`.

- [ ] **Step 5: Update `manifest.jl` to delegate for backward compatibility**

In `packages/HimalayaUI/src/manifest.jl`, replace the body of the no-config `parse_manifest(path)` (keeping `ManifestSample` struct unchanged):
```julia
# Backward-compat wrapper — uses simple.toml defaults
function parse_manifest(path::AbstractString)::Vector{ManifestSample}
    cfg = load_builtin_config("simple")
    parse_manifest(cfg, path)
end
```

- [ ] **Step 6: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass, including existing `test_manifest.jl` (which uses the no-config overload).

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/src/config.jl \
        packages/HimalayaUI/src/manifest.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/test/test_config.jl
git commit -m "feat(config): config-driven parse_manifest with positional+named column resolution"
```

---

## Task 5: DB schema migration — add config columns to experiments

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl`
- Modify: `packages/HimalayaUI/test/test_db.jl`

- [ ] **Step 1: Write failing test**

In `packages/HimalayaUI/test/test_db.jl`, add to the existing `@testset`:
```julia
@testset "experiments table has config columns" begin
    db = HimalayaUI.open_db(":memory:")
    row = Tables.rowtable(DBInterface.execute(db,
        "PRAGMA table_info(experiments)"))[1]
    col_names = [r.name for r in Tables.rowtable(DBInterface.execute(db,
        "PRAGMA table_info(experiments)"))]
    @test "config"          in col_names
    @test "experiment_type" in col_names
    @test "energy_kev"      in col_names
    @test "flight_path_m"   in col_names
end
```

- [ ] **Step 2: Run to confirm failure**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: assertion failure — columns don't exist yet.

- [ ] **Step 3: Update schema in db.jl**

In `packages/HimalayaUI/src/db.jl`, find the `CREATE TABLE IF NOT EXISTS experiments` statement and add the four new columns:
```sql
CREATE TABLE IF NOT EXISTS experiments (
    id             INTEGER PRIMARY KEY,
    name           TEXT,
    path           TEXT NOT NULL,
    data_dir       TEXT NOT NULL,
    analysis_dir   TEXT NOT NULL,
    manifest_path  TEXT,
    config         TEXT,
    experiment_type TEXT,
    energy_kev     REAL,
    flight_path_m  REAL,
    created_at     DATETIME DEFAULT CURRENT_TIMESTAMP
);
```

Also add a migration block in `open_db` (or wherever schema setup runs) to handle existing DBs:
```julia
DBInterface.execute(db, "ALTER TABLE experiments ADD COLUMN config TEXT")
DBInterface.execute(db, "ALTER TABLE experiments ADD COLUMN experiment_type TEXT")
DBInterface.execute(db, "ALTER TABLE experiments ADD COLUMN energy_kev REAL")
DBInterface.execute(db, "ALTER TABLE experiments ADD COLUMN flight_path_m REAL")
```
Wrap each `ALTER TABLE` in a try/catch so it silently skips if the column already exists (SQLite returns an error on duplicate column add, not a no-op):
```julia
for stmt in [
    "ALTER TABLE experiments ADD COLUMN config TEXT",
    "ALTER TABLE experiments ADD COLUMN experiment_type TEXT",
    "ALTER TABLE experiments ADD COLUMN energy_kev REAL",
    "ALTER TABLE experiments ADD COLUMN flight_path_m REAL",
]
    try DBInterface.execute(db, stmt) catch end
end
```

- [ ] **Step 4: Update `create_experiment` in db.jl**

Find the `create_experiment` function and add the new parameters:
```julia
function create_experiment(db::SQLite.DB;
    name          ::String,
    path          ::String,
    data_dir      ::String,
    analysis_dir  ::String,
    manifest_path ::Union{String,Nothing} = nothing,
    config        ::Union{String,Nothing} = nothing,
    experiment_type::Union{String,Nothing} = nothing,
    energy_kev    ::Union{Float64,Nothing} = nothing,
    flight_path_m ::Union{Float64,Nothing} = nothing,
)::Int
    res = DBInterface.execute(db,
        """INSERT INTO experiments
           (name, path, data_dir, analysis_dir, manifest_path,
            config, experiment_type, energy_kev, flight_path_m)
           VALUES (?,?,?,?,?,?,?,?,?)""",
        [name, path, data_dir, analysis_dir, manifest_path,
         config, experiment_type, energy_kev, flight_path_m])
    Int(DBInterface.lastrowid(res))
end
```

- [ ] **Step 5: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/db.jl \
        packages/HimalayaUI/test/test_db.jl
git commit -m "feat(db): add config, experiment_type, energy_kev, flight_path_m to experiments"
```

---

## Task 6: `config_from_db` — store and retrieve config blob

**Files:**
- Modify: `packages/HimalayaUI/src/config.jl`
- Modify: `packages/HimalayaUI/test/test_config.jl`

- [ ] **Step 1: Write failing test**

Add to `packages/HimalayaUI/test/test_config.jl`:
```julia
@testset "config_from_db roundtrip" begin
    import SQLite, DBInterface
    db = HimalayaUI.open_db(":memory:")
    cfg_orig = HimalayaUI.load_builtin_config("simple")
    blob = HimalayaUI.config_to_toml(cfg_orig)
    @test blob isa String
    @test contains(blob, "[experiment]")

    exp_id = HimalayaUI.create_experiment(db;
        name="Test/Exp", path="/tmp", data_dir="data",
        analysis_dir="analysis/automatic_analysis",
        config=blob, experiment_type="simple")

    cfg2 = HimalayaUI.config_from_db(db, exp_id)
    @test cfg2.data_dir == cfg_orig.data_dir
    @test cfg2.integration_pattern == cfg_orig.integration_pattern
    @test cfg2.delimiter == cfg_orig.delimiter
end

@testset "config_from_db falls back to simple when config is NULL" begin
    db = HimalayaUI.open_db(":memory:")
    exp_id = HimalayaUI.create_experiment(db;
        name="Legacy", path="/tmp", data_dir="data",
        analysis_dir="analysis/automatic_analysis")
    cfg = HimalayaUI.config_from_db(db, exp_id)
    @test cfg.data_dir == "data"
    @test cfg.integration_pattern == "{name}.dat"
end
```

- [ ] **Step 2: Run to confirm failure**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: `UndefVarError: config_from_db not defined`.

- [ ] **Step 3: Implement `config_to_toml` and `config_from_db` in config.jl**

Add inside `module ConfigModule` — `SQLite`, `DBInterface`, and `Tables` are already available in the parent HimalayaUI module, no additional `using` needed inside the submodule:
```julia
function config_to_toml(cfg::ExperimentConfig)::String
    # Serialize ExperimentConfig back to TOML text for DB storage
    col_or_str(v) = v isa String ? "\"$v\"" : string(v)
    """
    [experiment]
    name        = "$(cfg.name)"
    description = "$(cfg.description)"
    manifest    = "$(cfg.manifest_file)"

    [beamline]
    energy_kev    = $(something(cfg.energy_kev, 0.0))
    flight_path_m = $(something(cfg.flight_path_m, 0.0))

    [manifest]
    delimiter      = "$(escape_string(cfg.delimiter))"
    skip_rows      = $(cfg.skip_rows)
    header_row     = $(cfg.header_row)
    sample_id      = $(col_or_str(cfg.col_sample_id))
    label          = $(col_or_str(cfg.col_label))
    name           = $(col_or_str(cfg.col_name))
    filenames      = $(col_or_str(cfg.col_filenames))
    notes_sample   = $(col_or_str(cfg.col_notes_sample))
    notes_exposure = $(col_or_str(cfg.col_notes_exposure))

    [layout]
    data_dir      = "$(cfg.data_dir)"
    analysis_dir  = "$(cfg.analysis_dir)"
    exposure_type = "$(cfg.exposure_type)"

    [files]
    integration = "$(cfg.integration_pattern)"
    image       = "$(cfg.image_pattern)"
    """
end

function config_from_db(db::SQLite.DB, experiment_id::Int)::ExperimentConfig
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT config FROM experiments WHERE id = ?", [experiment_id]))
    isempty(rows) && error("Experiment $experiment_id not found")
    blob = rows[1].config
    if blob === nothing || blob === missing
        return load_builtin_config("simple")
    end
    mktempdir() do dir
        path = joinpath(dir, "experiment.toml")
        write(path, blob)
        load_config(path)
    end
end
```

Export `config_to_toml, config_from_db` from `HimalayaUI.jl`.

- [ ] **Step 4: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/config.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/test/test_config.jl
git commit -m "feat(config): config_to_toml + config_from_db with NULL fallback to simple"
```

---

## Task 7: Update `analyze_exposure!` to use config for path resolution

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl`
- Modify: `packages/HimalayaUI/test/test_pipeline.jl`

Currently `analyze_exposure!(db, exposure_id, analysis_dir)` constructs `dat_path = joinpath(analysis_dir, filename * ".dat")`. The new version reads `analysis_dir` and the integration pattern from the config stored in the DB.

- [ ] **Step 1: Write failing test**

Read `packages/HimalayaUI/test/test_pipeline.jl` to understand its setup. Then add a test that verifies the config pattern is used:
```julia
@testset "analyze_exposure! uses config integration pattern" begin
    db = HimalayaUI.open_db(":memory:")
    # Create experiment with a non-default integration pattern
    cfg = HimalayaUI.load_builtin_config("simple")
    blob = HimalayaUI.config_to_toml(cfg)

    mktempdir() do dir
        analysis_dir = joinpath(dir, "analysis", "automatic_analysis")
        mkpath(analysis_dir)
        write(joinpath(analysis_dir, "EX001.dat"),
            "0.01 100.0 5.0\n0.02 90.0 4.5\n0.03 80.0 4.0\n")

        exp_id = HimalayaUI.create_experiment(db;
            name="T/E", path=dir, data_dir=joinpath(dir, "data"),
            analysis_dir=analysis_dir, config=blob, experiment_type="simple")
        sample_id = HimalayaUI.create_sample(db;
            experiment_id=exp_id, label="S1", name="S1")
        exp_row_id = HimalayaUI.create_exposure(db;
            sample_id=sample_id, filename="EX001")

        # Should not throw — resolves via config
        HimalayaUI.analyze_exposure!(db, exp_row_id)
        peaks = HimalayaUI.get_peaks(db, exp_row_id)
        @test peaks isa Vector
    end
end
```

- [ ] **Step 2: Run to confirm failure**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: `MethodError: analyze_exposure! called with wrong arguments` or similar.

- [ ] **Step 3: Update `analyze_exposure!` signature in pipeline.jl**

Change the function to read `analysis_dir` and pattern from the DB config instead of accepting them as parameters:
```julia
function analyze_exposure!(db::SQLite.DB, exposure_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT e.filename, x.analysis_dir, x.id AS experiment_id
           FROM exposures e
           JOIN samples s ON s.id = e.sample_id
           JOIN experiments x ON x.id = s.experiment_id
           WHERE e.id = ?""", [exposure_id]))
    isempty(rows) && error("Exposure $exposure_id not found")
    row = rows[1]

    cfg = config_from_db(db, row.experiment_id)
    pattern = cfg.integration_pattern   # e.g. "{name}.dat"
    filename = replace(pattern, "{name}" => row.filename)
    dat_path = joinpath(row.analysis_dir, filename)

    isfile(dat_path) || error("dat file not found: $dat_path")
    # Everything below this line is unchanged from the existing analyze_exposure! body:
    # q, I, σ = load_dat(dat_path)
    # peaks_result = Himalaya.findpeaks(q, I, σ)
    # candidates = Himalaya.indexpeaks(q, peaks_result)
    # group_indices = auto_group(candidates)
    # persist_analysis!(db, exposure_id, peaks_result, candidates, group_indices)
    # Keep the existing implementation — only the dat_path construction above changes.
```

> **Important:** The old signature `analyze_exposure!(db, exposure_id, analysis_dir)` is called from `cli.jl` and route handlers. After this change, remove the `analysis_dir` parameter from all call sites (they can look it up from the DB via the experiment row — the config already does this).

- [ ] **Step 4: Update call sites**

Search for all calls to `analyze_exposure!`:
```bash
grep -rn "analyze_exposure!" packages/HimalayaUI/src/
```
Remove the `analysis_dir` argument from each call site — the function now fetches it from the DB.

- [ ] **Step 5: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass. Verify `test_pipeline.jl` and `test_routes_analysis.jl` still pass — these exercise the full analysis path.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/pipeline.jl \
        packages/HimalayaUI/src/cli.jl \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/test_pipeline.jl
git commit -m "feat(pipeline): analyze_exposure! reads path config from DB"
```

---

## Task 8: `himalaya config` subcommand (`new` + `list`)

**Files:**
- Modify: `packages/HimalayaUI/src/cli.jl`

- [ ] **Step 1: Write failing test**

In `packages/HimalayaUI/test/` create a minimal CLI invocation test (or extend an existing one). Check that `himalaya config list` prints `simple` and `himalaya config new` creates a file:
```julia
@testset "cli config list" begin
    buf = IOBuffer()
    HimalayaUI.cli_config_list(buf)
    output = String(take!(buf))
    @test contains(output, "simple")
end

@testset "cli config new creates experiment.toml" begin
    mktempdir() do dir
        HimalayaUI.cli_config_new(; type_name="simple", dir=dir)
        @test isfile(joinpath(dir, "experiment.toml"))
        # Must not overwrite an existing file
        @test_throws ErrorException HimalayaUI.cli_config_new(; type_name="simple", dir=dir)
    end
end
```

- [ ] **Step 2: Implement `cli_config_new` and `cli_config_list` in cli.jl**

Add these functions to `packages/HimalayaUI/src/cli.jl`:
```julia
function cli_config_list(io::IO=stdout)
    types = list_config_types()
    isempty(types) ? println(io, "(no built-in config types found)") :
                     println(io, join(types, "\n"))
end

function cli_config_new(; type_name::String="simple", dir::String)
    dest = joinpath(dir, "experiment.toml")
    isfile(dest) && error("experiment.toml already exists at $dest — will not overwrite")
    src = joinpath(configs_dir(), type_name * ".toml")
    isfile(src) || error("Unknown config type '$type_name'. Run 'himalaya config list' to see options.")
    cp(src, dest)
    println("Created $dest from template '$type_name'")
    println("Edit it to set your experiment name, beamline parameters, and manifest column mappings.")
end
```

- [ ] **Step 3: Wire into the CLI dispatch in `main`**

In the CLI argument parsing section of `cli.jl`, add a `config` subcommand:
```julia
# In the top-level command dispatch:
elseif command == "config"
    subcommand = length(args) > 1 ? args[2] : ""
    if subcommand == "list"
        cli_config_list()
    elseif subcommand == "new"
        s = ArgParseSettings()
        @add_arg_table! s begin
            "--type"
                default = "simple"
            "--dir"
                required = true
        end
        p = parse_args(args[3:end], s)
        cli_config_new(; type_name=p["type"], dir=p["dir"])
    else
        println(stderr, "Usage: himalaya config <list|new> [--type TYPE] [--dir DIR]")
        exit(1)
    end
```

- [ ] **Step 4: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/cli.jl \
        packages/HimalayaUI/test/test_config.jl
git commit -m "feat(cli): himalaya config new/list subcommands"
```

---

## Task 9: Update `himalaya init` to read `experiment.toml`

**Files:**
- Modify: `packages/HimalayaUI/src/cli.jl`
- Modify: `packages/HimalayaUI/test/test_pipeline.jl` (or new integration test)

`cli_init` currently takes a manifest path and uses hardcoded `data_dir`/`analysis_dir`. The new version reads `experiment.toml` from the experiment directory, then uses `resolve_files` for filesystem-based exposure discovery.

- [ ] **Step 1: Write the failing integration test**

Add to `packages/HimalayaUI/test/test_pipeline.jl`:
```julia
@testset "cli_init reads experiment.toml" begin
    db = HimalayaUI.open_db(":memory:")
    mktempdir() do dir
        # Set up experiment directory
        analysis_dir = joinpath(dir, "analysis", "automatic_analysis")
        data_dir     = joinpath(dir, "data")
        mkpath(analysis_dir)
        mkpath(data_dir)

        # Write fake .dat files
        write(joinpath(analysis_dir, "JC001.dat"),
            "0.01 100.0 1.0\n0.02 90.0 1.0\n0.03 80.0 1.0\n")
        write(joinpath(analysis_dir, "JC002.dat"),
            "0.01 100.0 1.0\n0.02 90.0 1.0\n0.03 80.0 1.0\n")

        # Write manifest
        manifest = joinpath(dir, "manifest.csv")
        write(manifest, join([
            "skip row",
            "1\tD1\tUX1\t\t\t\t\t\tJC001-002\tnote\t",
        ], "\n"))

        # Write experiment.toml
        write(joinpath(dir, "experiment.toml"), """
        [experiment]
        name = "Run/Exp"
        description = ""
        manifest = "manifest.csv"
        [beamline]
        energy_kev = 12.0
        flight_path_m = 2.5
        [manifest]
        delimiter = "\\t"
        skip_rows = 1
        header_row = 0
        sample_id = 1
        label = 2
        name = 3
        filenames = 9
        notes_sample = 10
        notes_exposure = 11
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "{name}.dat"
        image = "{name}.tiff"
        """)

        HimalayaUI.cli_init_with_db(db, dir)

        exps = HimalayaUI.list_experiments(db)
        @test length(exps) == 1
        @test exps[1].name == "Run/Exp"
        @test exps[1].energy_kev == 12.0

        samples = HimalayaUI.list_samples(db, exps[1].id)
        @test length(samples) == 1

        exposures = HimalayaUI.list_exposures(db, samples[1].id)
        @test length(exposures) == 2
        filenames = sort([e.filename for e in exposures])
        @test filenames == ["JC001", "JC002"]
    end
end
```

- [ ] **Step 2: Run to confirm failure**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: `UndefVarError: cli_init_with_db` or similar.

- [ ] **Step 3: Refactor `cli_init` in cli.jl**

Extract the DB-level logic into a testable `cli_init_with_db(db, exp_dir)` function and update `cli_init` to open the DB then call it:

```julia
function cli_init_with_db(db::SQLite.DB, exp_dir::String)
    exp_dir = abspath(exp_dir)
    toml_path = joinpath(exp_dir, "experiment.toml")
    isfile(toml_path) || error("experiment.toml not found in $exp_dir. Run 'himalaya config new --dir $exp_dir' first.")

    cfg = load_config(toml_path)
    blob = config_to_toml(cfg)

    # Resolve absolute directories
    data_dir     = isabspath(cfg.data_dir)     ? cfg.data_dir     : joinpath(exp_dir, cfg.data_dir)
    analysis_dir = isabspath(cfg.analysis_dir) ? cfg.analysis_dir : joinpath(exp_dir, cfg.analysis_dir)

    manifest_path = isabspath(cfg.manifest_file) ? cfg.manifest_file :
                    joinpath(exp_dir, cfg.manifest_file)

    exp_id = create_experiment(db;
        name           = cfg.name,
        path           = exp_dir,
        data_dir       = data_dir,
        analysis_dir   = analysis_dir,
        manifest_path  = isfile(manifest_path) ? manifest_path : nothing,
        config         = blob,
        experiment_type = cfg.exposure_type,
        energy_kev     = cfg.energy_kev,
        flight_path_m  = cfg.flight_path_m,
    )

    if isfile(manifest_path)
        samples = parse_manifest(cfg, manifest_path)
        for ms in samples
            s_id = create_sample(db;
                experiment_id = exp_id,
                label = ms.label,
                name  = ms.name)

            # Add sample note
            if !isempty(ms.notes_sample)
                DBInterface.execute(db,
                    "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, 'note', ?, 'manifest')",
                    [s_id, ms.notes_sample])
            end

            for prefix in ms.filenames
                stems = resolve_files(cfg, analysis_dir, prefix, cfg.integration_pattern)
                isempty(stems) && @warn "No files found for prefix '$prefix' in $analysis_dir"

                for stem in stems
                    image_candidate = joinpath(data_dir,
                        replace(cfg.image_pattern, "{name}" => stem))
                    image_path = isfile(image_candidate) ? image_candidate : nothing

                    e_id = create_exposure(db;
                        sample_id  = s_id,
                        filename   = stem,
                        image_path = image_path)

                    if !isempty(ms.notes_exposure)
                        DBInterface.execute(db,
                            "INSERT INTO exposure_tags (exposure_id, key, value, source) VALUES (?, 'note', ?, 'manifest')",
                            [e_id, ms.notes_exposure])
                    end
                end
            end
        end
        println("Imported $(length(samples)) samples from $(basename(manifest_path)).")
    end

    println("Initialized experiment '$(cfg.name)' (id=$exp_id).")
    exp_id
end

function cli_init(args)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "dir"
            help     = "experiment directory containing experiment.toml"
            required = true
    end
    p = parse_args(args, s)
    exp_dir = p["dir"]

    db_path = get(ENV, "HIMALAYA_DB", joinpath(dirname(exp_dir), "himalaya.db"))
    db = open_db(db_path)
    cli_init_with_db(db, exp_dir)
end
```

> **Read-only check:** Verify this function never calls `write`, `open(..., "w")`, `mkdir`, `cp`, or `mv` on any path derived from `exp_dir`, `data_dir`, or `analysis_dir`. All writes go to `db`.

> **`create_sample` / `create_exposure`:** Check `db.jl` — if these helper functions don't exist (only `create_experiment` was confirmed), replace the calls with raw `DBInterface.execute(db, "INSERT INTO samples ...", [...])` following the same pattern as `create_experiment`. Use `Int(DBInterface.lastrowid(res))` to get the inserted id.

- [ ] **Step 4: Remove the old `--beamline` flag**

The `--beamline` flag in the old `cli_init` is unused. It is not present in the new `cli_init_with_db` signature — confirm it's gone.

- [ ] **Step 5: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/cli.jl \
        packages/HimalayaUI/test/test_pipeline.jl
git commit -m "feat(cli): himalaya init reads experiment.toml; removes --beamline flag"
```

---

## Task 10: `himalaya reingest` command and REST route

**Files:**
- Modify: `packages/HimalayaUI/src/cli.jl`
- Modify: `packages/HimalayaUI/src/server.jl`
- Modify (or create): `packages/HimalayaUI/test/test_routes_experiments.jl`

`reingest` re-reads `experiment.toml` + `manifest.csv` for an already-registered experiment and updates the DB. It preserves curated exposures (those with `status IS NOT NULL` or manual peaks).

- [ ] **Step 1: Write failing test for reingest logic**

Add to `packages/HimalayaUI/test/test_pipeline.jl`:
```julia
@testset "reingest updates notes, inserts new exposures, preserves curated" begin
    db = HimalayaUI.open_db(":memory:")
    mktempdir() do dir
        analysis_dir = joinpath(dir, "analysis", "automatic_analysis")
        data_dir     = joinpath(dir, "data")
        mkpath(analysis_dir)
        mkpath(data_dir)
        for name in ["JC001", "JC002", "JC003"]
            write(joinpath(analysis_dir, name * ".dat"),
                "0.01 100.0 1.0\n0.02 90.0 1.0\n0.03 80.0 1.0\n")
        end

        # Initial manifest: only JC001-002
        manifest = joinpath(dir, "manifest.csv")
        write(manifest, join(["skip", "1\tD1\tUX1\t\t\t\t\t\tJC001-002\told_note\t"], "\n"))
        write(joinpath(dir, "experiment.toml"), """
        [experiment]
        name = "R/E"
        description = ""
        manifest = "manifest.csv"
        [beamline]
        energy_kev = 0.0
        flight_path_m = 0.0
        [manifest]
        delimiter = "\\t"
        skip_rows = 1
        header_row = 0
        sample_id = 1
        label = 2
        name = 3
        filenames = 9
        notes_sample = 10
        notes_exposure = 11
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "{name}.dat"
        image = "{name}.tiff"
        """)
        exp_id = HimalayaUI.cli_init_with_db(db, dir)

        # Curate JC001
        DBInterface.execute(db,
            "UPDATE exposures SET status = 'accepted' WHERE filename = 'JC001'")

        # Update manifest to add JC003 and change notes
        write(manifest, join(["skip", "1\tD1\tUX1\t\t\t\t\t\tJC001-003\tnew_note\t"], "\n"))
        HimalayaUI.reingest!(db, exp_id, dir)

        # JC003 should now exist
        all_exps = Tables.rowtable(DBInterface.execute(db,
            "SELECT filename, status FROM exposures ORDER BY filename"))
        filenames = [r.filename for r in all_exps]
        @test "JC003" in filenames

        # JC001 curation preserved
        jc001 = first(filter(r -> r.filename == "JC001", all_exps))
        @test jc001.status == "accepted"
    end
end
```

- [ ] **Step 2: Run to confirm failure**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: `UndefVarError: reingest! not defined`.

- [ ] **Step 3: Implement `reingest!` in cli.jl**

```julia
function reingest!(db::SQLite.DB, experiment_id::Int, exp_dir::String)
    exp_dir  = abspath(exp_dir)
    toml_path = joinpath(exp_dir, "experiment.toml")
    isfile(toml_path) || error("experiment.toml not found in $exp_dir")

    cfg  = load_config(toml_path)
    blob = config_to_toml(cfg)

    data_dir     = isabspath(cfg.data_dir)     ? cfg.data_dir     : joinpath(exp_dir, cfg.data_dir)
    analysis_dir = isabspath(cfg.analysis_dir) ? cfg.analysis_dir : joinpath(exp_dir, cfg.analysis_dir)
    manifest_path = isabspath(cfg.manifest_file) ? cfg.manifest_file :
                    joinpath(exp_dir, cfg.manifest_file)

    # Update experiment metadata
    DBInterface.execute(db,
        """UPDATE experiments
           SET name=?, config=?, experiment_type=?, energy_kev=?, flight_path_m=?,
               data_dir=?, analysis_dir=?, manifest_path=?
           WHERE id=?""",
        [cfg.name, blob, cfg.exposure_type, cfg.energy_kev, cfg.flight_path_m,
         data_dir, analysis_dir,
         isfile(manifest_path) ? manifest_path : nothing,
         experiment_id])

    isfile(manifest_path) || return

    samples = parse_manifest(cfg, manifest_path)
    for ms in samples
        # Upsert sample (match by name within experiment)
        existing_sample = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM samples WHERE experiment_id=? AND name=?",
            [experiment_id, ms.name]))
        s_id = if isempty(existing_sample)
            create_sample(db; experiment_id=experiment_id, label=ms.label, name=ms.name)
        else
            # Update label in case it changed
            DBInterface.execute(db, "UPDATE samples SET label=? WHERE id=?",
                [ms.label, existing_sample[1].id])
            existing_sample[1].id
        end

        for prefix in ms.filenames
            stems = resolve_files(cfg, analysis_dir, prefix, cfg.integration_pattern)
            for stem in stems
                existing_exp = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, status FROM exposures WHERE sample_id=? AND filename=?",
                    [s_id, stem]))
                if isempty(existing_exp)
                    image_candidate = joinpath(data_dir,
                        replace(cfg.image_pattern, "{name}" => stem))
                    image_path = isfile(image_candidate) ? image_candidate : nothing
                    create_exposure(db; sample_id=s_id, filename=stem, image_path=image_path)
                end
                # Existing exposures with status or manual peaks are untouched
            end
        end
    end
    println("Reingested experiment $experiment_id from $exp_dir.")
end

function cli_reingest(args)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "dir"
            help     = "experiment directory"
            required = true
    end
    p = parse_args(args, s)
    exp_dir = abspath(p["dir"])

    db_path = get(ENV, "HIMALAYA_DB", joinpath(dirname(exp_dir), "himalaya.db"))
    db = open_db(db_path)

    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM experiments WHERE path=?", [exp_dir]))
    isempty(rows) && error("No experiment registered at $exp_dir. Run 'himalaya init' first.")
    reingest!(db, rows[1].id, exp_dir)
end
```

Wire `cli_reingest` into the top-level dispatch in `main`:
```julia
elseif command == "reingest"
    cli_reingest(args[2:end])
```

- [ ] **Step 4: Add REST route in server.jl**

In `packages/HimalayaUI/src/server.jl`, add:
```julia
@post "/api/experiments/{id}/reingest" function(req::HTTP.Request, id::Int)
    exp = get_experiment(db[], id)
    exp === nothing && return HTTP.Response(404, "Experiment not found")
    try
        reingest!(db[], id, exp.path)
        return HTTP.Response(200, JSON3.write(Dict("ok" => true)))
    catch e
        return HTTP.Response(500, JSON3.write(Dict("error" => sprint(showerror, e))))
    end
end
```

- [ ] **Step 5: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass including new reingest test.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/cli.jl \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/test_pipeline.jl
git commit -m "feat(cli): himalaya reingest command + POST /api/experiments/:id/reingest"
```

---

## Task 11: Final verification

- [ ] **Step 1: Run the full test suite**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass.

- [ ] **Step 2: Run frontend tests**

```bash
cd packages/HimalayaUI/frontend
npm test
npm run build
```
Expected: all pass, build succeeds.

- [ ] **Step 3: Verify experiment directory is read-only at runtime**

After running `himalaya init` against a test experiment directory:
```bash
# Check no new files were created in the experiment directory other than experiment.toml
ls -la /tmp/test-exp/
```
Only `experiment.toml` (written by `config new`), `manifest.csv`, and the `data/` and `analysis/` directories should be present. No `himalaya.db`, no lock files, no generated files.

- [ ] **Step 4: Smoke-test the full CLI flow**

```bash
# Create a test experiment directory
mkdir -p /tmp/smoke-exp/data
mkdir -p /tmp/smoke-exp/analysis/automatic_analysis

# Generate config
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  config new --type simple --dir /tmp/smoke-exp

# Edit experiment.toml to set name, then init
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  init /tmp/smoke-exp

# List available config types
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  config list
```

- [ ] **Step 5: Final commit**

```bash
git add .
git commit -m "feat(config): experiment config system complete — read-only adapter-driven I/O"
```
