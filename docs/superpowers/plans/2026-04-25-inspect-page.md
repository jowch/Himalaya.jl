# Inspect Page Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an "Inspect" tab to HimalayaUI that lets users review 2D SAXS detector images, accept/reject exposures, mark one for indexing, and edit sample metadata — with a responsive three-breakpoint layout.

**Architecture:** Backend adds two new columns to `exposures` (status, image_path), discovers TIFF files at ingestion, and serves a new `GET /api/exposures/:id/image` route that returns a log-normalized grayscale PNG. The frontend renders this PNG on a `<canvas>` and applies a theme-aware colormap (CSS var → linear LUT). The Inspect page uses the same three-card grid pattern as IndexPage, reflowing across three breakpoints.

**Tech Stack:** Julia/Oxygen.jl (backend), FileIO + TiffImages + ImageIO + ImageTransformations + ImageCore (image processing), React 18 + TypeScript strict + Tailwind v4 + TanStack Query + Zustand (frontend), Vitest + Playwright (tests).

---

## File Map

**Backend — modify:**
- `packages/HimalayaUI/src/db.jl` — add `migrate_schema!`, update `create_exposure!`, `get_exposures`
- `packages/HimalayaUI/src/pipeline.jl` — add `find_tiff_for_dat`, call from ingestion
- `packages/HimalayaUI/src/cli.jl` — wire TIFF discovery in `cli_init`; add auto-fallback in `cli_analyze`
- `packages/HimalayaUI/src/routes_exposures.jl` — add image route, status route, exclude_rejected param
- `packages/HimalayaUI/src/json.jl` — no change needed (status/image_path pass through `row_to_json` automatically)
- `packages/HimalayaUI/Project.toml` — add FileIO, TiffImages, ImageIO, ImageTransformations, ImageCore

**Frontend — modify:**
- `packages/HimalayaUI/frontend/src/api.ts` — update `Exposure` type, add `setExposureStatus`, `selectExposure`, `addExposureTag`, `removeExposureTag`, `listExposures` params
- `packages/HimalayaUI/frontend/src/state.ts` — add `"inspect"` to `PageId`
- `packages/HimalayaUI/frontend/src/queries.ts` — update `useExposures`, add `useSetExposureStatus`, `useSelectExposure`, `useAddExposureTag`, `useRemoveExposureTag`
- `packages/HimalayaUI/frontend/src/components/TabRocker.tsx` — add Inspect tab before Index
- `packages/HimalayaUI/frontend/src/components/AppShell.tsx` — render InspectPage

**Frontend — create:**
- `packages/HimalayaUI/frontend/src/pages/InspectPage.tsx`
- `packages/HimalayaUI/frontend/src/components/DetectorImage.tsx`
- `packages/HimalayaUI/frontend/src/components/ThumbnailGallery.tsx`
- `packages/HimalayaUI/frontend/src/components/DetectorImageCard.tsx`
- `packages/HimalayaUI/frontend/src/components/SampleMetadataCard.tsx`

---

## Task 1: DB schema migration

Add `status` and `image_path` columns to the existing `exposures` table via a `migrate_schema!` function called from `open_db`. SQLite's `ALTER TABLE` has no `IF NOT EXISTS`, so wrap each statement in try/catch.

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl`
- Test: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Write the failing test**

In `packages/HimalayaUI/test/runtests.jl`, add inside the appropriate `@testset`:

```julia
@testset "exposures schema migration" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.migrate_schema!(db)
    # Verify columns exist by inserting and reading back
    DBInterface.execute(db,
        "INSERT INTO exposures (sample_id, filename, kind, selected, status, image_path)
         VALUES (NULL, 'test.dat', 'file', 0, 'accepted', '/tmp/test.tiff')")
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT status, image_path FROM exposures"))
    @test length(rows) == 1
    @test rows[1].status == "accepted"
    @test rows[1].image_path == "/tmp/test.tiff"
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: ERROR — `HimalayaUI.migrate_schema!` not defined.

- [ ] **Step 3: Implement `migrate_schema!` in `db.jl`**

Add after `create_schema!`:

```julia
function migrate_schema!(db::SQLite.DB)
    stmts = [
        "ALTER TABLE exposures ADD COLUMN status TEXT CHECK (status IN ('accepted', 'rejected'))",
        "ALTER TABLE exposures ADD COLUMN image_path TEXT",
    ]
    for stmt in stmts
        try
            DBInterface.execute(db, stmt)
        catch
            # column already exists — safe to ignore
        end
    end
end
```

Update `open_db` to call it:

```julia
function open_db(experiment_path::String)::SQLite.DB
    db_path = joinpath(experiment_path, "himalaya.db")
    db = SQLite.DB(db_path)
    create_schema!(db)
    migrate_schema!(db)
    DBInterface.execute(db, "PRAGMA foreign_keys = ON")
    db
end
```

Update `create_exposure!` signature to accept the new fields (default nothing):

```julia
function create_exposure!(db::SQLite.DB;
        sample_id::Int,
        filename::Union{String,Nothing}   = nothing,
        kind::String                       = "file",
        selected::Bool                     = false,
        status::Union{String,Nothing}      = nothing,
        image_path::Union{String,Nothing}  = nothing)
    result = DBInterface.execute(db,
        "INSERT INTO exposures (sample_id, filename, kind, selected, status, image_path)
         VALUES (?, ?, ?, ?, ?, ?)",
        [sample_id, filename, kind, Int(selected), status, image_path])
    Int(DBInterface.lastrowid(result))
end
```

Update `get_exposures` (no change needed — `SELECT *` picks up new columns automatically).

- [ ] **Step 4: Run test to verify it passes**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/runtests.jl
git commit -m "feat(db): add status and image_path columns to exposures via migrate_schema!"
```

---

## Task 2: Add image processing packages

**Files:**
- Modify: `packages/HimalayaUI/Project.toml`

- [ ] **Step 1: Add packages**

```bash
cd packages/HimalayaUI
julia --project=. -e '
using Pkg
Pkg.add(["FileIO", "TiffImages", "ImageIO", "ImageTransformations", "ImageCore", "ColorTypes"])
'
```

- [ ] **Step 2: Verify they load**

```bash
julia --project=packages/HimalayaUI -e '
using FileIO, TiffImages, ImageIO, ImageTransformations, ImageCore, ColorTypes
println("image packages OK")
'
```

Expected: `image packages OK`.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/Project.toml packages/HimalayaUI/Manifest.toml
git commit -m "deps: add image processing packages (FileIO, TiffImages, ImageIO, ImageTransformations, ImageCore)"
```

---

## Task 3: Image serving helpers and `GET /api/exposures/:id/image`

Add a helper module for image loading/normalization/encoding, then the route.

**Files:**
- Create: `packages/HimalayaUI/src/image.jl`
- Modify: `packages/HimalayaUI/src/routes_exposures.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl` (include new file)
- Test: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "image route" begin
    # Create a small test TIFF
    using FileIO, ImageCore, TiffImages
    tiff_path = tempname() * ".tiff"
    test_img = Gray.(rand(Float32, 32, 32))
    save(tiff_path, test_img)

    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.migrate_schema!(db)
    exp_id  = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
    exp_id2 = HimalayaUI.create_exposure!(db; sample_id=samp_id, image_path=tiff_path)

    server = HimalayaUI.start_test_server!(db, 8791)
    try
        r = HTTP.get("http://127.0.0.1:8791/api/exposures/$exp_id2/image")
        @test r.status == 200
        @test r.headers["Content-Type"] == "image/png"
        @test length(r.body) > 100

        # thumb variant
        rt = HTTP.get("http://127.0.0.1:8791/api/exposures/$exp_id2/image?thumb=1")
        @test rt.status == 200
        @test length(rt.body) < length(r.body)  # thumbnail is smaller

        # missing image_path → 404
        exp_id3 = HimalayaUI.create_exposure!(db; sample_id=samp_id)
        r404 = HTTP.get("http://127.0.0.1:8791/api/exposures/$exp_id3/image";
                        status_exception=false)
        @test r404.status == 404
    finally
        HimalayaUI.stop_test_server!()
        rm(tiff_path; force=true)
    end
end
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: ERROR — route not defined.

- [ ] **Step 3: Create `src/image.jl`**

```julia
using FileIO, TiffImages, ImageIO, ImageTransformations, ImageCore, ColorTypes

"""
    load_and_lognormalize(path) -> Matrix{Gray{Float32}}

Load a TIFF, convert to grayscale, apply log1p normalization to [0,1].
"""
function load_and_lognormalize(path::String)
    raw   = load(path)
    gray  = Gray.(raw)
    vals  = Float32.(channelview(gray))
    lv    = log1p.(vals)
    m     = maximum(lv)
    normed = m > 0 ? lv ./ m : lv
    colorview(Gray, normed)
end

"""
    encode_png(img) -> Vector{UInt8}

Encode any grayscale image as PNG bytes.
"""
function encode_png(img)
    buf = IOBuffer()
    save(Stream{format"PNG"}(buf), img)
    take!(buf)
end

"""
    resize_to_fit(img, max_px) -> image

Downscale so the longest side is ≤ max_px; no-op if already smaller.
"""
function resize_to_fit(img, max_px::Int)
    h, w = size(img)
    max(h, w) <= max_px && return img
    scale = max_px / max(h, w)
    imresize(img; ratio=scale)
end
```

- [ ] **Step 4: Add include to `HimalayaUI.jl`**

Find the `include(...)` block in `packages/HimalayaUI/src/HimalayaUI.jl` and add:

```julia
include("image.jl")
```

- [ ] **Step 5: Add the route in `routes_exposures.jl`**

Add inside `register_exposures_routes!()`:

```julia
@get "/api/exposures/{id}/image" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT image_path FROM exposures WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "exposure not found")))

    ip = rows[1].image_path
    (ip === nothing || ip isa Missing) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "no image for this exposure")))

    params   = HTTP.queryparams(req)
    is_thumb = get(params, "thumb", "0") == "1"

    img = load_and_lognormalize(String(ip))
    if is_thumb
        img = resize_to_fit(img, 128)
    end
    bytes = encode_png(img)

    HTTP.Response(200,
        ["Content-Type"  => "image/png",
         "Cache-Control" => "max-age=3600"],
        bytes)
end
```

- [ ] **Step 6: Run test to verify it passes**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/src/image.jl packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/src/routes_exposures.jl packages/HimalayaUI/test/runtests.jl
git commit -m "feat(api): GET /api/exposures/:id/image — log-normalized TIFF → PNG"
```

---

## Task 4: `PATCH /api/exposures/:id/status` route

**Files:**
- Modify: `packages/HimalayaUI/src/routes_exposures.jl`
- Test: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "PATCH /api/exposures/:id/status" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.migrate_schema!(db)
    exp_id  = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
    eid     = HimalayaUI.create_exposure!(db; sample_id=samp_id)

    server = HimalayaUI.start_test_server!(db, 8792)
    try
        # accept
        r = HTTP.patch("http://127.0.0.1:8792/api/exposures/$eid/status";
            body=JSON3.write(Dict(:status => "accepted")),
            headers=["Content-Type" => "application/json"])
        @test r.status == 200
        row = Tables.rowtable(DBInterface.execute(db,
            "SELECT status FROM exposures WHERE id = ?", [eid]))[1]
        @test row.status == "accepted"

        # reject
        HTTP.patch("http://127.0.0.1:8792/api/exposures/$eid/status";
            body=JSON3.write(Dict(:status => "rejected")),
            headers=["Content-Type" => "application/json"])
        row2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT status FROM exposures WHERE id = ?", [eid]))[1]
        @test row2.status == "rejected"

        # null (clear)
        HTTP.patch("http://127.0.0.1:8792/api/exposures/$eid/status";
            body=JSON3.write(Dict(:status => nothing)),
            headers=["Content-Type" => "application/json"])
        row3 = Tables.rowtable(DBInterface.execute(db,
            "SELECT status FROM exposures WHERE id = ?", [eid]))[1]
        @test row3.status === missing || row3.status === nothing

        # invalid value → 422
        r422 = HTTP.patch("http://127.0.0.1:8792/api/exposures/$eid/status";
            body=JSON3.write(Dict(:status => "garbage")),
            headers=["Content-Type" => "application/json"],
            status_exception=false)
        @test r422.status == 422
    finally
        HimalayaUI.stop_test_server!()
    end
end
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: connection refused or 404.

- [ ] **Step 3: Add route in `routes_exposures.jl`**

Add inside `register_exposures_routes!()`:

```julia
@patch "/api/exposures/{id}/status" function(req::HTTP.Request, id::Int)
    db   = current_db()
    body = json(req)

    raw_status = get(body, :status, nothing)
    status = raw_status === nothing || raw_status isa JSON3.Null ? nothing : String(raw_status)

    if status !== nothing && status ∉ ("accepted", "rejected")
        return HTTP.Response(422,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "status must be 'accepted', 'rejected', or null")))
    end

    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM exposures WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "exposure not found")))

    DBInterface.execute(db,
        "UPDATE exposures SET status = ? WHERE id = ?", [status, id])

    log_action!(db, req; action = "set_status",
        entity_type = "exposure", entity_id = id,
        note = something(status, "null"))

    HTTP.Response(200, ["Content-Type" => "application/json"],
        JSON3.write(Dict(:id => id, :status => status)))
end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_exposures.jl packages/HimalayaUI/test/runtests.jl
git commit -m "feat(api): PATCH /api/exposures/:id/status — accept/reject exposures"
```

---

## Task 5: `exclude_rejected` filter on exposures list

**Files:**
- Modify: `packages/HimalayaUI/src/routes_exposures.jl`
- Test: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "GET /api/samples/:id/exposures exclude_rejected" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.migrate_schema!(db)
    exp_id  = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
    HimalayaUI.create_exposure!(db; sample_id=samp_id, filename="good.dat")
    bad_id = HimalayaUI.create_exposure!(db; sample_id=samp_id, filename="bad.dat",
                                          status="rejected")

    server = HimalayaUI.start_test_server!(db, 8793)
    try
        # without filter — both returned
        r_all = HTTP.get("http://127.0.0.1:8793/api/samples/$samp_id/exposures")
        all_exps = JSON3.read(r_all.body)
        @test length(all_exps) == 2

        # with filter — rejected excluded
        r_fil = HTTP.get("http://127.0.0.1:8793/api/samples/$samp_id/exposures?exclude_rejected=true")
        fil_exps = JSON3.read(r_fil.body)
        @test length(fil_exps) == 1
        @test fil_exps[1][:filename] == "good.dat"
    finally
        HimalayaUI.stop_test_server!()
    end
end
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: FAIL — both exposures returned even with filter.

- [ ] **Step 3: Update the list route in `routes_exposures.jl`**

Replace the existing `@get "/api/samples/{id}/exposures"` handler:

```julia
@get "/api/samples/{id}/exposures" function(req::HTTP.Request, id::Int)
    db     = current_db()
    params = HTTP.queryparams(req)
    exclude_rejected = get(params, "exclude_rejected", "false") == "true"

    sql = exclude_rejected ?
        "SELECT * FROM exposures WHERE sample_id = ? AND (status IS NULL OR status != 'rejected') ORDER BY id" :
        "SELECT * FROM exposures WHERE sample_id = ? ORDER BY id"

    exs = Tables.rowtable(DBInterface.execute(db, sql, [id]))
    out = map(exs) do e
        tags = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, key, value, source FROM exposure_tags
             WHERE exposure_id = ? ORDER BY id", [Int(e.id)]))
        srcs = Tables.rowtable(DBInterface.execute(db,
            "SELECT source_exposure_id, role FROM exposure_sources
             WHERE averaged_exposure_id = ?", [Int(e.id)]))
        d = row_to_json(e; bool_keys = (:selected,))
        d[:tags]    = rows_to_json(tags)
        d[:sources] = rows_to_json(srcs)
        d
    end
    HTTP.Response(200, ["Content-Type" => "application/json"],
        JSON3.write(out))
end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_exposures.jl packages/HimalayaUI/test/runtests.jl
git commit -m "feat(api): exclude_rejected filter on GET /api/samples/:id/exposures"
```

---

## Task 6: TIFF discovery at ingestion + auto-fallback in `cli_analyze`

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl`
- Modify: `packages/HimalayaUI/src/cli.jl`
- Test: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "find_tiff_for_dat" begin
    dir = mktempdir()
    dat_path  = joinpath(dir, "sample_pos1.dat")
    tiff_path = joinpath(dir, "sample_pos1.tiff")
    touch(dat_path)
    touch(tiff_path)

    @test HimalayaUI.find_tiff_for_dat(dat_path) == tiff_path

    # .tif extension fallback
    tif_path = joinpath(dir, "sample_pos2.tif")
    dat_path2 = joinpath(dir, "sample_pos2.dat")
    touch(dat_path2); touch(tif_path)
    @test HimalayaUI.find_tiff_for_dat(dat_path2) == tif_path

    # no match → nothing
    dat_path3 = joinpath(dir, "sample_pos3.dat")
    touch(dat_path3)
    @test HimalayaUI.find_tiff_for_dat(dat_path3) === nothing
end
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: `HimalayaUI.find_tiff_for_dat` not defined.

- [ ] **Step 3: Add `find_tiff_for_dat` to `pipeline.jl`**

Add at the top of `pipeline.jl` (after existing using statements):

```julia
"""
    find_tiff_for_dat(dat_path) -> Union{String, Nothing}

Search for a TIFF image co-located with a .dat file, by replacing
the extension with .tiff or .tif. Returns the absolute path if found,
otherwise nothing.
"""
function find_tiff_for_dat(dat_path::String)::Union{String,Nothing}
    base = splitext(dat_path)[1]
    for ext in (".tiff", ".tif")
        candidate = base * ext
        isfile(candidate) && return candidate
    end
    nothing
end
```

- [ ] **Step 4: Wire into `cli_init` in `cli.jl`**

In the `for filename in ms.filenames` loop within `cli_init`, after `create_exposure!`, add TIFF discovery:

```julia
for filename in ms.filenames
    # Resolve absolute path for TIFF discovery
    abs_filename = isabspath(filename) ? filename :
                   joinpath(data_dir, filename)
    image_path = isfile(abs_filename) ? find_tiff_for_dat(abs_filename) : nothing
    create_exposure!(db; sample_id = s_id, filename = filename,
                     image_path = image_path)
    # ... existing exposure_tag insertion for notes_exposure stays unchanged
end
```

- [ ] **Step 5: Add auto-fallback in `cli_analyze`**

In `cli_analyze`, before calling `analyze_exposure!`, add:

```julia
for sample in samples
    exposures = get_exposures(db, Int(sample.id))
    # Auto-fallback: if no exposure is explicitly selected, use first accepted one
    has_selected = any(e -> Int(e.selected) == 1, exposures)
    if !has_selected
        first_accepted = findfirst(
            e -> !ismissing(e.status) && e.status == "accepted", exposures)
        if first_accepted !== nothing
            fallback_id = Int(exposures[first_accepted].id)
            DBInterface.execute(db,
                "UPDATE exposures SET selected = 0 WHERE sample_id = ?",
                [Int(sample.id)])
            DBInterface.execute(db,
                "UPDATE exposures SET selected = 1 WHERE id = ?",
                [fallback_id])
        end
    end
    for exp_row in exposures
        # ... existing analyze loop unchanged
    end
end
```

- [ ] **Step 6: Run test to verify it passes**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/src/pipeline.jl packages/HimalayaUI/src/cli.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(pipeline): discover TIFF at ingestion, auto-fallback to first accepted in analyze"
```

---

## Task 7: Frontend types and API functions

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts`

- [ ] **Step 1: Update `Exposure` type and add new fetchers**

Replace the `// Exposures` section in `api.ts`:

```typescript
// Exposures
export interface ExposureTag {
  id: number;
  key: string;
  value: string;
  source: string;
}

export interface Exposure {
  id: number;
  sample_id: number;
  filename: string | null;
  kind: "file" | "averaged" | "background_subtracted";
  selected: boolean;
  status: "accepted" | "rejected" | null;
  image_path: string | null;
  tags: ExposureTag[];
  sources: unknown[];
}

export const listExposures = (
  sample_id: number,
  opts?: { excludeRejected?: boolean },
) => {
  const qs = opts?.excludeRejected ? "?exclude_rejected=true" : "";
  return request<Exposure[]>("GET", `/api/samples/${sample_id}/exposures${qs}`);
};

export const setExposureStatus = (
  id: number,
  status: "accepted" | "rejected" | null,
  opts?: AuthOpts,
) => request<{ id: number; status: string | null }>(
  "PATCH", `/api/exposures/${id}/status`, { status }, opts);

export const selectExposure = (id: number, opts?: AuthOpts) =>
  request<{ id: number; selected: boolean }>(
    "PATCH", `/api/exposures/${id}/select`, {}, opts);

export const addExposureTag = (
  id: number, key: string, value: string, opts?: AuthOpts,
) => request<ExposureTag>(
  "POST", `/api/exposures/${id}/tags`, { key, value }, opts);

export const removeExposureTag = (
  id: number, tag_id: number, opts?: AuthOpts,
) => request<void>(
  "DELETE", `/api/exposures/${id}/tags/${tag_id}`, undefined, opts);
```

- [ ] **Step 2: Run TypeScript check**

```bash
cd packages/HimalayaUI/frontend
npm run build -- --noEmit 2>&1 | head -30
```

Expected: no errors in api.ts.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts
git commit -m "feat(frontend): update Exposure type, add setExposureStatus/selectExposure/tag mutations to api.ts"
```

---

## Task 8: State and queries updates

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/state.ts`
- Modify: `packages/HimalayaUI/frontend/src/queries.ts`

- [ ] **Step 1: Add `"inspect"` to `PageId` in `state.ts`**

```typescript
export type PageId = "index" | "compare" | "inspect";
```

Also update the Zustand `partialize` version bump to avoid stale persisted state:

```typescript
version: 3,   // was 2
```

- [ ] **Step 2: Update `useExposures` and add new mutations in `queries.ts`**

Replace `useExposures`:

```typescript
export function useExposures(
  sampleId: number | undefined,
  opts?: { excludeRejected?: boolean },
) {
  const excludeRejected = opts?.excludeRejected ?? false;
  return useQuery({
    queryKey: ["sample", sampleId ?? "none", "exposures", { excludeRejected }] as const,
    queryFn: () => api.listExposures(sampleId as number, { excludeRejected }),
    enabled: sampleId !== undefined,
  });
}
```

Add new mutations after the existing ones:

```typescript
export function useSetExposureStatus(sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: ({ exposureId, status }: {
      exposureId: number;
      status: "accepted" | "rejected" | null;
    }) => api.setExposureStatus(exposureId, status, authOpts(username)),
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) });
      qc.invalidateQueries({ queryKey: ["sample", sampleId, "exposures",
        { excludeRejected: true }] });
    },
  });
}

export function useSelectExposure(sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (exposureId: number) =>
      api.selectExposure(exposureId, authOpts(username)),
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) });
      qc.invalidateQueries({ queryKey: ["sample", sampleId, "exposures",
        { excludeRejected: true }] });
    },
  });
}

export function useAddExposureTag(sampleId: number, exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: ({ key, value }: { key: string; value: string }) =>
      api.addExposureTag(exposureId, key, value, authOpts(username)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) }),
  });
}

export function useRemoveExposureTag(sampleId: number, exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (tagId: number) =>
      api.removeExposureTag(exposureId, tagId, authOpts(username)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) }),
  });
}
```

- [ ] **Step 3: Run TypeScript check**

```bash
cd packages/HimalayaUI/frontend
npm run build -- --noEmit 2>&1 | head -30
```

Expected: no errors.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/state.ts \
        packages/HimalayaUI/frontend/src/queries.ts
git commit -m "feat(frontend): add inspect page state, exposure status/tag mutations"
```

---

## Task 9: TabRocker + AppShell wiring

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/TabRocker.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/AppShell.tsx`

- [ ] **Step 1: Write failing Vitest test**

Create `packages/HimalayaUI/frontend/test/TabRocker.test.tsx`:

```tsx
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { TabRocker } from "../src/components/TabRocker";
import { useAppState } from "../src/state";

vi.mock("../src/state", () => ({
  useAppState: vi.fn(),
}));

test("renders Inspect tab before Index tab", () => {
  (useAppState as vi.Mock).mockImplementation((sel: (s: unknown) => unknown) =>
    sel({ activePage: "index", setActivePage: vi.fn() }));

  render(<TabRocker />);
  const tabs = screen.getAllByRole("tab");
  expect(tabs[0]).toHaveTextContent("Inspect");
  expect(tabs[1]).toHaveTextContent("Index");
  expect(tabs[2]).toHaveTextContent("Compare");
});
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd packages/HimalayaUI/frontend
npm test -- TabRocker
```

Expected: FAIL — only 2 tabs found, "Inspect" not present.

- [ ] **Step 3: Update `TabRocker.tsx`**

```tsx
import { useAppState, type PageId } from "../state";

const TABS: readonly { id: PageId; label: string }[] = [
  { id: "inspect", label: "Inspect" },
  { id: "index",   label: "Index"   },
  { id: "compare", label: "Compare" },
];

export function TabRocker(): JSX.Element {
  const activePage = useAppState((s) => s.activePage);
  const setPage    = useAppState((s) => s.setActivePage);

  return (
    <div
      role="tablist"
      data-testid="tab-rocker"
      className="inline-flex items-center gap-0.5 p-0.5
                 bg-bg-elevated border border-border rounded-full"
    >
      {TABS.map((t) => {
        const active = t.id === activePage;
        return (
          <button
            key={t.id}
            role="tab"
            aria-selected={active}
            data-testid={`tab-${t.id}`}
            data-active={active || undefined}
            onClick={() => setPage(t.id)}
            className={
              "px-3.5 py-1 rounded-full font-sans text-[12px] font-medium " +
              (active
                ? "bg-bg-hover text-fg shadow-inner"
                : "text-fg-muted hover:text-fg")
            }
          >
            {t.label}
          </button>
        );
      })}
    </div>
  );
}
```

- [ ] **Step 4: Update `AppShell.tsx`**

Add import and conditional render:

```tsx
import { InspectPage } from "../pages/InspectPage";
```

Replace the page render line:

```tsx
{activePage === "index"   && <IndexPage />}
{activePage === "inspect" && <InspectPage />}
{activePage === "compare" && <ComparePage />}
```

- [ ] **Step 5: Create stub `InspectPage.tsx`** (so TypeScript compiles)

```tsx
export function InspectPage(): JSX.Element {
  return <div data-testid="inspect-page">Inspect</div>;
}
```

- [ ] **Step 6: Run Vitest to verify it passes**

```bash
cd packages/HimalayaUI/frontend
npm test -- TabRocker
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/TabRocker.tsx \
        packages/HimalayaUI/frontend/src/components/AppShell.tsx \
        packages/HimalayaUI/frontend/src/pages/InspectPage.tsx \
        packages/HimalayaUI/frontend/test/TabRocker.test.tsx
git commit -m "feat(nav): add Inspect tab to TabRocker, wire AppShell"
```

---

## Task 10: `DetectorImage` component

Fetches a grayscale PNG from the server and colorizes it with a theme-aware linear LUT on a `<canvas>`.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/DetectorImage.tsx`
- Test: `packages/HimalayaUI/frontend/test/DetectorImage.test.tsx`

- [ ] **Step 1: Write failing test**

Create `packages/HimalayaUI/frontend/test/DetectorImage.test.tsx`:

```tsx
import { render, screen, waitFor } from "@testing-library/react";
import { DetectorImage } from "../src/components/DetectorImage";

// Mock fetch to return a minimal 1x1 white PNG (base64)
const TINY_PNG = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwADhQGAWjR9awAAAABJRU5ErkJggg==";

beforeEach(() => {
  global.fetch = vi.fn().mockResolvedValue({
    ok: true,
    blob: () => Promise.resolve(
      new Blob([Uint8Array.from(atob(TINY_PNG), c => c.charCodeAt(0))],
               { type: "image/png" })),
  } as Response);
  // Mock createImageBitmap
  global.createImageBitmap = vi.fn().mockResolvedValue({
    width: 1, height: 1, close: vi.fn(),
  } as unknown as ImageBitmap);
});

test("renders a canvas element", async () => {
  render(<DetectorImage exposureId={1} imagePath="/tmp/foo.tiff" size="full" />);
  await waitFor(() => expect(screen.getByRole("img", { hidden: true })).toBeInTheDocument());
});

test("shows placeholder when imagePath is null", () => {
  render(<DetectorImage exposureId={1} imagePath={null} size="full" />);
  expect(screen.getByTestId("detector-image-placeholder")).toBeInTheDocument();
});
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd packages/HimalayaUI/frontend
npm test -- DetectorImage
```

Expected: FAIL — component not found.

- [ ] **Step 3: Implement `DetectorImage.tsx`**

```tsx
import { useCallback, useEffect, useRef } from "react";

interface Props {
  exposureId: number;
  imagePath: string | null;
  size: "thumb" | "full";
  className?: string;
}

function getCssColor(varName: string): [number, number, number] {
  const raw = getComputedStyle(document.documentElement)
    .getPropertyValue(varName).trim();
  const c = document.createElement("canvas");
  c.width = c.height = 1;
  const ctx = c.getContext("2d")!;
  ctx.fillStyle = raw;
  ctx.fillRect(0, 0, 1, 1);
  const d = ctx.getImageData(0, 0, 1, 1).data;
  return [d[0], d[1], d[2]];
}

export function DetectorImage({ exposureId, imagePath, size, className }: Props): JSX.Element {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const renderImage = useCallback(async () => {
    const canvas = canvasRef.current;
    if (!canvas || !imagePath) return;

    const url = `/api/exposures/${exposureId}/image${size === "thumb" ? "?thumb=1" : ""}`;
    const res = await fetch(url);
    if (!res.ok) return;

    const blob = await res.blob();
    const bitmap = await createImageBitmap(blob);

    // Draw grayscale to offscreen canvas to read pixel data
    const off = new OffscreenCanvas(bitmap.width, bitmap.height);
    const offCtx = off.getContext("2d")!;
    offCtx.drawImage(bitmap, 0, 0);
    bitmap.close();
    const imageData = offCtx.getImageData(0, 0, bitmap.width, bitmap.height);

    // Build LUT: bg (intensity=0) → fg (intensity=255)
    const [br, bg, bb] = getCssColor("--color-bg");
    const [fr, fg, fb] = getCssColor("--color-fg");
    const data = imageData.data;
    for (let i = 0; i < data.length; i += 4) {
      const t = data[i] / 255;
      data[i]   = Math.round(br + t * (fr - br));
      data[i+1] = Math.round(bg + t * (fg - bg));
      data[i+2] = Math.round(bb + t * (fb - bb));
      data[i+3] = 255;
    }

    canvas.width  = bitmap.width;
    canvas.height = bitmap.height;
    canvas.getContext("2d")!.putImageData(imageData, 0, 0);
  }, [exposureId, imagePath, size]);

  useEffect(() => {
    renderImage();
  }, [renderImage]);

  // Re-apply colormap when theme changes (AppShell toggles class on <html>)
  useEffect(() => {
    const observer = new MutationObserver(() => renderImage());
    observer.observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["class"],
    });
    return () => observer.disconnect();
  }, [renderImage]);

  if (!imagePath) {
    return (
      <div
        data-testid="detector-image-placeholder"
        className={`flex items-center justify-center text-fg-muted text-xs ${className ?? ""}`}
      >
        No image
      </div>
    );
  }

  return (
    <canvas
      ref={canvasRef}
      role="img"
      aria-label="Detector image"
      className={`object-contain ${className ?? ""}`}
      style={{ imageRendering: "pixelated" }}
    />
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd packages/HimalayaUI/frontend
npm test -- DetectorImage
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/DetectorImage.tsx \
        packages/HimalayaUI/frontend/test/DetectorImage.test.tsx
git commit -m "feat(ui): DetectorImage canvas component with theme-aware LUT colormap"
```

---

## Task 11: `ThumbnailGallery` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ThumbnailGallery.tsx`
- Test: `packages/HimalayaUI/frontend/test/ThumbnailGallery.test.tsx`

- [ ] **Step 1: Write failing test**

Create `packages/HimalayaUI/frontend/test/ThumbnailGallery.test.tsx`:

```tsx
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { ThumbnailGallery } from "../src/components/ThumbnailGallery";
import type { Exposure } from "../src/api";

const makeExposure = (overrides: Partial<Exposure>): Exposure => ({
  id: 1, sample_id: 1, filename: "pos1.dat", kind: "file",
  selected: false, status: null, image_path: null, tags: [], sources: [],
  ...overrides,
});

test("dims rejected exposures", () => {
  const exposures = [
    makeExposure({ id: 1, filename: "good.dat" }),
    makeExposure({ id: 2, filename: "bad.dat", status: "rejected" }),
  ];
  render(
    <ThumbnailGallery
      exposures={exposures}
      selectedId={1}
      onSelect={vi.fn()}
      columns={2}
    />
  );
  const rejected = screen.getByTestId("thumb-cell-2");
  expect(rejected).toHaveClass("opacity-40");
});

test("shows indexing chip on selected=true exposure", () => {
  const exposures = [makeExposure({ id: 1, selected: true })];
  render(
    <ThumbnailGallery
      exposures={exposures}
      selectedId={1}
      onSelect={vi.fn()}
      columns={1}
    />
  );
  expect(screen.getByText("⊙ Indexing")).toBeInTheDocument();
});

test("calls onSelect when thumbnail clicked", async () => {
  const onSelect = vi.fn();
  const exposures = [makeExposure({ id: 5, filename: "pos5.dat" })];
  render(
    <ThumbnailGallery
      exposures={exposures}
      selectedId={undefined}
      onSelect={onSelect}
      columns={1}
    />
  );
  await userEvent.click(screen.getByTestId("thumb-cell-5"));
  expect(onSelect).toHaveBeenCalledWith(5);
});
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd packages/HimalayaUI/frontend
npm test -- ThumbnailGallery
```

Expected: FAIL — component not found.

- [ ] **Step 3: Implement `ThumbnailGallery.tsx`**

```tsx
import type { Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";

interface Props {
  exposures: Exposure[];
  selectedId: number | undefined;
  onSelect: (id: number) => void;
  /** Number of columns: 1 for horizontal strip, 2 for three-col layout */
  columns: 1 | 2 | "auto";
  className?: string;
}

export function ThumbnailGallery({
  exposures, selectedId, onSelect, columns, className,
}: Props): JSX.Element {
  const gridClass =
    columns === 1
      ? "flex flex-row gap-2 overflow-x-auto"
      : columns === 2
      ? "grid grid-cols-2 gap-2"
      : "grid grid-cols-[repeat(auto-fill,minmax(72px,1fr))] gap-2";

  return (
    <div className={`${gridClass} ${className ?? ""}`}>
      {exposures.map((e) => {
        const isViewing  = e.id === selectedId;
        const isRejected = e.status === "rejected";
        const isIndexing = e.selected;

        /* border-radius: thumb is rounded-md (8px). chip offset is 5px.
           chip border-radius = 8 - 5 = 3px → rounded-[3px].
           Keep in sync if rounded-md or offset changes. */
        return (
          <div
            key={e.id}
            data-testid={`thumb-cell-${e.id}`}
            onClick={() => onSelect(e.id)}
            className={[
              "relative flex flex-col items-center gap-1 cursor-pointer",
              "shrink-0",
              isRejected ? "opacity-40" : "",
            ].join(" ")}
          >
            <div
              className={[
                "relative w-full overflow-hidden rounded-md",
                "aspect-[3/4]",
                isViewing ? "ring-2 ring-accent" : "",
              ].join(" ")}
            >
              <DetectorImage
                exposureId={e.id}
                imagePath={e.image_path}
                size="thumb"
                className="w-full h-full"
              />
              {isIndexing && (
                <span
                  className="absolute top-[5px] left-[5px] rounded-[3px]
                             bg-accent/85 text-bg text-[9px] font-semibold
                             px-1.5 py-0.5 leading-snug backdrop-blur-sm"
                >
                  ⊙ Indexing
                </span>
              )}
            </div>
            <span className="text-[9px] text-fg-muted truncate w-full text-center">
              {e.filename?.replace(/\.dat$/i, "") ?? `#${e.id}`}
            </span>
          </div>
        );
      })}
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd packages/HimalayaUI/frontend
npm test -- ThumbnailGallery
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ThumbnailGallery.tsx \
        packages/HimalayaUI/frontend/test/ThumbnailGallery.test.tsx
git commit -m "feat(ui): ThumbnailGallery with indexing chip and reject dimming"
```

---

## Task 12: `DetectorImageCard` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/DetectorImageCard.tsx`
- Test: `packages/HimalayaUI/frontend/test/DetectorImageCard.test.tsx`

- [ ] **Step 1: Write failing test**

Create `packages/HimalayaUI/frontend/test/DetectorImageCard.test.tsx`:

```tsx
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { DetectorImageCard } from "../src/components/DetectorImageCard";
import type { Exposure } from "../src/api";

vi.mock("../src/components/DetectorImage", () => ({
  DetectorImage: () => <canvas data-testid="mock-detector-image" />,
}));

const makeExposure = (overrides: Partial<Exposure>): Exposure => ({
  id: 1, sample_id: 1, filename: "pos1.dat", kind: "file",
  selected: false, status: null, image_path: "/tmp/foo.tiff",
  tags: [], sources: [], ...overrides,
});

test("shows Accept and Reject buttons", () => {
  render(
    <DetectorImageCard
      exposure={makeExposure({})}
      onSetStatus={vi.fn()}
      onSetIndexing={vi.fn()}
      onAddTag={vi.fn()}
    />
  );
  expect(screen.getByRole("button", { name: /accept/i })).toBeInTheDocument();
  expect(screen.getByRole("button", { name: /reject/i })).toBeInTheDocument();
});

test("Reject shows rejection note input after click", async () => {
  render(
    <DetectorImageCard
      exposure={makeExposure({})}
      onSetStatus={vi.fn()}
      onSetIndexing={vi.fn()}
      onAddTag={vi.fn()}
    />
  );
  await userEvent.click(screen.getByRole("button", { name: /reject/i }));
  expect(screen.getByPlaceholderText(/reason/i)).toBeInTheDocument();
});

test("Use for indexing is disabled when rejected", () => {
  render(
    <DetectorImageCard
      exposure={makeExposure({ status: "rejected" })}
      onSetStatus={vi.fn()}
      onSetIndexing={vi.fn()}
      onAddTag={vi.fn()}
    />
  );
  expect(screen.getByRole("button", { name: /indexing/i })).toBeDisabled();
});
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd packages/HimalayaUI/frontend
npm test -- DetectorImageCard
```

Expected: FAIL — component not found.

- [ ] **Step 3: Implement `DetectorImageCard.tsx`**

```tsx
import { useState } from "react";
import type { Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";

interface Props {
  exposure: Exposure;
  onSetStatus: (status: "accepted" | "rejected" | null) => void;
  onSetIndexing: () => void;
  onAddTag: (key: string, value: string) => void;
}

export function DetectorImageCard({
  exposure, onSetStatus, onSetIndexing, onAddTag,
}: Props): JSX.Element {
  const [rejectMode, setRejectMode] = useState(false);
  const [rejectNote, setRejectNote] = useState("");

  const isRejected  = exposure.status === "rejected";
  const isAccepted  = exposure.status === "accepted";
  const isIndexing  = exposure.selected;

  function handleRejectConfirm() {
    onSetStatus("rejected");
    if (rejectNote.trim()) {
      onAddTag("rejection_reason", rejectNote.trim());
    }
    setRejectMode(false);
    setRejectNote("");
  }

  const existingNote = exposure.tags.find((t) => t.key === "rejection_reason")?.value;

  return (
    <div className="card flex flex-col h-full min-h-0 p-3 gap-3">
      <div className="text-xs text-fg-muted truncate shrink-0">
        {exposure.filename ?? `Exposure #${exposure.id}`}
      </div>

      {/* Portrait image — takes available height */}
      <div className="flex-1 min-h-0 flex items-center justify-center">
        <DetectorImage
          exposureId={exposure.id}
          imagePath={exposure.image_path}
          size="full"
          className="max-h-full max-w-full object-contain"
        />
      </div>

      {/* Rejection note (read mode) */}
      {isRejected && existingNote && !rejectMode && (
        <p className="text-[10px] text-fg-muted italic shrink-0">
          {existingNote}{" "}
          <button
            className="underline"
            onClick={() => { setRejectNote(existingNote); setRejectMode(true); }}
          >
            edit
          </button>
        </p>
      )}

      {/* Rejection note input */}
      {rejectMode && (
        <div className="flex gap-1 shrink-0">
          <input
            className="flex-1 text-xs bg-bg border border-border rounded px-2 py-1"
            placeholder="Rejection reason (flare, missed sample…)"
            value={rejectNote}
            onChange={(e) => setRejectNote(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && handleRejectConfirm()}
            autoFocus
          />
          <button
            className="text-xs px-2 py-1 border border-red-400 text-red-400 rounded"
            onClick={handleRejectConfirm}
          >
            Confirm
          </button>
        </div>
      )}

      {/* Controls */}
      <div className="flex flex-col gap-1.5 shrink-0">
        <button
          onClick={() => onSetStatus(isAccepted ? null : "accepted")}
          className={[
            "w-full text-xs py-1.5 rounded border transition-colors",
            isAccepted
              ? "border-green-400 text-green-400 bg-green-400/10"
              : "border-border text-fg-muted hover:border-green-400 hover:text-green-400",
          ].join(" ")}
        >
          ✓ Accept
        </button>

        <button
          onClick={() => {
            if (isRejected) { onSetStatus(null); }
            else { setRejectMode(true); }
          }}
          className={[
            "w-full text-xs py-1.5 rounded border transition-colors",
            isRejected
              ? "border-red-400 text-red-400 bg-red-400/10"
              : "border-border text-fg-muted hover:border-red-400 hover:text-red-400",
          ].join(" ")}
        >
          ✗ Reject
        </button>

        <button
          disabled={isRejected}
          onClick={onSetIndexing}
          className={[
            "w-full text-xs py-1.5 rounded border transition-colors",
            "disabled:opacity-40 disabled:cursor-not-allowed",
            isIndexing
              ? "border-accent text-accent bg-accent/10"
              : "border-border text-fg-muted hover:border-accent hover:text-accent",
          ].join(" ")}
        >
          ⊙ Use for indexing
        </button>
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd packages/HimalayaUI/frontend
npm test -- DetectorImageCard
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/DetectorImageCard.tsx \
        packages/HimalayaUI/frontend/test/DetectorImageCard.test.tsx
git commit -m "feat(ui): DetectorImageCard with accept/reject/indexing controls"
```

---

## Task 13: `SampleMetadataCard` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/SampleMetadataCard.tsx`
- Test: `packages/HimalayaUI/frontend/test/SampleMetadataCard.test.tsx`

- [ ] **Step 1: Write failing test**

Create `packages/HimalayaUI/frontend/test/SampleMetadataCard.test.tsx`:

```tsx
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { SampleMetadataCard } from "../src/components/SampleMetadataCard";
import type { Sample } from "../src/api";

const makeSample = (overrides: Partial<Sample> = {}): Sample => ({
  id: 1, experiment_id: 1, label: "D1",
  name: "DOPC 37C", notes: "run 2",
  tags: [
    { id: 10, key: "concentration", value: "10mM", source: "manifest" },
    { id: 11, key: "temp", value: "37C", source: "manual" },
  ],
  ...overrides,
});

test("renders name, notes, label, and tags", () => {
  render(
    <SampleMetadataCard
      sample={makeSample()}
      exposureSummary={{ total: 5, accepted: 3, rejected: 2 }}
      onUpdateSample={vi.fn()}
      onAddTag={vi.fn()}
      onRemoveTag={vi.fn()}
    />
  );
  expect(screen.getByDisplayValue("DOPC 37C")).toBeInTheDocument();
  expect(screen.getByText("D1")).toBeInTheDocument();
  expect(screen.getByText("concentration: 10mM")).toBeInTheDocument();
  expect(screen.getByText("5 exposures · 3 accepted · 2 rejected")).toBeInTheDocument();
});

test("manifest tags have no delete button", () => {
  render(
    <SampleMetadataCard
      sample={makeSample()}
      exposureSummary={{ total: 5, accepted: 3, rejected: 2 }}
      onUpdateSample={vi.fn()}
      onAddTag={vi.fn()}
      onRemoveTag={vi.fn()}
    />
  );
  // manifest tag chip should not have a delete button
  const manifestTag = screen.getByText("concentration: 10mM").closest("span")!;
  expect(manifestTag.querySelector("button")).toBeNull();
});

test("calls onUpdateSample on name blur", async () => {
  const onUpdate = vi.fn();
  render(
    <SampleMetadataCard
      sample={makeSample()}
      exposureSummary={{ total: 1, accepted: 1, rejected: 0 }}
      onUpdateSample={onUpdate}
      onAddTag={vi.fn()}
      onRemoveTag={vi.fn()}
    />
  );
  const nameInput = screen.getByDisplayValue("DOPC 37C");
  await userEvent.clear(nameInput);
  await userEvent.type(nameInput, "New Name");
  await userEvent.tab(); // blur
  expect(onUpdate).toHaveBeenCalledWith({ name: "New Name" });
});
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd packages/HimalayaUI/frontend
npm test -- SampleMetadataCard
```

Expected: FAIL — component not found.

- [ ] **Step 3: Implement `SampleMetadataCard.tsx`**

```tsx
import { useState } from "react";
import type { Sample } from "../api";

interface ExposureSummary { total: number; accepted: number; rejected: number; }

interface Props {
  sample: Sample;
  exposureSummary: ExposureSummary;
  onUpdateSample: (patch: { name?: string; notes?: string }) => void;
  onAddTag: (key: string, value: string) => void;
  onRemoveTag: (tagId: number) => void;
}

export function SampleMetadataCard({
  sample, exposureSummary, onUpdateSample, onAddTag, onRemoveTag,
}: Props): JSX.Element {
  const [name,  setName]  = useState(sample.name  ?? "");
  const [notes, setNotes] = useState(sample.notes ?? "");
  const [newTagKey, setNewTagKey]   = useState("");
  const [newTagVal, setNewTagVal]   = useState("");
  const [addingTag, setAddingTag]   = useState(false);

  function handleAddTag() {
    const k = newTagKey.trim(); const v = newTagVal.trim();
    if (k && v) { onAddTag(k, v); setNewTagKey(""); setNewTagVal(""); setAddingTag(false); }
  }

  return (
    <div className="card flex flex-col gap-3 p-3 overflow-y-auto">
      <p className="text-[10px] text-fg-muted">
        {exposureSummary.total} exposures
        {" · "}{exposureSummary.accepted} accepted
        {" · "}{exposureSummary.rejected} rejected
      </p>

      {/* Name */}
      <div className="flex flex-col gap-0.5">
        <label className="text-[10px] text-fg-muted uppercase tracking-wide">Name</label>
        <input
          className="w-full bg-bg border border-border rounded px-2 py-1 text-sm text-fg"
          value={name}
          onChange={(e) => setName(e.target.value)}
          onBlur={() => onUpdateSample({ name })}
        />
      </div>

      {/* Notes */}
      <div className="flex flex-col gap-0.5">
        <label className="text-[10px] text-fg-muted uppercase tracking-wide">Notes</label>
        <textarea
          rows={2}
          className="w-full bg-bg border border-border rounded px-2 py-1 text-sm text-fg resize-none"
          value={notes}
          onChange={(e) => setNotes(e.target.value)}
          onBlur={() => onUpdateSample({ notes })}
        />
      </div>

      {/* Label (read-only) */}
      {sample.label && (
        <div className="flex flex-col gap-0.5">
          <label className="text-[10px] text-fg-muted uppercase tracking-wide">Label</label>
          <span className="text-xs text-fg-muted bg-bg-subtle px-2 py-1 rounded">
            {sample.label}
          </span>
        </div>
      )}

      {/* Tags */}
      <div className="flex flex-col gap-1.5">
        <label className="text-[10px] text-fg-muted uppercase tracking-wide">Tags</label>
        <div className="flex flex-wrap gap-1">
          {sample.tags.map((tag) => (
            <span
              key={tag.id}
              className="inline-flex items-center gap-1 text-[10px] px-2 py-0.5 rounded-full
                         bg-bg-subtle border border-border text-fg-muted"
            >
              {tag.key}: {tag.value}
              {tag.source !== "manifest" && (
                <button
                  onClick={() => onRemoveTag(tag.id)}
                  className="text-fg-muted hover:text-red-400 leading-none"
                  aria-label={`Remove ${tag.key} tag`}
                >
                  ×
                </button>
              )}
            </span>
          ))}
          {!addingTag && (
            <button
              onClick={() => setAddingTag(true)}
              className="text-[10px] text-fg-muted hover:text-fg px-2 py-0.5 rounded-full
                         border border-dashed border-border"
            >
              + tag
            </button>
          )}
        </div>
        {addingTag && (
          <div className="flex gap-1">
            <input
              className="flex-1 text-xs bg-bg border border-border rounded px-1.5 py-0.5"
              placeholder="key"
              value={newTagKey}
              onChange={(e) => setNewTagKey(e.target.value)}
              autoFocus
            />
            <input
              className="flex-1 text-xs bg-bg border border-border rounded px-1.5 py-0.5"
              placeholder="value"
              value={newTagVal}
              onChange={(e) => setNewTagVal(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && handleAddTag()}
            />
            <button
              onClick={handleAddTag}
              className="text-xs px-2 border border-accent text-accent rounded"
            >
              Add
            </button>
            <button
              onClick={() => setAddingTag(false)}
              className="text-xs px-2 text-fg-muted"
            >
              ×
            </button>
          </div>
        )}
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd packages/HimalayaUI/frontend
npm test -- SampleMetadataCard
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/SampleMetadataCard.tsx \
        packages/HimalayaUI/frontend/test/SampleMetadataCard.test.tsx
git commit -m "feat(ui): SampleMetadataCard with inline name/notes editing and tag management"
```

---

## Task 14: `InspectPage` — full composition

Replace the stub with the real page.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/InspectPage.tsx`

- [ ] **Step 1: Implement `InspectPage.tsx`**

```tsx
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useAppState } from "../state";
import { useExposures, useSamples, useSetExposureStatus,
         useSelectExposure, useAddExposureTag,
         useAddSampleTag, useRemoveSampleTag, useUpdateSample } from "../queries";
import { ThumbnailGallery } from "../components/ThumbnailGallery";
import { DetectorImageCard } from "../components/DetectorImageCard";
import { SampleMetadataCard } from "../components/SampleMetadataCard";
import type { Exposure } from "../api";

export function InspectPage(): JSX.Element {
  const username     = useAppState((s) => s.username);
  const experimentId = useAppState((s) => s.activeExperimentId);
  const sampleId     = useAppState((s) => s.activeSampleId);
  const openModal    = useAppState((s) => s.openNavModal);

  // Auto-open nav modal if no sample selected (mirrors IndexPage behaviour)
  const autoOpenedRef = useRef(false);
  useEffect(() => {
    if (autoOpenedRef.current) return;
    if (username === undefined) return;
    if (experimentId === undefined) {
      autoOpenedRef.current = true;
      openModal("experiment");
    } else if (sampleId === undefined) {
      autoOpenedRef.current = true;
      openModal("sample");
    }
  }, [username, experimentId, sampleId, openModal]);

  const exposuresQ = useExposures(sampleId);
  const samplesQ   = useSamples(experimentId ?? 0);
  const exposures  = exposuresQ.data ?? [];
  const sample     = samplesQ.data?.find((s) => s.id === sampleId);

  // Local selected exposure (which one is shown in image card)
  const defaultSelectedId = useMemo(() => {
    const indexing = exposures.find((e) => e.selected);
    if (indexing) return indexing.id;
    const firstAccepted = exposures.find((e) => e.status === "accepted");
    if (firstAccepted) return firstAccepted.id;
    return exposures[0]?.id;
  }, [exposures]);

  const [viewingId, setViewingId] = useState<number | undefined>(undefined);
  useEffect(() => {
    if (viewingId === undefined && defaultSelectedId !== undefined) {
      setViewingId(defaultSelectedId);
    }
  }, [defaultSelectedId, viewingId]);

  const viewingExposure = exposures.find((e) => e.id === viewingId);

  const setStatus     = useSetExposureStatus(sampleId ?? 0);
  const setIndexing   = useSelectExposure(sampleId ?? 0);
  const addExpTag     = useAddExposureTag(sampleId ?? 0, viewingId ?? 0);
  const updateSample  = useUpdateSample(experimentId ?? 0, sampleId ?? 0);
  const addSampleTag  = useAddSampleTag(experimentId ?? 0, sampleId ?? 0);
  const rmSampleTag   = useRemoveSampleTag(experimentId ?? 0, sampleId ?? 0);

  const exposureSummary = useMemo(() => ({
    total:    exposures.length,
    accepted: exposures.filter((e) => e.status === "accepted").length,
    rejected: exposures.filter((e) => e.status === "rejected").length,
  }), [exposures]);

  const handleSetStatus = useCallback((status: "accepted" | "rejected" | null) => {
    if (!viewingId) return;
    setStatus.mutate({ exposureId: viewingId, status });
  }, [viewingId, setStatus]);

  const handleSetIndexing = useCallback(() => {
    if (!viewingId) return;
    setIndexing.mutate(viewingId);
  }, [viewingId, setIndexing]);

  const handleAddTag = useCallback((key: string, value: string) => {
    addExpTag.mutate({ key, value });
  }, [addExpTag]);

  if (!sample) return <div className="flex-1 min-h-0" />;

  return (
    <div
      data-testid="inspect-page"
      className="flex-1 min-h-0 flex flex-col gap-3 px-4 pb-4 pt-2"
    >
      {/*
        Breakpoints:
          < 1100px  → single column stacked
          1100–1400 → two columns 3fr/2fr (60/40): left=(metadata+gallery) right=image
          ≥ 1400px  → three columns 28fr/22fr/50fr: metadata | gallery | image
      */}
      <div
        className="
          min-h-0 grid gap-3 flex-1
          grid-cols-1
          min-[1100px]:grid-cols-[3fr_2fr]
          min-[1400px]:grid-cols-[28fr_22fr_50fr]
          h-auto min-[1100px]:h-[700px] min-[1100px]:max-h-[calc(100vh-120px)]
        "
      >
        {/* Column 1 — metadata (top on medium+, top on small) */}
        <section
          className="
            card min-h-[200px] min-[1100px]:min-h-0 overflow-hidden
            order-1
            min-[1100px]:row-span-1 min-[1400px]:col-start-1
          "
        >
          <SampleMetadataCard
            sample={sample}
            exposureSummary={exposureSummary}
            onUpdateSample={(patch) => updateSample.mutate(patch)}
            onAddTag={(k, v) => addSampleTag.mutate({ key: k, value: v })}
            onRemoveTag={(id) => rmSampleTag.mutate(id)}
          />
        </section>

        {/* Column 2 — thumbnail gallery */}
        <section
          className="
            card min-h-[180px] min-[1100px]:min-h-0 overflow-hidden p-2
            order-2
            min-[1100px]:row-start-2 min-[1100px]:col-start-1
            min-[1400px]:row-start-1 min-[1400px]:col-start-2
          "
        >
          {/* small: horizontal strip; medium: auto-fill grid; large: 2-col grid */}
          <ThumbnailGallery
            exposures={exposures}
            selectedId={viewingId}
            onSelect={setViewingId}
            columns={
              /* We can't read media queries in JS without a hook; pass "auto"
                 and let the gallery use auto-fill. For the large breakpoint the
                 parent grid constrains the column width enough that auto-fill
                 naturally produces 2 columns at 22fr. */
              "auto"
            }
            className="h-full"
          />
        </section>

        {/* Column 3 — detector image (full height, right column) */}
        <section
          className="
            card min-h-[300px] min-[1100px]:min-h-0 overflow-hidden
            order-3
            min-[1100px]:row-span-2 min-[1100px]:col-start-2
            min-[1400px]:row-span-1 min-[1400px]:col-start-3
          "
        >
          {viewingExposure ? (
            <DetectorImageCard
              exposure={viewingExposure}
              onSetStatus={handleSetStatus}
              onSetIndexing={handleSetIndexing}
              onAddTag={handleAddTag}
            />
          ) : (
            <div className="flex items-center justify-center h-full text-fg-muted text-sm">
              Select an exposure
            </div>
          )}
        </section>
      </div>
    </div>
  );
}
```

- [ ] **Step 2: Run TypeScript check**

```bash
cd packages/HimalayaUI/frontend
npm run build -- --noEmit 2>&1 | head -40
```

Expected: no errors.

- [ ] **Step 3: Run full Vitest suite**

```bash
cd packages/HimalayaUI/frontend
npm test
```

Expected: all tests pass (new and existing).

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/pages/InspectPage.tsx
git commit -m "feat(ui): InspectPage — responsive three-breakpoint layout with metadata, gallery, image view"
```

---

## Task 15: Index page filtering

Pass `excludeRejected: true` to `useExposures` in IndexPage so rejected exposures are hidden from the trace selector.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/IndexPage.tsx` (or wherever `useExposures` is called for Index)

- [ ] **Step 1: Find where `useExposures` is called for the Index workflow**

```bash
grep -rn "useExposures" packages/HimalayaUI/frontend/src/
```

- [ ] **Step 2: Add `excludeRejected: true`**

In whichever component calls `useExposures` for the active exposure selection (likely `PlotCard.tsx` or `App.tsx`), update the call:

```typescript
// Before:
const exposuresQ = useExposures(sampleId);

// After:
const exposuresQ = useExposures(sampleId, { excludeRejected: true });
```

- [ ] **Step 3: Run TypeScript check and Vitest**

```bash
cd packages/HimalayaUI/frontend
npm run build -- --noEmit 2>&1 | head -20
npm test
```

Expected: all pass.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/
git commit -m "feat(index): hide rejected exposures from Index page exposure selector"
```

---

## Task 16: E2E tests

**Files:**
- Create: `packages/HimalayaUI/frontend/e2e/inspect.spec.ts`

- [ ] **Step 1: Write E2E tests**

Create `packages/HimalayaUI/frontend/e2e/inspect.spec.ts`:

```typescript
import { test, expect } from "@playwright/test";

test.describe("Inspect page", () => {
  test.beforeEach(async ({ page }) => {
    // Mock exposures list
    await page.route("/api/samples/*/exposures", (route) =>
      route.fulfill({
        json: [
          {
            id: 1, sample_id: 1, filename: "pos1.dat", kind: "file",
            selected: true, status: "accepted", image_path: "/tmp/pos1.tiff",
            tags: [], sources: [],
          },
          {
            id: 2, sample_id: 1, filename: "pos2.dat", kind: "file",
            selected: false, status: "rejected", image_path: null,
            tags: [{ id: 99, key: "rejection_reason", value: "flare", source: "manual" }],
            sources: [],
          },
          {
            id: 3, sample_id: 1, filename: "pos3.dat", kind: "file",
            selected: false, status: null, image_path: null,
            tags: [], sources: [],
          },
        ],
      })
    );
    await page.route("/api/exposures/*/image*", (route) =>
      route.fulfill({ contentType: "image/png", body: Buffer.alloc(100) })
    );
    await page.route("/api/exposures/*/status", (route) => route.fulfill({ json: { id: 1, status: "accepted" } }));
    await page.route("/api/exposures/*/select", (route) => route.fulfill({ json: { id: 1, selected: true } }));
    await page.route("/api/exposures/*/tags", (route) => route.fulfill({ status: 201, json: { id: 100, key: "rejection_reason", value: "test", source: "manual" } }));
  });

  test("Inspect tab appears before Index tab", async ({ page }) => {
    await page.goto("/");
    const tabs = page.getByRole("tab");
    await expect(tabs.first()).toHaveText("Inspect");
    await expect(tabs.nth(1)).toHaveText("Index");
  });

  test("navigating to Inspect shows the page", async ({ page }) => {
    await page.goto("/");
    await page.click('[data-testid="tab-inspect"]');
    await expect(page.getByTestId("inspect-page")).toBeVisible();
  });

  test("rejected exposure cell is dimmed", async ({ page }) => {
    await page.goto("/");
    await page.click('[data-testid="tab-inspect"]');
    const rejectedCell = page.getByTestId("thumb-cell-2");
    await expect(rejectedCell).toHaveClass(/opacity-40/);
  });

  test("clicking Reject button shows note input", async ({ page }) => {
    await page.goto("/");
    await page.click('[data-testid="tab-inspect"]');
    await page.getByRole("button", { name: /reject/i }).click();
    await expect(page.getByPlaceholder(/reason/i)).toBeVisible();
  });

  test("clicking thumbnail loads it in the image card", async ({ page }) => {
    await page.goto("/");
    await page.click('[data-testid="tab-inspect"]');
    await page.getByTestId("thumb-cell-3").click();
    // The image card header should update to pos3.dat
    await expect(page.getByText("pos3.dat")).toBeVisible();
  });
});
```

- [ ] **Step 2: Run E2E tests**

```bash
cd packages/HimalayaUI/frontend
npm run e2e -- --grep "Inspect page"
```

Expected: all 5 tests pass.

- [ ] **Step 3: Run full test suite**

```bash
npm test && npm run e2e
```

Expected: all pass.

- [ ] **Step 4: Final build check**

```bash
npm run build
```

Expected: `tsc --noEmit` clean, Vite build completes.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/e2e/inspect.spec.ts
git commit -m "test(e2e): Inspect page — tab order, reject dimming, image card, note input"
```

---

## Verification

**Backend:**
```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

All tests pass including: migration, image route (200/404/thumb), status route (accept/reject/null/422), expose filter, TIFF discovery.

**Frontend unit:**
```bash
cd packages/HimalayaUI/frontend
npm test
```

All existing + new Vitest tests pass.

**Frontend E2E:**
```bash
cd packages/HimalayaUI/frontend
npm run e2e
```

Inspect page E2E + existing E2E pass.

**Manual:**
1. Start server, open Inspect tab — gallery loads with TIFF thumbnails
2. Click a thumbnail — image card updates
3. Click Reject — note input appears; save → cell dims
4. Switch to Index tab — rejected exposure absent from selector
5. Toggle theme — detector image colormap updates immediately without page reload
