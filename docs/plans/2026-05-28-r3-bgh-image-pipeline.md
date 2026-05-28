# R3 B + G + H — Detector image pipeline (ONE branch / ONE PR)

Branch: `r3-bgh-image-pipeline` (off `main` @ `111ebb8`)
Bundles **#255 (B)**, **#260 (G)**, **#261 (H)** + folds **R3-S06**.
Rationale for collapsing: all three touch `packages/HimalayaUI/src/image.jl` and all
three require the SAME `IMAGE_PROCESSING_VERSION` bump — which must happen **exactly
once** (a deliberate one-responsibility exception: a shared backend cache-version bump
cannot be split across three PRs without two of them shipping stale-byte bugs).

---

## The single version bump

- **Site:** `packages/HimalayaUI/src/image.jl:7` — `const IMAGE_PROCESSING_VERSION = "v2"`.
- **Grep confirms ONE definition site.** Other matches in `image.jl` (lines 11, 20, 31)
  are docstring mentions + the `image_version_token` interpolation; they read the const.
- **New value:** `"v2"` -> **`"v3"`** (the value G's issue sketch already names).
- **Why once covers all three:** B (warm LUT — see finding below), G (full-path 1536px
  cap), H (thumb downscale-before-lognorm produces slightly different percentile clips)
  all change rendered bytes / display. All invalidate via the same `?v=<token>` key;
  `image_version_token` = `"v3-<mtime>"` after the bump, so every cached PNG (browser AND
  the new disk cache) misses and re-renders. There is one const; it can only bump once.

---

## CRITICAL design finding — U-1 LUT lives in the FRONTEND, not (only) the backend

**This is the load-bearing decision in the bundle and the reason both project reviewers
are attached at the plan gate. Read before approving.**

The U-1 issue body says "warm the LUT in `image.jl`; return `Matrix{RGB}` instead of
`Gray`." Verified against live source, that change **alone has zero visible effect**:

`DetectorImage.tsx:138-152` already runs its OWN look-up table in the browser. It reads
the decoded PNG's **R channel only** (`const t = data[i] / 255`) as a scalar intensity,
then rebuilds every pixel by interpolating between two CSS variables:

```
t=0  -> --color-bg  (paper,  oklch 0.978 0.006 85 — near-white)
t=1  -> --color-fg  (ink,    oklch 0.265 0.013 68 — dark)
```

So today the detector renders as **dark ink on light paper** (bright rings -> dark ink,
background -> paper). Any chroma the backend bakes into the PNG is discarded — only the R
channel survives, immediately remapped to the CSS-var ramp. A backend `Gray`->`RGB`
change round-trips through `data[i]` unchanged and repaints identically.

U-1's actual intent (DESIGN.md §5: "a dark window: `frame-edge` background, set into the
paper") is that the detector should read as **a dark warm window with bright signal**, not
ink-on-paper. That is a **frontend LUT-endpoint change**, not a backend pixel change.

### Decision (primary; flagged for reviewer sign-off)

**Warm the LUT at the frontend endpoints** — the lower-risk "warm-tinted" variant the
issue calls option (a):

- `t=0` (background / window backing) -> **`--color-frame-edge`** (oklch 0.150 0.010 55 —
  the warm near-black DESIGN.md already ships for the detector window).
- `t=1` (bright signal) -> a **warm-bright** endpoint (warm off-white, hue ~80, e.g. a new
  `--color-frame-signal` token ~ oklch 0.92 0.02 80, or reuse `frame-tag` 0.82 0.01 80 as
  the bright end if a separate token is judged unwarranted). Bright rings on warm
  near-black = "ink on dark paper" / a Print restatement of the window.

This makes the detector image a warm dark window that **bridges to** the surrounding
`frame-edge` border (`ContactSheetRow`, `LoupeFrame`, `FocusDetectorPanel` all already
wrap the canvas in `border-frame-edge bg-frame-edge`), instead of an inverted ink-on-paper
patch sitting cold inside a dark frame.

Endpoints come from CSS vars via the existing `getCssColor()` helper — no hard-coded hex,
theme-token-driven, consistent with the milestone's token discipline.

### Why still bump the backend version + still touch `image.jl`

Even though U-1's *visible* fix is frontend, the bundle still bumps
`IMAGE_PROCESSING_VERSION` because **G and H change the actual PNG bytes** and the LUT
endpoint change means already-cached client bitmaps should re-decode. The bump is the
shared invalidator. `image.jl` is also where G's full-path cap and H's thumb variant +
disk cache live, so the file is touched regardless.

### Reviewer questions to resolve at this gate

1. **Frontend-LUT vs backend-RGB**: confirm warming the frontend endpoints is the faithful
   U-1 fix (vs. backend RGB, which the frontend would discard). saxs-physics-reviewer:
   does a `frame-edge -> warm-bright` ramp preserve ring legibility / contrast as well as
   the current paper->ink ramp? (The percentile clip math is unchanged; only the display
   colour ramp moves.)
2. **New token vs reuse**: introduce `--color-frame-signal` (bright endpoint) or reuse
   `frame-tag`? DESIGN.md currently has no "bright signal" token; adding one is a small §2
   addition. himalaya-reviewer to weigh DESIGN.md drift vs token reuse.
3. **Backend RGB at all?** Recommendation: **do NOT** change `load_and_lognormalize`'s
   return type to RGB — it's dead chroma under the frontend LUT and would break
   `test_image.jl`'s `channelview` grayscale assertions for no visible gain. Keep the
   backend returning `Gray`; warm purely in the frontend. (A backend approach for some
   future no-frontend-LUT consumer would be a separate follow-up.)

---

## B / #255 — Detector window warmth (LUT + skeleton + placeholder + R3-S06)

### B.1 — Warm the detector LUT (U-1) — FRONTEND
- **File:** `packages/HimalayaUI/frontend/src/components/DetectorImage.tsx:138-148`.
- Change LUT endpoints from `--color-bg`/`--color-fg` to `--color-frame-edge` (t=0) /
  warm-bright (t=1, per reviewer token decision). Keep the R-channel-as-intensity read
  (`data[i]/255`); backend stays `Gray`.
- **Test (`DetectorImage.test.tsx`):** existing tests stub `getImageData` ->
  `Uint8ClampedArray(4)` (zeros) so they don't assert colour. Add a focused test mocking a
  non-trivial `getImageData` (a mid + a max pixel) and a `getComputedStyle` stub returning
  known frame-edge / signal values, asserting the written buffer interpolates between the
  warm endpoints (t=0 -> frame-edge RGB, t=255 -> signal RGB). Behaviour assertion on the
  LUT math, not a class-string test.

### B.2 — Boneyard skeletons fit The Print (U-2, non-detector)
- **Root cause (verified):** `main.tsx:17-22` `configureBoneyard({ color:"rgba(30,31,38,1)", ... })`
  — a dark near-black global bone fill = the "black blocks on paper" dark-era leftover.
  `boneyard.config.json` mirrors it (capture-time only).
- **Fix:** change the global `color` to the **`paper-sunk`** equivalent (oklch 0.951 0.008
  84). boneyard takes a CSS color string and `configureBoneyard` runs before any element
  exists, so a `var()` may not resolve — use an `oklch()`/computed literal matching
  `--color-paper-sunk`. Update `boneyard.config.json` in lockstep (per the main.tsx
  comment). Keep `animate:"pulse"` (tuned via `.bone{animation-duration:1.8s}`).

### B.3 — Detector-specific skeleton shimmer (U-2, detector)
- **Files:** `FocusDetectorPanel.tsx:159` `<Skeleton name="focus-detector">` and any other
  detector-window skeleton. Use a **`frame-edge`** bone fill (the per-component `color`
  prop overrides the global default per boneyard's API) — a detector window is dark even
  while loading, so its skeleton shimmers dark, matching the live `bg-frame-edge` window.
- **Test:** detector skeleton renders bones with the frame-edge fill — assert via the bone
  element's inline style / `data-*`, never a Tailwind class string.

### B.4 — Missing-image placeholder + R3-S06
- **File:** `DetectorImage.tsx:204-213` (the `if (!imagePath)` branch; R3-S06 = the
  `text-fg-muted text-xs` on line 208).
- **Fix:** render the placeholder **inside a `frame-edge` window** with **`frame-tag`**
  caption text (`text-frame-tag`): `bg-frame-edge` window + `text-frame-tag` "No image"
  caption, matching the live detector treatment. Removes the last `text-fg-muted` survivor
  (R3-S06) — **no surviving `text-fg-*`/`bg-bg-*` class in the file** (#255 Done-when).
- **Test (`DetectorImage.test.tsx`):** existing `shows placeholder when imagePath is null`
  stays green (keeps `data-testid="detector-image-placeholder"`). Add an assertion the
  placeholder carries the frame-edge window + frame-tag caption via `data-*`
  (e.g. `data-variant="frame-window"`), and assert the rendered class list excludes
  `text-fg-muted`.

---

## G / #260 — Cap full-image PNG at 1536px (backend)
- **File:** `routes_exposures.jl:48-55`.
- **Fix:** in the `else` (full) branch, `img = resize_to_fit(img, 1536)` **after**
  `load_and_lognormalize(path)` (percentile-clip math runs on full res — science view
  stays accurate; only the display raster is capped). `resize_to_fit` no-ops when
  `max(h,w)<=max_px`, so sub-1536 detectors are untouched.
- `Cache-Control: private, max-age=31536000, immutable` stays; `?v=v3-...` invalidates.
- **Test (`test_routes_image.jl`):** request `size=full` (no `thumb`) on a >1536px fixture
  TIFF, assert returned PNG max side <= 1536; regression that a <=1536 fixture returns
  unresized. Synthesize TIFFs in-test (the `test_image.jl` Q0f31 pattern, scaled up).

---

## H / #261 — Thumbnail render perf (backend)

### H.1 — `load_and_lognormalize_thumb(path, max_px=128)` — downscale-first
- **File:** `image.jl`. Load TIFF, recover Q0f31 counts
  (`reinterpret.(Int32, channelview(raw))`, `max(.,0)`), build `Gray` view,
  **`resize_to_fit` to `max_px` FIRST**, then `log1p`+percentile-clip+clamp/normalize on
  the ~16K-pixel downscaled matrix. Order-of-magnitude per-thumb speedup. Clip points
  differ slightly (downsampled data) — perceptually identical for 128px; documented in the
  docstring. Keep `load_and_lognormalize` unchanged for G's capped-full path.
- **Test (`test_image.jl`):** `load_and_lognormalize_thumb` on the Q0f31 fixture returns a
  <=128px image, values in `[0,1]`, bright -> ~1, dead/background -> 0. Floor-style
  assertions (`~0`,`~1`,`0<mid<1`), not exact clip points.

### H.2 — Disk-backed PNG cache
- **Cache dir:** `<db_dir>/cache/thumb-128/`, `db_dir = dirname(current_db().file)`.
  **Verified:** `SQLite.DB` has a `.file` field — a real path for disk DBs, `":memory:"`
  for in-memory test DBs. Helper **must guard `:memory:`**: when the file is `":memory:"`
  (or empty), skip the disk cache and always recompute. Keeps the `:memory:`/unit fixtures
  correct; cache-hit/miss tests use an explicit temp-file DB
  (`open_db(joinpath(tmp,"himalaya.db"))`).
- **Key:** `thumb-128/{exposure_id}-{image_version_token(path)}.png`. Token already encodes
  `IMAGE_PROCESSING_VERSION` ("v3") + TIFF mtime -> invalidates on TIFF rewrite or version
  bump (stale `v2-*` files never read again).
- **Helpers in `image.jl`:** `thumb_cache_dir(db)->Union{String,Nothing}` (nothing for
  `:memory:`); `thumb_cache_path(db, exposure_id, token)`;
  `ensure_thumb_cached(db, exposure_id, path)->Vector{UInt8}` — read-if-present, else
  render (`load_and_lognormalize_thumb`+`encode_png`), `mkpath`, **atomic write** (temp
  name + `mv`, so prewarm threads / concurrent requests never read a half-written file),
  return bytes. On `:memory:`/no-dir: render + return, no write.
- **Route wiring:** `is_thumb` branch -> `ensure_thumb_cached`; full branch ->
  `load_and_lognormalize` + `resize_to_fit(.,1536)` (G).
- **Tests (temp-file DB):** (a) miss writes the file; (b) second call serves from disk
  (delete the source TIFF after the first call; assert the second still returns bytes);
  (c) token change writes a NEW file, leaves the stale one; (d) missing-TIFF skip via H.3.

### H.3 — `prewarm_thumbnails!(db; threads=true)` wired into init + reingest
- **File:** `image.jl` (fn) + `cli.jl` + `pipeline.jl`/`cli.jl` wiring.
- `SELECT id, image_path FROM exposures WHERE image_path IS NOT NULL`; for each present
  TIFF `ensure_thumb_cached`; **skip-with-`@info` on missing TIFF** (matches the route's
  graceful 404). Thread-parallel `Threads.@threads`; `Atomic{Int}` counters.
- **Wire AFTER ingest completes** (data-validation failures abort before slow image work):
  - `cli_init_with_db!`: call after `_cli_init_inner!` returns, **outside** the init tx
    (a TIFF-render failure must not roll back a good ingest).
  - reingest: call from `reingest!` **after** the `lock/SQLite.transaction` block commits
    (`_reingest_inner!` runs inside the tx; the prewarm goes after). The route
    `POST /api/experiments/:id/reingest` calls `reingest!`, so CLI + HTTP both prewarm.
  - **Concurrency:** prewarm only READS the DB (one SELECT) then writes the FS cache — no
    `_DB_WRITE_LOCK` needed. The SELECT materializes via `Tables.rowtable` BEFORE the
    threaded loop, so no statement-cache race across threads.
- **Tests (temp-file DB + synthetic TIFFs):** init/reingest populates `cache/thumb-128/`
  (#files == #exposures-with-image); missing-TIFF exposure skipped, ingest still succeeds;
  the **read-only-experiment-dir snapshot test in `test_pipeline.jl` stays green** — the
  cache dir is under `db_dir`, NOT `exp_dir`, so the invariant holds (**confirm
  `db_dir != exp_dir` in the fixture**). Run `threads=false` for deterministic asserts +
  one `threads=true` smoke run.

### H.4 — gitignore the cache
- Cache at `<db_dir>/cache/`. Default `db_dir` is `~/.himalaya` or `HIMALAYA_DB_PATH`'s
  dir; the repo's deployment `data/` is already gitignored. Add `cache/` under the `data/`
  block in root `.gitignore` defensively; tests use `mktempdir` so nothing leaks into the
  tracked tree.

---

## R3-S06 — folded into B.4 (no separate work)
Last legacy `text-fg-muted` survivor (`DetectorImage.tsx:208`; issue lists `:197`) is
removed when B.4 redesigns the placeholder as a `frame-edge` window with `frame-tag` text.

---

## TDD task list (one PR, grouped by concern, landing together)
Order: version bump + cache key in place before byte-changing tests.

1. **Bump `IMAGE_PROCESSING_VERSION` "v2"->"v3"** (image.jl:7). Commit alone (shared
   invalidator; everything keys off it).
2. **G:** test full-path <=1536 cap (red) -> `resize_to_fit(img,1536)` in route else-branch.
3. **H.1:** test `load_and_lognormalize_thumb` (red) -> downscale-first variant (green);
   route thumb branch keeps existing path until H.2.
4. **H.2:** disk-cache miss-write / hit-read / token-rotate tests (red) ->
   `thumb_cache_dir`/`thumb_cache_path`/`ensure_thumb_cached` (+ `:memory:` guard) + route
   thumb rewire (green).
5. **H.3:** prewarm populates cache + skips missing TIFF + ingest ok (red) ->
   `prewarm_thumbnails!` + wire into `cli_init_with_db!` and `reingest!` (green).
6. **H.4:** `.gitignore` `cache/`.
7. **B.1:** DetectorImage LUT-endpoint warmth test (red) -> swap endpoints to
   frame-edge/signal (green). (Pending reviewer token decision — may add styles.css token.)
8. **B.4 / R3-S06:** placeholder frame-edge window + frame-tag text test (red) -> redesign,
   remove `text-fg-muted` (green).
9. **B.2:** flip `configureBoneyard` global `color` to paper-sunk (main.tsx +
   boneyard.config.json). Covered by existing skeleton tests staying green + a targeted
   bone-fill assertion if feasible.
10. **B.3:** detector `<Skeleton color=frame-edge>` override + test.
11. Update `packages/HimalayaUI/src/AGENTS.md` "Image route" section: it currently still
    claims `Cache-Control: no-store` (STALE — live code is `immutable, max-age`). Document
    the new thumb disk cache + version key + prewarm. (Doc accuracy, no behaviour.)

Each numbered item = its own commit (TDD red->green where a test exists).

---

## Verify gate (single capture each — slow-suite discipline)
- **Julia (capture ONCE, grep the file — never re-run per grep):**
  ```
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' \
    > /Users/me/.claude/jobs/67cdbfb1/tmp/jl-test.out 2>&1
  grep -E "Test Summary|did not pass|fail|Error" .../jl-test.out ; tail -50 .../jl-test.out
  ```
  Core `Pkg.test()` is unaffected (only HimalayaUI image paths change) — skip unless a
  cross-dep surfaces.
- **Frontend (the DetectorImage / skeleton half):**
  ```
  (cd packages/HimalayaUI/frontend && npm test)        # Vitest one-shot
  (cd packages/HimalayaUI/frontend && npm run build)   # tsc --noEmit + vite build
  ```
- **Image-byte staleness guard:** no committed golden/cached PNG fixtures in-tree
  (verified — tests synthesize TIFFs via `mktempdir`). The disk cache keys on "v3", so no
  `v2-*` artifact can leak. Nothing to regenerate.

## Visual Done-when verification
1. **Unit (deterministic):** B.1 LUT-math test asserts the canvas buffer interpolates
   frame-edge -> signal; B.4 test asserts the placeholder's frame-edge/frame-tag tokens via
   `data-*`. Durable regression guards, no screenshot flake.
2. **Build:** `npm run build` proves the token refs + `getCssColor` calls type-check and
   the styles.css token (if added) resolves.
3. **Operator eyeball (out-of-band, noted in PR):** U-1/U-2 were operator-raised from the
   running app (screenshots couldn't surface them); the PR body calls out a running-app
   eyeball on contact-sheet + loupe as final visual sign-off, with the unit tests as the
   standing regression floor. Visual won't be claimed "done" from tests alone.

## Risks / watch-items
- **U-1 frontend-vs-backend** (the big one — reviewer gate). If reviewers insist on
  backend-RGB, B.1 changes scope and I'll re-plan that branch before implementing.
- **`:memory:` cache guard** — most likely green-locally/red-in-CI split if forgotten.
  Pinned as an explicit H.2 requirement + test.
- **Read-only-experiment-dir snapshot** (`test_pipeline.jl`) — cache writes go to `db_dir`,
  not `exp_dir`; fixtures keep them disjoint.
- **Prewarm placement** — OUTSIDE the ingest transaction so a TIFF-render error can't roll
  back a committed ingest.
