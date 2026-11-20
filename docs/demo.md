# Demo: run the real app on synthetic data

How to run HimalayaUI end to end (Julia backend + Vite frontend) on a synthetic
SAXS experiment, locally or in a Claude Code cloud session, without any lab data.
Used for UI walkthroughs and Playwright screenshots.

## Quick start

```bash
make setup        # once per checkout: npm ci + Pkg.develop/instantiate (≈5 min cold)
make demo         # generate data if missing, start backend :8080, ingest, start Vite :5173
# open http://127.0.0.1:5173/   (Ctrl-C stops both servers)

make demo-empty   # same, but on a FRESH empty DB with nothing ingested, so the
                  # "+ New experiment" funnel can be walked from the start; point it
                  # at the absolute path of scratch/demo-experiment
make demo-data    # (re)generate the synthetic experiment only
make demo-clean   # delete the generated data and demo DBs
```

Overrides: `DEMO_DIR=/tmp/demo make demo` (data location), `DEMO_DB_DIR` (demo DB + logs, default `scratch/demo`), `DEMO_API_PORT` /
`DEMO_UI_PORT` (ports, for running beside another backend; see
[frontend-dev-loop.md](frontend-dev-loop.md)).

`make demo` keeps its DB in `scratch/demo/` between runs (edits persist; ingest is
skipped when the experiment is already there, and a scan that was interrupted or
failed is re-run). `make demo-empty` wipes
`scratch/demo-empty/` on every start. Backend logs go next to each DB.

## The synthetic data

`packages/HimalayaUI/scripts/make_demo_data.jl` writes an SSRL-style tree that the
funnel's layout auto-detect (`resolve_experiment_layout`, routes_fs.jl) recognises:

```
scratch/demo-experiment/
  data/<label>_S<scan>_0_001.tif     Pilatus 1M frame, signed Int32 counts, -1 in module gaps
  data/<label>_S<scan>_0_001.prp     PRP sidecar (energy, pipe length, horizontal position, time…)
  analysis/setup_info_*.txt          beam center + AgBe distance (geometry source)
  analysis/automatic_analysis/tot_files/<label>_S<scan>_tot.dat   q  I  sigma
  TRUTH.csv                          per sample: scans, generated phases + lattice, expected
                                     auto call (`auto`), phase-tagged reflection q values,
                                     scans to drop (`drop`)
```

- **Not real data, never committed.** It lives under gitignored `scratch/` (or
  wherever `DEMO_DIR` points). Only the generator is in git.
- **Peak positions come from `Himalaya.phaseratios`**, so the demo tracks any
  change to the core ratio tables (e.g. the √11 removal in #304).
- 12 samples in two capillary holders (two loads, 75 min apart), one sample per
  slot (3.95 mm pitch, so the slot clustering in grouping.jl splits them), 1–3
  scans each (28 exposures): clean Pn3m ×2, Pn3m+Im3m, Im3m, Ia3d, HII, Lα,
  Lα+HII, weak Pn3m, Pn3m+HII, a buffer blank, a ULV form-factor-only trace, and
  a leaked capillary. Each scan is one frame (`_0_001`, the frame the funnel's
  layout detection keys on).
- **Some scans are bad on purpose**, so the drop and representative rules have
  something to act on. TRUTH's `drop` column lists them (`BAD_SCANS` in the
  generator):
  - `MO_H2O`: the newest scan is a beam-off shot: a flat dark floor of about
    0.01 counts per pixel plus hot pixels, so Loupe shows a black frame.
  - `DOPE_H2O`: the oldest scan missed the capillary: steep q⁻⁴ air and slit
    scatter over a faint air floor, so Loupe shows a glow at the beamstop and
    darkness beyond it.
  - `leaked_capillary`: both scans missed the capillary in the same way, so the
    sample has no usable frame.

  Nothing drops them at ingest; a triager does. `bin/verify triage`
  drops them through the endpoint Loupe's Drop uses, so `bin/verify check-truth` can check the corpus
  after triage.
- **Expected calls are not always the generated truth**, by design. TRUTH's `auto`
  column is what Himalaya's own pipeline calls. Two samples differ:
  - `MO_PEG5`: Pn3m √6 overlaps Im3m √10, so `auto_group` keeps Im3m. Pn3m stays
    a candidate, as the second hypothesis to promote by hand.
  - `MO_AMP_weak`: too weak to index automatically.

  `test/test_demo_data.jl` holds a regression floor on these calls over three
  seeds. If a change to peak finding, scoring or `phaseratios` moves them, fix
  the change or retune the sample deliberately.
- Seeded (`--seed`, default 7) and deterministic within a Julia release; ~100 MB,
  ~10 s. Traces and detector frames draw from separate random streams, so
  `--no-images` (empty `.tif` placeholders, ~1 MB) writes the same `.dat` files.
  `make demo` refuses such a tree; regenerate with `make demo-data`.
- Safety: the generator refuses to write into a non-empty directory that has no
  `TRUTH.csv`, and `make demo-clean` only deletes `DEMO_DIR` when it has one.
- TIFFs are written as `Gray{Q0f31}`, which TiffImages stores as a signed 32-bit
  integer TIFF, the same on-disk format as the Pilatus. Checked both ways: a
  plain int32 TIFF from Python `tifffile` loads in Julia as `Gray{Q0f31}` with
  the counts intact, and the generator's files read back in `tifffile` as int32.

`make demo` ingests through `scripts/demo_ingest.jl`, which makes the same three
calls the funnel does (`/api/fs/resolve` → `/api/fs/manifest` → `POST /api/experiments`)
and waits for the first scan. The first browser visit still shows the
"Who are you?" onboarding; any name works.

## Cloud sessions (claude.ai/code)

The app needs Julia (compat says ≥ 1.9; tested on 1.13), Node ≥ 20 and npm. A cloud environment that runs it needs:

1. **Network:** unrestricted, or at least `julialang.org`, `*.julialang.org`
   (installer, package server, registry), `github.com` and the npm registry.
   The stock "Default" environment blocks the Julia hosts.
2. **Setup script** that installs Julia, e.g.:

   ```bash
   #!/usr/bin/env bash
   set -uo pipefail
   if [ ! -x /root/.juliaup/bin/juliaup ]; then
     curl -fsSL https://install.julialang.org | sh -s -- --yes || exit 1
   fi
   ln -sf /root/.juliaup/bin/julia   /usr/local/bin/julia
   ln -sf /root/.juliaup/bin/juliaup /usr/local/bin/juliaup
   julia --version
   ```

   Node, npm and Playwright's Chromium are already in the cloud image. Keep the
   script machine-level; run `make setup` in the session (it needs the repo).

Then, in the session: `make setup && make demo`. Playwright can drive
`http://127.0.0.1:5173` from the same container (`require('playwright')` works
from any directory). To skip onboarding in scripted runs, seed localStorage
before navigating:

```js
await page.addInitScript(() => localStorage.setItem('himalaya-ui:state',
  JSON.stringify({ state: { username: 'demo', tutorialSeen: true }, version: 5 })));
```

### Driving it from an agent

`bin/verify` wraps all of this for agents: `bin/verify seed` brings the demo up in
the background, `open`/`snapshot`/`click`/`press`/`screenshot` drive a shared
headless browser (with onboarding already skipped), and `bin/verify check-truth`
compares the app's phase calls and frame counts with `TRUTH.csv`. See the
`verify-himalaya` skill (`.claude/skills/verify-himalaya/`) and its Feature Map.

### Known noise

- `npm run e2e` (mocked Playwright) also binds :5173; stop the demo first, or it
  fails fast on the port (`--strictPort`).

- Vite logs `[boneyard] capture failed: … Executable doesn't exist at /opt/pw-browsers/…`
  in cloud sessions: the boneyard skeleton-capture plugin uses the frontend's own
  Playwright version, whose browser build isn't in the image. Harmless; skeletons
  still load from the committed `src/bones/` captures.
- On Julia ≥ 1.11, `Pkg.develop(path=".")` rewrites `packages/HimalayaUI/Project.toml`
  (adds `[sources]`, strips comments). `make setup` restores the tracked file; the
  gitignored Manifest keeps the path dependency.
