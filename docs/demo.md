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

Overrides: `DEMO_DIR=/tmp/demo make demo` (data location), `DEMO_API_PORT` /
`DEMO_UI_PORT` (ports, for running beside another backend; see
[frontend-dev-loop.md](frontend-dev-loop.md)).

`make demo` keeps its DB in `scratch/demo/` between runs (edits persist; ingest is
skipped when the experiment is already there). `make demo-empty` wipes
`scratch/demo-empty/` on every start. Backend logs go next to each DB.

## The synthetic data

`packages/HimalayaUI/scripts/make_demo_data.jl` writes an SSRL-style tree that the
funnel's layout auto-detect (`resolve_experiment_layout`, routes_fs.jl) recognises:

```
scratch/demo-experiment/
  data/<label>_S<scan>_0_<NNN>.tif   Pilatus 1M frame, signed Int32 counts, -1 in module gaps
  data/<label>_S<scan>_0_<NNN>.prp   PRP sidecar (energy, pipe length, horizontal position…)
  analysis/setup_info_*.txt          beam center + AgBe distance (geometry source)
  analysis/automatic_analysis/tot_files/<label>_S<scan>_tot.dat   q  I  sigma
  TRUTH.csv                          generated phase, lattice and peak q for each sample
```

- **Not real data, never committed.** It lives under gitignored `scratch/` (or
  wherever `DEMO_DIR` points). Only the generator is in git.
- **Peak positions come from `Himalaya.phaseratios`**, so the demo tracks any
  change to the core ratio tables (e.g. the √11 removal in #304).
- 12 samples, one per capillary slot (3.95 mm pitch, so the slot clustering in
  grouping.jl splits them): clean Pn3m ×2, Pn3m+Im3m, Im3m, Ia3d, HII, Lα,
  Lα+HII, weak Pn3m, Pn3m+HII, a buffer blank, and a ULV form-factor-only trace.
- Seeded (`--seed`, default 7) and deterministic within a Julia release; ~118 MB,
  ~10 s. `--no-images` writes empty `.tif` placeholders (~1 MB) when only traces matter.
- TIFFs are written as `Gray{Q0f31}`, which TiffImages stores as a signed 32-bit
  integer TIFF, the same on-disk format as the Pilatus. Checked both ways: a
  plain int32 TIFF from Python `tifffile` loads in Julia as `Gray{Q0f31}` with
  the counts intact, and the generator's files read back in `tifffile` as int32.

`make demo` ingests through `scripts/demo_ingest.jl`, which makes the same three
calls the funnel does (`/api/fs/resolve` → `/api/fs/manifest` → `POST /api/experiments`)
and waits for the first scan. The first browser visit still shows the
"Who are you?" onboarding; any name works.

## Cloud sessions (claude.ai/code)

The app needs Julia ≥ 1.10, Node ≥ 20 and npm. A cloud environment that runs it needs:

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

### Known noise

- Vite logs `[boneyard] capture failed: … Executable doesn't exist at /opt/pw-browsers/…`
  in cloud sessions: the boneyard skeleton-capture plugin uses the frontend's own
  Playwright version, whose browser build isn't in the image. Harmless; skeletons
  still load from the committed `src/bones/` captures.
- On Julia ≥ 1.11, `Pkg.develop(path=".")` rewrites `packages/HimalayaUI/Project.toml`
  (adds `[sources]`, strips comments). `make setup` restores the tracked file; the
  gitignored Manifest keeps the path dependency.
