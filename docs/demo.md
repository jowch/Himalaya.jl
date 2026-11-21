# Demo: run the real app on synthetic data

How to run HimalayaUI end to end (Julia backend + Vite frontend) on a synthetic
SAXS experiment, locally or in a Claude Code cloud session, without any lab data.
Used for UI walkthroughs and Playwright screenshots.

## Quick start

```bash
make dev-setup    # once per checkout: npm ci + Pkg.develop/instantiate (≈5 min cold), with
                  # the ~100 s route precompile workload off for dev (CLAUDE.md)
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
[frontend-dev-loop.md](frontend-dev-loop.md)), `DEMO_UI_MODE=build` (serve a `vite build`
through `vite preview` instead of the dev server).

`make verify-truth` sets `DEMO_UI_MODE=build`. It loads the corpus and the workspace dozens of times
in one Chrome tab, and every load of the dev server's ~280 unbundled modules leaves the
renderer holding shared-memory handles until garbage collection. In CI that piles up past
2,000 handles and Chrome starts failing loads with `ERR_INSUFFICIENT_RESOURCES`. The built
app loads a few chunks, stays flat, and is what people actually run.

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
  data/<label>_S<scan>_0_001.prp     PRP sidecar (energy, pipe length, horizontal position, time, beam monitors…)
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
  generator), and its `flare` column says which of those still hold the sample.
  Flares are the commonest reason to drop a frame at the beamline, so they are
  the commonest here:
  - Flares, in `MO_H2O_37C` (newest scan), `MO_F127` (oldest), `DOPC_H2O`
    (newest) and `MO_AMP_10` (newest): the beam clipped the capillary wall. The detector image
    shows a bright horizontal streak through the beamstop, perpendicular to the
    capillary. The sample's rings are still there at about 60% strength, and the
    1D trace gains a low-q excess. Across seeds 1–12, 24 of 48 flares kept the
    sample's auto call and 24 lost a phase or every phase, so a flare can mislead
    the auto call without looking wrong in 1D.
    A flare must be dropped, and its auto call may lose phases but must not add
    one.
  - `MO_H2O`: the newest scan is a beam-off shot: a flat dark floor of about
    0.01 counts per pixel plus hot pixels, so the detector image is black. Its
    `.prp` I₀ is the ion chamber's negative electrometer offset (about −0.07 to
    −0.12 V, as on real SSRL beam-off frames), not 0. The workspace's "No counts
    on the detector" caption shows on it because its trace is dark, and the
    negative I₀ adds a supporting hint.
  - `DOPE_H2O`: the oldest scan missed the capillary: steep q⁻⁴ air and slit
    scatter over a faint air floor, so the detector image shows a glow at the beamstop and
    darkness beyond it.
  - `leaked_capillary`: both scans missed the capillary in the same way, so the
    sample has no usable frame.

  Nothing drops them at ingest; a triager does. `bin/verify triage`
  drops them through the endpoint the workspace's Drop (X) uses, so `bin/verify check-truth` can check the corpus
  after triage. `make verify-truth` (also a CI job) runs the whole sequence on a
  freshly regenerated demo, all with `--ui`: check-truth expecting every frame
  kept, triage, check-truth expecting every bad frame dropped, confirm
  `MO_H2O`'s Pn3m call as a person would, then clear the analysis hashes and
  re-run `himalaya analyze` so every kept frame is re-indexed from scratch. The
  last step checks that nothing a person sees changed, the confirmed call
  included. `analyze` exits non-zero if it skipped any exposure.
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
   script machine-level; run `make dev-setup` in the session (it needs the repo).

Then, in the session: `make dev-setup && make demo`. Playwright can drive
`http://127.0.0.1:5173` from the same container (`require('playwright')` works
from any directory). To skip onboarding in scripted runs, seed localStorage
before navigating:

```js
await page.addInitScript(() => localStorage.setItem('himalaya-ui:state',
  JSON.stringify({ state: { username: 'demo', tutorialSeen: true }, version: 5 })));
```

### How close it is to real data

The generator is tuned to aggregate figures from SSRL Pilatus 1M lipid
beamtimes, so images and traces load and index like real ones:

| | Demo | Beamline |
|---|---|---|
| Counts per pixel, far from the beam | ~100 | ~100 |
| Beamstop shadow radius | 10 px | ~10 px |
| Trace scatter at high q | ~0.3% | ~0.35% |
| Low-q rise (I(0.006) ÷ high-q floor) | ~100–200× | ~110× (median) |
| Bilayer hump near 0.11 Å⁻¹ | ~2× the floor | ~2× (median) |
| Strongest peak ÷ local background | 0.1–8 | median ~0.4, p90 ~8 |
| Full detector PNG the image route sends | ~700 KB | ~850 KB |

Two differences remain on purpose. The curated samples sit at the clear end of
the real range, because each has a known answer `test_demo_data.jl` checks; most
real frames are weaker, and over a third have no automatic peaks. And there are
13 samples with 1–3 frames each, where a real experiment has hundreds of samples
with a median of 4 frames. Performance work needs real scale and a real network
(the app is usually reached through an SSH tunnel), so test it with the
`--large` option below and a throttled browser.

### Large mode

`make demo-large-data` (generate only, into `scratch/demo-large/`) or
`make demo-large` (generate if missing, ingest every experiment into one DB at
`scratch/demo-large-db/`, start the app) writes realistic scale for performance
work. Override with `DEMO_LARGE_DIR` and `DEMO_LARGE_DB_DIR`. Directly:
`make_demo_data.jl [out_dir] --large [--experiments N] [--samples M] [--seed S] [--no-images]`
(defaults 10 experiments of 300 samples).

- One complete beamtime tree per experiment, `out_dir/large_01/` onward, each with
  its own `TRUTH.csv` (same columns; `auto` is empty because nothing is pinned),
  start date and flight path (0.8 m or 1.7 m setups). Sample labels are `L01_N0001`.
- Per sample: 4 frames about 82% of the time, with short, long and very long
  series in the tail; about 40% have no Bragg peaks; the rest have one phase
  (85%) or two (15%), weighted towards Lamellar, Im3m and Square; peak height and
  width are log-normal at the real medians; about 6% of frames are bad (flare,
  missed, beam off). Holders of 10-16 samples sit more than an hour apart, so
  grouping recovers the loads and samples.
- **Images come from a pool.** Each experiment renders about 24 frames (no-peak,
  weak and strong examples of each phase, flare, missed, beam off) and every
  exposure's `.tif` is a hardlink to one of them. A frame's image therefore does
  not match its own trace; use large mode for scale, not for checking images
  against traces. `--no-images` writes empty `.tif` files instead.
- **Geometry is approximate.** Pool frames use each experiment's flight path,
  but every trace keeps the curated demo's q grid and noise model, and
  `setup_info` keeps its 1.8 m calibration rings. Fine for scale and loading
  tests; don't use large mode to check q calibration.
- A finished run writes `.demo-large-complete`. An interrupted run has none, so
  `make demo-large` regenerates instead of ingesting a partial tree, and the
  generator only writes into an empty folder or an existing large tree.
- Cost, measured with the defaults: about 35 s, 1.4 GB on disk (12,000 `.tif`
  names, 240 distinct files), 3,000 samples and 12,000 exposures. Ingesting all
  ten through `demo_ingest.jl` takes about 3 minutes.
- The defaults are a worst case: each experiment is as large as the largest real
  one, and the total is about 1.7× a typical lab database (about 1,800 samples
  over 11 experiments). Use `--samples 180` to match that total. Costs that grow
  with the whole database show the difference: `GET /api/experiments` took 6.3 s
  on the defaults and 1.6 s on the real database (#348).

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
  (adds `[sources]`, strips comments). `make setup` (and `make dev-setup`) restores the tracked file; the
  gitignored Manifest keeps the path dependency.
