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
make demo-clean   # delete the generated data and demo DBs, large demo included
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

`packages/HimalayaUI/scripts/make_demo_data.jl` writes an SSRL-style tree with two
setup folders, which the funnel's layout auto-detect (`resolve_experiment_layout`,
routes_fs.jl) recognises as one beamtime of two setups. The capillary setup has to
live in `capillary/`: `setup_roots` stops at a folder that has its own `data/`, so a
second setup beside root-level `data/` would never be found.

```
scratch/demo-experiment/
  capillary/data/<label>_S<scan>_0_001.tif     Pilatus 1M frame, signed Int32 counts, -1 in module gaps
  capillary/data/<label>_S<scan>_0_001.prp     PRP sidecar (energy, pipe length, horizontal position, time, beam monitors…)
  capillary/analysis/setup_info_*.txt          beam center + AgBe distance (geometry source)
  capillary/analysis/automatic_analysis/tot_files/<label>_S<scan>_tot.dat   q  I  sigma
  autosampler/                                 the solution SAXS setup, see "The autosampler setup"
  TRUTH.csv                                    per capillary sample: scans, generated phases + lattice, expected
                                               auto call (`auto`), phase-tagged reflection q values,
                                               scans to drop (`drop`)
  TRUTH_autosampler.csv                        per autosampler acquisition, see below
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
- **Some scans are bad on purpose**, so the drop and selection rules have
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
- **`make upgrade-drill`** (also a CI job) checks that a DB written by older code
  still reads right. It ingests the demo with the base commit's backend
  (`UPGRADE_BASE`, default `origin/main`; CI uses the PR's base), reviews a few
  samples through the base's routes with the base's own copy of
  `upgrade_history.mjs` (a call, then only an image without that phase
  selected; a call on one image and a peak edit on another; a call on an image
  that is then dropped; calls on two images; no phase), then opens that DB with
  this checkout. It fails if any sample's call changed, if any check-truth
  invariant broke or went unchecked, or if the first sample doesn't show its
  call ticked and marked "Not found". Auto calls
  that differ from TRUTH are only reported, since the base wrote those analyses.
  A fresh demo never runs a migration, so this is the check that does.
- **Expected calls are not always the generated truth**, by design. TRUTH's `auto`
  column is what Himalaya's own pipeline calls. One sample differs:
  - `MO_AMP_weak`: too weak to index automatically.

  `MO_PEG5` (Pn3m + Im3m) used to be the other one. Its Pn3m √6 sits about 1% from
  Im3m √10, and `auto_group` now gives each such peak to the lattice that fits it
  best, so both phases are called (`docs/scoring.md`).

  `test/test_demo_data.jl` holds a regression floor on these calls over three
  seeds. If a change to peak finding, scoring or `phaseratios` moves them, fix
  the change or retune the sample deliberately.
- Seeded (`--seed`, default 7) and deterministic within a Julia release; ~180 MB,
  ~18 s (about 110 MB and 15 s of that is the capillary setup, which is byte-identical to
  the single-setup demo). Traces and detector frames draw from separate random streams, so
  `--no-images` (empty `.tif` placeholders, ~1 MB) writes the same `.dat` files.
  `make demo` refuses such a tree; regenerate with `make demo-data`.
- Safety: the generator refuses to write into a non-empty directory that has no
  `TRUTH.csv`, and `make demo-clean` only deletes `DEMO_DIR` when it has one.
- TIFFs are written as `Gray{Q0f31}`, which TiffImages stores as a signed 32-bit
  integer TIFF, the same on-disk format as the Pilatus. Checked both ways: a
  plain int32 TIFF from Python `tifffile` loads in Julia as `Gray{Q0f31}` with
  the counts intact, and the generator's files read back in `tifffile` as int32.

`make demo` ingests through `scripts/demo_ingest.jl`, which makes the same two
calls the funnel does (`/api/fs/preview` → `POST /api/experiments`, with both
setups, capillary first so its sample ids stay 1–13) and waits for the first scan.
On a cold backend the ingest takes about 40 s, against 30 s for the capillary setup alone. The first browser visit still shows the
"Who are you?" onboarding; any name works.

## The autosampler setup

`autosampler/` is a solution SAXS setup in the format of a real autosampler
folder (every value invented):

```
autosampler/
  data/<label>_<PP><W>_<S|B><scan>_0_<frame>.tif   10 frames per acquisition, 001–010; frame 001 is rendered, 002–010 are hardlinks to it
  data/<label>_<PP><W>_<S|B><scan>_0_<frame>.prp   one per frame; same horizontal position in every file, 1 s exposure
  data/AT_SolSAXS_<label>_<scan>.txt               one per queue; only its name is read
  data/agbe_S0100_0_001.{tif,prp}                  calibration shot, no well token, no averaged trace
  analysis/integ.mpp                               `-i yes i2 2.5e7`, so traces are divided by I_2
  analysis/setup_info_*.txt
  analysis/automatic_analysis/tot_files/<stem>_tot.dat      the average of the kept frames: q  I  sigma
  analysis/automatic_analysis/sastool/<stem>_0_001.log      the averaging log the app reads frame counts from
```

Three runs in two queues: `lysA` on plate 01 (first buffer A, samples B–F, closing
buffer H), `micB` on plate 02 (first buffer A, samples B–E, no closing buffer) and
`lysA` again after a second `AT_SolSAXS_lysA_0113.txt` (first buffer A, samples C
and G, closing buffer H). That is 16 acquisitions, 160 frames and the calibration
shot; the app ingests 17 samples with one image each, in loads 101, 108, 113 and 0.

A buffer is a flat solvent level with a low-q upturn. A sample is the same buffer
plus a polydisperse sphere form factor with its own radius and amplitude, so
sample minus its run's first buffer has a minimum at 4.4934 / radius
(`minimum_q`). Noise is Poisson-like for the kept frames and sigma matches it.
Traces keep the capillary demo's scale (about 17 at high q).

Cases kept on purpose, each flagged in `TRUTH_autosampler.csv`:

| Flag | Where | What it exercises |
|---|---|---|
| `several_rejected` | `lysA 01D` | 4 of 10 frames rejected |
| `one_frame_left` | `micB 02C` | 9 of 10 rejected; its trace is `sqrt(10)` noisier |
| `i2_collapsed` | `lysA 01F` | `.prp` I_2 at 4% of normal; the trace is 25 times too high, as SASTool's division makes it |
| `no_closing_buffer` | every acquisition of `micB` | a run with no buffer after its samples |
| `repeated_well` | `lysA 01C (S0114)` | same label and well as `lysA 01C (S0103)`, new scan, in a later run |
| `first_buffer` | the first buffer of each run | the default background of the run's samples; its log has no summary lines, so its frame count is read from the per-frame lines, which needs the log to be older than a minute (the generator backdates it by an hour) |
| `calibration` | `agbe_S0100` | a shot with no well token: load 0, no default background, no frame counts |

The autosampler streams are separate from the capillary ones, so adding or changing
this setup never moves a capillary trace or `TRUTH.csv` row.

`TRUTH_autosampler.csv` has one row per acquisition and one for the calibration shot:

| Column | Meaning |
|---|---|
| `sample` | the sample name the app shows, `<label> <PP><W> (<S\|B><scan>)`; the calibration shot's is its stem |
| `stem` | the exposure's filename, without the frame suffix |
| `run` | series label (empty for the calibration shot) |
| `load` | the run's load index: the scan of its first acquisition, `0` for the calibration shot |
| `well` | plate row and well, like `01B` |
| `kind` | `sample`, `buffer` or `calibration` |
| `default_background` | `sample` of the run's first buffer; empty for a buffer and the calibration shot |
| `frames_total`, `frames_rejected` | the counts the log holds; empty for the calibration shot |
| `flags` | space-separated flags from the table above |
| `repeat_of` | for `repeated_well`, the earlier sample's name |
| `radius_A`, `amplitude`, `minimum_q` | the form factor that was put in (empty for buffers); `amplitude` is relative to the solvent level |
| `note` | free text |

`bin/verify check-truth` reads this file into its own `autosampler` block and leaves
these samples out of everything else. Per row it checks, from the file and not from
reading the app back:

- ingestion: one sample with one image whose filename is `stem`, in the autosampler
  setup, and no autosampler sample is unlisted;
- `GET /api/experiments/:id/compares`: the effective background is `default_background`
  (the first buffer of the sample's own run, so the run without a closing buffer and the
  repeated-label run resolve to their own), its source is `default`, and `is_buffer` matches
  `kind`; a buffer and the calibration shot have none;
- `frames_total` and `frames_rejected` on the exposure (empty for the calibration shot);
- no automatic analysis: no peaks, no indices, no call, no ticked indices, no analysis error;
- with `--ui`, the sample's workspace: the Background module reads `This is a buffer` for a
  buffer, `<default_background>run's buffer` for a sample and `No buffer in this run` for the
  calibration shot; a sample is drawn in Minus with its difference line and the Minus button
  pressed; no call candidates, peak markers or custom index panel are shown.

The collapsed I_2 and the form factors are in the file for a later check;
`test/test_demo_autosampler.jl` checks them against the ingested DB.

## Cloud sessions (claude.ai/code)

The app needs Julia (compat says ≥ 1.9; tested on 1.13), Node 22 (22.22.2+), 24 (24.15+) or 26+ and npm. A cloud environment that runs it needs:

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

### Late traces

At SSRL the integration pipeline writes a `_tot.dat` some time after its image,
so a pass can see an exposure before its trace (see docs/experiment-config.md,
"Scan and grouping"). `make_demo_data.jl --late-trace` (`make demo-data
DEMO_DATA_FLAGS=--late-trace`, or `bin/verify seed --reset --late-trace`) holds
`MO_lowH2O`'s traces back in `<demo>/late/`, outside the setup folders. The
workspace notice should then report them as waiting for a trace, and
`check-truth` expects no auto call on those frames. `bin/verify deliver-late`
moves them into `tot_files` with a fresh mtime. The next pass has to analyze
them, and `check-truth` fails while any delivered trace is still reported as
waiting. With `--ui` it also checks the plot: a held frame says "Waiting for
this image's trace." and every other frame draws its trace. `make verify-truth`
(CI's "Demo check-truth" job) runs this sequence first, and delivers the traces
while the workspace is open on `MO_lowH2O`. It then waits, without a reload,
for the plot to draw the trace and its peak markers. That step is the only guard for a trace
that lands while its sample is open: `plot_matches_trace` cannot see a stale
plot, because each sample click remounts the plot and refetches a failed trace.

### Large mode

`make demo-large-data` (generate only, into `scratch/demo-large/`) or
`make demo-large` (generate if missing, ingest every experiment into one DB at
`scratch/demo-large-db/`, start the app) writes realistic scale for performance
work. Override with `DEMO_LARGE_DIR` and `DEMO_LARGE_DB_DIR`. Directly:
`make_demo_data.jl [out_dir] --large [--experiments N] [--samples M] [--seed S] [--no-images]`
(defaults 10 experiments of 300 samples).

- Large mode has no autosampler setup: it exists for scale and loading tests, each
  `large_NN/` is a single capillary setup (its own `data/`), and the autosampler's
  per-acquisition files and logs would not add to what it measures. The curated demo
  is the one place to exercise the mode.
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
- `make demo-clean` also deletes `DEMO_LARGE_DIR` and `DEMO_LARGE_DB_DIR`. It
  deletes the large tree only when everything in it is a `large_NN` folder or
  the `.demo-large-complete` marker, so an interrupted run is cleaned up too.
  Any other file, a hidden one, or a link makes it leave the folder alone and
  say why, and so does a folder with no `large_NN` folder in it. The check runs
  in Julia (`--clean-large` in `make_demo_data.jl`) and only when
  `DEMO_LARGE_DIR` exists, so it needs the `make setup` environment only then.
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

### Slow-network checks

`VERIFY_THROTTLE=<mbit>mbit,<ms>ms bin/verify …` applies Chrome's network
emulation to every command that drives the browser, for example `VERIFY_THROTTLE=20mbit,30ms`. The latency
is a minimum round trip that overlaps server time, so a slow response is not
delayed further.

### Known noise

- `npm run e2e` (mocked Playwright) uses its own port (:5191) and build folder
  (`dist-e2e/`), so the demo can keep running.

- On Julia ≥ 1.11, `Pkg.develop(path=".")` rewrites `packages/HimalayaUI/Project.toml`
  (adds `[sources]`, strips comments). `make setup` (and `make dev-setup`) restores the tracked file; the
  gitignored Manifest keeps the path dependency.
