# Himalaya Feature Map

How a user reaches each feature, the `bin/verify` commands that drive it, what to check, and its gotchas. Every recipe here was run against `make demo` data on 2026-11-19. The ids below (experiment 1, sample 3 = `MO_PEG5`) hold for a fresh `bin/verify seed`; on another DB, read them from `bin/verify api GET /api/samples`.

Names in `--name` come from `bin/verify snapshot --raw`. When a recipe stops matching, the UI changed: fix the entry in the same PR.

Keys are global unless noted; `?` opens the shortcut overlay on any page.

## Experiments home

- **Route:** `/experiments` (`/` redirects here). The top nav link "Experiments" returns here from anywhere.
- **Shows:** a timeline of experiments, each a button named `<name> <status> <data dir> <n> samples · <n> exp · <n> loads`.
- **Drive:**
  ```bash
  bin/verify open /experiments
  bin/verify click --role button --name "Demo (synthetic)"     # opens the corpus
  ```
- **Check:** `bin/verify api GET /api/experiments` (`ingest_status`, `stats`).

## Ingest funnel (new experiment)

- **Route:** `/experiments/new`, then `/experiments/new/config`, then `/experiments/<id>/grouping` after the scan. Reached from the "+ New experiment" button on the home page.
- **Needs an empty DB.** Start a second instance so the main demo is untouched:
  ```bash
  export DEMO_API_PORT=8081 DEMO_UI_PORT=5174 VERIFY_CDP_PORT=9334 VERIFY_STATE_DIR=scratch/verify-empty
  bin/verify seed --empty
  bin/verify open /experiments/new
  bin/verify fill --role textbox --name "Data directory" --value "$PWD/scratch/demo-experiment"
  bin/verify click --role button --name "Review →"
  bin/verify wait-settle --timeout 60000               # waits out the "Indexing directory" progressbar
  bin/verify click --role button --name "Approve"
  bin/verify wait-settle --url '/experiments/\d+/grouping' --timeout 120000
  bin/verify check-truth
  bin/verify close && bin/verify stop                  # then unset the four variables
  ```
- **Check:** the config page reads "Indexed 26 exposures found", and `check-truth` passes on the new experiment.
- **Gotchas:**
  - The three checks under the path field (`status: ✓ directory exists`, and so on) are the funnel's validation. Snapshot them before clicking "Review →".
  - The config step keeps its state in memory. Reloading `/experiments/new/config` (or a Vite full reload, which can follow an edit under `frontend/`) shows "No experiment selected"; start again from `/experiments/new`.
  - The scan after Approve runs on the server and reports over SSE, so a plain `wait-settle` returns before it ends. Use `--url` as above, or poll `bin/verify api GET /api/experiments/<id>` for `ingest_status: complete`.

## Corpus (contact sheet)

- **Route:** `/experiments/<id>/corpus`. Reached by clicking an experiment on the home page, or Esc from Focus or Loupe.
- **Shows:** stats (loads, samples, exposures, span, sessions), a grouping banner with a "Review grouping →" link, and a `table "Samples"` with a row per sample: checkbox "Select sample", the sample name button, frame thumbnails ("Frame 1"…), kept count and phase ("Not indexed" until a human assigns one).
- **Drive:**
  ```bash
  bin/verify open /experiments/1/corpus
  bin/verify press ArrowDown                          # move the sample cursor
  bin/verify press Enter                              # Focus on the cursor sample
  bin/verify press l                                  # Loupe on the cursor sample
  bin/verify press x                                  # drop the cursor frame; x again restores it (writes to the DB)
  bin/verify click --role checkbox --name "Select sample" --nth 2
  ```
- **Check:** `bin/verify api GET "/api/samples?experiment_id=1"` (`assignment_state`, `phase`).
- **Gotchas:**
  - Returning from Loupe with Esc resets the cursor to the first sample (seen 2026-11-19: ↓↓, `l` opens sample 3's Loupe, Esc, then Enter opens sample 1). Click the sample name instead of relying on the cursor after a round trip.
  - The bottom dock ("Drop X", "Set representative R", "Loupe L", "Focus ↵") acts on the cursor sample, not the checkbox selection.

## Loupe (frame review)

- **Route:** `/sample/<sampleId>/loupe`. Reached with `l` from the corpus or Focus, or the "Loupe L" dock button.
- **Shows:** the detector image, "frame i of n · kept|dropped", frame buttons, "Drop" and "Set as representative", sample tags.
- **Drive:**
  ```bash
  bin/verify open /sample/3/loupe
  bin/verify press ArrowRight                         # next frame
  bin/verify click --role button --name "Set as representative"
  bin/verify press Escape                             # back to the corpus
  ```
- **Check:** `bin/verify api GET /api/samples/3/exposures` (`status` is `rejected` when dropped, `selected` marks the representative).
- **Gotchas:**
  - Drop writes `status: rejected`; undo with the same button or `bin/verify api PATCH /api/exposures/<id>/status --data '{"status":null}'`.
  - Set as representative can be moved but not cleared (the button disappears and the frame reads "Representative for indexing"). Only `bin/verify seed --reset` returns a sample to "no frame selected".

## Focus (indexing workspace)

- **Route:** `/sample/<sampleId>`. Reached with Enter from the corpus, a click on the sample name, or "Focus ↵" in Loupe.
- **Shows:** the trace plot (`img "Integration trace: intensity vs q"`) with a toolbar (log/linear q, Auto-fit, + Peak, Copy, Download formats), detector frames, the reflection comb, and the assignment rail (`complementary`) with candidate buttons such as `button "Im3m"` ("explains 5 peaks 0.35"), "Form factor, no Bragg peaks" and "+ custom index…".
- **Drive:**
  ```bash
  bin/verify open /sample/3
  bin/verify click --role button --name Pn3m --within aside        # add a phase to the call
  bin/verify click --role button --name Im3m --exact               # exact: "Im3m, in assignment" also matches once checked
  bin/verify snapshot --selector aside --raw                       # "Assignment 2 Coexistence · 2 phases"
  bin/verify click --role button --name "Remove Pn3m"              # undo
  bin/verify press p                                               # add-peak mode; peaks show as "Auto peak at q = …"
  bin/verify press ArrowDown                                       # next sample
  ```
- **Check:** `bin/verify check-truth` (`assigned`, `assigned_matches_truth`, `assigned_frames`), or `bin/verify api GET /api/exposures/<id>/assignment`.
- **Gotchas:**
  - Focus opens frame 1, but the agreed representative rule picks the newest kept frame. An assignment made here therefore shows up in `check-truth` under `assigned_off_representative`. This is the known representative-frame inconsistency (audit PR 2), not a CLI bug.
  - `MO_PEG5` is the coexistence test case: the auto call is Im3m only, and a human must add Pn3m.

## Grouping review

- **Route:** `/experiments/<id>/grouping`. Reached from "Review grouping →" on the corpus, and automatically after an ingest.
- **Shows:** "0 flagged of 12 samples" for the demo, a "Needs review" / "All loads" filter, "Split", "Confirm groups".
- **Drive:**
  ```bash
  bin/verify open /experiments/1/grouping
  bin/verify click --role button --name "All loads"
  ```
- **Check:** `check-truth` reports `strays` (an app sample that mixes two truth samples) and `split into n app samples`.

## Series

- **Folio route:** `/series` (top nav "Series"). Lists saved series with search, filters (All, Has transition, Cross-experiment) and sort. Empty on a fresh demo ("No series yet").
- **New series:** select two or more samples on the corpus, then "+ New series", which opens `/series/new` (scoping: name the ordering variable and give each sample a value).
  ```bash
  bin/verify open /experiments/1/corpus
  bin/verify click --role checkbox --name "Select sample" --nth 0
  bin/verify click --role checkbox --name "Select sample" --nth 1
  bin/verify click --role button --name "+ New series"
  bin/verify fill --role textbox --name "Ordering variable name" --value "temperature"
  bin/verify fill --role textbox --name "Value for MO_H2O_S2401 (S01P01)" --value "25"
  bin/verify fill --role textbox --name "Value for MO_H2O_37C_S2404 (S01P02)" --value "37"
  bin/verify click --role button --name "Confirm & build →"      # opens /series/<id>, the builder
  bin/verify screenshot scratch/verify/shots/series-builder.png
  ```
- **Builder** (`/series/<id>`): the overlay plot with a log/linear q toolbar, Copy and Download formats, "Peak ticks" / "Peak labels" toggles, a "Trace offset" slider, a Traces list, and "Add sample A" in the dock.
- **Check:** `bin/verify api GET /api/series`, `bin/verify api GET /api/series/<id>/traces`.
- **Gotchas:**
  - "Confirm & build" saves a series and records the variable as a `scoping` tag on every sample. `bin/verify api DELETE /api/series/<id>` removes the series but not the tags; `seed --reset` removes both. "Discard" leaves nothing behind.
  - Each trace is the sample's newest frame ("exp 3" for `MO_H2O`), and a sample with no assignment is listed as "Form factor" (seen 2026-11-19).

## Not mapped yet

Experiment configuration (`/experiments/<id>/config`), the custom-index modal, sample tags, figure export (Copy / Download formats), and the onboarding flow (`bin/verify open / --onboarding`). Add each the first time a PR touches it.
