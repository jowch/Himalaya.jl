# Himalaya Feature Map

How a user reaches each feature, the `bin/verify` commands that drive it, what to check, and its gotchas. Every recipe here was run against `make demo` data on 2026-11-19, and the corpus, Loupe and Focus recipes again on 2026-11-20 after the demo gained bad frames. The ids below (experiment 1, sample 1 = `MO_H2O`, sample 3 = `MO_PEG5`) hold for a fresh `bin/verify seed`; on another DB, read them from `bin/verify api GET /api/samples`.

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
- **Check:** the config page reads "Indexed 28 exposures found", and `check-truth` passes on the new experiment.
- **Gotchas:**
  - The three checks under the path field (`status: ✓ directory exists`, and so on) are the funnel's validation. Snapshot them before clicking "Review →".
  - The page probes `GET /api/fs/resolve` after each path change, and the server answers 400 for anything that is not a directory (a partial path, a file, a missing path). To drive that validation, pass `--expect-http 400:/api/fs/resolve` to the `fill`; without it the command fails on the 400, and with it the command fails if the 400 never comes.
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
- **Check:** `bin/verify check-truth` (`shown`, and the `invariants` block), or `bin/verify api GET "/api/samples?experiment_id=1"` (`assignment_state`, `phase`, `representative_exposure_id`: the frame the status cell describes). `check-truth --ui` reads the status cell text itself.
- **Gotchas:**
  - Returning from Loupe with Esc resets the cursor to the first sample (seen 2026-11-19: ↓↓, `l` opens sample 3's Loupe, Esc, then Enter opens sample 1). After a round trip, `open /sample/<id>` directly instead of relying on the cursor.
  - The bottom dock ("Drop X", "Set representative R", "Loupe L", "Focus ↵") acts on the cursor sample, not the checkbox selection.
  - Bad frames to drop (TRUTH's `drop` column): flares in `MO_H2O_37C` frame 2, `MO_F127` frame 1, `DOPC_H2O` frame 2 and `MO_AMP_10` frame 3; `MO_H2O` frame 3 (beam-off); `DOPE_H2O` frame 1 (missed capillary); and both frames of `leaked_capillary`. `bin/verify triage` drops them; `bin/verify triage --undo` keeps them again. A dropped thumbnail reads "Frame 3, dropped" and the Kept cell "2 / 3". In Loupe a flare shows the sample's rings crossed by a bright horizontal streak through the beamstop, the beam-off frame is black with sparse speckle, and a missed capillary is dark apart from a glow around the beamstop.

## Loupe (frame review)

- **Route:** `/sample/<sampleId>/loupe`. Reached by clicking a sample name on the corpus, with `l` from the corpus or Focus, or the "Loupe L" dock button.
- **Shows:** the detector image, "frame i of n · kept|dropped", frame buttons, "Drop" and "Set as representative", sample tags.
- **Drive:**
  ```bash
  bin/verify open /sample/3/loupe
  bin/verify press ArrowRight                         # next frame
  bin/verify click --role button --name "Set as representative"
  bin/verify press Escape                             # back to the corpus
  ```
- **Check:** `bin/verify api GET /api/samples/3/exposures` (`status` is `rejected` when dropped, `selected` marks the frame picked as representative; a dropped pick is ignored and the newest kept frame is used). The pressed `Frame N` thumbnail is the frame on screen; `check-truth --ui` reads it as `views.loupe`.
- **Gotchas:**
  - Drop writes `status: rejected`; undo with the same button or `bin/verify api PATCH /api/exposures/<id>/status --data '{"status":null}'`.
  - Set as representative can be moved but not cleared (the button disappears and the frame reads "Representative for indexing"). Only `bin/verify seed --reset` returns a sample to "no frame selected".

## Focus (indexing workspace)

- **Route:** `/sample/<sampleId>`. Reached with Enter from the corpus or "Focus ↵" in Loupe. (Clicking a sample name on the corpus opens its Loupe, `/sample/<id>/loupe?experiment=<eid>`, not Focus.)
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
- **Check:** `bin/verify check-truth` (`assigned_frames`, and `shown` for what the corpus then says), or `bin/verify api GET /api/exposures/<id>/assignment`. `check-truth --ui` reads the pressed `Frame N` thumbnail as `views.focus`.
- **Gotchas:**
  - Focus opens the representative, the frame the corpus status describes. If you switch to another frame and assign there, the corpus still describes the representative, which `check-truth` reports under `status_matches_assignment`.
  - Adding or excluding a peak re-analyses the frame and can change its auto call, which `check-truth` then reports as a regression. Undo the edit, or `bin/verify seed --reset`.
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

Hidden pending the redesign (decided 2026-11-19: the Series page is dropped for now and its routes and nav are hidden in the first redesign PR). Not mapped; check-truth does not depend on it.

## Experiment configuration

- **Route:** `/experiments/<id>/config`.
- **Shows:** the experiment name with "Rename experiment", the description button "Edit experiment description" (it reads "Add a description..." when empty), the Geometry card with an "Override <field>" button per field, the Acquisition timeline, and the Sources card with the three patterns and "Rescan now".
- **Drive:** an edit, then a clear, each followed by a reload so the check proves the value persisted on the server.
  ```bash
  bin/verify open /experiments/1/config
  bin/verify click --role button --name "Edit experiment description"
  bin/verify fill --role textbox --name "Experiment description" --value "April beamtime"
  bin/verify press Control+Enter
  bin/verify open /experiments/1/config          # reload: the button reads "April beamtime"
  bin/verify click --role button --name "Edit experiment description"
  bin/verify fill --role textbox --name "Experiment description" --value ""
  bin/verify press Control+Enter
  bin/verify open /experiments/1/config          # reload: the button reads "Add a description..."
  ```
- **Check:** `bin/verify api GET /api/experiments/1` shows `"description": null` after the clear.
- **Gotchas:**
  - Clearing sends `{"description": null}`, not an empty string. Exercise the clear, not only an edit: the server once answered it with 400 and the page showed no error. `press` now fails when a request it triggered returns 4xx/5xx, so that case reads `"ok": false` with the failing request under `settle.httpErrors`. `make verify-truth` runs this edit-then-clear in CI (step 9).
  - "Rescan now" changes the rescan schedule (it goes back to hourly checks when new files are found), which no page shows. Check the schedule in the backend tests, not here.

## Not mapped yet

The custom-index modal, sample tags, figure export (Copy / Download formats), and the onboarding flow (`bin/verify open / --onboarding`). Add each the first time a PR touches it.
