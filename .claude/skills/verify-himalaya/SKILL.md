---
name: verify-himalaya
description: Prove a HimalayaUI change works in the real app. Bring up the synthetic demo, drive the browser with bin/verify (snapshot, click, press, screenshot), check results against the demo TRUTH.csv, and attach the evidence to the PR. Use after any change a user could see, and before asking for review.
---

# verify-himalaya

Unit tests show that a function returns what it returns. This skill shows that the app does what it should: it runs the real backend and frontend on the synthetic demo experiment, drives them the way a user would, and checks the outcome against known answers.

Everything goes through one CLI, `bin/verify`. Every command prints one JSON object and exits 0 on success or 1 on failure. A failure carries `error` and a `hint` that says what to do next. Run `bin/verify help` or `bin/verify <command> --help` for the full option list.

## The loop

```bash
bin/verify doctor                 # what is up, what is missing, and the fix for each failing check
bin/verify seed                   # demo up in the background (generate, backend :8080, ingest, Vite :5173); no-op when already up
bin/verify open /experiments/1/corpus
bin/verify snapshot --raw         # accessibility tree: the roles and names you can click
bin/verify click --role button --name "MO_PEG5_S2406 (S01P03)"   # opens it in the workspace
bin/verify press Enter            # opens the workspace; keyboard shortcuts are first-class, see the Feature Map
bin/verify screenshot scratch/verify/shots/focus-peg5.png
bin/verify check-truth            # phase calls, frame counts and the corpus status against TRUTH.csv
bin/verify triage                 # drop the bad frames TRUTH lists, as a triager would
bin/verify check-truth --ui       # also reads the corpus and the workspace in the browser
bin/verify stop                   # stop the servers seed started (bin/verify close stops the browser)
```

Read [references/features/README.md](references/features/README.md) first. It is the Feature Map: every feature, how a user reaches it, the commands that drive it, and its gotchas.

## Commands

| Command | Use it to |
|---|---|
| `doctor` | Check node, julia, frontend deps, Chromium, demo data, backend, DB, ingest, Vite and Vite's `/api` proxy. |
| `seed [--empty] [--reset]` | Start the demo detached and wait until it is ready. `--empty` starts on a fresh DB to walk the ingest funnel. `--reset` deletes the demo DBs and data first, and closes the browser and deletes its profile, since an open tab would keep showing the old DB (`--dry-run` lists what it would delete). |
| `open <route>` | Load a route in the shared headless browser, starting it if needed, and wait for it to settle. A user is seeded so onboarding is skipped; `--onboarding` clears the user so it shows. Fails when the page does not settle, or when any request made while loading returns 4xx/5xx (a missing thumbnail image included). |
| `snapshot [--selector css] [--raw]` | Read the page as an aria tree. Use it before every click to get exact names. |
| `click` / `fill` / `press` | Act on exactly one element, chosen by `--role` + `--name`, `--label`, `--text`, `--testid` or `--selector`, narrowed by `--within <css>`, `--exact` or `--nth <i>`. Each then waits for the page to settle, including requests the action started, and fails if it does not settle within `--timeout`, or if one of those requests returned 4xx/5xx (`settle.httpErrors`). When an error status is the app's normal answer, allow exactly that one with `--expect-http <status>:<path regex>` (repeatable, `x` wildcards like `4xx`); the command then fails if no matching response arrives, so an error-path check can't pass vacuously. |
| `wait-settle [--url re]` | Wait for a server-triggered navigation (with `--url`), then for the page to settle. Fails on 4xx/5xx like the actions; its result has `settled` and `httpErrors` at the top level, not under `settle`. |
| `screenshot [file] [--full] [--selector css]` | Save evidence. Default path `scratch/verify/shots/`. |
| `api <METHOD> </api/path> [--data json] [--dry-run]` | Read or change backend state directly. Writes are real; use `--dry-run` to check the request first. |
| `check-truth [--strict] [--ui]` | Compare what the app shows with the demo's `TRUTH.csv`, and report the invariants. `--ui` reads the corpus status, and the status and frame the workspace shows, in the browser. |
| `triage [--undo] [--dry-run]` | Drop the frames in TRUTH's `drop` column through the same endpoint the workspace's Drop (X) uses; `--undo` keeps them again, including any a person dropped by hand. |

## check-truth

`TRUTH.csv` (written by `make demo-data`) records two answers per sample: `phases`, the phases the generator put in, and `auto`, what Himalaya's own pipeline is expected to call. They differ on purpose for samples a human must finish by hand (see `docs/demo.md`).

- Default gate: every frame's auto call must equal `auto`, and every generated frame must be present in the right sample. A failure here is a regression in peak finding, indexing, grouping or ingest.
- `known_gaps` lists samples where `auto` misses `phases`. Today that is `MO_PEG5` (Pn3m + Im3m coexistence, auto keeps only Im3m) and `MO_AMP_weak` (too weak to index). `--strict` fails on these too. Use it when working on the physics.
- Bad frames: TRUTH's `drop` column lists frames a triager should drop. Its `flare` column marks the ones that still hold the sample (the beam clipped the capillary wall); their auto call may lose phases, because the streak can hide one, but must not add one (`expected: "within Pn3m"`). The rest hold no sample (a beam-off shot, missed capillaries) and must call `none`. `frames.bad_still_kept` and `frames.good_dropped` report triage progress. Neither fails the gate unless you pass `--expect fresh` (every frame still kept) or `--expect triaged` (every bad frame dropped). With `--ui`, check-truth also reads the corpus Kept column and fails if it disagrees with the API.
- `make verify-truth` runs the whole live check on a regenerated demo, as CI's "Demo check-truth" job does: seed, check-truth `--expect fresh`, triage, check-truth `--expect triaged`, confirm one sample's call, a forced re-index (`analyze` after clearing the analysis hashes), check-truth, then a comparison that nothing a person sees changed, the description edit-and-clear, a custom index added and removed in the workspace, and clicks on an auto and a manual peak that share an id, excluding and restoring the auto peak while its manual twin is still drawn and then removing the manual one (the real DB numbers the two sources from separate sequences, so ids collide; the step adds the manual peak where the next manual id lands on an auto peak of a kept frame, and fails rather than passing vacuously when it cannot make them collide). Each step's output lands in `scratch/verify/truth/`. It serves the built app (`DEMO_UI_MODE=build`), not the Vite dev server; see docs/demo.md for why.
- `shown` is the sample's status: `unreviewed`, `phase` (with its `call`) or `form_factor`. It is read from `GET /api/samples`, the data the corpus renders, together with the frame it describes (`representative_exposure_id`, reported as `corpus_api` because the page never marks that frame). check-truth does not decide which frame is the representative: it reads the app's choice. With `--ui` it also reads the corpus status cell's text and the current `Frame N` thumbnail in the workspace. It fails if the text disagrees with the API, or if a view shows no single pressed frame, because then check-truth itself is misreading the app. (The workspace showing no frame is expected only when every frame is dropped.)
- `assigned_frames` lists every frame carrying a human assignment, wherever it sits.

### Invariants

The `invariants` block states properties the app should hold and lists each violation by sample. A broken invariant fails `ok` and the exit code, like a TRUTH regression. One that could not be checked (`checked: false`, e.g. `views_agree` without `--ui`) reports `ok: null` and does not gate. `summary.invariants_failed` names the broken ones.

| Invariant | Breaks when |
|---|---|
| `status_matches_assignment` | A sample with no confirmed call (no assignment on a kept frame) is not shown as unreviewed, or a sample with one is shown as unreviewed or with a call it does not carry. An auto call alone is not confirmed (Jonathan, 2026-11-20: an unconfirmed sample reads "Unreviewed", whether or not Himalaya has an auto call). |
| `representative_kept` | The frame the sample status describes (`corpus_api`), or the frame the workspace opens on (`--ui`), is dropped while the sample has a kept frame. |
| `views_agree` | `corpus_api` and the workspace are on different frames of one sample. Needs `--ui`. |

`checked` is false for `views_agree` without `--ui`, and for both view invariants when a workspace read failed (that read is also a gate problem).

The representative is one rule everywhere: the selected frame if it is kept, else the newest kept frame, else none (`frontend/src/lib/exposure.ts`, `src/representative.jl`, held together by `frontend/test/fixtures/representative-exposure.json`). The workspace opens on it. Switching frames in the workspace and assigning there can still break `status_matches_assignment`, since the corpus status describes the representative.

## Evidence in a PR

Attach to every PR that touches something a user can see:

1. `bin/verify doctor` and `bin/verify check-truth` output (the `summary` block is enough when both pass; name any invariant your change breaks or fixes).
2. Screenshots of each flow you touched, taken with `bin/verify screenshot` after driving it with the CLI.
3. One line per Feature Map entry you used, saying whether it was still accurate.

## Gotchas

- The demo DB in `scratch/demo/` persists between runs, so anything you click or `api`-write stays. Undo what you changed, or `bin/verify seed --reset` (slow: it regenerates the data). Never commit `scratch/`.
- The browser is one shared headless Chromium (CDP port 9333, profile in `scratch/verify/profile`). Each command reconnects to the same tab, so commands compose across calls. `bin/verify close` ends it. If something else already listens on the CDP port, the CLI refuses to attach; set `VERIFY_CDP_PORT`.
- Settling means 500 ms with no request in flight (SSE, fonts and in-page `blob:`/`data:` URLs excluded), no visible boneyard overlay, no `aria-busy="true"` and no `role="progressbar"`. Actions track the requests they start. A standalone `wait-settle` only sees requests that start after it connects, so prefer letting the action settle. A failed settle lists the requests still in flight. `open`, the actions and `wait-settle` also fail when a request they saw returned 4xx/5xx, because the app can show nothing while the server rejected the change. Limits of that check:
  - A request that starts more than 500 ms after the action goes quiet is missed by that action and shows up on whichever command runs next.
  - `--no-settle` skips it; the result says `httpChecked: false`, and a later `wait-settle` only sees requests that start after it connects.
  - A 5xx fails the command even if the app retries it to success. That is deliberate: a 5xx is a server bug, not a normal answer.
- On a name collision the CLI refuses rather than guessing. Use `--exact` (for example `--name Im3m --exact` against "Im3m, in assignment"), `--within aside` or `--nth`.
- Ports: `DEMO_API_PORT`, `DEMO_UI_PORT` and `VERIFY_CDP_PORT` move everything, including `seed`. Mocked Playwright (`npm run e2e`) also wants :5173, so `bin/verify stop` before running it.
- The browser runs with the back/forward cache off. With it on, each page left by a full navigation stayed cached with its `/api/events` stream open, and the sixth `open` hung on Chromium's six-connections-per-host limit. In-app navigation (clicks, keys) never hit this.
- Speed checks: the default demo is small and local, so it hides costs that grow with sample count or network delay. Use the large demo (docs/demo.md, "A beamtime-sized demo") and set `VERIFY_THROTTLE=<mbit>mbit,<ms>ms` (for example `20mbit,30ms`), which every command applies to the page through CDP. Object results echo `throttle` so the evidence says which conditions it ran under. Chrome's latency is a minimum round trip that overlaps server time, so it is not added on top of a slow response.
- In cloud sessions the frontend's Playwright has no matching browser, so the CLI falls back to the preinstalled `/opt/pw-browsers/chromium-*`. Set `VERIFY_CHROMIUM` to force a binary.

## Keeping it true

The CLI has its own tests: `npm run test:verify` in `packages/HimalayaUI/frontend` (unit tests plus a real-browser test against a stub page). Change the CLI only with a test.

The Feature Map is maintained, not frozen. When a PR changes how a feature is reached or driven, update its entry in the same PR, re-running the recipe to confirm it. When a recipe fails because the UI changed, fix the entry rather than working around it.
