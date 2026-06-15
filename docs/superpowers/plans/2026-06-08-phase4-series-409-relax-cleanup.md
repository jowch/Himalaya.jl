# Series 409-Relax + Conflict/Dead-Hook Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Relax the series-commit optimistic-concurrency 409 gate to last-write-wins (LWW), and remove the now-dead conflict UI + experiment-scoped picker / recently-used hooks — fulfilling the "no conflict UI" redesign decision (`docs/redesign-notes.md` 2026-06-03).

**Architecture:** This is **Plan 6a** — the self-contained backend + cleanup slice of the Series-builder cutover (Phase-4 strategy step 6). It lands on the *existing legacy* builder page (which keeps working on LWW) and is independent of the greenfield builder page (Plan 6b). Three coordinated changes that must land on one branch (atomic so no transition-window `ConflictError` goes unhandled — but note the 409 only ever fired when a hash was *both* sent and mismatched, so the changes are order-independent within the branch): (1) Julia backend stops gating the commit; (2) frontend stops sending `expected_content_hash` and deletes the conflict modal; (3) frontend drops three orphaned query hooks.

**Tech Stack:** Julia/Oxygen.jl + SQLite backend; React 18 + TypeScript (strict, `exactOptionalPropertyTypes`), TanStack Query, Zustand; Julia stdlib `Test`; Vitest.

---

## Provenance ledger

**MODIFY (backend):** `packages/HimalayaUI/src/routes_series.jl` (commit handler) + the backend test asserting the 409.

**MODIFY (frontend):** `src/lib/series/buildSeriesCommitBody.ts` (stop sending hash), `src/api.ts` (`CommitSeriesPlateBody` type; `ConflictError` + the 409-parse if orphaned), `src/components/App.tsx` (drop the modal mount), `src/lib/queue/conflictBridge.ts` (drop the series arm), `src/state.ts` (`pendingConflict`/`setPendingConflict` if orphaned), `src/lib/series/seriesDraft.ts` + `seriesDraftFactories.ts` (`baseHash` if orphaned), `src/queries.ts` (delete 3 hooks).

**DELETE (frontend):** `src/components/SeriesCommitConflictModal.tsx` + `src/components/ConflictModalShell.tsx` + their tests; the 3 dead-hook tests.

**KEEP (do NOT touch):** `compute_series_content_hash` / `current_series_content_hash` (still write `content_hash` for the `series_plate_committed` post_state + future fork/stale checks); the legacy `SeriesBuilderPage` / `SeriesRecipeEditor` (replaced later in Plan 6b — they keep working on LWW); the queue/idempotency infra; `useCorpusPickerSamples` / `useCorpusSampleTags` (carry — used by scoping).

---

## Cross-cutting notes (read before starting)

- **The Julia backend suite is slow (5–10 min) and has 2 pre-existing failures** from the core-cache dev-link gotcha ([[feedback_himalayaui_core_devlink]]) — those are NOT regressions. Per [[feedback_teammate_bg_jobs_reap_on_idle]], the **orchestrator runs the full Julia gate**, not the implementer subagent. The backend implementer runs only the **targeted** series-route test file (fast) to prove its change; the orchestrator runs the full suite as the final gate, capturing to a file and grepping (ignore the 2 known dev-link errors).
- **Orphan-gated removal:** several removals (`ConflictError`, `pendingConflict`, `baseHash`) are "delete IF grep shows no remaining consumer after the primary deletion." Always grep first; if a consumer remains, leave the symbol and note it. Never `git add -A` — stage named paths only. Every commit's exact last line: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

## Task 1: Backend — relax the commit 409 gate to LWW

**Files:**
- Modify: `packages/HimalayaUI/src/routes_series.jl` (the `@post "/api/series/{id}/commit"` handler, currently lines ~228–256)
- Modify: the backend test asserting the commit 409 (locate by grep — see Step 1)

- [ ] **Step 1: Locate + update the failing test (TDD, inverted)**

Find the backend test(s) asserting the commit conflict:
```bash
cd packages/HimalayaUI
grep -rn "expected_content_hash\|409\|conflict\|current_hash" test/ | grep -i series
```
The current test asserts: *commit with a stale `expected_content_hash` → 409*. Invert it: assert that a commit with a **mismatched/stale `expected_content_hash` now succeeds (200) and applies the members (LWW)**, and that a commit with **no** `expected_content_hash` still succeeds. Delete any assertion that a 409 body / `:error => "conflict"` is returned from this route. Keep the 404-on-missing-series assertion. (The exact testset is in a `routes_series`/`series` test file — read it, rewrite the conflict testset, don't guess.)

- [ ] **Step 2: Run the targeted test; verify it fails**

```bash
# Run ONLY the series-routes test file (fast — do not run the full suite here):
julia --project=. -e 'using Pkg; Pkg.test()' # WRONG — too slow.
```
Instead run the single test file directly (locate it from Step 1), e.g.:
```bash
julia --project=packages/HimalayaUI packages/HimalayaUI/test/<series_routes_test_file>.jl 2>&1 | tail -30
```
Expected: FAIL — the stale-hash commit currently returns 409, not 200. (If the test file isn't directly runnable standalone due to harness setup, note that and rely on the orchestrator's full-suite gate to confirm red→green; still write the test.)

- [ ] **Step 3: Relax the handler**

In `routes_series.jl`, the commit handler currently reads:
```julia
        expected_hash = haskey(body, :expected_content_hash) &&
                        body.expected_content_hash !== nothing ?
                        String(body.expected_content_hash) : nothing

        return with_idempotency(db, req) do
            # Existence (404) before the conflict check (409) — HTTP semantics.
            # No author gate (architecture decision 3).
            if !series_exists(db, id)
                return _json_error(404, "series not found")
            end
            # Optimistic-concurrency check (NOT the author gate): the stored
            # hash must match the client's expected_content_hash, else 409.
            current_hash = current_series_content_hash(db, id)
            if expected_hash !== nothing && current_hash !== expected_hash
                current_state = fetch_series_with_plate(db, id)
                return HTTP.Response(409, ["Content-Type" => "application/json"],
                    JSON3.write(Dict(
                        :error         => "conflict",
                        :current_hash  => current_hash,
                        :current_state => current_state,
                    )))
            end

            members_payload = [_series_member_payload(db, m) for m in body.members]
```
Replace it with (delete the `expected_hash` extraction + the `current_hash` read + the 409 branch; keep the 404 + the members payload):
```julia
        return with_idempotency(db, req) do
            # Existence (404) only — the optimistic-concurrency 409 gate was
            # relaxed to last-write-wins (no conflict UI; docs/redesign-notes.md
            # 2026-06-03). `content_hash` is still written by the committed event
            # (post_state + future fork/stale checks), so compute_series_content_hash
            # / current_series_content_hash stay defined and used in the event layer.
            if !series_exists(db, id)
                return _json_error(404, "series not found")
            end

            members_payload = [_series_member_payload(db, m) for m in body.members]
```
Do NOT remove the `body.members` validation above this block, the `apply_event!(... "series_plate_committed" ...)`, the `out = fetch_series_with_plate(db, id)`, or the `_enqueue_broadcast_from_result!(... post_state = out)` that follow. Do NOT touch `compute_series_content_hash` / `current_series_content_hash` definitions.

- [ ] **Step 4: Run the targeted test; verify it passes** — same command as Step 2 → PASS (stale-hash commit now 200; members applied).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_series.jl packages/HimalayaUI/test/<series_routes_test_file>.jl
git commit -m "$(cat <<'EOF'
Series commit: relax the optimistic-concurrency 409 gate to last-write-wins

The no-conflict-UI decision (docs/redesign-notes.md 2026-06-03): the commit no
longer compares the client's expected_content_hash against the stored hash.
404-on-missing-series stays; content_hash is still written by the committed
event for post_state + future fork/stale checks.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

> ORCHESTRATOR GATE (not the implementer): after this task, run the full Julia suite once, capture + grep, and confirm only the 2 known dev-link failures remain:
> `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-6a.out 2>&1` then grep for `Test Summary`/`Error`/`Fail`.

---

## Task 2: Frontend — stop sending `expected_content_hash` + delete the conflict modal

**Files:**
- Modify: `src/lib/series/buildSeriesCommitBody.ts`, `src/api.ts`, `src/components/App.tsx`, `src/lib/queue/conflictBridge.ts`, `src/state.ts`, `src/lib/series/seriesDraft.ts`, `src/lib/series/seriesDraftFactories.ts`
- Delete: `src/components/SeriesCommitConflictModal.tsx`, `src/components/ConflictModalShell.tsx`, and their tests
- Modify/Delete tests: `test/buildSeriesCommitBody.test.ts` (or wherever it lives), `test/SeriesCommitConflictModal.test.tsx`, `test/ConflictModalShell.test.tsx`, `test/conflictBridge.test.ts`

- [ ] **Step 1: Stop sending the hash (test-first)**

`src/lib/series/buildSeriesCommitBody.ts` currently sets `expected_content_hash` from `draft.baseHash` (line ~42). Update its test first: assert the returned body has **no** `expected_content_hash` key (use `expect(body).not.toHaveProperty("expected_content_hash")`). Run it → FAIL. Then edit `buildSeriesCommitBody.ts` to stop adding the field (return just `{ members }`). Re-run → PASS. Remove `expected_content_hash` from the `CommitSeriesPlateBody` interface in `src/api.ts`.

- [ ] **Step 2: Delete the conflict modal + shell + App mount**

```bash
git rm src/components/SeriesCommitConflictModal.tsx src/components/ConflictModalShell.tsx
git rm test/SeriesCommitConflictModal.test.tsx test/ConflictModalShell.test.tsx
```
In `src/components/App.tsx`, remove the `<SeriesCommitConflictModal />` mount + its import (the Explore located it near App.tsx:100 — verify by grep).

- [ ] **Step 3: Remove the conflict bridge's series arm**

In `src/lib/queue/conflictBridge.ts`, remove the `series_commit` / `ConflictError → setPendingConflict` dispatch arm (the Explore: lines ~6–50). Update `test/conflictBridge.test.ts` accordingly (remove the series-commit-dispatch case; if the bridge has no other arm, the file may delete entirely — grep for other `conflictBridge` consumers first).

- [ ] **Step 4: Orphan-gated removals (grep before each)**

After Steps 1–3, grep for remaining consumers and remove each symbol ONLY if orphaned:
```bash
grep -rn "pendingConflict\|setPendingConflict" src test        # if no consumers → remove the slot + setter from state.ts
grep -rn "ConflictError" src test                               # if only the api.ts 409-parse references it → remove ConflictError + the 409→throw parse in the commit api fn
grep -rn "baseHash" src test                                    # if only seriesDraft type/factory reference it → remove baseHash from SeriesDraft + fromSeries
```
For each: if a consumer remains (e.g. the queue uses `pendingConflict` generically, or another mutation throws `ConflictError`), LEAVE it and note the surviving consumer in your report. Update the corresponding tests for whatever you remove.

- [ ] **Step 5: Gate**

```bash
npx tsc --noEmit -p tsconfig.build.json     # 0
npm run lint:design                          # 0
# guard: no dangling references to deleted symbols
grep -rn "SeriesCommitConflictModal\|ConflictModalShell" src test && echo DANGLING || echo clean
npm test > /tmp/vitest-6a-t2.out 2>&1; tail -15 /tmp/vitest-6a-t2.out   # whole suite green
```

- [ ] **Step 6: Commit**

```bash
git add src/lib/series/buildSeriesCommitBody.ts src/api.ts src/components/App.tsx \
        src/lib/queue/conflictBridge.ts src/state.ts src/lib/series/seriesDraft.ts \
        src/lib/series/seriesDraftFactories.ts test/  # name the exact modified test files
git commit -m "$(cat <<'EOF'
Series commit: stop sending expected_content_hash, delete the conflict modal

Frontend half of the LWW relax — buildSeriesCommitBody no longer sends the hash,
and the SeriesCommitConflictModal / ConflictModalShell / conflict-bridge series
arm (and any now-orphaned pendingConflict / ConflictError / baseHash) are removed.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Frontend — drop the dead picker / recently-used hooks

The experiment-scoped picker + recently-used hooks have NO non-test consumers (confirmed: the builder picker decision is DROP; corpus-wide `useCorpusPickerSamples`/`useCorpusSampleTags` are the carry). Remove them.

**Files:**
- Modify: `src/queries.ts` (remove `usePickerSamples` ~937–943, `useSampleTags` ~921–927, `useRecentlyPickedExposures` ~907–915), `src/api.ts` (remove the api functions they call, if orphaned)
- Delete/modify their tests

- [ ] **Step 1: Confirm zero non-test consumers (grep)**

```bash
for h in usePickerSamples useSampleTags useRecentlyPickedExposures; do
  echo "== $h =="; grep -rn "\b$h\b" src --include=*.tsx --include=*.ts | grep -v "src/queries.ts"
done
```
Each must show only test files (or nothing). If a `src/**` non-test consumer appears, STOP and report it — do not delete a hook still in use. (Note: `useSampleTags` is the EXPERIMENT-SCOPED variant — do NOT confuse it with `useCorpusSampleTags`, which scoping uses and is CARRY.)

- [ ] **Step 2: Remove the hooks + their api functions**

Delete the three hook definitions from `src/queries.ts`. Grep the api functions they call (e.g. `pickerSamples(experimentId)`, `sampleTags(experimentId)`, `recentlyPickedExposures(...)`) in `src/api.ts`; remove any now-orphaned (grep-confirm no other caller). Delete or trim their tests:
```bash
grep -rln "usePickerSamples\|useSampleTags\b\|useRecentlyPickedExposures" test/
git rm <each test file that exclusively tests these hooks>   # or trim the relevant cases from a shared file
```

- [ ] **Step 3: Gate**

```bash
npx tsc --noEmit -p tsconfig.build.json     # 0
npm run lint:design                          # 0
grep -rn "usePickerSamples\|useRecentlyPickedExposures" src test && echo DANGLING || echo clean
npm test > /tmp/vitest-6a-t3.out 2>&1; tail -15 /tmp/vitest-6a-t3.out   # whole suite green
```

- [ ] **Step 4: Commit**

```bash
git add src/queries.ts src/api.ts test/   # name exact modified/removed test paths
git commit -m "$(cat <<'EOF'
Series: drop the dead experiment-scoped picker + recently-used hooks

usePickerSamples / useSampleTags (experiment-scoped) / useRecentlyPickedExposures
have no consumers after the Compare retirement; the builder picker uses the
corpus-wide useCorpusPickerSamples. Removed with their orphaned api fns + tests.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Final gate (orchestrator)

- [ ] Full frontend suite: `npm test` (whole vitest — green), `npx tsc --noEmit -p tsconfig.build.json` (0), `npm run lint:design` (0), `npm run build` (passes).
- [ ] **Full Julia suite (orchestrator runs this — slow, 5–10 min):** `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-6a-final.out 2>&1`; grep `Test Summary` / `Fail` / `Error`; confirm only the 2 known dev-link failures remain (no new failures from the relax).
- [ ] Final whole-change review (the project `himalaya-reviewer` for the Julia/SQLite half + `frontend-reviewer` for the TS half).
- [ ] Do NOT run `finishing-a-development-branch` — the greenfield branch stays unmerged until the whole rebuild lands.

## Deferred to Plan 6b (the greenfield builder page)

- The greenfield unified-"Compose" `src/print/pages/SeriesBuilderPage.tsx` (wiring `BuilderRail` + `SeriesPlate` + `MemberList` + the recipe editor + add-sample picker + Save/Commit), the route repoint, deleting the legacy `SeriesBuilderPage` / `SeriesRecipeEditor` / `SeriesBuilderRail` / `MultiTracePlot` / etc., and the e2e migration. **Open architecture decision for 6b:** the mode model (eager auto-draft vs lazy-draft-on-first-edit vs preserve read/edit gate) — to be settled before drafting 6b.
