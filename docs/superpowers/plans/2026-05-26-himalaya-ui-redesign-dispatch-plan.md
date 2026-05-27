# HimalayaUI redesign — agent-team dispatch plan

**Status:** v1, 2026-05-26. The *execution* layer on top of the issue map ([`2026-05-17-himalaya-ui-redesign-issue-map.md`](2026-05-17-himalaya-ui-redesign-issue-map.md)). The issue map sliced six phases into 31 issues with a dependency DAG and a file-contention graph (its §3). This doc assigns the **remaining open issues** to a concrete agent team, defines the per-issue lifecycle, and pins the land-order rules that keep parallel worktrees from colliding.

**Scope note:** ~⅓ of the plan has already shipped. This doc plans only what is open as of 2026-05-26.

---

## 1. Where the work actually stands (2026-05-26)

**11 of 31 issues are merged**, and — critically — the issue map's worst contention point is already behind us.

**Merged:** I0.1, I0.2, I0.3 (Phase 0 corpus routes) · I1.1, I1.2, I1.3, I1.4, I1.5 (Phase 1 through the loupe) · I2.1, I2.2 (series schema + routes) · I3.2 (`SeriesMember` type migration).

**Landed-but-issue-open:** the series event-kind cluster **I2.3 / I2.4 / I2.5** shipped via PR #198 — all six `series_*` dispatcher branches in `events.jl`, the three mutators (`saveSeries` / `commitSeriesPlate` / `deleteSeries`), the `applyRemoteToCache` cases, and `mutatorRegistry` wiring are on `main`. The GitHub issues remain open (trailing contract-test layers). **Treat I2.7 as the verification that closes them** (§4, Lane A).

> **Why this reshapes the parallelism story.** Issue-map §3 named three switch files — `events.jl`, `applyRemoteToCache.ts`, `mutatorRegistry.ts` — as the high-contention hotspot that forced single-ownership of the I1.3/I2.3/I2.4/I2.5 cluster and dropped Wave C's *effective* concurrency to ~4–5. **All of that has merged.** The remaining open work has far less file contention; only two collision points survive, both late (§5).

### 1.1 Audit-verified status (2026-05-26 — four-agent fan-out + spot-checks)

GitHub open/closed state is **unreliable** — five issues are code-complete but never closed. Ground truth from the codebase:

| Issue | GH state | TRUE status | Evidence |
|---|---|---|---|
| I2.3 #166 (created/recipe/deleted) | open | **CODE DONE** (unclosed) | dispatcher + mutators + `applyRemoteToCache` + `resolveMutatorForEvent` + round-trip tests all present (PR #198) |
| I2.4 #167 (plate_committed) | open | **CODE DONE** (unclosed) | `post_state` handler + `synthesizeFromSse` present |
| I2.5 #168 (pins) | open | **CODE DONE** (unclosed) | five-layer wiring, `entity_type='user'` |
| I2.6 #169 (batch-tags route) | open | **CODE DONE** (unclosed) | `routes_samples.jl:215` |
| I2.7 #170 (native round-trip test) | open | **CODE DONE** (unclosed) | round-trips for all kinds in `test_routes_series.jl` (created L499, recipe L562, committed L642, pins L697) |
| I1.6 #162 (culling) | open | **NOT STARTED** | `ContactSheetRow.tsx:69` comment: "selection (#162) … wire in separately"; thumbnails inert |
| I1.7 #163 (Inspect cutover) | open | **NOT STARTED** | `InspectPage.tsx` + `/inspect` routes + `activePage:"inspect"` all live |
| **I3.1 #171 (migration)** | open | **NOT STARTED** | `migrate_comparisons_to_series!` only referenced in `db.jl` comments — never implemented |
| I3.3 #173 (folio) | open | **NOT STARTED** | no `useSeriesList` hook, no `/series` page (key shape stubbed `queries.ts:94`) |
| I3.4 #174 (scoping) | open | **NOT STARTED** | no `/series/new`, no scoping modal |
| I3.5a #175 (builder surface) | open | **NOT STARTED** | no `/series/:id` page |
| I3.5b #176 (builder mutations) | open | **PARTIAL** | `saveSeries`/`commitSeriesPlate` mutator scaffolding present (#198); UI wiring, `nextOptimisticId` placeholders, permalink resolution unbuilt — gated behind I3.5a |
| I3.6 #177 (Compare cutover) | open | **NOT STARTED** | `Compare.tsx`, `routes_comparisons.jl`, `comparisons.jl` all live |
| I4.1–I4.4 #178–181 | open | **NOT STARTED** | no `/sample/:sampleId` route; `hoveredQ` absent; `IndexPage` + three-card + `/index` routes live |
| I5.1–I5.3 #182–184 | open | **NOT STARTED** | `AppShell`/`TabRocker`/`WorkspaceGrid` live; persist `version:3` **no `migrate`**; queue `schema_version:2` |

**Done & closed (trusted):** Phase 0 (I0.1–I0.3) · Phase 1 foundation (I1.1, I1.2, I1.3, I1.4, I1.5) · Phase 2 schema/routes (I2.1, I2.2) · I3.2 (SeriesMember type).

> **Green-signal caveat (Lane A gate).** I confirmed the series code *exists and is wired*, but my standalone run of `test_routes_series.jl` errored on 14/27 testsets — `UndefVarError: with_test_server` (defined in `test_http.jl`, which `runtests.jl` includes first). A **harness artifact of running one file in isolation, not a regression** (pure-logic testsets passed; route testsets need the shared server helper). **A definitive green requires a full-harness run** — folded into Round 0 below.

### 1.2 Corrected remaining critical path (4 issues)

The two predecessors the issue map drew (I2.7 → I3.1) are **already satisfied** — I2.7's fold is proven by existing round-trip tests, and I3.1 is unblocked *now*. The real remaining critical path collapses to:

```
I3.1 → I3.6 → I5.1 → I5.2 → I5.3   (I5.x serial; I3.6 also waits on I3.3/I3.4/I3.5b)
```

**Front-load review on I3.1** (`migrate_comparisons_to_series!` — the synthesize-events migration is the most-reviewed mechanism in the spec; master plan §6.1, risk register). It is now Lane A's *first* real task, not its third.

### 1.3 PR-history reconciliation (three-agent fan-out)

The redesign is exactly **PRs #185–#200, no open/draft PRs** (everything ≤#184 is pre-redesign). Reconciling PR ↔ issue ↔ code:

- **The five open Phase-2 issues are pure bookkeeping debt.** #198 wired I2.3/I2.4/I2.5 and #191 delivered I2.6, but neither used GitHub `Closes #N` syntax (`closingIssuesReferences` empty) — so #166–169 stayed open despite merging. **I2.7 (#170) is genuinely done too**: its exact deliverable (series via `POST`, empty view rows, re-fold) is `test_routes_series.jl:499`, added by #198 — it just never had its own PR or close link. → **Action (Round 0): after a green suite, manually close #166, #167, #168, #169, #170.**
- **#198 quietly landed I3.3/I3.5b *scaffolding* the issue map never listed under them** — confirming the `[[series-event-cluster-frontend-scaffolding]]` memory. Already merged: `api.ts` Series/SeriesSample/SeriesMember types, `queries.ts` key shapes (`queryKeys.series`/`seriesList`/`seriesPins`), `mutatorRegistry` series cases, the `saveSeries`/`deleteSeries`/`commitSeriesPlate` mutators, `applyRemoteToCache` series handlers, + frontend tests (~550 lines). **This shrinks two lanes' work:**
  - **I3.3 (folio)** — only the `useSeriesList`/`useSeries` *read hooks* + the `/series` page remain (key shapes + api types already exist).
  - **I3.5b (builder mutations)** — the mutator/cache layer is done; only the *UI-facing wiring* (negative `nextOptimisticId` placeholders, spinner-commit hookup, slug-permalink series resolution) remains, on top of I3.5a.
- **No frozen machinery was touched.** `comparison*` tables, `comparisons.jl`, `routes_comparisons.jl`, the `comparison_*` dispatcher branches are untouched across #185–200 (#188 only read `comparisons.jl`). Clean.
- **PR #200 is incidental core-Himalaya churn**, not redesign — `src/sharpness.jl` SG edge handling (closed #199). Lanes branching from `main` simply carry it; no interaction with redesign files.

---

## 2. Team composition

**5 agents + the orchestrator (me) + you (final merge authority).**

| Agent | Role | Issues (in order) | Domain |
|---|---|---|---|
| **Orchestrator** | me (this session's lead) | — assigns, gates plans, presents to you | — |
| **Lane A** | Backend / critical path | I2.6 → I2.7 → **I3.1** | Julia / SQLite / events |
| **Lane B** | Phase 1 closeout | I1.6 → I1.7 | Frontend (sample table) |
| **Lane C** | Phase 3 series UI | I3.3, I3.5a → I3.5b, I3.4, **I3.6** | Frontend (series) |
| **Lane D** | Phase 4 focus workspace | I4.1 → I4.2 → I4.3 → I4.4 | Frontend (focus) |
| **Reviewer** | Dedicated PR-review agent | runs the `review-pr` loop on every PR, every lane | all domains |
| **(Phase 5)** | folded into Lane A once it frees up | I5.1 → I5.2 → I5.3 | Frontend (serial) |

**Why 4 implementation lanes, not 6.** The unblocked issues partition into four disjoint file-ownership groups. A 5th implementation agent could split Lane C's folio (I3.3) from the builder (I3.5a/b), but both append to `queries.ts` / `api.ts` — you'd buy concurrency with rebase cost. **Hold Lane C as one agent unless I3.5a's converged plan shows the builder is heavy enough to justify the split** (then spawn Lane C2, with an explicit `queries.ts` land-order against Lane C1).

**Lane A lands I2.6 first.** I2.6 (batch-tags route) is what unblocks Lane C's I3.4 (scoping). Sequencing it ahead of I2.7 keeps Lane C from stalling.

**Phase 5 is serial — one owner.** I5.1→I5.2→I5.3 cannot parallelize. Lane A (its critical-path work done after I3.1) inherits them, or the Reviewer agent promotes to closeout. No new agent needed.

---

## 3. Per-issue lifecycle (the loop every lane runs)

For **each** issue, the owning lane runs this loop. One issue = one worktree = one PR.

```
1. WORKTREE      git worktree add .claude/worktrees/<issue> -b <issue>   (superpowers:using-git-worktrees)
                 first-time setup in the new tree: /worktree-setup
2. PLAN          write the file-level TDD plan            (superpowers:writing-plans)
                 — only for issues marked "detailed plan: needed" (§3.1); "direct"
                   issues execute straight from the issue card + master-plan §11.
3. PLAN REVIEW   fan out the relevant project reviewers ON THE PLAN (§3.2)
                 iterate plan ↔ reviewers until a CLEAN pass (converged)
4. PLAN APPROVAL bring the converged plan to the orchestrator → then to YOU
                 ── no code is written before your approval ──
5. EXECUTE       TDD: failing test → minimal impl → green → commit, per step
                 (superpowers:test-driven-development, subagent-driven-development)
6. PR            open the PR                              (superpowers:request-pr-review author half)
7. PR REVIEW     the dedicated Reviewer agent runs the review-pr loop to convergence
8. MERGE         YOU do the final merge                  (superpowers:finishing-a-development-branch)
9. CLEANUP       remove the worktree; lane advances to its next issue
```

**Plans are reviewer-converged *before* you see them** (your answer 2): the fan-out + iterate in step 3 happens inside the lane; you receive a plan that already survived its domain reviewers, not a first draft.

### 3.1 Which open issues need a detailed plan

**Detailed plan: needed** — I3.1, I3.3, I3.4, I3.5a, I3.5b, I4.1, I4.2, I4.3.
**Direct** (execute from card + master-plan §11) — I1.6, I1.7, I2.6, I2.7, I3.6, I4.4, I5.1, I5.3.
**Borderline** — I5.2 (persist-version bumps): the issue map marks it "needed" because a bump without `migrate` silently wipes user state — keep it as **needed**.

### 3.2 Reviewer → issue mapping (for both the plan review in step 3 and the PR review in step 7)

| Issue(s) | Reviewers to fan out |
|---|---|
| I2.6, I2.7 | `himalaya-reviewer` |
| **I3.1** (migration, event replay) | `himalaya-reviewer` **+** `queue-reviewer` |
| I1.6, I1.7, I3.3, I3.4, I3.5a, I4.1, I4.2, I4.4, I3.6 | `frontend-reviewer` |
| I3.5b (builder mutations), I4.3 (q-link channel) | `frontend-reviewer` **+** `queue-reviewer` |
| I5.1, I5.2, I5.3 | `frontend-reviewer` (+ `queue-reviewer` for I5.2's queue `schema_version` bump) |

`saxs-physics-reviewer` is **out of scope** for this batch — none of the open issues touch peak-finding / scoring / phase physics.

---

## 4. The four lanes in detail

### Lane A — Backend / critical path  (Julia)
- **I2.6** `POST /api/samples/tags/batch` — *direct*. One `with_idempotency` tx, N inserts, atomic. **Land this first** (unblocks I3.4). Touches `routes_samples.jl` + `server.jl` (append-only).
- **I2.7** native-series `rebuild_views_from_log!` round-trip — *direct*. New Julia test file (no contention). **This is the verification that effectively closes #166/#167/#168** — it proves the pure-replace branches fold from empty.
- **I3.1** `migrate_comparisons_to_series!` — *detailed plan*. **Critical path; highest review priority.** Synthesize-events copy: construct payloads → raw-`INSERT` `user_actions` (no broadcast) → fold through the I2.3/I2.4 dispatcher branches → sentinel row written last in-transaction. Placed after `migrate_compare_view_choices!` / `_relax_nullability!`. Edits the `migrate_schema!` sequence body (wave-separated from I2.1, already merged).

### Lane B — Phase 1 closeout  (frontend)
- **I1.6** culling wiring — *direct*. Multi-select batch reject + representative pick via existing `useSetExposureStatus` / `useSelectExposure` / `useAddExposureTag`. Largest "direct" issue (optimistic-update implications on the batch path).
- **I1.7** Inspect cutover — *direct*. Delete `InspectPage` + Inspect E2E spec; `/inspect*`→`/samples` redirect; **retire the `activePage:"inspect"` branch in `state.ts`**. Playwright mocked spec for cull/batch-reject/loupe-flip/tag.

### Lane C — Phase 3 series UI  (frontend)  ← heaviest lane
- **I3.3** folio — *detailed plan*. `/series` masonry, `useSeriesList()`. Appends to `queries.ts`.
- **I3.5a** builder surface — *detailed plan*. `/series/:id` visual surface; composes `MultiTracePlot` render core, `FigureExportControls`, `QNumInput`. **Step one: re-audit I2.2's route shapes against the builder's needs** (master plan §6.3) — drift spawns a fast-follow patch owned by I2.2's author. Composes the render core; does **not** re-edit `MultiTracePlot`/`MemberTraceLayer` (those land as I3.2 follow-ups).
- **I3.5b** builder mutations — *detailed plan*, **needs I3.5a**. Optimistic recipe edits (negative `series_samples` placeholder ids) vs spinner plate-commit; slug-permalink `series` resolution; `ConflictModal` reuse. Intentionally decoupled from I3.4 (no tag-criteria predicate in v1).
- **I3.4** scoping — *detailed plan*, **needs I2.6** (Lane A). `/series/new`; writes `sample_tags` via the batch route with `source='scoping'`, stored as structured `(key,value)`. `useFocusTrap` on the modal.
- **I3.6** Compare cutover — *direct*. Delete `Compare.tsx`, `routes_comparisons.jl`, `comparisons.jl`, Compare E2E specs; `/compare*`→`/series`; **retire `activePage:"compare"` in `state.ts`**; migrated-DB round-trip test. `comparison*` tables + dispatcher branches **stay forever**.

### Lane D — Phase 4 focus workspace  (frontend)
- **I4.1** Zustand wiring shim — *detailed plan*. `/sample/:sampleId`→`activeSampleId` (route-param→Zustand sync, or prop-drill — pick one). Adds field to `state.ts`.
- **I4.2** layout — *detailed plan*, **needs I4.1**. Trace hero / co-resident detector / phase rail / Notes drawer. **Carried-over trace/index interaction tests must pass unchanged** (regression floor).
- **I4.3** q-link — *detailed plan*, **needs I4.2**. Ephemeral `hoveredQ` Zustand field + action (excluded from `partialize`); rotation-aware `DetectorImage` ring overlay; must not bypass the `QNumInput` focus-gate.
- **I4.4** Index cutover — *direct*. Delete `IndexPage` + three-card composition + Index E2E specs; `/index*`→`/sample/...`; **retire `activePage:"index"` in `state.ts`**.

---

## 5. The two surviving file-collision constraints

Almost all remaining open work is new-file or disjoint. Only two collisions survive, both late and both manageable because **each issue rebases on latest `main` at worktree-creation and merges serially** (one PR at a time, you merge):

1. **`state.ts` `activePage` narrowing — I3.6 (Lane C) vs I4.4 (Lane D).** Both retire an `activePage` branch, on different tracks, in the same window. **Rule: land-order — whichever cutover is ready second rebases onto the first.** The Reviewer agent enforces (won't approve the second cutover PR until it's rebased on the first). Field-adds (I4.1, I4.3) and other-branch retirements (I1.7's `inspect`) are append/disjoint — no special handling.
2. **`queries.ts` / `api.ts` within Lane C** (I3.3, I3.4, I3.5a all append hooks/keys). **Solved by single ownership of Lane C** — one agent serializes them. (Only relevant if Lane C is split into C1/C2; then C2 rebases.)

Everything else flagged in issue-map §3 (`events.jl`, `applyRemoteToCache.ts`, `mutatorRegistry.ts`, `routes_samples.jl`) is either already merged or append-only with the last PR rebasing.

---

## 6. Execution order (what to dispatch, when) — revised for audit-verified status

Because I2.6 and I2.7 are already done, Lane A skips straight to the migration, and **every lane's first real issue is unblocked now**.

**Round 0 — bookkeeping gate (orchestrator, before/parallel to Round 1):**
- Run the full HimalayaUI Julia suite once (`Pkg.test`, capture to file, grep — slow, ~5–10 min). It defines `with_test_server`, so it exercises the series route + round-trip tests properly.
- **If green:** close #166, #167, #168, #169, #170 (I2.3–I2.7 — code done, never `Closes`-linked; §1.3). This unblocks I3.1 *with confidence*.
- **If red:** triage before Lane A touches I3.1 (the migration folds through those exact dispatcher branches).

**Round 1 — dispatch immediately (4 lanes, all deps satisfied):**
- Lane A → **I3.1** (the migration; critical path; gate on Round 0 green).
- Lane B → **I1.6**.
- Lane C → **I3.3** (folio — lighter warm-up) and/or **I3.5a**; **I3.4 is also available now** (I2.6 done).
- Lane D → **I4.1**.

**Round 2 — unlocks as Round 1 merges:**
- Lane B: I1.6 → **I1.7** (Phase 1 done).
- Lane C: **I3.4**, then I3.5a → **I3.5b** (UI wiring on top of the existing mutator scaffolding).
- Lane D: I4.1 → **I4.2** → **I4.3**.

**Round 3 — cutovers (apply the §5 land-order):**
- Lane C: **I3.6** (needs I3.1, I3.3, I3.4, I3.5b).
- Lane D: **I4.4** (needs I4.2, I4.3).

**Round 4 — serial tail (single owner):**
- **I5.1** (needs I1.7, I3.6, I4.4) → **I5.2** → **I5.3**. No parallelism possible.

---

## 7. Open items for the orchestrator before Round 1

- [x] ~~Confirm I2.3–I2.7 status~~ — **RESOLVED (§1.1):** all code-complete, issues unclosed. Round 0 confirms green + closes them.
- [ ] **Run Round 0** (full suite) — the only remaining unknown is the green signal; standalone runs can't produce it (`with_test_server` lives in the harness).
- [ ] Decide whether to keep Lane C as one agent or split C1/C2 — defer until I3.5a's converged plan is in hand.
- [ ] Each lane runs `/worktree-setup` in its fresh worktree before first execution (npm install + Pkg.instantiate).
- [ ] Reviewer agent is briefed on the §5 land-order so it gates the second cutover PR.

---

## 8. Tail-state addendum (2026-05-27 — Rounds 1–3 complete)

**Merged:** all of Phases 0–2, plus I1.6, I1.7, I3.1, I3.2, I3.3, I3.4, I3.5a, I4.1, I4.2, I4.3. `main` @ `b516671`. **In flight:** I3.5b (#176) + I4.4 (#181), both executing under fresh lanes (review topology: lanes run `request-pr-review`, the human runs `review-pr`).

**Remaining after the in-flight pair:** I3.6 (#177) → I5.1 (#182) → I5.2 (#183) → I5.3 (#184), then the deferred follow-ups.

Three corrections to the body above, found in the 2026-05-27 tail survey:

1. **The `state.ts` `activePage` collision is 2-way (and the series stage was never in it).** Originally thought 3-way, but **I3.5b's execution (PR #213) proved the series stage is CorpusShell — URL-owned by react-router — and `useStateFromUrl`/`useUrlFromState` mount only in AppShell**, so `"series"` is NOT a `PageId` member (I3.5b dropped its planned union-widening as dead code; its `state.ts` change is append-only — `seriesDraft` slot + `ConflictError` type). The `activePage` union only ever held the three legacy AppShell pages: `"inspect"` (cut by I1.7 ✅), `"index"` (cut by I4.4), `"compare"` (cut by I3.6). So only **I4.4 removes `"index"`** and **I3.6 removes `"compare"`** — sequential rebases, clean. **Endgame:** after I4.4, `PageId = "compare"` (one member, the default/fallback); I3.6 removing it **empties the union**, at which point the `activePage` model is fully vestigial. A zero-member `PageId` breaks `coerceActivePage`/default — so I3.6 must either leave the model for **I5.1** (#182, scoped to exactly "retire the `activePage` model") to delete (cleaner — keeps I3.6 a cutover) or pull I5.1's model-removal forward. This is I3.6's headline open decision (relayed to lane-i36).

2. **I3.6's reviewer set should be himalaya-reviewer + queue-reviewer, not just frontend-reviewer (§3.2 understates it).** I3.6 is *direct* but deletes `routes_comparisons.jl` + `comparisons.jl` (Julia) and adds a `rebuild_views_from_log!` round-trip test on a **real migrated comparison-era DB**. That backend deletion + replay test needs the Julia/queue reviewers, not only frontend. `comparison*` tables + dispatcher branches stay forever.

3. **Follow-ups #207/#208/#209 are deferred until after I5.3** (human decision, 2026-05-27). They block nothing (all off critical path); deferring keeps the single human reviewer's load sequential, avoids #207's `e2e/live` spec-repoint colliding with the I4.4/I3.6 cutover spec churn, and gives I5.3's dead-code sweep a stable tree. Batch all three after the redesign closes. (#208 edits the shared render core the series builder composes; #209 completes I4.3's 2-of-3 q-link triple; #207 is corpus tag/rename UI + orphaned-spec repoint.)

**Phase 5 = one fresh lane, serial.** I5.1 → I5.2 → I5.3. I5.2 (#183) is *detailed-plan: needed* — a persist-version bump without a real `migrate` silently wipes user prefs; its queue `schema_version` bump's rationale is I3.6 (stale `comparison_*` ops must drop as a toast, not throw in `mutatorRegistry`). I5.3 is the build-green gate that closes the redesign.

## 9. Completion addendum (2026-05-27 — redesign COMPLETE)

**The redesign is fully merged.** `main` @ `1e0b274`. Final merge sequence (all cleared 2-round `review-pr` loops; lanes ran `request-pr-review`, the human ran `review-pr`, the orchestrator ran the `.claude/agents/` reviewer fan-out at plan stage):

| PR | Issue | What landed |
|---|---|---|
| #214 | I4.4 (#181) | Retire the Index page — Phase 4 cutover |
| #215 | I3.6 (#177) | Retire the Compare page — Phase 3 cutover |
| #216 | I5.1 (#182) | Retire the dual-navigation scaffolding + the `activePage` model |
| #218 | I5.2 (#183) | Persisted-state version bumps + migrations |
| #219 | I5.3 (#184) | Final dead-code sweep + build verification |

**The §8 union endgame played out as predicted:** the 2-way collision was correct, I3.6 emptied the `PageId` union to a `"none"` sentinel, and I5.1 deleted the whole `activePage` model. Final state: single `CorpusShell`, `persist` version 4 + a real `migrate` (no prefs wipe), queue `SCHEMA_VERSION` 3, dead code swept, `npm run build` green. The `comparison*` tables + `comparison_*` dispatcher/replay machinery are kept forever (verified intact by queue-reviewer at I5.3).

**Three execution-time deviations worth recording (all caught by the review loop / grep-gate):**
1. **I5.1 — `UtilityCluster` (theme toggle + multiplayer switch-user).** The frontend-reviewer endorsed re-homing it into `CorpusTopbar`, but the corpus-topbar mockup (`docs/superpowers/specs/2026-05-17-corpus-app-shell-design.md:121`) is wordmark + 3 stage-tabs + Beamtime chip only — no utilities. Design authority chose "match mockup, defer relocation" → `UtilityCluster` deleted; the switch-user affordance relocation is **#217** (filed, deferred).
2. **I5.3 — `peakCycle.ts` mis-classified as KEEP.** Plan claimed `MemberTraceLayer` kept it alive; frontend-reviewer proved its only consumer was the deleted `cyclePeakDisplayForMember` action → would have broken the build + left an orphan. Reclassified to DELETE; a report-only `ts-prune` discovery pass was added as the mechanical backstop.
3. **`Closes #N` backticks gotcha** — a backtick-wrapped `` `Closes #182` `` on PR #216 parsed as inline code and did NOT auto-close #182 (manual close needed). Subsequent PRs wrote it un-backticked. (Recorded in memory `pr-closes-link-convention`.)

**Remaining (deferred to the user's post-redesign review + feature-restoration passes — not started):** #217 (switch-user/theme home), #207 (corpus tag/rename UI + orphaned-spec repoint), #208 (heatmap + cross-trace peak-tracking render-core), #209 (reflections-table + q-link rows).
