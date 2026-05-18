# HimalayaUI redesign — issue map (parallelizable work breakdown)

> **For agentic workers:** this is the *decomposition layer*. The spec is the master plan ([`2026-05-17-himalaya-ui-redesign.md`](2026-05-17-himalaya-ui-redesign.md), v5 review-converged); this file slices its six phases into 31 individually-dispatchable issues with an explicit dependency DAG. To *execute* an issue, use `superpowers:subagent-driven-development` — each issue marked **detailed plan: needed** gets a per-issue file-level TDD plan written when it is picked up; issues marked **direct** can be executed from the master plan + the issue card.

**Goal:** turn the phase-sequenced master plan into a dependency graph of bite-sized issues so multiple agents can work concurrently.

**Status:** v2, 2026-05-17 — revised after a four-reviewer pass (backend · queue · frontend · decomposition-logic). v2 fixes: the critical path was miscomputed (two phantom edges, claimed 11 issues — actually 9); the `I2.7→I3.1` edge was missing; §3 missed `routes_samples.jl`, the Wave-E `state.ts`/`queries.ts` contention; I3.5 (builder) was oversized and is split into I3.5a/I3.5b. Reviewer history in §6.

**How to read an issue card:** `Depends on` must be *merged* before the issue starts — these edges are the **only hard constraints**. `Unblocks` is the reverse edge. Waves (§1) are *recommended scheduling cohorts*, not constraints. `Shared files` flags textual-conflict contention with other issues (see §3). `Detailed plan` says whether the issue needs its own file-level TDD plan before an agent executes it.

---

## 1. The dependency graph

```
                          ┌─────────────────────────────────────────┐
  WAVE A   I0.1 ─────────┐ │  I0.2   I0.3   (Phase 0 — all parallel) │
           (samples)     │ └─────────────────────────────────────────┘
                         │   I1.1 (shell)   I2.1 (schema)   I2.6 (batch-tags route)   ← no deps
                         │
  WAVE B   I1.2 (corpus query) ◄── I0.1        I2.2 (routes) ◄── I2.1
                         │
  WAVE C   I1.3 (add_tag) ◄── I1.2    I1.4 (contact sheet)  I1.5 (loupe) ◄── I1.1,I1.2
                         │   I2.3  I2.4  I2.5  (series event kinds) ◄── I2.1,I2.2
                         │   I3.2 (type-migration) ◄── I2.2
  WAVE D   I1.6 (cull) ─► I1.7 (P1 cutover)         I2.7 (native round-trip test) ◄── I2.2..I2.5
                         │
  WAVE E   ── TRACK 3 ──────────────────────────  ── TRACK 4 ──────────────────────
           I3.1 migration ◄── I2.3,I2.4,I2.7      I4.1 zustand-shim ◄── I1.1
           I3.3 folio  I3.4 scoping ◄── I3.2       I4.2 layout  I4.3 q-link
           I3.5a builder-surface  I3.5b builder-mut  I4.4 (P4 cutover)
           I3.6 (P3 cutover) ◄── I3.1,I3.3,I3.4,I3.5b
                         │
  WAVE F   I5.1 (retire dual-nav) ─► I5.2 (persist bumps) ─► I5.3 (sweep + build)
```

**Two long-lived parallel tracks.** After Wave A, **Phase 1** (sample table) and **Phase 2** (series backend) are fully independent — different files, different data model — and run as two tracks. After they land, **Phase 3** (depends on Phase 2) and **Phase 4** (depends on Phase 1) run as two more independent tracks. Phase 5 is the only strictly-serial tail.

### 1.1 Waves are scheduling cohorts, not constraints

The **only hard constraints are the per-issue `Depends on` edges.** Waves are a suggested batching for a dispatcher; an issue may start the moment its dependencies are merged regardless of wave label. In particular:

- **I1.1, I2.1, I2.6 have no dependencies** — they can begin in Wave A alongside Phase 0. They are drawn at Wave A above. The master plan's narrative "Phase 2 depends on Phase 0" is a *phase-ordering convention*; it does **not** manifest as any issue-level edge — no Phase 2 issue consumes a Phase 0 route. Phase 2 is genuinely independent of Phase 0.
- **Within Wave C, I2.2 must merge before I2.3/I2.4/I2.5 start** (they list it in `Depends on`). Wave C is therefore staged: {I1.3, I1.4, I1.5, I2.2, I3.2} then the {I2.3, I2.4, I2.5} cluster.
- **Within Wave D, I1.6 must merge before I1.7 starts** — that pair is a serial chain, not two parallel issues.

### 1.2 Critical path

The longest dependency chain — **9 issues** — recomputed from the edge list (not the phase order):

```
I2.1 → I2.2 → I2.3 → I2.7 → I3.1 → I3.6 → I5.1 → I5.2 → I5.3
```

(A co-equal 9-issue chain runs `I2.1→I2.2→I3.2→I3.5a→I3.5b→I3.6→I5.1→I5.2→I5.3`.) The critical path **starts at I2.1**, the series schema — not at any Phase 0 issue (§1.1). Phase 1 and Phase 4 have slack against it. Front-load reviewer attention on the Phase 2 event-kind issues (I2.3–I2.5), the I2.7 round-trip test, and the I3.1 migration — a defect there stalls the whole chain.

### 1.3 Wave / parallelism table

| Wave | Issues | Effective concurrency | Gate to next wave |
|---|---|---|---|
| A | I0.1, I0.2, I0.3, I1.1, I2.1, I2.6 | 6 | I0.1, I1.1, I2.1 merged |
| B | I1.2, I2.2 | 2 | I1.2, I2.2 merged |
| C | I1.3, I1.4, I1.5, I3.2 + the I2.3/I2.4/I2.5 cluster | nominal 7 → **~4–5** (I2.2 intra-wave gate; §3 single-owns the event-kind cluster) | per-track (see cards) |
| D | I1.6→I1.7 (serial), I2.7 | 2 | I1.7, I2.7 merged |
| E | Track 3 (I3.1, I3.3, I3.4, I3.5a, I3.5b, I3.6) ∥ Track 4 (I4.1, I4.2, I4.3, I4.4) | ~5–6 | all of E merged |
| F | I5.1 → I5.2 → I5.3 | 1 (serial) | — |

---

## 2. Issues

### Phase 0 — corpus foundation

#### I0.1 — `GET /api/samples` corpus route (+ `q_units`)
- **Wave A** · **Depends on:** — · **Unblocks:** I1.2
- **Files:** modify `packages/HimalayaUI/src/routes_samples.jl` (sibling of `routes_samples.jl:4`); register in `server.jl`.
- **Shared files:** `routes_samples.jl` with I1.3, I2.6 — append-only, wave-serialized (see §3).
- **Deliverable:** `GET /api/samples?experiment_id=` returns the whole corpus; projection **includes `q_units` per sample** (from the owning experiment's config). No schema change, no events.
- **Acceptance:** route returns the corpus with `q_units`; existing `GET /api/experiments/{id}/samples` untouched and green; Julia route test added.
- **Detailed plan:** direct (master plan §3.1 is sufficient).

#### I0.2 — `GET /api/sample-tags` corpus route
- **Wave A** · **Depends on:** — · **Unblocks:** I3.4 (scoping reads tags)
- **Files:** modify `packages/HimalayaUI/src/routes_picker.jl` (sibling of `routes_picker.jl:55`); register in `server.jl`.
- **Shared files:** `routes_picker.jl` with I0.3 — append-only, low conflict; land I0.2 then rebase I0.3.
- **Deliverable:** corpus-wide `GET /api/sample-tags`.
- **Acceptance:** route returns corpus tags; experiment-scoped sibling green; Julia route test added.
- **Detailed plan:** direct.

#### I0.3 — `GET /api/picker-samples` corpus route
- **Wave A** · **Depends on:** — · **Unblocks:** I3.4
- **Files:** modify `packages/HimalayaUI/src/routes_picker.jl` (sibling of `routes_picker.jl:71`); register in `server.jl`.
- **Shared files:** `routes_picker.jl` with I0.2 — rebase on I0.2.
- **Deliverable:** corpus-wide `GET /api/picker-samples`.
- **Acceptance:** route returns the corpus picker projection; experiment-scoped sibling green; Julia route test added.
- **Detailed plan:** direct.

### Phase 1 — sample table

#### I1.1 — App shell + hoisted `<Routes>` + nav-bridge
- **Wave A** · **Depends on:** — · **Unblocks:** I1.4, I1.5, I4.1
- **Files:** create the new shell component (topbar: corpus wordmark, three stage-tabs, Beamtime facet chip); modify `AppShell.tsx`, `TabRocker.tsx`, the app router; modify `state.ts` (nav-bridge + stale-`activePage` coercion).
- **Shared files:** `state.ts` with I4.1, I4.3, I3.6, I4.4, I5.1, I5.2 — wave-separated from most; see §3.
- **Deliverable:** one hoisted top-level `<Routes>` — old routes → old `AppShell` body, new routes → new shell body. `useGlobalShortcuts` + theme effect hoisted. URL↔`activePage` nav-bridge that **coerces** a stale persisted `activePage` (master plan §2.3, §4.1). The route *table* lives here; later issues register their route slot into it.
- **Acceptance:** both shells render under one router; nav-bridge tested (Vitest); no persist version bump (additive field reads `undefined`).
- **Detailed plan:** **needed** — the dual-router hoist and nav-bridge are the structural spine of Phases 1–4.

#### I1.2 — Corpus query layer
- **Wave B** · **Depends on:** I0.1 · **Unblocks:** I1.3, I1.4, I1.5
- **Files:** modify `queries.ts` (`queries.ts:119` has the per-experiment `useSamples`), `api.ts`, the query-key module.
- **Shared files:** `queries.ts` / `api.ts` / query-key module with I3.3, I3.5a — wave-separated; see §3.
- **Deliverable:** `useCorpusSamples()` with a **distinct** corpus query-key shape (`["corpus","samples"]`) — not the per-experiment `queryKeys.samples`.
- **Acceptance:** `useCorpusSamples()` fetches `GET /api/samples`; Vitest covers the hook + key.
- **Detailed plan:** direct.

#### I1.3 — Corpus `add_tag` six-layer change
- **Wave C** · **Depends on:** I1.2 · **Unblocks:** I1.7
- **Files:** `mutatorRegistry.ts` (`resolveMutator` + `resolveMutatorForEvent`), `applyRemoteToCache.ts` (`add_tag`/`remove_tag` branch), `lib/queue/mutators/trivial.ts` (mutator `onMutate`/`onSuccess`), `routes_samples.jl:71` (parameterize `source='manual'`).
- **Shared files:** `mutatorRegistry.ts`, `applyRemoteToCache.ts` — **high contention** with I2.3/I2.4/I2.5; `routes_samples.jl` with I0.1/I2.6 (see §3). I1.3's edit is a *modification* of the existing `add_tag` branch (binary→tri-scope), not an append — rebasing it onto the I2.3 cluster is not purely mechanical.
- **Deliverable:** corpus sample-tags round-trip — tri-scope discriminator keyed on **`experiment_id` payload-field presence** (corpus and experiment tags both arrive `entity_type='sample'`); the `applyRemoteToCache` `add_tag`/`remove_tag` branch made **tri-scope** (binary today); corpus query-key optimistic patch (master plan §4.2).
- **Acceptance:** corpus tag add/remove round-trips through the queue *and* a foreign-tab SSE frame; six-layer contract tests green.
- **Detailed plan:** **needed** — six-layer change, binary→tri-scope conversion touches two `resolve*` functions.

#### I1.4 — Contact sheet UI
- **Wave C** · **Depends on:** I1.1, I1.2 · **Unblocks:** I1.6, I1.7
- **Files:** create the contact-sheet view + row component; reuse `ThumbnailGallery` / `DetectorImage size="thumb"`.
- **Deliverable:** `/samples` corpus table — each row: identity, exposure-thumbnail strip, Kept count, **Tags** column (light free-form chips), Status (mockup `sample-table.html`).
- **Acceptance:** `/samples` renders the corpus; boneyard skeleton gates on `query.isLoading`; Vitest for the table.
- **Detailed plan:** **needed** — new surface, mockup translation.

#### I1.5 — Loupe UI
- **Wave C** · **Depends on:** I1.1, I1.2 · **Unblocks:** I1.7
- **Files:** create the loupe view (`/samples/loupe/:sampleId`); absorb `DetectorImageCard` / `SampleMetadataCard` internals.
- **Deliverable:** `DetectorImage size="full"` + exposure strip + metadata sidecar with a Sample-tags section. **Must not** hard-code file-per-exposure (deferred derived-exposure feature — master plan §4.2, risk register).
- **Acceptance:** loupe flips frames; Vitest; boneyard skeleton.
- **Detailed plan:** **needed** — new surface, derived-exposure constraint.

#### I1.6 — Culling wiring
- **Wave D** · **Depends on:** I1.4 · **Unblocks:** I1.7
- **Files:** wire into the contact sheet; use existing `useSetExposureStatus` / `useSelectExposure` / `useAddExposureTag`.
- **Deliverable:** reject-only culling, multi-select batch reject, representative pick.
- **Acceptance:** cull + batch-reject + representative-pick round-trip through the queue; Vitest.
- **Detailed plan:** direct (the largest "direct" issue — the multi-select batch path has optimistic-update implications; the hooks already exist).

#### I1.7 — Phase 1 E2E + Inspect cutover
- **Wave D** · **Depends on:** I1.3, I1.4, I1.5, I1.6 · **Unblocks:** I5.1
- **Files:** delete `InspectPage` + Inspect-only components + the Inspect E2E spec; add `/inspect*` → `/samples` redirect; retire the `activePage:"inspect"` branch (`state.ts`).
- **Shared files:** `state.ts` `activePage` union — see §3.
- **Deliverable:** Playwright mocked spec (cull / batch-reject / loupe-flip / tag); Inspect fully retired.
- **Acceptance:** master plan §4.3 exit criteria; full suite green.
- **Detailed plan:** direct.

### Phase 2 — series backend

#### I2.1 — Schema: `migrate_series!` + `schema_migrations` sentinel
- **Wave A** · **Depends on:** — · **Unblocks:** I2.2, I2.3, I2.4, I2.5, I3.1
- **Files:** modify the migration module — add `migrate_series!` (`CREATE TABLE IF NOT EXISTS`), register it in the `migrate_schema!` sequence, and add the `schema_migrations(name TEXT PRIMARY KEY, applied_at TEXT)` sentinel.
- **Shared files:** the `migrate_schema!` sequence body with I3.1 — wave-separated (A vs E); see §3.
- **Deliverable:** `series`, `series_members`, `series_samples`, `series_messages`, `series_pins` per master plan §5.1 — all `CHECK` constraints, `UNIQUE(series_id, position)`, all child FKs `ON DELETE CASCADE`, `series_samples.sample_id` `ON DELETE CASCADE`. `comparison*` tables stay in `SCHEMA` permanently.
- **Acceptance:** fresh DB and `migrate_schema!` on an existing DB both produce the tables; idempotent on re-run; Julia schema test.
- **Detailed plan:** **needed** — exact DDL, constraint placement.

#### I2.2 — `routes_series.jl` + `series.jl` skeleton
- **Wave B** · **Depends on:** I2.1 · **Unblocks:** I2.3, I2.4, I2.5, I2.7, I3.2, I3.3, I3.5a
- **Files:** create `routes_series.jl` (mirror `routes_comparisons.jl`, corpus-wide) and `series.jl` (logic adapted from `comparisons.jl`); register in `server.jl`.
- **Deliverable:** `GET/POST /api/series`, `GET /api/series/{id}`, `GET /api/series/{id}/forks`, `PATCH /api/series/{id}`, `POST /api/series/{id}/commit`, `DELETE /api/series/{id}`, messages + pins. **Fix** the `last_event_at` mixed-timestamp sort bug (`comparisons.jl:669`, issue #76). **Drop** the `is_author` gate (architecture decision 3). The event-emitting routes are stubbed until I2.3–I2.5 land. **This route contract is the frozen Phase 2 deliverable** I3.5a re-audits — a post-merge drift finding spawns a fast-follow patch issue owned by this issue's author.
- **Acceptance:** GET routes round-trip on empty `series` tables; `last_event_at` sort test; Julia route tests.
- **Detailed plan:** **needed** — route-by-route contract is the frozen Phase 2 deliverable Phase 3 audits against.

#### I2.3 — Event kinds: `series_created` + `series_recipe_updated` + `series_deleted`
- **Wave C** · **Depends on:** I2.1, I2.2 · **Unblocks:** I2.7, I3.1, I3.5b
- **Files:** `events.jl` (dispatcher branches), `applyRemoteToCache.ts`, `mutatorRegistry.ts` (both `resolve*`), the relevant *new* mutator file(s), paired contract tests.
- **Shared files:** `events.jl`, `applyRemoteToCache.ts`, `mutatorRegistry.ts` — **high contention** with I1.3/I2.4/I2.5 (see §3).
- **Deliverable:** three kinds. `series_created` — **upsert** the parent `series` row (route mints the AUTOINCREMENT id), then `DELETE`+`INSERT` `series_samples` from the full snapshot. `series_recipe_updated` — pure-replace `series_samples`. `series_deleted` — one-line `DELETE FROM series` (four-table cascade). Wire `entity_type='series'`. All six layers each (master plan §5.2).
- **Acceptance:** six-layer contract tests per kind; dispatcher folds correctly from empty.
- **Detailed plan:** **needed** — covers all three kinds; pure-replace-vs-`comparison_submitted`-shaped is the load-bearing distinction. May be implemented as three commits under one owner.

#### I2.4 — Event kind: `series_plate_committed`
- **Wave C** · **Depends on:** I2.1, I2.2 · **Unblocks:** I2.7, I3.1, I3.5b
- **Files:** `events.jl`, `applyRemoteToCache.ts` (genuinely new `post_state` handler), `mutatorRegistry.ts` (case + `synthesizeFromSse`), a new mutator file, contract tests.
- **Shared files:** `events.jl`, `applyRemoteToCache.ts`, `mutatorRegistry.ts` — **high contention** (see §3).
- **Deliverable:** `series_plate_committed` — `DELETE`+`INSERT` `series_members` (members carry no ids — dispatcher mints them), set `state='committed'`, compute `content_hash` from the plate only. First comparison-shaped event to carry a `post_state` envelope (`fetch_series_with_plate` projection); `synthesizeFromSse` so the deferred resolves with the route response shape (master plan §5.2).
- **Acceptance:** six-layer contract tests; `post_state` round-trips; `content_hash` NULL while draft.
- **Detailed plan:** **needed** — the hardest event kind (`post_state` + id-less payload).

#### I2.5 — Event kinds: `series_pinned` / `series_unpinned`
- **Wave C** · **Depends on:** I2.1, I2.2 · **Unblocks:** I2.7
- **Files:** `events.jl`, `applyRemoteToCache.ts`, `mutatorRegistry.ts`, contract tests.
- **Shared files:** `events.jl`, `applyRemoteToCache.ts`, `mutatorRegistry.ts` — **high contention** (see §3).
- **Deliverable:** five-layer pair (no `post_state`, no mutator). Stored `entity_type='user'`, `entity_id=user_id` (the `comparison_pinned` precedent) — **outside** the `entity_type="series"` round-trip fold.
- **Acceptance:** five-layer contract tests; pins replay correctly under `entity_type='user'`.
- **Detailed plan:** direct (the `comparison_pinned` precedent is a close template).

#### I2.6 — `POST /api/samples/tags/batch` route
- **Wave A** · **Depends on:** — · **Unblocks:** I3.4
- **Files:** modify `routes_samples.jl`; register in `server.jl`.
- **Shared files:** `routes_samples.jl` with I0.1, I1.3 — append-only, wave-serialized (see §3).
- **Deliverable:** one `with_idempotency` tx, N tag inserts — non-optional atomic boundary. Reuses the existing `add_tag` kind; scoping (Phase 3) passes `source='scoping'`.
- **Acceptance:** batch insert is atomic (mid-batch failure rolls back fully); Julia route + idempotency test.
- **Detailed plan:** direct.

#### I2.7 — Native-series `rebuild_views_from_log!` round-trip test
- **Wave D** · **Depends on:** I2.2, I2.3, I2.4, I2.5 · **Unblocks:** I3.1
- **Files:** Julia test file for the round-trip.
- **Deliverable:** create a series natively via `POST /api/series`, **empty** the target series' view rows (`series`/`series_members`/`series_samples`), re-fold the log, assert the rebuilt state matches — proves the pure-replace branches fold from empty.
- **Acceptance:** round-trip test green; `/api/series*` round-trips through queue + SSE; `series` tables ship empty (master plan §5.3).
- **Detailed plan:** direct.

### Phase 3 — series stage UI

#### I3.1 — `migrate_comparisons_to_series!` data migration
- **Wave E (Track 3)** · **Depends on:** I2.1, I2.3, I2.4, **I2.7** · **Unblocks:** I3.6
- **Files:** modify the migration module — add `migrate_comparisons_to_series!` and **insert it into the `migrate_schema!` sequence body after** `migrate_compare_view_choices!` and `migrate_compare_relax_nullability!` (this edits the dispatch sequence, not a pure append).
- **Shared files:** the `migrate_schema!` sequence body with I2.1 — wave-separated; see §3.
- **Why I2.7:** the migration's correctness criterion (master plan §6.1 step 4) is that synthesized payloads fold through the *pure-replace* dispatcher branches and survive a `rebuild_views_from_log!` re-fold. I2.7 is the issue that *proves the pure-replace branches fold from empty* — I3.1 must not be built on an unverified fold.
- **Deliverable:** the event-sourced copy — sentinel gate; construct `series_created` + `series_plate_committed` payloads; raw-`INSERT` the `user_actions` rows (`client_op_id`/`client_id` NULL, no broadcast); fold payloads through the I2.3/I2.4 dispatcher branches; sentinel row written last in the same transaction; copy messages + pins (master plan §6.1).
- **Acceptance:** migration is idempotent, runs before `serve` accepts connections; messages/pins copied; no broadcast fan-out.
- **Detailed plan:** **needed** — on the critical path; the synthesize-events mechanism is the most-reviewed part of the spec.

#### I3.2 — `ComparisonMember` → `SeriesMember` type migration
- **Wave C** · **Depends on:** I2.2 · **Unblocks:** I3.3, I3.4, I3.5a
- **Files:** `MultiTracePlot`, `MemberTraceLayer`, `multiTraceAdapter`, `lib/figure-export/marks/multiTraceExportMarks.ts`, `MemberMetaGutter`, `MemberMetaRow`, `lib/comparison/{coloring,labels,draftFactories}`. **Not** `lib/comparison/yBands.ts` (pure numeric, no type).
- **Deliverable:** introduce `SeriesMember`; migrate the render pipeline's input type. Optional fields as `T | null` (mirror `ComparisonMember` in `api.ts`) for `exactOptionalPropertyTypes`.
- **Acceptance:** `tsc --noEmit` green; the render pipeline's existing Vitest passes on the new type.
- **Detailed plan:** **needed** — wide type change, foundational for the Phase 3 UI.

#### I3.3 — Folio UI
- **Wave E (Track 3)** · **Depends on:** I3.2, I2.2 · **Unblocks:** I3.6
- **Files:** create `/series` folio view; `useSeriesList()` hook (`queries.ts`).
- **Shared files:** `queries.ts` / query-key module with I1.2, I3.5a — see §3.
- **Deliverable:** corpus masonry of saved series (mockup `series-folio.html`).
- **Acceptance:** `/series` lists series corpus-wide; boneyard skeleton; Vitest + Playwright mocked spec.
- **Detailed plan:** **needed** — new surface.

#### I3.4 — Scoping UI
- **Wave E (Track 3)** · **Depends on:** I3.2, I2.6 · **Unblocks:** I3.6
- **Files:** create `/series/new` scoping view; `useFocusTrap` on the scoping modal.
- **Deliverable:** machine-proposed ordering variable + confirm-and-build gate; confirmations write `sample_tags` via the batch route with `source='scoping'`, stored as a structured `(key,value)` pair (master plan §6.2).
- **Acceptance:** scoping proposes → confirms → writes structured tags; Vitest + Playwright mocked spec.
- **Detailed plan:** **needed** — new surface, structured-tag contract.

#### I3.5a — Builder surface (layout + render core)
- **Wave E (Track 3)** · **Depends on:** I3.2, I2.2 · **Unblocks:** I3.5b
- **Files:** create `/series/:id` builder view; compose `MultiTracePlot` render core, `FigureExportControls`, `QNumInput` (offset dock).
- **Deliverable:** the visual surface of `series-builder.html` — offset waterfall default, heatmap toggle, collapsible rail, peak-tracking annotation layer. **Step one:** re-audit the Phase 2 series route response shapes against the builder's needs (master plan §6.3) — treat I2.2's contract as frozen; a drift finding spawns a fast-follow patch on I2.2. If the peak-tracking layer needs `MultiTracePlot`/`MemberTraceLayer` *internal* changes, those land as an I3.2 follow-up, not inside this issue (the builder *composes* the render core, it does not re-edit it).
- **Acceptance:** the builder renders a series in all three representations; Vitest + figure-export specs.
- **Detailed plan:** **needed** — the largest visual surface.

#### I3.5b — Builder mutations (recipe edits + plate commit + permalink)
- **Wave E (Track 3)** · **Depends on:** I3.5a, I2.3, I2.4 · **Unblocks:** I3.6
- **Files:** the builder view's mutation wiring; slug-permalink hooks gain `series` resolution; `ConflictModal` reuse.
- **Deliverable:** optimistic recipe edits (`series_recipe_updated` — negative `series_samples` placeholder ids via `nextOptimisticId()`); spinner plate-commit (`series_plate_committed` — no optimistic write); figure export round-trip. The builder is **intentionally decoupled from I3.4** — recipe edits do not read scoping's tags (the tag-criteria predicate is deferred, master-plan open decision 2), so there is no I3.4→I3.5b edge.
- **Acceptance:** scope → build → edit → commit → export round-trips; Vitest; Playwright mocked spec.
- **Detailed plan:** **needed** — optimistic/spinner split, negative-id placeholder machinery.

#### I3.6 — Phase 3 E2E + migrated-DB round-trip test + Compare cutover
- **Wave E (Track 3)** · **Depends on:** I3.1, I3.3, I3.4, I3.5b · **Unblocks:** I5.1
- **Files:** delete `Compare.tsx` + Compare-only components, `routes_comparisons.jl`, `comparisons.jl`, the Compare E2E specs; `/compare*` → `/series` redirect; retire the `activePage:"compare"` branch (`state.ts`).
- **Shared files:** `state.ts` `activePage` union with I4.4 (same wave, different track) — see §3.
- **Deliverable:** a `rebuild_views_from_log!` round-trip for `entity_type="series"` on a row migrated from a **real comparison-era DB** (empties view rows, re-folds, asserts equality); folio→scoping→builder Playwright specs.
- **Acceptance:** master plan §6.3 exit criteria; `comparison*` tables + dispatcher branches **kept**.
- **Detailed plan:** direct.

### Phase 4 — focus workspace

#### I4.1 — Zustand wiring shim
- **Wave E (Track 4)** · **Depends on:** I1.1 · **Unblocks:** I4.2
- **Files:** modify the router / a sync shim — `/sample/:sampleId` → `activeSampleId`; `state.ts` if implemented as a Zustand field + action.
- **Shared files:** `state.ts` with I4.3, I3.6, I4.4 (Wave E) — see §3.
- **Deliverable:** route-param→Zustand sync (or the prop-drill alternative — pick one, master plan §7) so `PlotCard`/`IndicesCard`/`PhasePanel` (which read `active*Id` directly from Zustand) work under a URL-routed surface.
- **Acceptance:** the route param drives `activeSampleId`; Vitest.
- **Detailed plan:** **needed** — the Zustand-direct wiring decision affects three reused components.

#### I4.2 — Focus workspace layout
- **Wave E (Track 4)** · **Depends on:** I4.1 · **Unblocks:** I4.3, I4.4
- **Files:** create `/sample/:sampleId` view (registers its route slot into I1.1's hoisted `<Routes>`); reuse `TraceViewer`, `PlotCard`, `PhasePanel`, `MillerPlot` unchanged.
- **Deliverable:** trace as hero plate, co-resident detector image, phase call as rail, Notes as margin/drawer (mockup `focus-workspace.html`).
- **Acceptance:** indexing one sample works as before in the new layout; the **carried-over trace/index interaction tests pass unchanged** (regression floor); boneyard skeletons.
- **Detailed plan:** **needed** — new surface, regression-floor constraint.

#### I4.3 — q-link
- **Wave E (Track 4)** · **Depends on:** I4.2 · **Unblocks:** I4.4
- **Files:** add an ephemeral `hoveredQ` Zustand field + named action (`state.ts`, excluded from `partialize`); make the `DetectorImage` ring overlay rotation-aware.
- **Shared files:** `state.ts` with I4.1, I3.6, I4.4 (Wave E) — additive field; see §3.
- **Deliverable:** hovering a trace peak / detector ring / reflection row lights all three; the hover channel must not bypass the `QNumInput` focus-gate (master plan §7).
- **Acceptance:** q-link lights all three surfaces; Vitest + Playwright mocked spec.
- **Detailed plan:** **needed** — rotation-aware overlay + focus-gate interaction.

#### I4.4 — Phase 4 E2E + Index cutover
- **Wave E (Track 4)** · **Depends on:** I4.2, I4.3 · **Unblocks:** I5.1
- **Files:** delete `IndexPage` + the three-card composition + the Index-workspace E2E specs; `/index*` → `/sample/...` redirect; retire the `activePage:"index"` branch (`state.ts`).
- **Shared files:** `state.ts` `activePage` union with I3.6 (same wave, different track) — see §3.
- **Deliverable:** the Index workspace fully retired.
- **Acceptance:** master plan §7 exit criteria; full suite green.
- **Detailed plan:** direct.

### Phase 5 — final cutover (serial)

#### I5.1 — Retire dual-nav scaffolding
- **Wave F** · **Depends on:** I1.7, I3.6, I4.4 · **Unblocks:** I5.2
- **Files:** delete the old `AppShell` shell body, `WorkspaceGrid`, `TabRocker`, the Zustand `activePage` model, the nav-bridge.
- **Deliverable:** single-shell app; dual-nav model gone.
- **Acceptance:** all four new surfaces reachable; no `activePage` references remain.
- **Detailed plan:** direct.

#### I5.2 — Persist-version bumps
- **Wave F** · **Depends on:** I5.1 · **Unblocks:** I5.3
- **Files:** `state.ts` (remove `activePage` from `partialize`, bump `persist` version **with a real `migrate`** preserving `username`/`theme`/`tutorialSeen`); `lib/queue/persistence.ts` (bump queue `schema_version` — independent counter).
- **Deliverable:** both version counters bumped deliberately; the Zustand bump ships a `migrate`.
- **Behavioral note:** the queue `schema_version` bump's *correctness* rationale traces to **I3.6** (which retired `routes_comparisons.jl` + the comparison mutators) — after that, a stale `sessionStorage` `comparison_*` op must drop as a toast, not throw in `mutatorRegistry`. The `Depends on: I5.1` edge is wave-ordering; the logical reason is I3.6.
- **Acceptance:** a pre-cutover persisted state survives the migrate (keeps user prefs); a stale queued `comparison_*` op drops as a toast, not a `mutatorRegistry` throw; Vitest.
- **Detailed plan:** **needed** — a bump without `migrate` silently wipes every user's state (risk register).

#### I5.3 — Dead-code sweep + build green
- **Wave F** · **Depends on:** I5.2 · **Unblocks:** —
- **Files:** repo-wide sweep.
- **Deliverable:** final dead-code removal; `npm run build` (`tsc --noEmit` + vite) green. `comparison*` tables + dispatcher branches **never deleted**.
- **Acceptance:** build green; full suite green; master plan §8.
- **Detailed plan:** direct.

---

## 3. Shared-file contention (read before dispatching a wave)

True logical independence does not mean conflict-free merges. These files are edited by multiple "parallel" issues:

| File | Issues that touch it | Conflict kind | Handling |
|---|---|---|---|
| `events.jl` (`update_view_for_event!` dispatcher) | I2.3, I2.4, I2.5 | adjacent `case` arms | **Single owner** for I2.3–I2.5, or land in order I2.3→I2.4→I2.5, each rebasing. |
| `applyRemoteToCache.ts` | I1.3, I2.3, I2.4, I2.5 | I2.3–I2.5 *append* new branches; I1.3 *modifies* the existing `add_tag`/`remove_tag` branch | Cluster I2.3–I2.5 single-owner. I1.3 is a different track — rebase it last; its rebase is not purely mechanical (see I1.3 card). |
| `mutatorRegistry.ts` (`resolveMutator` + `resolveMutatorForEvent`) | I1.3, I2.3, I2.4, I2.5 | adjacent `case` arms in two functions | Same as `events.jl`. |
| `state.ts` | I1.1, I1.7, I3.6, I4.1, I4.3, I4.4, I5.1, I5.2 | I1.7/I3.6/I4.4 each narrow the `activePage` union; I4.1/I4.3 add fields; I5.1/I5.2 delete the model | **I3.6 and I4.4 are the same wave on different tracks and both narrow `activePage`** — give them a land-order (last rebases), same as the event-kind cluster. I4.1/I4.3 field-adds are append-only. I1.7 (Wave D) and I5.1/I5.2 (Wave F) are wave-separated. |
| `queries.ts` / `api.ts` / query-key module | I1.2, I3.3, I3.5a | append-only hook/key additions | I1.2 (Wave B) wave-separated; I3.3 and I3.5a (Wave E Track 3) — last-merged rebases, low risk. |
| `routes_samples.jl` | I0.1, I1.3, I2.6 | I0.1/I2.6 append routes; I1.3 edits line ~71 | Append-only / disjoint edit; wave-serialized (A→A→C); last-merged rebases. |
| `routes_picker.jl` | I0.2, I0.3 | append-only | Land I0.2, rebase I0.3. |
| `server.jl` (`register_routes!`) | I0.1, I0.2, I0.3, I2.2, I2.6 | flat append-only list | Low risk; last-merged rebases. |
| the `migrate_schema!` sequence body | I2.1, I3.1 | ordered insert into the dispatch sequence (not pure append) | Wave-separated (A vs E); I3.1 inserts after the `migrate_compare_*` steps — see I3.1 card. |
| series mutator files (`lib/queue/mutators/series*.ts`) | I2.3, I2.4 | **new files** — no contention | Each event-kind issue adds its own mutator file. Only `trivial.ts` is edited in-place, by I1.3 alone. |

**Recommendation for the Wave C event-kind cluster:** assign I2.3, I2.4, I2.5 to **one agent in sequence**, or to three agents with an explicit land-order. They are *logically* three issues — kept separate for reviewability and detailed-planning — but they collide on three switch-statement files. Six agents rebasing three files costs more than one agent doing three vertical slices. I1.3 lives on the Phase 1 track and rebases onto whatever the Phase 2 cluster produced.

**Recommendation for Wave E `state.ts`:** I3.6 (Track 3) and I4.4 (Track 4) both narrow the `activePage` union — give them a land-order. The per-surface `activePage`-branch retirement is mandated by master plan §2.3 (a stale persisted `activePage` otherwise renders an empty `PageBody`); it is kept, not deferred to I5.1.

`★ Insight ─────────────────────────────────────`
This is the difference between a *task graph* and a *file graph*. The dependency DAG in §1 is the logical/data graph — what must be *correct* before what. §3 is the file graph — what must be *merged* before what. They disagree: I2.3/I2.4/I2.5 are siblings with no logical edge between them, yet they serialize on `events.jl`. A decomposition that only draws the logical graph hands agents a guaranteed three-way merge conflict and calls it "parallel."
`─────────────────────────────────────────────────`

**Cross-cutting (master plan §11):** no sample-delete route exists today, so there is no issue for it. If one is ever added, it must emit `series_recipe_updated` for affected series — the `series_samples.sample_id` cascade is otherwise invisible to event replay. Any future sample-delete issue inherits this constraint.

---

## 4. Dispatch summary

- **31 issues**, 6 waves. Nominal peak concurrency ~7 (Wave C); **effective ~4–5** once the I2.2 intra-wave gate and §3 single-ownership of the event-kind cluster are applied.
- **Issue grain is deliberately uneven:** small additive routes (I0.2, one GET route) sit beside irreducible surfaces (I3.5a, the builder). I3.5 was split into I3.5a/I3.5b to bring the largest surface closer to the median; I2.3 bundles three event kinds under one owner *because* §3 single-owns the cluster — splitting it would only add switch-file contention.
- **18 issues need a per-issue detailed file-level TDD plan** before an agent executes them (marked **detailed plan: needed**); **13 are direct** from the master plan + their issue card.
- **Critical path** is 9 issues — `I2.1→I2.2→I2.3→I2.7→I3.1→I3.6→I5.1→I5.2→I5.3` (§1.2). It runs Phase 2 → 3 → 5 and starts at the series schema, not Phase 0. Phase 1 and Phase 4 have slack — schedule reviewer attention on the critical path first.
- The master plan's §11 per-phase constraint ledger is the **acceptance backstop** for every issue in that phase — a detailed plan that violates a §11 constraint is wrong even if this card looks satisfied.

## 5. Next step

Write the **detailed file-level plans** for the two widest fan-out unblockers — both have zero dependencies and start in Wave A: **I1.1** (the app-shell/dual-router spine, unblocks the Phase 1 + Phase 4 tracks) and **I2.1** (the series schema, head of the critical path). Each detailed plan follows `superpowers:writing-plans` proper: bite-sized TDD steps, exact code, exact commands.

---

## 6. Reviewer history

- **v1 → v2:** four reviewers (himalaya-reviewer · queue-reviewer · frontend-reviewer · a decomposition-logic reviewer) reviewed the issue map as a decomposition. Findings resolved in v2:
  - **Critical path miscomputed** (decomposition-logic blocker): v1's path had two phantom edges (`I0.1→I2.1`, `I2.7→I3.1`) and claimed 11 issues. v2 recomputes from the edge list — 9 issues, starting at I2.1 (§1.2); §1.1 states Phase 2 has no issue-level edge to Phase 0.
  - **Missing `I2.7→I3.1` edge** (himalaya blocker; decomposition-logic B2; queue): I3.1's migration reuses the pure-replace fold that I2.7 proves. v2 adds I2.7 to I3.1's `Depends on` and makes I2.7's `Unblocks` concrete.
  - **`state.ts` Wave-E collision** (frontend blocker): I3.6 and I4.4 both narrow the `activePage` union in the same wave on different tracks. v2's §3 `state.ts` row captures all eight writers and prescribes an I3.6/I4.4 land-order.
  - **§3 missed `routes_samples.jl`** (all three domain reviewers, independently): added — I0.1/I1.3/I2.6.
  - **I3.5 oversized** (frontend; decomposition-logic): split into I3.5a (surface) + I3.5b (mutations).
  - **§3 missed `queries.ts`/`api.ts` Wave-E contention** (frontend; decomposition-logic): added.
  - Smaller fixes: wave table now shows effective concurrency (I2.2 intra-wave gate, I1.6→I1.7 serial); `migrate_schema!` sequence-edit noted on I3.1; series mutators noted as new files; I5.2's behavioral dependency on I3.6 noted; I3.5b's intentional non-edge to I3.4 stated; the cross-cutting sample-delete constraint recorded in §3.
- **v2 confirmation round:** the decomposition-logic reviewer independently re-derived the full edge list, confirmed the graph acyclic, recomputed the longest chain at exactly 9 (both claimed critical paths real and co-equal), and verified the `I2.7→I3.1` edge, the I3.5 split, whole-graph edge/reverse-edge symmetry, the §3 additions, and the 18/13/31 detailed-plan count. **v2 confirmed clean** — ready to dispatch.
