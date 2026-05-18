# HimalayaUI redesign — master implementation plan

**Status:** v5, 2026-05-17 — **review-converged.** Hardened across five four-reviewer passes (backend/SQLite · mutation queue · frontend · architecture); the fifth came back clean from all four. The sequencing spine for translating the redesign into shippable engineering work.
**Design record:** [`docs/redesign-notes.md`](../../redesign-notes.md) — the workflow model, the converged three-surface model, the three architecture decisions, and the five canonical mockups (`docs/redesign-mockups/{sample-table,focus-workspace,series-builder,series-folio,series-scoping}.html`).
**Companion design docs:** [event-log](../../event-log.md), [mutation-queue](../../mutation-queue.md), [contract-testing](../../contract-testing.md).

This is a **master plan**: the architecture target, the incremental-delivery principle, and six execution phases. Phases 0–2 are specified to actionable detail; Phases 3–5 at blueprint level. Each phase gets its own detailed plan in this directory when it is picked up — §11 lists the constraints those detailed plans must honour.

**v5 changes (fourth reviewer pass):** architecture and frontend reviews came back **clean**. The backend review found one blocker — §5.2's `series_created` dispatcher branch must **upsert** the parent `series` row (the route mints it to capture the AUTOINCREMENT id; a plain `INSERT` collides on the live path), with only the child tables pure-replace. The queue review pinned the wire `entity_type='series'` for the four non-pin `series_*` events and identified `resolveMutatorForEvent` as a *third* per-kind switch site. v5 corrects §5.2, §2.2, §5.1 (`content_hash` folds the plate only), and §11. Full reviewer history in §12.

---

## 1. What the redesign is

The app becomes **three surfaces, one per workflow stage**, replacing today's three-card Index workspace + Inspect page + Compare page:

| New surface | Replaces | Mockup |
|---|---|---|
| Sample table (contact sheet + loupe) | Inspect page | `sample-table.html` |
| Focus workspace (index one sample) | three-card Index | `focus-workspace.html` |
| Series stage — folio + scoping + builder | Compare page | `series-{folio,scoping,builder}.html` |

Three architecture decisions (redesign-notes §5, "architecture revisited") reshape the data model:

1. **Experiment → provenance facet.** The top-level scope becomes the whole corpus of samples; experiment becomes a filter facet. Additive route work, no migration of existing scoping.
2. **Saved comparison → series = recipe + plate.** A series gains a living *recipe* (criteria) layer above today's frozen *plate* (snapshot).
3. **Series editing un-gated.** The `is_author` gate is dropped; the event log already carries attribution.

---

## 2. The shape of the work

### 2.1 Incremental-delivery principle — and what "additive" actually means

**The app is never broken between phases.** Every new surface ships at a *new route*. Each phase ends by retiring the surface it replaces, *after* the replacement is accepted — so within a phase the old surface is the rollback path, and across phases nothing co-equal competes.

"Additive" holds for **new routes** (`server.jl`'s `register_routes!` is a flat list) and **new tables / new columns** (the migration system — `SCHEMA` const + idempotent `migrate_schema!` on every `open_db`, with `ALTER … ADD COLUMN` in duplicate-tolerant try/catch and `CREATE … IF NOT EXISTS`).

It does **not** hold for **renaming or dropping a table**. A rename breaks every route still reading the old name, and — critically — **event-log replay**: `rebuild_views_from_log!` re-folds historical events through dispatcher branches that `INSERT INTO` named tables. Therefore:

> **The `comparison*` tables (`comparisons`, `comparison_members`, `comparison_messages`, `comparison_pins`) are permanent.** They are event-replay machinery. The `series` model is built as **new tables**, data **copied** forward. The historical `comparison_*` dispatcher branches stay frozen and bit-faithful — never aliased, never repointed.

There is **no feature flag** — route coexistence plus per-surface cutover is the mechanism. "Old route still reachable by URL" is the accepted rollback story; the schema migration is forward-only (a conscious trade-off, §10).

### 2.2 What carries over — and how the claim must be qualified

A large validated core is reused. The honest split:

- **Behaviourally validated, reused as-is:** the trace-editing interaction model and the figure-export pipeline. `TraceViewer.tsx`, `PlotCard.tsx` (`computeFit`, `QNumInput`, `useAutoPickExposure`), `PhasePanel.tsx`, `MillerPlot.tsx`, all of `lib/figure-export/`, the peak/index queue mutators. **Do not rewrite these.**
- **Reused, with new wiring:** `PlotCard`, `IndicesCard`, `PhasePanel` read `active*Id` *directly from Zustand* — no props. Phase 4 must either keep a route-param→Zustand sync or prop-drill `sampleId`/`exposureId`. They carry over behaviourally, not literally.
- **Reused, with new feature surface:** the q-link (Phase 4) adds a hovered-q highlight channel that `TraceViewer` and `DetectorImage` must emit/consume; `DetectorImage`'s ring overlay must be rotation-aware (the `ResizeObserver` auto-rotate trap).
- **Reused, with a typed-contract migration:** `MultiTracePlot`, `MemberTraceLayer`, `multiTraceAdapter`, `lib/figure-export/marks/multiTraceExportMarks.ts`, `lib/comparison/{coloring,labels,draftFactories}`, `MemberMetaGutter`, `MemberMetaRow` are typed on `ComparisonMember`. Phase 3 migrates them to a `SeriesMember` type — the *render* pipeline carries over, its *input type* does not. (`lib/comparison/yBands.ts` is a pure numeric module — no `ComparisonMember` type — and carries over unchanged.)
- **Frameworks extended, not modified:** the queue / SSE / event-log core (`useQueueMutation`, `replayCoordinator`, `with_idempotency`, `InTransaction`, post-commit broadcast, optimistic-id negativity). New work adds new event kinds + mutators *and extends three per-kind switch sites* (`update_view_for_event!`, `applyRemoteToCache.ts`, and `resolveMutatorForEvent`) plus the mutator-scope discriminator (`mutatorRegistry`).
- **Unchanged:** `DetectorImage`/`ThumbnailGallery` rendering; ingestion (per-experiment, stays so).

### 2.3 Cross-cutting: navigation unification

Today's navigation is **dual**: Compare is URL-routed; Index/Inspect are driven by Zustand `activePage` (`AppShell.tsx`, `TabRocker.tsx`). Every new surface is **URL-routed from the start**, owning its URL via `useParams`/`useNavigate` — *not* by extending the fragile `parseLocation` 3-page union. The dual model runs for Phases 1–4; Phase 1 owns a URL↔`activePage` nav-bridge. Each per-surface cutover retires not only the old *URL routes* but the old surface's `activePage` *branch* (a stale persisted `activePage:"inspect"` at URL `/` otherwise renders an empty `PageBody` — the issue-#77 class). Phase 5 retires the dual model itself.

**Persist-version caution.** Adding persisted Zustand fields needs **no** `persist` version bump (missing keys read `undefined`). A bump *without* a `migrate` callback (the store has none, `state.ts`) **discards every user's persisted state**. Bump only when an existing field's shape changes; Phase 5's removal of `activePage` from `partialize` is the one warranted bump + `migrate`.

### 2.4 Recommended execution order

Six phases, numbered by recommended execution order:

```
Phase 0  Corpus foundation     (corpus listing routes — additive, invisible)
Phase 1  Sample table          (replaces Inspect)   ── depends on 0
Phase 2  Series backend        (new tables/routes/events — invisible) ── depends on 0
Phase 3  Series stage UI       (replaces Compare)   ── depends on 2; soft-depends on 1
Phase 4  Focus workspace       (replaces Index)     ── depends on 1
Phase 5  Final cutover         (retire dual-nav)    ── depends on 1–4
```

Phases 1 and 2 are independent and may run in parallel. Phases 1–4 each end by retiring the surface they replace (per-surface cutover); Phase 2 retires nothing (it ships invisibly). Phase 5 is only the dual-nav retirement and the final sweep. Phase 3 *soft*-depends on Phase 1 (scoping is normally entered from sample-table multi-select, but the folio's `+ New` is an independent entry point).

---

## 3. Phase 0 — Corpus foundation

**Goal:** corpus-level listing routes. Genuinely additive, genuinely invisible.

**Depends on:** nothing.

### 3.1 Corpus listing routes

Three experiment-scoped *listing* routes get corpus siblings; the existing experiment-scoped routes stay:

| New route | Sibling of | File |
|---|---|---|
| `GET /api/samples` (optional `?experiment_id=`) | `GET /api/experiments/{id}/samples` (`routes_samples.jl:4`) | `routes_samples.jl` |
| `GET /api/sample-tags` | `GET /api/experiments/{eid}/sample-tags` (`routes_picker.jl:55`) | `routes_picker.jl` |
| `GET /api/picker-samples` | `GET /api/experiments/{eid}/picker-samples` (`routes_picker.jl:71`) | `routes_picker.jl` |

The `GET /api/samples` projection **includes `q_units` per sample** (from the owning experiment's config) — so the corpus surfaces and the future cross-experiment normalization have it without a later route-body change. No schema change, no event kinds, no migration.

### 3.2 Build, test, exit criteria

- **Tests:** Julia route tests for the three endpoints. The slow suite — capture once, grep.
- **Exit:** `GET /api/samples` returns the whole corpus with `q_units`; existing experiment-scoped routes untouched and green; the running frontend genuinely unaffected.

---

## 4. Phase 1 — Sample table (contact sheet + loupe)

**Goal:** ship the sample table at `/samples`, replacing Inspect. Introduces the new app shell and the URL-routing model.

**Depends on:** Phase 0.

**Mockup:** `sample-table.html` — contact-sheet view + loupe view.

### 4.1 New shell + routing model

- New shell: the topbar from the mockups (corpus wordmark, three stage-tabs, Beamtime facet chip). **One hoisted top-level `<Routes>`**: old routes dispatch into the old `AppShell` body, new routes into the new shell body — two page bodies under one router. Hoist `useGlobalShortcuts` and the theme effect.
- New URL routes: `/samples` (`?beamtime=<id>` as URL query state), `/samples/loupe/:sampleId`. The Samples stage-tab points here.
- Phase 1 owns the URL↔`activePage` **nav-bridge**, which must also **coerce a stale persisted `activePage`** (not only redirect URLs).

### 4.2 Surface

- **Contact sheet:** corpus-wide via a **distinctly-named** corpus query (`useCorpusSamples()` — `useSamples(experimentId)` already exists, `queries.ts:119`). A new corpus query key shape (e.g. `["corpus","samples"]` — `queryKeys.samples` is per-experiment-keyed today). Each row: identity, exposure-thumbnail strip (`ThumbnailGallery`/`DetectorImage size="thumb"`), Kept count, **Tags** column (light free-form chips), Status.
- **Loupe:** `DetectorImage size="full"` + exposure strip + metadata sidecar with a Sample-tags section. Absorbs `DetectorImageCard`/`SampleMetadataCard` internals. Must **not** hard-code a file-per-exposure assumption (the deferred derived-exposure feature).
- **Culling:** reject-only, multi-select batch reject, representative pick — wired to existing `useSetExposureStatus`/`useSelectExposure`/`useAddExposureTag`.
- **Corpus `add_tag` — a full six-layer change, not a one-liner.** The contact sheet is the first reader of a corpus `samples` query, so this lands here. It requires: (a) the new corpus query key; (b) a **tri-scope discriminator** in `resolveMutator`/`resolveMutatorForEvent` (`mutatorRegistry.ts`) — today it discriminates sample-vs-exposure tags on `experimentId !== undefined`; a corpus sample-tag op carries no `experimentId` and would misroute to the exposure mutator; (c) the **`applyRemoteToCache.ts` `add_tag`/`remove_tag` branch made tri-scope too** — it is itself binary today (`samples(experiment_id)` vs `exposures(sample_id)`), so a corpus tag SSE frame would invalidate a garbage key and the contact sheet would never refresh from a foreign tab; (d) the mutator `onMutate`/`onSuccess` (`trivial.ts`) patching the corpus key optimistically. Parameterize the hard-coded `source='manual'` at `routes_samples.jl:71`.

### 4.3 Build, test, exit criteria

- **Tests:** Vitest for table/loupe; boneyard skeletons (gate on `query.isLoading`); a Playwright mocked spec for cull / batch-reject / loupe-flip / tag.
- **Cutover (end of phase):** delete `InspectPage` + Inspect-only components **and the Inspect E2E spec**; redirect `/inspect*` → `/samples` **and** retire the `activePage:"inspect"` branch.
- **Exit:** `/samples` renders the corpus, culling + tagging round-trip through the queue (corpus cache included), the loupe flips frames; Inspect is gone; full suite green.

---

## 5. Phase 2 — Series backend

**Goal:** the entire series data model — new tables, routes, event kinds — shipped invisibly (empty `series` tables, unused routes). No UI, no data migration yet. Profile identical to Phase 0: zero user-facing change, independently verifiable, independently mergeable.

**Depends on:** Phase 0. Independent of Phase 1.

### 5.1 Schema — new tables, never a rename

Add via a `migrate_series!` (`CREATE TABLE IF NOT EXISTS`); the `comparison*` tables stay in `SCHEMA` permanently.

- `series` — the comparison columns (incl. `forked_from_id` as a `series(id)` self-ref, view-choice columns) **plus** recipe columns: `ordering_variable TEXT`, `order_rule TEXT NOT NULL DEFAULT 'manual' CHECK (order_rule IN ('ascending','descending','manual'))`, `state TEXT NOT NULL DEFAULT 'committed' CHECK (state IN ('draft','committed'))`. `content_hash` is **NULL while `state='draft'`**, computed on `series_plate_committed` from the **plate only** (`series_members`) — never the recipe, so `series_recipe_updated` never touches it.
- `series_members` — the plate: the `comparison_members` shape (per-exposure frozen snapshot rows, `json_valid` guards intact), `series_id` FK **`ON DELETE CASCADE`**.
- `series_samples` — the recipe membership: `(id, series_id REFERENCES series ON DELETE CASCADE, sample_id REFERENCES samples ON DELETE CASCADE, position INTEGER NOT NULL, pinned INTEGER NOT NULL DEFAULT 0 CHECK (pinned IN (0,1)), excluded INTEGER NOT NULL DEFAULT 0 CHECK (excluded IN (0,1)), UNIQUE(series_id, position))`. `sample_id` is **`CASCADE`**, not `SET NULL`: a `series_samples` row is a pure pointer with no snapshot, so an orphan (NULL `sample_id`) is unrenderable. (Sample deletion has no route today; if one is added it must emit `series_recipe_updated` for affected series — §11.)
- `series_messages`, `series_pins` — the `comparison_messages`/`comparison_pins` shapes, `series_id` FK **`ON DELETE CASCADE`**. With all four child tables cascading, the `series_deleted` dispatcher branch stays a one-line `DELETE FROM series`.
- A minimal `schema_migrations(name TEXT PRIMARY KEY, applied_at TEXT)` sentinel table — there is no migration-version table today, and Phase 3's copy needs a real sentinel (§6.1).

### 5.2 Routes, event kinds, gate

New `routes_series.jl` mirroring `routes_comparisons.jl` (corpus-wide): `GET/POST /api/series`, `GET /api/series/{id}`, `GET /api/series/{id}/forks`, `PATCH /api/series/{id}` (recipe edit), `POST /api/series/{id}/commit`, `DELETE /api/series/{id}`, messages + pins. New `series.jl` for business logic adapted from `comparisons.jl` — **fix** the `last_event_at` mixed-timestamp sort bug (`comparisons.jl:669`, issue #76) rather than copy it. Drop the `is_author` gate (architecture decision 3; lives only at `routes_comparisons.jl:255,312`).

`POST /api/samples/tags/batch` — one `with_idempotency` tx, N tag inserts — **non-optional** (N single calls have no atomic boundary; a mid-batch reload leaves a half-confirmed recipe). Scoping-step confirmations (Phase 3) reuse the existing `add_tag` kind with `source='scoping'`.

**Event kinds — the six layers per kind:** dispatcher branch · SSE `broadcast_event!` (+ `post_state` envelope where noted) · `applyRemoteToCache` handler · mutator · paired contract test · `rebuild_views_from_log!` round-trip.

| Kind | Payload | `post_state`? | Notes |
|---|---|---|---|
| `series_created` | **full `series_samples` snapshot** + recipe fields; **zero `series_members`** | no | A draft. New lifecycle — `comparison_created` rejects empty members. |
| `series_recipe_updated` | **full `series_samples` snapshot** (never a delta) | no | Mutable view table. |
| `series_plate_committed` | **full member list** (the plate) — members carry **no ids** | **yes** — `fetch_series_with_plate` projection | The old `submit`. First comparison-shaped event to carry `post_state` — the `applyRemoteToCache` handler is genuinely new code. |
| `series_deleted` | — | no | One-line `DELETE FROM series` (four-table cascade). |
| `series_pinned` / `series_unpinned` | — | no | **Five-layer** (no `post_state`, no mutator). Stored with **`entity_type='user'`, `entity_id=user_id`** (the `comparison_pinned` precedent) — therefore *outside* the `entity_type="series"` round-trip fold. |

**The view-producing `series_*` dispatcher branches are pure-replace on the child tables and upsert on the parent — explicitly NOT `comparison_submitted`-shaped.** `_update_view_for_comparison_submitted!` discriminates each member on `m.id` (null → INSERT, integer → UPDATE, integer-not-found → `error()`); folding such an event from an empty view (as `rebuild_views_from_log!` and the copy-forward replay both do) would throw. Instead:

- `series_created` — **upsert** the `series` row (`SELECT` for an existing row at the event's `entity_id`, `UPDATE` if present else `INSERT`): the route mints the row to capture the AUTOINCREMENT id, exactly as `routes_comparisons.jl` does and as `_update_view_for_comparison_created!` upserts — a plain `INSERT` would collide on the live route path. Then `DELETE FROM series_samples WHERE series_id=?` + `INSERT` all from the payload snapshot.
- `series_recipe_updated` — `DELETE FROM series_samples WHERE series_id=?`, then `INSERT` all from the payload snapshot.
- `series_plate_committed` — `DELETE FROM series_members WHERE series_id=?`, then `INSERT` all from the payload member list (members carry no ids — the dispatcher mints them); set `state='committed'`, compute `content_hash` from the plate.

Parent-row upsert + child `DELETE`+`INSERT` makes every fold idempotent and order-independent given a full-snapshot payload — the live route path, the migration's synthesized events (§6.1), and `rebuild_views_from_log!` all replay correctly.

All four non-pin `series_*` events carry **`entity_type='series'`** on the wire. `resolveMutatorForEvent` (defined in `mutatorRegistry.ts`, called from `replayCoordinator.ts`) is a **third** per-kind switch site — beyond `update_view_for_event!` and `applyRemoteToCache` — and must gain a case for each kind; `series_plate_committed` needs a `synthesizeFromSse` so its deferred resolves with the route's response shape (its payload is an id-less member list, not that shape).

The historical `comparison_*` dispatcher branches (`events.jl:386-421`) stay **frozen, bit-faithful, targeting `comparison*` tables** — never aliased.

### 5.3 Build, test, exit criteria

- **Tests:** Julia route + dispatcher + idempotency tests; the six-layer paired contract tests per kind; a **`rebuild_views_from_log!` round-trip on series created natively via `POST /api/series`** — the test **empties the target series' view rows** (`series`/`series_members`/`series_samples`) before re-folding, then asserts the rebuilt state matches (proves the pure-replace branches fold correctly from empty). Vitest for the new mutators + cache handlers.
- **Exit:** `/api/series*` round-trips through the queue + SSE; the native-series round-trip test is green; `series` tables ship empty; the running frontend is unaffected.

---

## 6. Phase 3 — Series stage UI (folio + scoping + builder)

**Goal:** ship the Series stage at `/series`, replacing Compare. Owns the comparison→series **data migration** and the Compare retirement.

**Depends on:** Phase 2; soft-depends on Phase 1 (§2.4).

**Mockups:** `series-folio.html`, `series-scoping.html`, `series-builder.html`.

### 6.1 The comparison→series data migration (ships with this phase)

`migrate_comparisons_to_series!` is added to `migrate_schema!` **in this phase's PR** (so the copy captures the final comparison state — no dual-write window). It must be placed **after `migrate_compare_view_choices!` and `migrate_compare_relax_nullability!`** in the `migrate_schema!` sequence — it reads the view-choice columns and the relaxed-nullability columns those late migrations add. It runs at `open_db`, single-threaded, **before `serve` accepts connections**, wrapped in its own `SQLite.transaction`.

**The mechanism (resolves former open decision 3).** A direct table-seed would leave migrated `series` rows non-event-sourced — `rebuild_views_from_log!` rebuilds per `(entity_type, entity_id)`, so such a row would replay to *empty*. The copy is therefore event-sourced, and the log and the view are produced from **one source of truth**:

1. **Gate:** at the start of the transaction, check the `schema_migrations` sentinel (§5.1); if the marker row is present, skip. (Not a `series`-populated check — a live `POST /api/series` row would defeat that.)
2. **For each `comparisons` row, construct the synthetic event payloads first:** a `series_created` payload (full `series_samples` snapshot derived from the comparison's members' samples; recipe fields; zero `series_members`) and a `series_plate_committed` payload (the full member list, members carrying no ids).
3. **Raw-`INSERT` the two `user_actions` rows** — never via `apply_event!` (its public method broadcasts; a migration must not fan out N replay-as-reruns to every tab). `entity_type='series'`, `client_op_id = NULL`, `client_id = NULL`, carrying the comparison's original `created_by`/`created_at`; fresh autoincrement ids appended after all historical `comparison_*` events. Because `client_op_id` is NULL, these rows write **no `idempotent_responses` cache rows** — intentional; a migration is not a client op.
4. **Produce the view rows by folding those payloads through the same `update_view_for_event!` `series_*` dispatcher branches** that live replay uses (§5.2, pure-replace) — *not* by an independent direct write. Log and view agree by construction.
5. After the loop, write the `schema_migrations` marker row **last, inside the same transaction** — the gate flips only on a fully-committed copy.
6. Copy `comparison_messages`→`series_messages`, `comparison_pins`→`series_pins` the same way.

After Phase 3, `/api/comparisons` routes, `comparisons.jl`, `routes_comparisons.jl` are deleted (the Compare UI is gone); the `comparison*` tables and `comparison_*` dispatcher branches **stay forever**.

### 6.2 UI

- `/series` — **folio:** corpus masonry, `useSeriesList()`. `/series/new` — **scoping step:** machine-proposed ordering variable, confirm-and-build gate; confirmations write `sample_tags` via the batch route with `source='scoping'`, storing the ordering variable as a **structured `(key,value)`** (key = the ordering-variable name, value = the parsed value) so the deferred tag-criteria predicate (decision 2) can query it without a later data migration. `/series/:id` — **builder:** the consolidated `series-builder.html`.
- **Type migration:** `ComparisonMember`→`SeriesMember` across `MultiTracePlot`, `MemberTraceLayer`, `multiTraceAdapter`, `multiTraceExportMarks.ts`, `MemberMetaGutter`, `MemberMetaRow`, `lib/comparison/{coloring,labels,draftFactories}`. Model optional fields as `T | null` (mirror `ComparisonMember` in `api.ts`) for `exactOptionalPropertyTypes`.
- **Optimistic vs spinner:** recipe edits (`series_recipe_updated`) are **optimistic** — negative `series_samples` placeholder ids via `nextOptimisticId()`. Plate commit (`series_plate_committed`) is **spinner** — no optimistic write, the `saveComparison` precedent.
- **Reuse:** `lib/figure-export/`, `MultiTracePlot` render core, `FigureExportControls`, `ConflictModal`; slug-permalink hooks gain `series` resolution following Compare's URL-ownership pattern. New overlays (scoping modal) use `useFocusTrap`. Externally-updatable numeric inputs reuse `QNumInput`.

### 6.3 Build, test, exit criteria

- **First step:** re-audit the Phase 2 series route response shapes against the builder's needs — treat the §5.2 event-kind table and the route contract as the frozen Phase 2 deliverable; flag drift before building UI on it.
- **Tests:** Vitest for folio/scoping/builder; the six-layer contract tests; a **`rebuild_views_from_log!` round-trip for `entity_type="series"` on a row migrated from a real comparison-era DB** — empties the migrated series' view rows, re-folds, asserts equality (exercises the *migrated* path, not only native series); Playwright mocked specs for folio→scoping→builder; figure-export specs extended.
- **Cutover (end of phase):** delete `Compare.tsx` + Compare-only components, `routes_comparisons.jl`, `comparisons.jl`, **and the Compare E2E specs**; redirect `/compare*` → `/series` **and** retire the `activePage:"compare"` branch. `comparison*` tables + dispatcher branches stay.
- **Exit:** a series scopes → builds → edits → commits a plate → exports; the folio lists series corpus-wide; Compare is gone; the migrated-series round-trip test is green on a comparison-era DB.

---

## 7. Phase 4 — Focus workspace

**Goal:** ship the focus workspace at `/sample/:sampleId`, replacing the three-card Index. Layout + the q-link.

**Depends on:** Phase 1.

**Mockup:** `focus-workspace.html`.

- **Layout:** trace as hero plate; co-resident detector image; phase call as rail; Notes as a margin/drawer. The trace/index components are reused (§2.2) — Phase 4 resolves their **Zustand-direct wiring**: a route-param→Zustand sync shim (`/sample/:sampleId` → `activeSampleId`) or prop-drill into `PlotCard`/`IndicesCard`/`PhasePanel`.
- **The q-link:** hovering a trace peak / detector ring / reflection-row lights all three — an **ephemeral `hoveredQ` Zustand field** + named action, excluded from `partialize` (the `hoveredPeakId` pattern, not a context). The `DetectorImage` ring overlay is **rotation-aware**. The hover channel must not bypass the `QNumInput` focus-gate where those components host numeric inputs.
- New route `/sample/:sampleId`; the Index stage-tab re-points here.
- **Tests:** Vitest for the layout + q-link; **the carried-over trace/index interaction tests must pass unchanged** (the regression floor); boneyard skeletons for the data-driven panels; a Playwright mocked spec for the q-link.
- **Cutover (end of phase):** delete `IndexPage` + the three-card composition **and the Index-workspace E2E specs**; redirect `/index*` → `/sample/...` **and** retire the `activePage:"index"` branch.
- **Exit:** indexing one sample works exactly as before, in the new layout, q-link live; Index is gone.

---

## 8. Phase 5 — Final cutover

**Goal:** retire the transitional scaffolding. Small — the per-surface cutovers happened in Phases 1, 3, 4.

**Depends on:** Phases 1–4.

- Delete the old `AppShell` shell body, `WorkspaceGrid`, `TabRocker`, the Zustand `activePage` model and the nav-bridge.
- `state.ts`: remove `activePage` from `partialize`; bump the `persist` version **with a real `migrate`** that drops the dead key while preserving `username`/`theme`/`tutorialSeen`.
- `lib/queue/persistence.ts`: bump the queue persistence `schema_version` so any pre-cutover queued `comparison_*` op in `sessionStorage` drops cleanly (a toast, not a `mutatorRegistry` throw). This counter is independent of the Zustand `persist` version — bump both, deliberately.
- Final dead-code sweep; `npm run build` (`tsc --noEmit` + vite) green.
- **The `comparison*` tables and their dispatcher branches are never deleted** — event-replay machinery.

---

## 9. Open decisions

1. **Series schema shape — settled.** New `series`/`series_members`/`series_samples` tables; recipe = explicit ordered `series_samples` list + `ordering_variable` + `order_rule` + `state`. A real table (not a JSON blob) because recipe membership is mutated by `series_recipe_updated` and benefits from row-level event folding, and the deferred predicate (decision 2) wants it queryable.
2. **Tag-criteria predicate** — deferred to a post-Phase-3 follow-up. v1 series recipes are explicit sample lists; the folio's "+N new match" pill comes with the predicate. Phase 3 scoping writes `sample_tags` in a structured `(key,value)` shape (§6.2) so decision 2 inherits no tag-data migration. Confirm the staging is acceptable.
3. **Copy-forward event sourcing — RESOLVED (v3/v4).** The migration synthesizes `series_*` events: it constructs full-snapshot payloads, raw-`INSERT`s the `user_actions` rows (`client_op_id` NULL, no broadcast, `migrate_schema!`-time), and folds those payloads through the pure-replace dispatcher branches to produce the view rows — §6.1. A direct table-seed was rejected (per-entity replay would rebuild migrated rows to empty).

---

## 10. Risk register

| Risk | Phase | Mitigation |
|---|---|---|
| Treating replay-machinery tables/branches as renamable — breaks `rebuild_views_from_log!` | 2,3 | `series*` are new tables; `comparison*` tables + branches frozen forever (§2.1). |
| Migrated `series` rows not event-sourced → round-trip test passes vacuously | 3 | The copy synthesizes `series_*` events (§6.1); the Phase 3 round-trip test exercises `entity_type="series"` on a *migrated* row and empties view rows before re-folding. |
| Synthesized event payload and view rows disagree | 3 | The migration constructs payloads, then folds them through the dispatcher to produce the view rows — one source of truth (§6.1 step 4). |
| `comparison_submitted`-shaped dispatcher throws on fold-from-empty | 2,3 | `series_*` view branches are pure-replace (`DELETE by series_id` + `INSERT`), not member-id-discriminating (§5.2). |
| Migration synthesizing events via `apply_event!` fans out N broadcasts | 3 | Raw `INSERT` only, `client_op_id`/`client_id` NULL, runs in `migrate_schema!` before `serve` accepts connections (§6.1). |
| `migrate_comparisons_to_series!` non-idempotent, mis-ordered, or partially applied | 3 | `schema_migrations` sentinel written last in the same transaction; placed after `migrate_compare_view_choices!`/`_relax_nullability!`; round-trip test on a real comparison-era DB gates the merge. |
| New event kind wired in < 6 layers (or 5 for pin kinds) | 2 | Per-kind layer table (§5.2); `series_plate_committed` `post_state`; pin kinds explicitly five-layer, `entity_type='user'`. |
| Persist-version bump without `migrate` silently wipes user state | 1,5 | Don't bump to add fields; Phase 5's `activePage` removal ships a real `migrate` (§2.3). |
| Stale persisted `activePage` after a surface is retired → empty `PageBody` | 1,3,4 | Each per-surface cutover retires the `activePage` branch; Phase 1's nav-bridge coerces stale values. |
| Dual-nav coexistence (URL ↔ `activePage`) desync across Phases 1–4 | 1 | Phase 1 owns and tests the nav-bridge. |
| Corpus `add_tag` misroutes / fails to invalidate the corpus cache | 1 | Tri-scope discriminator in `resolveMutator`/`resolveMutatorForEvent` *and* the `applyRemoteToCache` branch; full six-layer corpus-key work (§4.2). |
| Focus-workspace reshell disturbs validated trace/index interactions | 4 | Carried-over interaction tests must pass *unchanged* — the regression floor. |
| Phase 2/3 contract drift (route shapes shift before the UI consumes them) | 3 | Phase 3's first step re-audits the Phase 2 route shapes against the §5.2 frozen contract (§6.3). |
| `comparison`→`series` test churn | 2,3 | Budget new `series` suites as a Phase 2/3 line item; slow Julia suite — capture once, grep. |
| Scoping parser cold-start: `sample_tags` starts empty | 3 | Accepted by design; state it in the Phase 3 acceptance criteria. |
| Cross-experiment series overlay mismatched `q_units` | 0,3 | `q_units` in the Phase 0 `/api/samples` projection; normalization + caption signal a Phase 3 line item. |
| Loupe hard-codes a file-per-exposure assumption | 1 | Phase 1 detailed plan carries the derived-exposure constraint forward. |
| No feature flag → no instant rollback | all | Accepted: the old route stays reachable until its phase's cutover; the migration is forward-only by design. |

---

## 11. Constraints for the detailed phase plans

**Phase 0:** `GET /api/samples` projection includes `q_units`.

**Phase 1:** distinct corpus-query name (`useCorpusSamples`) + a new corpus query key shape; corpus `add_tag` is full six-layer — corpus key, tri-scope `resolveMutator`/`resolveMutatorForEvent` discriminator, **tri-scope `applyRemoteToCache` `add_tag`/`remove_tag` branch** (binary today), mutator `onMutate`/`onSuccess` corpus patches; the tri-scope split keys on `experiment_id` **payload-field presence** (corpus and experiment-scoped sample tags both arrive as `entity_type='sample'` — the split cannot use a new `entity_type` value), and the binary→tri-scope conversion must touch `resolveMutator` **and** `resolveMutatorForEvent` consistently (the latter's `add_tag`/`remove_tag` case is `entityType`-keyed today); parameterize `routes_samples.jl:71` tag `source`; new surfaces own their URL via `useParams`/`useNavigate`, not the `parseLocation` union; `?beamtime=` as URL query state; one hoisted top-level `<Routes>`; `useGlobalShortcuts`/theme effect hoisted; the nav-bridge coerces stale `activePage`; loupe must not assume file-per-exposure; Tailwind `@theme` for "The Print" palette; boneyard skeletons; delete the Inspect E2E spec at cutover.

**Phase 2:** `series*` created as new tables, `comparison*` retained permanently; `CHECK` constraints on `order_rule`/`state`/`pinned`/`excluded`; `UNIQUE(series_id, position)` on `series_samples`; all four `series_*` child tables `ON DELETE CASCADE` on `series_id`; `series_samples.sample_id` `ON DELETE CASCADE`; the view-producing `series_*` dispatcher branches are **pure-replace** (`DELETE by series_id` + `INSERT from full-snapshot payload`), not `comparison_submitted`-shaped; `series_created` payload self-reconstructs the row from empty; `series_plate_committed` members carry no ids and the branch carries a `post_state` envelope (genuinely new frontend handler code); `series_pinned`/`unpinned` are five-layer, `entity_type='user'`, outside the series round-trip fold; `content_hash` NULL for drafts, folds the plate only; keep all `with_idempotency`-wrapped series routes JSON-only (`idempotent_responses` replay is JSON-only); fix the `last_event_at` sort bug in `series.jl`; `schema_migrations` sentinel table; the native-series round-trip test empties view rows before re-folding; the four non-pin `series_*` events carry `entity_type='series'` on the wire; `resolveMutatorForEvent` gains a case per kind and `series_plate_committed` a `synthesizeFromSse`; the `applyRemoteToCache` `series_recipe_updated` branch is a wholesale snapshot-replace (no id-splice); `series_samples.id` is replay-volatile — never key client state on it (use `(series_id, position)`).

**Phase 3:** `migrate_comparisons_to_series!` placed after `migrate_compare_view_choices!`/`_relax_nullability!`, in its own `SQLite.transaction`, sentinel row written last; the copy constructs payloads then folds them through the dispatcher (§6.1); the migrated-series round-trip test empties view rows and exercises `entity_type="series"`; first step re-audits the Phase 2 route shapes against the §5.2 frozen contract; `ComparisonMember`→`SeriesMember` migration covers `MultiTracePlot`, `MemberTraceLayer`, `multiTraceAdapter`, `multiTraceExportMarks.ts`, `MemberMetaGutter`, `MemberMetaRow`, `lib/comparison/{coloring,labels,draftFactories}` (not `yBands.ts` — pure numeric, no type); optional fields as `T | null`; optimistic recipe edits vs spinner plate-commit; scoping writes structured `(key,value)` `sample_tags`; `q_units` normalization + caption; `useFocusTrap` on the scoping modal; `QNumInput` for the offset dock; the Plot `.scale().invert()` cast; re-audit the `add_tag` cache fan-out for the scoping surface's query; delete the Compare E2E specs at cutover.

**Phase 4:** route-param→Zustand sync vs prop-drill for `PlotCard`/`IndicesCard`/`PhasePanel`; q-link via an ephemeral `hoveredQ` Zustand field + named action (not a context); rotation-aware `DetectorImage` ring overlay; the q-link channel must not bypass the `QNumInput` focus-gate; boneyard skeletons for the focus-workspace panels; carried-over interaction tests as the regression floor; delete the Index-workspace E2E specs at cutover.

**Phase 5:** `activePage` removed from `partialize` with a real `persist` `migrate`; queue persistence `schema_version` bump (independent counter); `comparison*` tables + dispatcher branches kept forever.

**Cross-cutting:** if a sample-delete path is ever added, it must emit `series_recipe_updated` for affected series (the `series_samples.sample_id` cascade is otherwise invisible to event replay).

---

## 12. Reviewer history

- **v1 → v2:** four reviewers (himalaya/SQLite · queue · frontend · architecture) all independently flagged the `comparisons`→`series` *rename* mislabelled additive. v2: new tables + copy, frozen `comparison*` machinery, series backend Phase 0→2, per-surface cutover, "carries over" qualified.
- **v2 → v3:** re-review confirmed v1 findings resolved and Phases 0–1 sound. Three reviewers converged: a direct table-seed makes `rebuild_views_from_log!` (per-entity) replay migrated rows to empty. v3 resolved open decision 3 (synthesize events), split the series backend into its own Phase 2, renumbered.
- **v3 → v4:** re-review found frontend and architecture **clean**; the queue raised one blocker and himalaya four should-fixes, all on the synthesize-events mechanism. v4: the migration constructs payloads then folds them through the dispatcher (§6.1); the `series_*` view branches are pure-replace, not `comparison_submitted`-shaped (§5.2); migration ordering, sentinel-in-transaction, the tri-scope `applyRemoteToCache` branch, `UNIQUE(series_id, position)`, `entity_type='user'` pins, the Phase 3 anti-drift checkpoint, and per-phase E2E-spec deletion are pinned.
- **v4 → v5:** re-review found architecture and frontend **clean**. The backend review caught one blocker — §5.2's `series_created` must **upsert** the parent `series` row (route-minted id; a plain `INSERT` collides), only the child tables pure-replace. The queue review pinned the wire `entity_type='series'` for the non-pin kinds and `resolveMutatorForEvent` as a third per-kind switch site. v5 corrects §5.2 (parent upsert), §2.2 (three switch sites), §5.1 (`content_hash` folds the plate only), and §11 (the `resolveMutatorForEvent` cases, the `series_recipe_updated` wholesale-replace cache branch, the payload-presence tri-scope key, the replay-volatile `series_samples.id` invariant).
- **v5 confirmation round:** all four reviewers returned **clean** — no blockers, no should-fixes. Two one-line factual corrections were applied (the `resolveMutatorForEvent` file citation — `mutatorRegistry.ts`, not `replayCoordinator.ts`; and an explicit note that the binary→tri-scope conversion touches both registry functions). The plan is **review-converged** and ready to build from.

---

## 13. Next step

Phases 0–1 are confirmed sound across all reviewer passes — implementation can begin. Spin up the detailed, file-level **Phase 0** plan (`2026-05-17-himalaya-ui-redesign-phase-0.md`); Phase 1's can follow immediately. Phase 2's detailed plan inherits the §11 Phase-2 constraints; Phase 3's resolves open decision 2.
