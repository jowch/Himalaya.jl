# HimalayaUI redesign — implementation gap report

**Date:** 2026-05-27 · `main` @ 1e0b274 (redesign fully merged)
**Validated against:** the authoritative planning spine — [`2026-05-17-himalaya-ui-redesign.md`](superpowers/plans/2026-05-17-himalaya-ui-redesign.md) (master plan: §11 constraint ledger, per-phase exit criteria, §10 risk register) and [`2026-05-17-himalaya-ui-redesign-issue-map.md`](superpowers/plans/2026-05-17-himalaya-ui-redesign-issue-map.md) (31-issue acceptance criteria) — plus the five canonical mockups in [`redesign-mockups/`](redesign-mockups/).
**Method:** six read-only phase-agents, one per phase, each traversing its issues' `Acceptance` lines × §11 constraints × exit criteria × tagged risk rows × its mockup, verified against shipped source (not the plans' own status claims). Findings are gaps/partials/divergences only; everything not listed verified as met.

---

## 1. Scorecard

**Issue ledger:** 31 issues — **26 met · 5 partial · 0 gap.**

| Phase | Issues | Met | Partial | Gap |
|---|---|---|---|---|
| 0 — Corpus foundation | I0.1–I0.3 | 3 | — | — |
| 1 — Sample table / loupe | I1.1–I1.7 | 5 | 2 (I1.3, I1.4) | — |
| 2 — Series backend | I2.1–I2.7 | 7 | — | — |
| 3 — Series stage UI | I3.1–I3.6 | 5 | 2 (I3.2, I3.5a) | — |
| 4 — Focus workspace | I4.1–I4.4 | 3 | 1 (I4.3) | — |
| 5 — Final cutover | I5.1–I5.3 | 3 | — | — |

**Findings by lens × severity** (13 total):

| Lens | High | Med | Low |
|---|---|---|---|
| Completeness | 1 | 1 | 1 |
| Visual | — | 3 | 2 |
| Workflow | — | 1 | 1 |
| Regression | — | — | 2 |

**Findings by status** — the more decision-relevant cut:

| Status | Count | Findings |
|---|---|---|
| **Known-deferral** (already in the parked batch #207/#208/#209) | 6 | F-1, F-2 (#207) · F-6, F-7 (#208) · F-8, F-10 (#209) |
| **Net-new / unintended** (not covered by the batch) | 5 | F-3, F-4, F-9, F-11, F-12 |
| Acceptable-as-is (noted only) | 1 | F-5 |

The load-bearing backend invariants — the pure-replace series dispatcher, event-sourced `comparison→series` migration, six-layer event wiring, frozen `comparison*` replay machinery, persist `migrate` callback — are **all correct in source**. The most important result of this report: **every Med-or-higher UI gap is already captured by the deferred batch.** The redesign shipped its load-bearing parts faithfully, and the user's instinct to park #207/#208/#209 for a restoration pass correctly bracketed the real holes. The net-new findings (§5 Theme B) are smaller — one latent bug and some polish.

---

## 2. Issue ledger (all 31)

| Issue | Status | Evidence |
|---|---|---|
| I0.1 `GET /api/samples` (+`q_units`) | met | route `routes_samples.jl:28`; per-sample `q_units` `:77-78`; tests `test_routes_samples.jl:316-398` |
| I0.2 `GET /api/sample-tags` | met | `routes_picker.jl:77`; sibling intact `:61`; tests `test_picker_routes.jl:150,195` |
| I0.3 `GET /api/picker-samples` | met | `routes_picker.jl:99`; sibling intact `:93`; tests `test_picker_samples_route.jl:199,226` |
| I1.1 shell + hoisted `<Routes>` + nav-bridge | met | `AppRoutes.tsx:62-135`; nav-bridge shipped #189, retired by I5.1 per §2.3/§8 |
| I1.2 corpus query layer | met | `queries.ts:122,157` key `["corpus","samples"]`; `api.ts:112` |
| I1.3 corpus `add_tag` six-layer | **partial** | all six layers present & queue-tested, but **no UI consumes the mutator** (F-1) |
| I1.4 contact sheet UI | **partial** | rows/thumbs/Kept/Tags/Status present; several mockup affordances absent (F-3) |
| I1.5 loupe UI | met | `LoupeFrame.tsx:22-67` iterates all exposures (no file-per-exposure assumption) |
| I1.6 culling wiring | met | `ContactSheetRow.tsx:145-182`; E2E `corpus-culling.spec.ts:63,82,106` |
| I1.7 Inspect cutover | met | InspectPage+components+E2E deleted; `/inspect/*`→`/samples` `AppRoutes.tsx:110` |
| I2.1 schema + sentinel | met | `db.jl:732-835` — all CHECKs, `UNIQUE(series_id,position)`, CASCADE FKs, sentinel |
| I2.2 routes + `series.jl` | met | `#76` sort bug fixed `series.jl:34-37`; `is_author` gate dropped `routes_series.jl:11` |
| I2.3 created/recipe/deleted | met | dispatcher `events.jl:428-439`; parent UPSERT `:757-822`; pure-replace `:836-862` |
| I2.4 plate_committed | met | `events.jl:905-921`; `post_state` `routes_series.jl:262`; `synthesizeFromSse` `commitSeriesPlate.ts:61` |
| I2.5 pinned/unpinned | met | `events.jl:449-470` `entity_type='user'`; five-layer |
| I2.6 batch-tags route | met | `routes_samples.jl:215-300` single `with_idempotency` tx, atomic |
| I2.7 native round-trip test | met | `test_routes_series.jl:476-712` empties view rows then re-folds every kind |
| I3.1 `migrate_comparisons_to_series!` | met | `db.jl:956-1016` gate→synthesize→raw-INSERT→fold→sentinel-last; ordered after compare migrations `:367-384` |
| I3.2 `ComparisonMember`→`SeriesMember` | **partial** | all 8 modules migrated, `yBands.ts` untouched, build green; but `SeriesMember` declared twice in `api.ts` (F-9) |
| I3.3 folio UI | met | `SeriesFolioPage.tsx`; `useSeriesList()` `queries.ts:628` |
| I3.4 scoping UI | met | `SeriesScopingPage.tsx`; structured `(key,value)` batch write `api.ts:693-707` |
| I3.5a builder surface | **partial** | core+export+offset-dock+collapsible rail shipped; heatmap disabled (F-6), peak-tracking layer absent (F-7) |
| I3.5b builder mutations | met | optimistic recipe `seriesDraftFactories.ts:47-50`; spinner commit; series permalink |
| I3.6 Compare cutover | met | `Compare.tsx`/`routes_comparisons.jl` deleted; `comparisons.jl` kept as helper; `/compare*`→`/series` |
| I4.1 Zustand wiring shim | met | `useSyncActiveSampleFromRoute.ts`; test `FocusWorkspacePage.routing.test.tsx:60` |
| I4.2 focus layout | met | `FocusWorkspaceLayout.tsx` grid `[1fr_348px_250px]`; reused components unchanged |
| I4.3 q-link | **partial** | trace-peak ↔ detector-ring wired + rotation-aware overlay; **reflection-row third surface absent** (F-10) |
| I4.4 Index cutover | met | `IndexPage` deleted; redirects + `IndexSlugRedirect.tsx`; `activePage:"index"` gone |
| I5.1 retire dual-nav | met | `AppShell`/`WorkspaceGrid`/`TabRocker` deleted; no live `activePage` ref |
| I5.2 persist bumps | met | Zustand `version:4` + real `migrate` `state.ts:493,509-522`; queue `SCHEMA_VERSION:3` |
| I5.3 dead-code sweep + build | met | `npm run build` exit 0; `comparison*` tables/branches intact |

---

## 3. Findings by lens

Severity: **High** = workflow-breaking or promised-and-missing · **Med** = degraded but functional · **Low** = cosmetic/polish.
Status: **unintended** = silent gap · **known-deferral** = documented (future-feature-ideas.md / deferred batch #217/#207/#208/#209).

### Completeness

**F-1 · High · I1.3 / §4.2 / §4.3 exit · known-deferral (#207) — Sample tagging has no UI on-ramp.**
The full six-layer corpus `add_tag` machinery exists and is queue-tested (`mutatorRegistry.ts:78-90,188-191`, `applyRemoteToCache.ts:321-347`, `trivial.ts:309-380`, `routes_samples.jl:142,160`), but **nothing invokes it**. The contact-sheet "+ tag" button is `disabled` with title "coming soon" (`ContactSheetRow.tsx:282-291`, asserted disabled in `contact-sheet.test.tsx:163`); the loupe Sample-tags section is comment-marked "intentionally read-only" (`LoupeSidebar.tsx:142-160`). The §4.3 exit criterion "culling **+ tagging** round-trip through the queue" is met for culling, not tagging. **This is #207's scope** (corpus tag/rename UI) — so it's a planned restoration, not a silent gap. Worth flagging as the **highest-value item in the deferred batch**: the hard part (six-layer plumbing) already shipped and is tested, and it restores a capability the contact-sheet header advertises.

**F-2 · Med · §4.3 ("Playwright spec for cull / batch-reject / loupe-flip / tag") · known-deferral (#207) — E2E omits the `tag` scenario.**
`corpus-culling.spec.ts` covers cull/batch-reject/representative/loupe-flip (`:63/:82/:106/:125`) but has no tag test — consistent with F-1 (nothing to drive). #207 includes repointing the tag specs.

**F-9 · Low · I3.2 / §2.2 · unintended — `SeriesMember` declared twice in `api.ts`.**
`api.ts:453` and `api.ts:781` are byte-identical `export interface SeriesMember`. The doc comment at `:443-452` ("twinned until ComparisonMember is retired") implies the second block should be `ComparisonMember` (which still exists separately at `:424`). Compiles via TS declaration-merging, so build is green — but editing one block silently merges with the other. Latent maintenance defect.

### Visual

**F-3 · Med · I1.4 / mockup `sample-table.html` · unintended — Contact-sheet affordances missing.**
Absent vs mockup: the per-sample "screened" checkmark dot (`ContactSheetRow.tsx:197-201` is "identity only (no screened mark; that is #162)" — but #162 shipped without it); the sheet-header progress block ("6/9 samples screened" + bar; mockup lines 470-485 vs `SamplesPage.tsx:60-68` kicker only); the floating bottom-center multi-select cull bar (mockup `.cullbar` 569-575 vs the per-row inline action bar `ContactSheetRow.tsx:226-249`). Palette match needs a live screenshot.

**F-4 · Low · I1.4 / mockup `.thumb` · unintended/needs-screenshot — Thumbnail aspect ratio.**
Contact-sheet thumbnails render 3:4 portrait (`ContactSheetRow.tsx:51` `aspect-[3/4]`); mockup frames are square (62×62). Loupe big-frame is correctly square. May be intentional given real detector aspect — flag for the live pass.

**F-6 · Med · I3.5a / §6.3 "all three representations" / mockup `series-builder.html:448-451` · known-deferral (#208) — Heatmap representation is a disabled stub.**
`RepresentationToggle.tsx:35-40` ships the Heatmap button `disabled`, title "coming soon (#208)"; `SeriesBuilderPage.tsx:194` "heatmap deferred to #208". The §6.3 exit "renders a series in all three representations" is therefore not satisfiable today (waterfall + full-bleed ship; heatmap does not). Documented deferral, but note it leaves an exit criterion formally unmet.

**F-7 · Med · I3.5a / mockup `series-builder.html:471,728,845,939` · known-deferral (#208) — Peak-tracking ("Track reflections") layer absent.**
The builder mockup's "Track reflections" toggle and the lines linking each reflection across traces as it migrates with dose are not shipped; only the global "Peak ticks"/"Peak labels" toggles exist (`AnnotationToggles.tsx:75-82`; no reflection-connecting layer in `MultiTracePlot.tsx`/`MemberTraceLayer.tsx`). **#208 explicitly bundles the cross-trace peak-tracking render-core** with the heatmap (it's the same I3.2 render-core follow-up). The validating agent flagged this as undocumented because it couldn't see #208's body — it is in fact planned. *(See §5 — the reflection theme.)*

**F-8 · Med · I4.3 / mockup `focus-workspace.html:487-512,853-888` · known-deferral (#209) — Reflections table panel missing from the focus workspace.**
The mockup's detector row is two-panel (detector image | reflections table `#refl-body`); `FocusWorkspaceLayout.tsx` renders only `FocusDetectorPanel` (`DetectorImage` + `DetectorRingOverlay`). The right-hand reflections table is structurally absent. **This is #209's scope** (reflections-table surface). *(See §5.)*

### Workflow

**F-10 · Med · I4.3 / §7 / §11 Phase-4 "lights all three" · known-deferral (#209) — q-link lights only 2 of 3 surfaces.**
`hoveredQ` wires trace-peak ↔ detector-ring (`state.ts:80,331`, `TraceViewer.tsx:438`, `PlotCard.tsx:66-67,286`; rotation-aware overlay + focus-gate both correct). But there is **no reflection-row surface** to be the third target (grep clean across `IndicesCard`/`PhasePanel`/`MillerPlot`). The mockup's hover comment (line 649) is literally "peak ↔ ring ↔ reflection row". Both `qLinkHover.test.tsx` and `qlink.spec.ts` assert only ring↔trace — so the §7 "lights all three" acceptance is unmet as shipped. **#209 explicitly bundles the q-link row wiring with the reflections-table** (F-8), so completing #209 closes both this and the dangling §7 criterion in one stroke.

**F-5 · Low · §6.1/§6.3 "migrated from a real comparison-era DB" · known-deferral — Migration round-trip uses an in-test fixture, not an on-disk legacy DB.**
`test_migrate_comparisons_to_series.jl:31-220` builds the comparison via the real schema + real tables, clears the sentinel, re-migrates, empties series view rows, and re-folds for `entity_type="series"` asserting `content_hash` equality (`:159-216`). This exercises the real migrated path and satisfies the invariant's intent; the only divergence from a strict reading is that no committed legacy `.db` artifact is checked in. Acceptable; noted for completeness.

### Regression

**F-11 · Low · §4 regression / focus topbar · unintended — In-page sample stepper + Notes drawer not shipped.**
Mockup focus topbar (lines 451-458) has a prev/next sample stepper ("sample 4 of 9") and a Notes toggle with unread badge that becomes a drawer below the breakpoint. `FocusWorkspacePage.tsx` renders neither; `FocusNotesMargin` is `hidden xl:block` with **no drawer fallback** below xl (mockup line 431 "the margin yields to the topbar Notes toggle" not honored — Notes is unreachable below xl). Mitigated app-wide by the `,`/`.` global keyboard shortcut (`useGlobalShortcuts`), so stepping isn't fully lost.

**F-12 · Low · §4 regression / ChatCard · known-deferral — ChatCard dropped, now dead code.**
`ChatCard.tsx` exists but has zero live mount sites (grep finds only the file + AGENTS.md). This matches the redesign rationale (ChatCard was an explicit de-emphasis target — 1.5% of activity for 22% of screen) and the mockups have no chat panel, so the *drop* is intentional. Flagged only because chat @-mention peak-highlighting is now unreachable dead code with no entry point.

---

## 4. §11 constraints & risk register — spot results

All §11 per-phase constraints and §10 risk-register mitigations were verified present **except** where a finding above notes otherwise. Highlights confirmed in source (not exhaustive — see agent detail):

- **Replay safety (§2.1, risk rows):** `comparison*` tables permanent (`db.jl:474,490,512,548`); `comparison_*` dispatcher branches frozen (`events.jl:386-421`); `series_*` branches pure-replace, parent UPSERT — never `comparison_submitted`-shaped. ✓
- **`content_hash` NULL while draft, plate-only on commit** (`series.jl:268-300`, excludes `series_samples`). ✓
- **`entity_type='series'`** on the four non-pin wire frames; **`='user'`** on pins, outside the series fold. ✓
- **`series_samples.id` replay-volatile** — client keys on `(series_id, position)`; cache branch wholesale-replaces, no id-splice (`applyRemoteToCache.ts:280-287`). ✓
- **Persist bump ships a real `migrate`** preserving `username`/`theme`/`tutorialSeen` + wipe-guard (`state.ts:509-522`, tested). Two counters bumped independently. ✓
- **q-link `hoveredQ` ephemeral**, excluded from `partialize`; rotation-aware ring overlay; does not bypass `QNumInput` focus-gate. ✓
- **Loupe makes no file-per-exposure assumption** (`LoupeFrame.tsx:22-67`). ✓

One documented plan-deviation, verified sound and not counted as a gap: `resolveMutatorForEvent`'s `add_tag`/`remove_tag` arm stays 2-arm rather than the literal tri-scope §11 wording, because corpus and experiment sample-tag SSE frames are byte-identical on the wire (a third arm would be unreachable). Backed by a design doc at the call site (`mutatorRegistry.ts:174-191`).

---

## 5. Cross-cutting synthesis

**Theme A — the deferred batch is well-targeted; it captures every Med+ UI gap.** Reconciling the findings against the *actual scope* of the parked batch (not just its issue numbers):

| Finding | Deferred issue | Coverage |
|---|---|---|
| F-1 / F-2 — sample-tag UI on-ramp + tag E2E | **#207** corpus tag/rename UI + repoint specs | full |
| F-6 — heatmap representation | **#208** heatmap + cross-trace peak-tracking | full |
| F-7 — "Track reflections" peak-tracking layer | **#208** (same render-core follow-up) | full |
| F-8 — reflections-table panel (focus) | **#209** reflections-table + q-link row | full |
| F-10 — q-link 2-of-3 (reflection row) | **#209** (same issue) | full |

So the "reflection surfaces" cluster the agents surfaced (F-7/F-8/F-10) is **already a planned two-issue arc** — #208 for the builder's cross-trace tracking, #209 for the focus reflections-table + the third q-link target. The agents flagged them as undocumented only because they were given issue *numbers* without bodies. **One note for the restoration pass:** #209 should also amend or close the §7 "lights all three" q-link acceptance, which is left formally dangling until the reflection row exists.

**Theme B — the genuinely net-new findings (not in the batch) are small.** Five items the parked issues don't cover:
- **F-9 (Low, latent bug)** — duplicate `SeriesMember` declaration in `api.ts:453/:781`; the second was meant to be `ComparisonMember`. Compiles via declaration-merging today but is a footgun. The only real *defect* in the report.
- **F-3 (Med, visual)** — contact-sheet affordances absent from the mockup: per-sample screened-mark, sheet-header progress block, floating cull bar. Not obviously inside #207's "tag/rename" scope — confirm whether it rides along or needs its own issue.
- **F-11 (Low)** — focus-workspace in-page sample stepper + Notes drawer below the `xl` breakpoint (Notes is unreachable on narrow viewports; the `,`/`.` global shortcut mitigates stepping).
- **F-4 (Low)** — contact-sheet thumbnail aspect ratio (3:4 vs square) — defer to the live-screenshot pass.
- **F-12 (Low)** — `ChatCard` is now unmounted dead code (intentional drop); decide delete vs document.

**Incremental-delivery invariant (§2.1) — clean.** Every per-surface cutover landed: Inspect (I1.7), Compare (I3.6), Index (I4.4) deleted with redirects in place; the dual-nav model and `activePage` fully retired (I5.1) with a guarded persist `migrate` (I5.2); `npm run build` green (I5.3). No surface orphaned, no replay machinery deleted.

---

## 6. Recommended triage order

The Med+ UI gaps are already owned by #207/#208/#209 — they belong to the restoration pass the user has parked, not to a new bug-fix cycle. What this report adds on top:

1. **F-9** (Low, but a real defect) — fix the duplicate `SeriesMember` in `api.ts` (rename `:781` → `ComparisonMember`). Trivial, independent of the restoration pass, prevents a silent declaration-merge footgun. The one item worth fixing now.
2. **Confirm batch scope covers F-3** — does #207 include the contact-sheet screened-mark / progress header / cull bar, or do those need their own issue? A scoping decision, not work.
3. **When #209 is picked up** — close or amend the §7 "lights all three" q-link acceptance alongside it.
4. **F-11 / F-4** — fold into the planned **live-screenshot pass** (the visual items this code-level report flagged for pixel confirmation).
5. **F-12** — decide `ChatCard`'s fate (delete dead code, or record its parking next to #217's switch-user note).
