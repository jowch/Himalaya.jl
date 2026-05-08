# Compare picker rethink — design spec

**Status:** Draft (2026-05-08, rev. after frontend-reviewer + himalaya-reviewer pass)
**Bundle:** GitHub issues #84, #85, #87 (§6–7 only)
**PR cut:** Two PRs, ordered #84 → polish

## Problem

The Compare page's "Add traces" picker (`ComparisonPicker.tsx`) is exposure-centric: every exposure across every sample is its own row. Users reason about samples, not raw integration files, so this gets the granularity wrong. The visual hierarchy mirrors the data wrong (filename emphasized over sample name), the "Recently used" section duplicates the top of "All exposures", and the modal eats the plot the user was just editing. The Compare edit mode wastes the right slot of the WorkspaceGrid on a hint card while the most-important CTA hides in the bottom-right corner of the plot area.

## Design decisions (locked)

1. **Drift policy: freeze at pick time.** When a user adds a sample to a comparison via sample-default (no explicit override), the picker resolves the sample's `selected = 1` exposure once at pick time and stores that resolved `exposure_id` on the comparison member. If Inspect later flips the selected flag for the sample, the comparison does *not* silently retarget. Schema unchanged.
2. **PR cut: two PRs, #84 → polish.** PR1 lands the sample-first data model and the visual cleanups from #85 that fall out for free. PR2 extracts the inline-panel variant for edit mode (#87 §6–7) by separating the picker body from its shell.
3. **Out of scope.** #87 §1–5 and §8–10 (action-button hierarchy, member-row label fixes, drag-handle affordances, draft-mode signal, undo/redo, editing-specific plot affordances) are deferred to a follow-on Compare-edit redesign track.

## PR1 — sample-first picker

### Backend

**New route:** `GET /api/experiments/:eid/picker-samples`

Returns one row per sample in the experiment, enriched with the resolved indexing exposure plus the full exposure list (so the override caret has its data without a second round-trip).

```json
[
  {
    "sample": { "id": …, "name": …, "label": …, "notes": …, "tags": [...] },
    "indexing_exposure_id": 4711,
    "all_exposures": [
      { "id": 4711, "filename": "JC068P_E257_S1418_tot", "selected": true },
      { "id": 4712, "filename": "JC068P_E258_S1418_tot", "selected": false }
    ]
  }
]
```

`indexing_exposure_id` resolution per spec §"Default exposure per sample":
1. The exposure with `exposures.selected = 1` for the sample. SQL: `… WHERE selected = 1 ORDER BY id DESC LIMIT 1`. The `selected` LWW invariant in `routes_exposures.jl:99–117` guarantees at most one in current data, but the explicit `LIMIT 1` + tiebreaker is defensive against legacy rows that pre-date the invariant.
2. Else the highest `exposures.id` for the sample. (`id` is autoincrementing; without a dedicated `analyzed_at` column this is the cleanest proxy for "most-recently-ingested-and-typically-most-recently-analyzed". Avoids a join through `user_actions` for what is a defensive fallback path — most production samples have `selected = 1` set on the Inspect page.)
3. Else `null` (sample has zero exposures — UI renders the row disabled with no checkbox interaction; see frontend §"Selection state model" below).

**Query shape (locked):** *two* `Tables.rowtable(DBInterface.execute(db, …))` calls, not one JOIN'd query. (1) `SELECT * FROM samples WHERE experiment_id = ?` for the sample list, plus its tags via existing `samples_with_tags` helper. (2) `SELECT id, sample_id, filename, selected FROM exposures WHERE sample_id IN (?, …) ORDER BY sample_id ASC, id ASC` for the exposures, then group by `sample_id` in Julia. Avoids a Cartesian flatten that would require deduping in-Julia. Per CLAUDE.md "SQLite.jl" gotchas: materialize raw rows via `Tables.rowtable` before fields are read.

**Serialization (locked):** when projecting each exposure into `all_exposures`, pass `bool_keys = (:selected,)` to `row_to_json` (`json.jl:13`) so the boolean serializes as `true`/`false`, not `1`/`0`. Mirrors the existing `routes_exposures.jl:21` pattern. `indexing_exposure_id` serializes as JSON `null` (not an absent key) when no exposure resolves — frontend TS strict mode reads it as `number | null`.

Read-only, no idempotency, no SSE, no `user_actions` row. Lives in `routes_picker.jl` next to the existing `recently-picked-exposures` and `sample-tags` routes. Helper `picker_samples(db, experiment_id)` lives in `comparisons.jl` next to `recently_used_exposures`.

**Backend tests** (`test/test_picker_samples_route.jl`):

- Sample with one exposure marked `selected = 1` → that exposure is `indexing_exposure_id`.
- Sample with multiple exposures, none marked selected → highest-id exposure.
- Sample with one exposure → that exposure regardless of selected flag.
- Sample with zero exposures → row included, `indexing_exposure_id = null`. (Decision: include rather than filter — UI can show a disabled row instead of silently dropping samples.)
- Unknown experiment id → empty list, HTTP 200 (matches existing `sample-tags` semantics).
- `all_exposures` ordered by `id ASC` (assert ordering, not just set-equality — SQLite returns insertion order without explicit `ORDER BY`).
- **Multi-experiment isolation:** sample in experiment A whose exposure ids are smaller than any in experiment B — confirms the helper filters on `samples.experiment_id`, not on a global `MAX(id)`. Two-experiment fixture, assert no cross-leakage.
- **`exposures.sample_id IS NULL` orphans:** seed an orphan exposure row, confirm it appears in no sample's `all_exposures`. (Schema has no `NOT NULL` on `sample_id`.)
- **Sample with `name IS NULL` and/or `label IS NULL`:** confirm the JSON output is `null` (not absent key, not `missing`). The `routes_users.jl` NULL-fill pattern uses `ismissing` to normalize — assert the picker route does the same.
- **Defensive multi-`selected` legacy:** seed two rows with `selected = 1` for the same sample (simulates pre-LWW data). Confirm `LIMIT 1 ORDER BY id DESC` resolves deterministically to the higher id.
- **`selected` JSON-shape regression:** assert exact JSON bytes for a sample row to lock in `"selected": true`/`false`, not `1`/`0`. This is the canonical contract-test layer per CLAUDE.md "Multi-layer contract testing."
- **`indexing_exposure_id`: null JSON, not absent key** — assert via `JSON3.read` that the key is present with value `nothing`/`null`.

### Frontend — components

Decomposition:

- **`ComparisonPickerBody.tsx`** (new) — extracted: filters (search, experiment chips, tag chips), Recents section, main list, Add-N-selected footer. Stateless on container choice; both modal and (PR2) inline shells consume it.
- **`ComparisonPicker.tsx`** (rewritten) — thin modal shell only: overlay, dialog role, focus trap, Esc-to-close, restores focus to trigger. Wraps `ComparisonPickerBody`.
- **`SamplePickerRow.tsx`** (new) — primary row. Sample name (font-medium, primary slot), sample label as secondary, sample notes as tertiary line-clamp, exposure-count badge (e.g. "3 exposures"), checkbox, override disclosure caret. The caret expands an inline `<ul>` of override rows. When `indexing_exposure_id === null` (zero-exposure sample), the row renders disabled: no checkbox, no caret, `data-disabled` set, no `data-exposure-id`.
- **`ExposureListRow.tsx`** (extended, not duplicated) — gains a `control: "checkbox" | "radio"` prop so the override leaf reuses the same component. Filename labelling, locked-state styling, and notes truncation stay shared. The picker's primary row is `SamplePickerRow` (new); `ExposureListRow` becomes the override leaf when caret is expanded. No new `ExposureOverrideRow` component — avoids the two-near-identical-rows footgun the reviewer flagged.

### Frontend — selection state model

`ComparisonPickerBody` is a pure render layer — it does not own selection state. It accepts a callback prop:

```ts
export type Pick = {
  sample_id: number;
  exposure_id: number;          // resolved at pick time (frozen)
  source: "default" | "override";
};
interface BodyProps {
  experimentId: number | undefined;
  picks: Pick[];                // controlled
  onPicksChange: (next: Pick[]) => void;  // batch shell
  onPick?: (pick: Pick) => void;          // immediate shell — fires per pick
  alreadyAddedExposureIds: Set<number>;
}
```

- The shell decides what to do with picks: modal accumulates into `picks` and flushes via "Add N selected"; inline panel ignores `picks` and uses `onPick` to call `addMember` immediately.
- A row whose `indexing_exposure_id === null` is rendered disabled — its checkbox is absent (not just disabled) so a `Pick` with `exposure_id: null` cannot be constructed. This makes line 112's "every pick resolves to an exposure id" claim a true invariant rather than a comment.
- `Pick.source` is required (not optional) so TS-strict's `exactOptionalPropertyTypes` rule doesn't bite. It's debug-only metadata, not persisted, but its presence on the type is guaranteed.

Replaces the earlier `commitMode: "batch" | "immediate"` proposal — the reviewer flagged that as a load-bearing branch inside the body. Pushing the commit decision to the shell collapses the body to "render rows, emit picks."

### Frontend — recents

Reuse existing `GET /api/users/:id/recently-picked-exposures` (returns exposure ids). Client maps each id → its sample id (via the picker-samples response), dedupes to one row per sample (the most recent pick wins), and renders as `SamplePickerRow` with the picked exposure pre-selected — as the default if the picked exposure is still the sample's `indexing_exposure_id`, else as an explicit override (caret expanded, radio set). Recents section dedupes against the main list by sample id (the same sample never appears twice in the visible list).

The dedup is a `useMemo([recentsQ.data, pickerSamplesQ.data])` inside `ComparisonPickerBody` — recents is server state and must not leak into Zustand. (Per CLAUDE.md state-split: TanStack Query owns server state, Zustand owns client state.)

**Future-bug seed (documented for follow-on, not addressed in PR1):** `recently_used_exposures` reads `comparison_members.exposure_id`. Freeze-at-pick + sample-default means two distinct user intents collapse to the same row: "user explicitly chose this exposure" vs "user picked sample, system resolved it." A future "recently picked samples" route would either require a new column on `comparison_members` (e.g. `picked_via TEXT CHECK IN ('default','override')`) or a parallel derivation through `user_actions.payload`. No `user_actions` invariant requires the source attribution today, so PR1 doesn't add one. Out of scope.

If recents granularity ends up noisy (e.g. lots of cross-sample picks dropping rows from view), follow-on work can add a `GET /api/users/:id/recently-picked-samples` route. Out of scope for PR1.

### Frontend — visual changes (#85 fall-out)

- Primary text is sample name (was filename).
- Sample notes appear as a line-clamp-2 secondary line.
- Filename only appears under the override caret, where it's actually meaningful.
- Section dividers: `<hr>`-style border, not just spacing. Headings as `text-xs text-fg-muted` (legibility upgrade from current `text-fg-dim`).
- Selection feedback beyond "N selected": picked rows get a subtle accent ring. Tailwind v4 slash-opacity (`ring-accent/30`) decomposes against `--color-accent` only if the token is in `oklch`/`rgb` form in `styles.css` `@theme`. Implementer step: verify the existing accent token shape works with the slash syntax; if not, define an explicit `--color-accent-30` (or use `bg-accent/15` which is already in use elsewhere in this file at `text-accent` chip styling). A full chip strip is deferred — it duplicates information the row's checked state already conveys.

### Frontend — tests (Vitest)

- **`test/SamplePickerRow.test.tsx`** (new) — sample-first rendering, caret toggle behavior, override radio behavior, exposure-count badge, **disabled-row branch** (zero-exposure samples render no checkbox + `data-disabled`).
- **`test/ComparisonPickerBody.test.tsx`** (new) — filter chip behavior, search filter, recents section dedup against main list, controlled-`picks` wiring (verifies `onPicksChange` is called with the right `Pick` shape), `onPick` immediate path (calls fired in pick order with resolved exposure id).
- **`test/ComparisonPicker.test.tsx`** (updated) — narrowed to modal-shell concerns: open/close, Esc, focus trap, focus-restore. Body-level assertions move to `ComparisonPickerBody.test.tsx`. **Add focus-trap regression test** — the modal shell uses `useFocusTrap`; pulling shared filter chips out of it shouldn't allow Tab to escape the dialog.
- **`test/ExposureListRow.test.tsx`** (updated) — new `control: "checkbox" | "radio"` prop branch.
- **Skeleton-gating test** in `ComparisonPickerBody.test.tsx` — confirms the picker doesn't flicker on background refetch by gating on `isLoading`, not `isPending` (per CLAUDE.md "Skeleton loading via boneyard-js"). The four queries (`useSamples`, the new `usePickerSamples`, `useRecentlyPickedExposures`, `useSampleTags`) all need this discipline.

### Frontend — E2E (Playwright, mocked)

- **`e2e/comparisonPicker.spec.ts`** (updated) — mocks for the new picker-samples route. Covers: open picker → see sample rows (not exposure rows), open caret → see exposures, override selection → resolved exposure id flows to addMember, recents section dedups against main list.
- **Headless verification I run myself** (per user request):
  - `cd packages/HimalayaUI/frontend && npm run e2e -- --grep "ComparisonPicker"` from the worktree, capture output to `/tmp/picker-e2e.out`, verify zero failures before opening the PR.
  - `npm test -- --run ComparisonPicker SamplePickerRow ComparisonPickerBody` for the Vitest slice.
  - `npm run build` for the TS-strict typecheck gate.
  - Backend slice: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI", test_args=["test_picker_samples_route.jl", "test_picker_routes.jl"])'` — captured to `/tmp/jl-picker.out` per the long-run-once-grep-many CLAUDE.md rule.

### Migration / rollout

- Existing `data-testid="picker-row"` on rows stays. New `data-sample-id` added. `data-exposure-id` is set to the resolved exposure id on every *selectable* row (the default unless the override is active). Disabled rows (zero-exposure samples) get `data-disabled` and **omit** `data-exposure-id`. Existing E2E selectors that target `data-exposure-id` continue to work for selectable rows; tests that need to count rows must filter on `:not([data-disabled])` if zero-exposure samples are in the fixture.
- No schema migration. No event-log change. No SSE change. No mutation-queue contact surface.
- `comparison_members.exposure_id` semantics unchanged; the picker is just smarter about which exposure id to send.

## PR2 — inline picker panel in edit mode (#87 §6–7)

### Components

- **`ComparisonPickerPanel.tsx`** (new) — inline shell wrapping the same `ComparisonPickerBody`. No overlay, no focus trap, no Esc-to-close. Different selection-commit semantics (see below).
- **`ComparePageEdit.tsx`** (updated) — right WorkspaceGrid slot swaps from the hint card to `<ComparisonPickerPanel experimentId={eid} searchInputRef={pickerSearchRef} />`. The bottom-right "+ Add traces" button (rendered when `plotMembers.length > 0`) is removed. The empty-state CTA inside the plot host stays but its `onClick` switches to `pickerSearchRef.current?.focus()` — focuses the inline panel's search input rather than opening a modal.
- **No `useImperativeHandle`.** The panel accepts a `searchInputRef: RefObject<HTMLInputElement>` prop and threads it down to the body's search input. Idiomatic for this codebase (avoids introducing a one-off imperative-handle pattern).

### Selection-commit semantics

The modal shell accumulates picks and flushes on "Add N selected" because the modal is dismissive. The inline panel is always-on, so batching is visual noise — each toggle should commit immediately. The body's `picks` / `onPicksChange` / `onPick` prop split (see PR1 §"Selection state model" above, revised per reviewer feedback) makes this a shell decision, not a body branch:

- **Modal shell (`ComparisonPicker`):** owns `picks` state, passes via `picks` + `onPicksChange`. Renders the "Add N selected" footer. On click, iterates `picks` and calls `addMember(p.exposure_id, qc)` per entry.
- **Inline shell (`ComparisonPickerPanel`):** does NOT own `picks` state — passes `picks={[]}` + `onPick={(p) => addMember(p.exposure_id, qc)}`. Each toggle fires `addMember` directly. The picked row visually transitions to the locked state (`alreadyAdded = true`) on the next render via the existing draft-membership lookup (`alreadyAddedExposureIds` prop). No footer rendered.
- **Removal in inline mode** stays with the existing member gutter's "remove" affordance — the panel is add-only. PR2 does not introduce a `removeMemberByExposure` Zustand action.

### Frontend — tests (Vitest)

- **`test/ComparisonPickerPanel.test.tsx`** (new) — inline shell renders without overlay/dialog role, immediate-commit fires `addMember` on first toggle, picked rows transition to locked.
- **`test/ComparePageEdit.test.tsx`** (updated) — right slot now hosts the panel (testid `comparison-picker-panel`); assert the hint card is gone. The bottom-right "+ Add traces" button assertion deletes; empty-state CTA assertions stay.

### Frontend — E2E

- **`e2e/live/comparePickerInline.spec.ts`** (new, live-mode) — real backend roundtrip. Open `/compare/new`, panel is visible in the right slot, click a sample row, see a band appear on the plot, click the override caret, switch exposure, see the band's resolved exposure id update. Per the live-mode runbook (`packages/HimalayaUI/frontend/e2e/live/README.md`): wait ~800ms after `page.goto("/")` before the first mutation that expects an SSE echo. If that wait time changes, update the runbook in lockstep.
- **Headless verification I run myself:**
  - `npm run e2e -- --grep "ComparisonPickerPanel|ComparePageEdit"` for the mocked slice.
  - For the live spec, the user (operator) brings up backend + Vite per the live-mode runbook; I run `npm run e2e:live -- --grep "comparePickerInline"` from there.
  - `npm test` (full Vitest), `npm run build`, and the same Julia slice command as PR1.

### Trade-offs flagged

- **Right slot width.** The 1400px breakpoint floor is `minmax(320px, 22fr)`. Sample-first rows are denser than exposure-first (one row per sample, not per file), so this should fit. The override caret expands inline, so it doesn't need extra width.
- **<1400px stacking.** WorkspaceGrid drops to a single column under 1400px; the inline panel lands below the plot. Acceptable per existing IndexPage IA — the plot stays the dominant pane.
- **Modal entry points.** `ComparisonPicker` (modal) stays mounted for any non-edit-mode entry point that opens it (e.g. future fork-target seeding). PR2 doesn't remove the modal, just unhooks it from the edit page's right slot.

## Component dependency graph

```
ComparisonPickerBody  ◄── ComparisonPicker        (modal shell, batch commit)
                      ◄── ComparisonPickerPanel   (inline shell, immediate commit, PR2)

SamplePickerRow       ◄── ComparisonPickerBody  (primary row)
ExposureListRow       ◄── SamplePickerRow       (caret-expanded leaf, control="radio")
                      (still standalone for any future Inspect-side use)

picker_samples(db, experiment_id)  ◄── routes_picker.jl
                                   ◄── test_picker_samples_route.jl
```

## Acceptance criteria mapping

**Issue #84:**
- [x] Modal title is "Add traces", not "Add exposures". (PR1 trivial copy fix bundled with the rewrite.)
- [x] Primary list is samples; one row per sample.
- [x] Each sample contributes its indexing exposure to the comparison by default.
- [x] Override path exists (caret + radio).
- [x] Recently used dedupes against the main list (by sample id).
- [x] Existing testids migrate cleanly (`picker-row` retained; `data-sample-id` added).

**Issue #85:**
- [x] No technical filename in the row's primary visual slot.
- [x] Recently used does not duplicate items rendered immediately below.
- [x] Selected state has visible representation beyond the counter (accent ring on picked rows).
- [x] Section headings are legible at first glance.

**Issue #87 §6–7:**
- [x] In edit mode, the right slot hosts an inline picker.
- [x] The standalone modal continues to work for non-edit entry points.

## Risks & mitigations

| Risk | Mitigation |
|------|-----------|
| Existing E2E selectors break when picker-row content changes | Retain `picker-row` testid + `data-exposure-id` on the resolved exposure id; add `data-sample-id` for new assertions. |
| Recents section becomes empty for users whose history is across sister samples that no longer match the picker's sample-set | Acceptable — recents is best-effort. Document in component header. |
| `picker_samples` cost on experiments with hundreds of samples | Two queries (samples list, then exposures bulk-grouped via `WHERE sample_id IN (?, …)`) — *not* one JOIN'd query (the JOIN'd shape would Cartesian-flatten samples × exposures and need in-Julia deduping, which the previous spec draft mistakenly recommended). Benchmark on a 200-sample fixture if needed; falls under the existing perf envelope. |
| Inline panel misses the "Esc to dismiss" muscle memory | Inline panel doesn't dismiss — the panel is always present in edit mode. Document in `ComparisonPickerPanel` header. The bottom-right "+ Add traces" button being gone is the discoverability cue. |
| Drift surprise — user expects the picker to retarget when Inspect's selected exposure changes | Drift policy locked to freeze-at-pick. Deliberate; documented in component + spec. Stale-banner work is a separate follow-on if user feedback demands it. |

## Implementation sequence (high-level)

PR1:
1. Backend: `picker_samples(db, experiment_id)` helper + route + test.
2. Frontend: extract `ComparisonPickerBody` from `ComparisonPicker` (no behavior change).
3. Frontend: `SamplePickerRow` + `ExposureOverrideRow` + sample-first row source. Update tests.
4. Frontend: visual cleanups (hr dividers, heading legibility, accent ring).
5. E2E: update mocks for picker-samples; verify spec passes.
6. Headless verification (Vitest + Playwright + Julia slice + `npm run build`).
7. Cut PR.

PR2:
1. Frontend: `ComparisonPickerBody` gains `commitMode: "batch" | "immediate"` prop.
2. Frontend: new `ComparisonPickerPanel`.
3. `ComparePageEdit` swaps right slot; remove buried "+ Add traces" button.
4. Tests: new `ComparisonPickerPanel.test.tsx`, update `ComparePageEdit.test.tsx`.
5. New live E2E `comparePickerInline.spec.ts`.
6. Headless verification (mocked + live + Vitest + `npm run build` + Julia slice).
7. Cut PR.

## Reviewer agents to run before each PR's writing-plans step

- `frontend-reviewer` — picker rewrite is in their wheelhouse (TS-strict, Zustand, TanStack, testid stability).
- `himalaya-reviewer` — backend route + helper for the new endpoint, plus the cross-cutting Plan/SQLite gotchas.

`saxs-physics-reviewer` and `queue-reviewer` are out of scope — no peakfinding/scoring/index changes, no mutation-queue contact surface.
