# Compare picker rethink — design spec

**Status:** Draft (2026-05-08)
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
1. The exposure with `exposures.selected = 1` for the sample, if any.
2. Else the highest `exposures.id` for the sample. (`id` is autoincrementing; without a dedicated `analyzed_at` column this is the cleanest proxy for "most-recently-ingested-and-typically-most-recently-analyzed". Avoids a join through `user_actions` for what is a defensive fallback path — most production samples have `selected = 1` set on the Inspect page.)
3. Else `null` (sample has zero exposures — no row should be selectable).

Read-only, no idempotency, no SSE, no `user_actions` row. Lives in `routes_picker.jl` next to the existing `recently-picked-exposures` and `sample-tags` routes. Helper `picker_samples(db, experiment_id)` lives in `comparisons.jl` next to `recently_used_exposures`.

**Backend tests** (`test/test_picker_samples_route.jl`):

- Sample with one exposure marked `selected = 1` → that exposure is `indexing_exposure_id`.
- Sample with multiple exposures, none marked selected → highest-id exposure.
- Sample with one exposure → that exposure regardless of selected flag.
- Sample with zero exposures → row included, `indexing_exposure_id = null`. (Decision: include rather than filter — UI can show a disabled row instead of silently dropping samples.)
- Unknown experiment id → empty list, HTTP 200 (matches existing `sample-tags` semantics).
- `all_exposures` ordered by `id ASC` so the override list is stable across renders.

### Frontend — components

Decomposition:

- **`ComparisonPickerBody.tsx`** (new) — extracted: filters (search, experiment chips, tag chips), Recents section, main list, Add-N-selected footer. Stateless on container choice; both modal and (PR2) inline shells consume it.
- **`ComparisonPicker.tsx`** (rewritten) — thin modal shell only: overlay, dialog role, focus trap, Esc-to-close, restores focus to trigger. Wraps `ComparisonPickerBody`.
- **`SamplePickerRow.tsx`** (new) — primary row. Sample name (font-medium, primary slot), sample label as secondary, sample notes as tertiary line-clamp, exposure-count badge (e.g. "3 exposures"), checkbox, override disclosure caret. The caret expands an inline `<ul>` of override rows.
- **`ExposureOverrideRow.tsx`** (new) — leaf row with a radio (one-active-per-sample). Filename as primary, "selected" badge on the indexing exposure. Hidden until parent caret is expanded.
- **`ExposureListRow.tsx`** (kept) — no longer the picker's primary row, but retained for any future Inspect-side use. No changes.

### Frontend — selection state model

Local state in `ComparisonPickerBody`:

```ts
type Pick = {
  sample_id: number;
  exposure_id: number;          // resolved at pick time (frozen)
  source: "default" | "override";
};
const [picks, setPicks] = useState<Pick[]>([]);
```

- Toggling a sample's checkbox adds `{sample_id, exposure_id: indexing_exposure_id, source: "default"}`.
- Opening the caret + selecting a different exposure replaces the entry with `{exposure_id: <override>, source: "override"}` while keeping the same `sample_id`.
- "Add N selected" iterates `picks` and calls existing `addMember(exposureId, qc)` for each. The Zustand action is unchanged. `source` is debug-only metadata, not persisted.

### Frontend — recents

Reuse existing `GET /api/users/:id/recently-picked-exposures` (returns exposure ids). Client maps each id → its sample id (via the picker-samples response), dedupes to one row per sample (the most recent pick wins), and renders as `SamplePickerRow` with the picked exposure pre-selected as the override (so re-adding from recents lands the same exposure the user picked last time). Recents section dedupes against the main list by sample id (the same sample never appears twice in the visible list).

If recents granularity ends up noisy (e.g. lots of cross-sample picks dropping rows from view), follow-on work can add a `GET /api/users/:id/recently-picked-samples` route. Out of scope for PR1.

### Frontend — visual changes (#85 fall-out)

- Primary text is sample name (was filename).
- Sample notes appear as a line-clamp-2 secondary line.
- Filename only appears under the override caret, where it's actually meaningful.
- Section dividers: `<hr>`-style border, not just spacing. Headings as `text-xs text-fg-muted` (legibility upgrade from current `text-fg-dim`).
- Selection feedback beyond "N selected": picked rows get a subtle accent ring (`ring-1 ring-accent/30`). A full chip strip is deferred — it duplicates information the row's checked state already conveys.

### Frontend — tests (Vitest)

- **`test/SamplePickerRow.test.tsx`** (new) — sample-first rendering, caret toggle behavior, override radio behavior, exposure-count badge.
- **`test/ComparisonPickerBody.test.tsx`** (new) — filter chip behavior, search filter, recents section, add-N-selected wiring (verifies `addMember` called once per pick with the resolved exposure id).
- **`test/ComparisonPicker.test.tsx`** (updated) — narrowed to modal-shell concerns: open/close, Esc, focus trap, focus-restore. Body-level assertions move to `ComparisonPickerBody.test.tsx`.

### Frontend — E2E (Playwright, mocked)

- **`e2e/comparisonPicker.spec.ts`** (updated) — mocks for the new picker-samples route. Covers: open picker → see sample rows (not exposure rows), open caret → see exposures, override selection → resolved exposure id flows to addMember, recents section dedups against main list.
- **Headless verification I run myself** (per user request):
  - `cd packages/HimalayaUI/frontend && npm run e2e -- --grep "ComparisonPicker"` from the worktree, capture output to `/tmp/picker-e2e.out`, verify zero failures before opening the PR.
  - `npm test -- --run ComparisonPicker SamplePickerRow ComparisonPickerBody` for the Vitest slice.
  - `npm run build` for the TS-strict typecheck gate.
  - Backend slice: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI", test_args=["test_picker_samples_route.jl", "test_picker_routes.jl"])'` — captured to `/tmp/jl-picker.out` per the long-run-once-grep-many CLAUDE.md rule.

### Migration / rollout

- Existing `data-testid="picker-row"` on rows stays. New `data-sample-id` added. `data-exposure-id` is set to the resolved exposure id (the default unless the override is active). Existing E2E selectors that target `data-exposure-id` continue to work because every pick resolves to an exposure id.
- No schema migration. No event-log change. No SSE change. No mutation-queue contact surface.
- `comparison_members.exposure_id` semantics unchanged; the picker is just smarter about which exposure id to send.

## PR2 — inline picker panel in edit mode (#87 §6–7)

### Components

- **`ComparisonPickerPanel.tsx`** (new) — inline shell wrapping the same `ComparisonPickerBody`. No overlay, no focus trap, no Esc-to-close. Different selection-commit semantics (see below).
- **`ComparePageEdit.tsx`** (updated) — right WorkspaceGrid slot swaps from the hint card to `<ComparisonPickerPanel experimentId={eid} />`. The bottom-right "+ Add traces" button (rendered when `plotMembers.length > 0`) is removed. The empty-state CTA inside the plot host stays but its onClick switches to `panelRef.current?.focusSearch()` (or similar) — focuses the inline panel's search input rather than opening a modal.

### Selection-commit semantics

The modal shell has a "Add N selected" footer because the modal is dismissive — the user commits a batch and the modal goes away. The inline panel is always-on, so batching is visual noise. New body-level prop:

```ts
type CommitMode = "batch" | "immediate";
```

- `"batch"` (modal): existing behavior. Picks accumulate; "Add N selected" footer fires `addMember` per pick on click.
- `"immediate"` (inline panel): each toggle/override fires `addMember` (or a future `removeMemberByExposure`) immediately. The picked row visually transitions to the locked state (`alreadyAdded = true`) on the next render via the existing draft-membership lookup. No footer rendered.

`ComparisonPickerBody` already does the lookup; the new mode flips the toggle handler.

### Frontend — tests (Vitest)

- **`test/ComparisonPickerPanel.test.tsx`** (new) — inline shell renders without overlay/dialog role, immediate-commit fires `addMember` on first toggle, picked rows transition to locked.
- **`test/ComparePageEdit.test.tsx`** (updated) — right slot now hosts the panel (testid `comparison-picker-panel`); assert the hint card is gone. The bottom-right "+ Add traces" button assertion deletes; empty-state CTA assertions stay.

### Frontend — E2E

- **`e2e/live/comparePickerInline.spec.ts`** (new, live-mode) — real backend roundtrip. Open `/compare/new`, panel is visible in the right slot, click a sample row, see a band appear on the plot, click the override caret, switch exposure, see the band's resolved exposure id update. Per the live-mode runbook: wait ~800ms after `page.goto("/")` before the first mutation that expects an SSE echo.
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
                      
SamplePickerRow       ◄── ComparisonPickerBody (primary row)
ExposureOverrideRow   ◄── SamplePickerRow      (caret-expanded leaf)

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
| `picker-samples` N+1 against `all_exposures` if experiments have hundreds of samples | One JOIN'd query in `picker_samples`; benchmark in backend test on a 200-sample fixture if needed. Falls under existing perf envelope. |
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
