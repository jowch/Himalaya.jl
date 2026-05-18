# `SeriesMember` type migration (I3.2 / #172) — design

> Issue: #172 — "Migrate the trace-render pipeline to the `SeriesMember` type (I3.2)".
> References: master plan §2.2 / §6.2 (`docs/superpowers/plans/2026-05-17-himalaya-ui-redesign.md`), issue map I3.2 (`docs/superpowers/plans/2026-05-17-himalaya-ui-redesign-issue-map.md`).

## Summary

The multi-trace render pipeline is typed on `ComparisonMember`. Phase 3 builds a
series builder on the same render pipeline, so the pipeline's **input type** must
stop being comparison-specific. This issue introduces a `SeriesMember` type and
retypes the render pipeline onto it. Render *behaviour* is untouched — only the
type changes.

Dependency #165 (I2.2 — series REST routes) has merged; `SeriesMember` is defined
against the real `fetch_series_with_plate` response shape.

## The type: a full standalone interface

`SeriesMember` is defined as a **complete, self-standing interface** in
`api.ts` — not derived from `ComparisonMember`. It mirrors the I2.2 series
member route response (`fetch_series_with_plate`, `packages/HimalayaUI/src/series.jl`):

```ts
/** Per-member shape returned by the series GET / POST endpoints. Mirrors `fetch_series_with_plate`. */
export interface SeriesMember {
  id: number;
  series_id: number;
  exposure_id: number | null;
  display_order: number;
  band_height: number;
  y_offset: number;
  normalization: string;
  color_override: string | null;
  label_override: string | null;
  q_window_min: number | null;
  q_window_max: number | null;
  peak_display: unknown;
  snapshot: MemberSnapshot | null;
  is_stale: boolean;
  created_by: number | null;
  created_at: string | null;
}
```

It is field-for-field identical to `ComparisonMember` except `comparison_id` →
`series_id`. `MemberSnapshot` is a shared, non-comparison-specific type and is
reused as-is. Optional fields are modelled `T | null`, mirroring `ComparisonMember`
for `exactOptionalPropertyTypes`.

### Rejected alternative: `type SeriesMember = Omit<ComparisonMember, 'comparison_id'>`

A derived definition was considered and rejected:

- **The end state retires `ComparisonMember`.** I3.6 (#177) deletes `Compare.tsx`
  and the Compare-only components — the only frontend consumers of
  `ComparisonMember`. The type then becomes dead frontend code, swept at I5.3.
  A definition derived via `Omit<ComparisonMember, …>` cannot survive its base
  type's deletion; the full standalone interface must be written eventually
  regardless. `Omit` only defers that work into an unrelated cleanup PR.
- **The derivation arrow points the wrong way.** `SeriesMember` is the permanent,
  going-forward type; `ComparisonMember` is frozen and scheduled for deletion.
  Deriving the survivor from the doomed type inverts the dependency and blocks
  `SeriesMember` from gaining series-specific fields without accreting
  `Omit<…> & { … }`.
- The `Omit` form also dropped `series_id`, so it failed the issue's own
  "mirrors the route response shape" requirement.

The standalone interface wins on permanence, honesty, independent evolution, and
literal route fidelity. Its only cost is the Compare-era bridge below.

## Scope: modules retyped

Retype the member-type annotations from `ComparisonMember` to `SeriesMember` in
the modules the issue lists:

- `components/MultiTracePlot.tsx`
- `components/MemberTraceLayer.tsx`
- `lib/figure-export/adapters/multiTraceAdapter.ts`
- `lib/figure-export/marks/multiTraceExportMarks.ts`
- `components/MemberMetaGutter.tsx`
- `components/MemberMetaRow.tsx`
- `lib/comparison/coloring.ts`
- `lib/comparison/labels.ts`
- `lib/comparison/draftFactories.ts` — partial; see below and the bridge note.

**Excluded:** `lib/comparison/yBands.ts` — a pure numeric module with no
`ComparisonMember` type; carries over unchanged.

The first eight modules are the render pipeline proper; it never reads the
member's foreign key (`comparison_id` / `series_id`) — verified by grep across
all eight — so the migration there is a pure annotation swap with no behavioural
change.

`draftFactories.ts` is listed by issue #172 and master plan §2.2/§6.2, and
belongs in scope for a forward-looking reason rather than a render one: its
`memberFromSaved` is the reusable saved-member → draft helper (it recomputes the
snapshot against the live cache), and the Phase 3 series builder (I3.5b) will
reuse it. Migrating its parameter to `SeriesMember` now readies it for that
reuse — the same rationale that migrates the render modules. **Only
`memberFromSaved` migrates**; the comparison-only factory wrappers
(`fromComparison` / `fromComparisonAsFork`) keep their `Comparison` parameter and
`memberFromNewExposure` has no member-type parameter — none of those three
migrate (see the bridge below). The detailed plan confirms, per module, that no
other comparison-specific field or import is touched.

### Test fixtures

The repo's *enforced* typecheck gate is `npm run build`, which runs
`tsc --noEmit -p tsconfig.build.json` — and `tsconfig.build.json` **excludes**
`test/` and `e2e/`. The bare `tsc --noEmit` (resolving to `tsconfig.json`,
whose `include` covers `test/`) is **not** a clean gate and never has been: on
`main` it reports ~183 pre-existing errors, almost all in `test/` — Vitest
globals (`test`, `expect`, `describe`) that `tsconfig.json` never wires up, plus
some stale pre-existing fixtures. The repo never typechecks `test/` with `tsc`;
Vitest type-erases. (Issue #172's acceptance bullet says "`tsc --noEmit`
passes"; read against the repo that means the enforced `npm run build`
typecheck, which is `src/`-only and green.)

The render-module Vitest suites build their own `ComparisonMember`-typed member
fixtures inline (no shared fixture module) and pass them to the migrated
components. Once a component is retyped to `SeriesMember`, a fixture left as
`ComparisonMember` would become a *new* bare-`tsc` error (missing `series_id`,
excess `comparison_id`). So the fixtures retype to `SeriesMember`
(`comparison_id` → `series_id`) — not to make bare `tsc` *pass* (it cannot,
given the pre-existing errors), but to keep this migration from *adding* to that
error count and to keep each fixture honest against the component it feeds. A
mechanical change, not a behavioural one. Test files that reference
`ComparisonMember` but do **not** feed a migrated component (`draftPersistence`,
`MentionChip`, `queue/saveComparison`) stay on `ComparisonMember`.

## The Compare-era bridge

`SeriesMember` carries a required `series_id`, so `ComparisonMember` is **not**
structurally assignable to it. Between this issue (Wave C) and I3.6 (Wave E),
`Compare.tsx` is still live and still drives comparison data into the retyped
render pipeline. Two bridge sites are needed; both are localized, both are
deleted wholesale at I3.6, and neither changes render behaviour:

1. **`Compare.tsx`** drives *two* member arrays into the retyped pipeline. In
   **edit mode** it builds `plotMembers` through its local
   `draftToMember(d): ComparisonMember` converter. In **review mode** it sorts
   `compQ.data.members` (a `ComparisonMember[]`) into a `members` memo. Both
   arrays feed `MultiTracePlot` / `MemberMetaGutter` / `resolveDisplayLabels` /
   `buildMultiTraceExportSpec`. The `draftToMember` return type, the review-mode
   `members` memo's element type, and the two `sampleIdFor` annotations all
   retarget to `SeriesMember`. `draftToMember` already fabricates a placeholder
   foreign key; the review-mode memo spreads `{ ...m, series_id: 0 }` per member
   — the same value-irrelevant bridge as `draftFactories.ts`, since the pipeline
   never reads the key.
2. **`draftFactories.ts`** — once `memberFromSaved` is retyped to `SeriesMember`
   (see Scope above), its comparison-only callers `fromComparison` /
   `fromComparisonAsFork` (parameter type `Comparison`) still map
   `Comparison.members` (`ComparisonMember[]`) into it. Each call site bridges by
   spreading in a placeholder `series_id` (`{ ...m, series_id: 0 }`) — the
   spread-introduced `comparison_id` is exempt from excess-property checking, and
   `memberFromSaved` never reads `series_id`. `memberFromNewExposure` has no
   member-type parameter and is untouched.

The bridge is throwaway by construction: every bridge site lives in a file I3.6
deletes. No `RenderMember` base type is introduced — the master plan describes a
1:1 `ComparisonMember` → `SeriesMember` swap, and a third named type would be
scope the plan does not sanction.

## Coordination notes

- **`api.ts` contention with the series event-kind cluster.** #166/#167/#168
  (I2.3–I2.5) add a `Series` parent type, mutating fetchers, and query keys to
  `api.ts` (per project memory `series-event-cluster-frontend-scaffolding`).
  This issue adds only `SeriesMember`. If the cluster lands first and has
  defined a member type for `Series.members`, reconcile to a single
  `SeriesMember` owned by this issue; otherwise this issue owns the type and the
  cluster / I3.3 rebase `Series.members` onto it. The additions are textually
  separate (distinct interface blocks) — low conflict, but flag it at rebase.
- **`Series` parent type is out of scope.** This issue defines only the member
  type. The `Series` parent interface belongs to the event-kind cluster / I3.3.

## Acceptance criteria

- `SeriesMember` is defined as a standalone interface in `api.ts`, mirroring the
  I2.2 route response; optional fields modelled `T | null`.
- The listed modules are retyped onto `SeriesMember`.
- `lib/comparison/yBands.ts` is unchanged.
- `npm run build` passes — the enforced typecheck gate
  (`tsc --noEmit -p tsconfig.build.json`, `src/`-only) plus `vite build`. The
  Compare-era bridge keeps `Compare.tsx` and `draftFactories.ts` compiling.
- Bare `tsc --noEmit` reports no *new* errors beyond the ~183 pre-existing ones
  (see "Test fixtures") — the retyped fixtures keep the render-module suites
  from adding to the count.
- The render pipeline's existing Vitest passes unchanged — no behavioural change.

## Out of scope

- The `Series` parent type, the folio / scoping / builder UIs (#173–#176).
- The comparison-to-series data migration (#171).
- Any render-behaviour change — the pipeline carries over as-is.
- Retiring `ComparisonMember` / `Compare.tsx` — that is I3.6 (#177); this issue
  only adds the bridge that keeps them compiling until then.
