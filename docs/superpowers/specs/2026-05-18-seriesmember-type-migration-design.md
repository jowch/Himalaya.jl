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
- `lib/comparison/draftFactories.ts` — see the bridge note below.

**Excluded:** `lib/comparison/yBands.ts` — a pure numeric module with no
`ComparisonMember` type; carries over unchanged.

The render pipeline never reads the member's foreign key (`comparison_id` /
`series_id`) — verified by grep across all listed modules. The migration is
therefore a pure annotation swap with no behavioural change. The detailed plan
confirms, per module, that no other comparison-specific field or import is
touched.

### Test fixtures

`tsc --noEmit` runs against `tsconfig.json`, whose `include` covers `test/`.
The render-module Vitest suites build their own `ComparisonMember`-typed member
fixtures inline (no shared fixture module) and pass them to the migrated
components. Those fixtures must therefore retype to `SeriesMember`
(`comparison_id` → `series_id`) for `tsc` to pass — a mechanical change, not a
behavioural one. Test files that reference `ComparisonMember` but do **not**
feed a migrated component (`draftPersistence`, `MentionChip`,
`queue/saveComparison`) stay on `ComparisonMember`.

## The Compare-era bridge

`SeriesMember` carries a required `series_id`, so `ComparisonMember` is **not**
structurally assignable to it. Between this issue (Wave C) and I3.6 (Wave E),
`Compare.tsx` is still live and still drives comparison data into the retyped
render pipeline. Two bridge sites are needed; both are localized, both are
deleted wholesale at I3.6, and neither changes render behaviour:

1. **`Compare.tsx`** builds `plotMembers` through its own local
   `draftToMember(d): ComparisonMember` converter and passes the result to
   `MultiTracePlot` / `MemberMetaGutter` / the export adapter. The converter's
   return type and the two `sampleIdFor` annotations retarget to `SeriesMember`;
   `draftToMember` supplies a placeholder `series_id` (consistent with the
   placeholder `id` it already fabricates for unsaved drafts).
2. **`draftFactories.ts`** — `fromComparison` / `fromComparisonAsFork` are
   comparison-only factories (parameter type `Comparison`) that map
   `Comparison.members` (`ComparisonMember[]`) into the now-`SeriesMember`-typed
   `memberFromSaved`. Those call sites get a localized cast. `memberFromSaved`
   itself migrates cleanly (it reads only shared fields); `memberFromNewExposure`
   has no member-type parameter and is untouched.

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
- `tsc --noEmit` passes (the Compare-era bridge keeps `Compare.tsx` and
  `draftFactories.ts` compiling).
- The render pipeline's existing Vitest passes unchanged — no behavioural change.

## Out of scope

- The `Series` parent type, the folio / scoping / builder UIs (#173–#176).
- The comparison-to-series data migration (#171).
- Any render-behaviour change — the pipeline carries over as-is.
- Retiring `ComparisonMember` / `Compare.tsx` — that is I3.6 (#177); this issue
  only adds the bridge that keeps them compiling until then.
