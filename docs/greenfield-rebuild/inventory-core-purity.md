# Core-purity ledger (greenfield rebuild)

## Shared core importable by src/print/
- `src/api.ts`, `src/queries.ts`, `src/state.ts`, `src/phases.ts`
- `src/lib/**` — all pure-logic files are clean; the only impure file is `src/lib/figure-export/**` (see below)
- `src/hooks/**` — all logic hooks are clean
- `src/ErrorBoundary.tsx` (shared)

## Forbidden from src/print/ (enforced by no-legacy-import guard)
- `src/components/**` (incl. old `src/components/ui/**`)
- `src/pages/**`

## Known coupling — lib/figure-export → renderer components (NOT pure core)

figure-export imports old renderer components (peakMark, MemberHeatmapLayer,
CrossTraceTrackingLayer, RepresentationToggle). Expected: it reuses renderer mark
logic. NOT a foundation blocker (print/ does not import figure-export). Resolved
when renderers are rebuilt in print/ and figure-export is repointed.

Exact Step-1 hits:

```
src/lib/figure-export/marks/traceExportMarks.ts:8:import { peakMark } from "../../../components/ui/peakMark";
src/lib/figure-export/adapters/multiTraceAdapter.ts:11:import type { Representation } from "../../../components/RepresentationToggle";
src/lib/figure-export/marks/multiTraceExportMarks.ts:14:import { buildMemberHeatmapMarks } from "../../../components/MemberHeatmapLayer";
src/lib/figure-export/marks/multiTraceExportMarks.ts:15:import { peakMark } from "../../../components/ui/peakMark";
src/lib/figure-export/marks/multiTraceExportMarks.ts:16:import { buildCrossTraceTrackingMarks } from "../../../components/CrossTraceTrackingLayer";
src/lib/figure-export/marks/multiTraceExportMarks.ts:17:import type { Representation } from "../../../components/RepresentationToggle";
```

## Surprise entanglements to extract (must be empty before component rebuilds)

none

## Renderer geometry functions to reuse (Option 2: rebuild visuals, reuse math)

Exact Step-2 output:

```
src/lib/plot/formatAxis.ts:10:export function formatAxis(d: number): string
src/lib/plot/invertQ.ts:38:export function invertQ(plot: unknown, px: number): number | null
src/lib/plot/invertQ.ts:50:export function applyQ(plot: unknown, q: number): number | null
src/lib/plot/sparkline.ts:3:export const SPARK_W = 76;
src/lib/plot/sparkline.ts:4:export const SPARK_H = 28;
src/lib/plot/sparkline.ts:15:export function sparklinePath(trace: Trace): string | null
src/lib/qRing.ts:11:export const RING_VIEWBOX = 100;
src/lib/qRing.ts:13:export const RING_R_MIN = 12;
src/lib/qRing.ts:15:export const RING_R_SPAN = 33;
src/lib/qRing.ts:18:export function qToRadius(q: number, qLo: number, qHi: number): number
src/lib/qRing.ts:25:export function nearestRingQ(
src/lib/detectorOrient.ts:6:export const ROTATE_THRESHOLD = 1.25;
src/lib/detectorOrient.ts:9:export const ROTATE_MIN_VIEWPORT = 1400;
src/lib/detectorOrient.ts:24:export function decideOrient(i: OrientInput): OrientResult
src/lib/plot/labelDodge.ts:72:export function layoutPeakLabels(
```

Grouped by file:

### src/lib/plot/formatAxis.ts
- `formatAxis(d: number): string` — axis tick label formatter

### src/lib/plot/invertQ.ts
- `invertQ(plot, px): number | null` — pixel → q-value inversion
- `applyQ(plot, q): number | null` — q-value → pixel projection

### src/lib/plot/sparkline.ts
- `SPARK_W = 76` — sparkline canvas width constant
- `SPARK_H = 28` — sparkline canvas height constant
- `sparklinePath(trace: Trace): string | null` — SVG path string for trace sparkline

### src/lib/plot/labelDodge.ts
- `layoutPeakLabels(...)` — peak label layout/dodging algorithm

### src/lib/qRing.ts
- `RING_VIEWBOX = 100` — q-ring SVG viewBox size
- `RING_R_MIN = 12` — minimum ring radius
- `RING_R_SPAN = 33` — ring radius span
- `qToRadius(q, qLo, qHi): number` — maps q-value to SVG ring radius
- `nearestRingQ(...)` — snaps a pointer position to the nearest ring q-value

### src/lib/detectorOrient.ts
- `ROTATE_THRESHOLD = 1.25` — aspect ratio threshold for detector rotation
- `ROTATE_MIN_VIEWPORT = 1400` — minimum viewport width for rotation
- `decideOrient(i: OrientInput): OrientResult` — decides detector panel orientation from viewport geometry

## Notes on src/lib geometry inventory

The following geometry/plot files exist under `src/lib/`:
- `src/lib/plot/` — `formatAxis.ts`, `invertQ.ts`, `labelDodge.ts`, `sparkline.ts`
- `src/lib/qRing.ts`
- `src/lib/detectorOrient.ts`

All other `src/lib/` subdirectories (`assignment.ts`, `authOpts.ts`, `clientId.ts`,
`clientOpId.ts`, `color-distance.ts`, `comparison/`, `customIndex.ts`,
`figure-export/`, `queue/`, `sample/`, `scoping/`, `series/`, `seriesRatio.ts`,
`toast.ts`, `units.ts`, `url/`) are application-logic or adapter utilities
(not SAXS/rendering geometry), and are importable by `src/print/` without
requiring a rebuild.
