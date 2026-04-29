# Loading Transitions & Skeleton Screens Design

## Context

HimalayaUI's current loading states are uniformly bare — every pending data fetch shows a plain italic `<HintText>` string ("Loading trace…", "Loading phase assignments…", "loading experiments…"). Switching samples or exposures triggers cascading TanStack Query fetches across multiple cards simultaneously, leaving large empty regions during what can be a frequently-repeated action in normal SAXS curation workflows.

The goal is subtle polish: make data arrival feel intentional rather than sudden, without adding visual noise that competes with the scientific content. Loud spinners or dramatic transitions are explicitly out of scope.

## Approach

Use **boneyard-js** throughout — a library that auto-generates pixel-perfect skeleton screens by capturing real component layout via a headless browser. The Vite plugin re-captures skeletons on every HMR update, so skeletons stay in sync with the UI automatically. Combined with boneyard's built-in `stagger` and `transition` props, this delivers pulse skeletons + stagger-on-arrival with zero per-card animation code.

Page tab transitions (Index → Inspect → Compare) are deferred.

## Scope

Seven targets:

| Target | File | Loading condition |
|---|---|---|
| `PlotCard` | `components/PlotCard.tsx` | `traceQ.isLoading \|\| peaksQ.isLoading` |
| `PhasePanel` | `components/PhasePanel.tsx` | `indicesQ.isLoading \|\| groupsQ.isLoading` |
| `ChatCard` | `components/ChatCard.tsx` | `messagesQ.isLoading` |
| `DetectorImageCard` | `components/DetectorImageCard.tsx` | local image loading state |
| `SampleMetadataCard` | `components/SampleMetadataCard.tsx` | `exposuresQ.isLoading` |
| NavModal — experiments list | `components/NavModal.tsx` | `experimentsQ.isLoading` |
| NavModal — samples list | `components/NavModal.tsx` | `samplesQ.isLoading` |

## Build Setup

### Installation

```bash
cd packages/HimalayaUI/frontend
npm install boneyard-js
```

### Vite plugin

```ts
// vite.config.ts
import { boneyardPlugin } from 'boneyard-js/vite'

export default defineConfig({
  plugins: [
    // ... existing plugins
    boneyardPlugin()
  ]
})
```

### Global config

```json
// packages/HimalayaUI/frontend/boneyard.config.json
{
  "breakpoints": [1024, 1440, 1920],
  "out": "./src/bones",
  "animate": "pulse",
  "boneClass": "bone",
  "color": "rgba(42, 44, 52, 1)",
  "darkColor": "rgba(50, 52, 62, 1)"
}
```

The color values are rgba equivalents of HimalayaUI's card background tokens (`oklch(0.20 0.01 250)` / `oklch(0.22 0.01 250)`). Boneyard accepts hex/rgba but not oklch.

Breakpoints cover both the stacked single-column layout (<1400px) and the three-column workspace grid (≥1400px) plus large monitor widths.

### Registry

```ts
// main.tsx — add once at the top
import './bones/registry'
```

The `./src/bones/` directory is gitignored (generated artefacts). Add `src/bones/` to `.gitignore`. `boneyard.config.json` is committed — it's config, not a generated artefact.

## Animation Parameters

Applied globally via `boneyard.config.json` and overridden per component only if needed:

| Parameter | Value | How set | Rationale |
|---|---|---|---|
| `animate` | `"pulse"` | `boneyard.config.json` | Calmer than shimmer; suits a focused scientific tool |
| Pulse duration | `1.8s` | CSS via `boneClass` | Slower than boneyard's default — less visually insistent |
| `stagger` | `50` (ms) | Per `<Skeleton>` prop | Barely perceptible delay between bones on reveal |
| `transition` | `200` (ms) | Per `<Skeleton>` prop | Gentle fade-out when real content replaces the skeleton |

`stagger` and `transition` are set on each `<Skeleton>` component (boneyard config doesn't expose them globally). The pulse duration is not a boneyard prop — override it via a global CSS rule targeting the class passed to `boneClass`:

```css
/* styles.css */
.bone { animation-duration: 1.8s; }
```

Pass `boneClass="bone"` on each `<Skeleton>` (or set it globally in the config as `"boneClass": "bone"`).

## Workspace Card Skeletons

The `<Skeleton>` wrapper goes around the **inner content area** of each card — not the outer card chrome (border, header label, section labels). The card frame stays visible during loading, keeping the layout stable. Only the data-driven region is replaced by bones.

### Loading condition

Use `isLoading` (`isPending && isFetching`) not `isPending` alone. Disabled queries (e.g. `useTrace` when `exposureId` is `undefined`) have `isPending: true` but `isFetching: false`, so `isLoading` is false. Skeletons only appear when actively fetching. The existing `HintText` empty states for the no-selection case are left untouched.

### Fixtures

All five workspace cards require a `fixture` prop because their queries are disabled until an experiment/sample/exposure is selected — with no real content, the headless CLI has nothing to measure. Fixtures are `dev`-only (boneyard strips them in production) and contain minimal plausible mock data:

| Card | Fixture content |
|---|---|
| `PlotCard` | Fake `q` / `I` arrays (20 points), 3 peak objects with `q` and `prominence` |
| `PhasePanel` | 3 index rows: `{ phase: "Pn3m", score: 0.91, npeaks: 4 }`, etc. |
| `ChatCard` | 2 short message objects with `username` and `body` |
| `DetectorImageCard` | 1×1 transparent PNG data URI (just to trigger render) |
| `SampleMetadataCard` | 2 exposure stubs with `id`, `filename`, `status: "accepted"` |

Fixtures are defined inline as constants in each component file, wrapped in `process.env.NODE_ENV === 'development'` guards or used directly as the `fixture` prop — boneyard ignores the prop in production builds.

### Example (PlotCard)

```tsx
import { Skeleton } from 'boneyard-js/react'

const FIXTURE = /* minimal fake trace + peaks */

function PlotCard() {
  const { traceQ, peaksQ } = /* existing hooks */
  const loading = traceQ.isLoading || peaksQ.isLoading

  return (
    <div className="card">
      <SectionLabel>Trace</SectionLabel>
      <Skeleton name="plot-card" loading={loading} stagger={50} transition={200} fixture={FIXTURE}>
        {/* existing content unchanged */}
      </Skeleton>
    </div>
  )
}
```

## NavModal Skeletons

The NavModal renders two conditional branches based on `navModalStep`: the experiments list and the samples list. Each branch gets its own `<Skeleton>` wrapping the list area. Skeletons are independent — only the active list shows a skeleton.

```tsx
// NavModal.tsx — experiments branch
<Skeleton name="nav-experiments" loading={experimentsQ.isLoading} stagger={50} transition={200}
  fixture={NAV_FIXTURE_EXPERIMENTS}>
  {/* existing experiments list */}
</Skeleton>

// NavModal.tsx — samples branch
<Skeleton name="nav-samples" loading={samplesQ.isLoading} stagger={50} transition={200}
  fixture={NAV_FIXTURE_SAMPLES}>
  {/* existing samples list */}
</Skeleton>
```

Fixtures are 4 hardcoded row objects each — enough to fill the typical list height.

```ts
const NAV_FIXTURE_EXPERIMENTS = [
  { id: 1, name: "Experiment A" },
  { id: 2, name: "Experiment B" },
  { id: 3, name: "Experiment C" },
  { id: 4, name: "Experiment D" },
]

const NAV_FIXTURE_SAMPLES = [
  { id: 1, name: "JC001" },
  { id: 2, name: "JC002" },
  { id: 3, name: "JC003" },
  { id: 4, name: "JC004" },
]
```

## Files to Create / Modify

| File | Change |
|---|---|
| `packages/HimalayaUI/frontend/package.json` | Add `boneyard-js` dependency |
| `packages/HimalayaUI/frontend/vite.config.ts` | Add `boneyardPlugin()` |
| `packages/HimalayaUI/frontend/boneyard.config.json` | Create — global animation config |
| `packages/HimalayaUI/frontend/.gitignore` | Add `src/bones/` |
| `packages/HimalayaUI/frontend/src/main.tsx` | Add `import './bones/registry'` |
| `packages/HimalayaUI/frontend/src/components/PlotCard.tsx` | Wrap content in `<Skeleton>`, add fixture, switch to `isLoading` |
| `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx` | Wrap content in `<Skeleton>`, add fixture, switch to `isLoading` |
| `packages/HimalayaUI/frontend/src/components/ChatCard.tsx` | Wrap content in `<Skeleton>`, add fixture, switch to `isLoading` |
| `packages/HimalayaUI/frontend/src/components/DetectorImageCard.tsx` | Wrap content in `<Skeleton>`, add fixture |
| `packages/HimalayaUI/frontend/src/components/SampleMetadataCard.tsx` | Wrap content in `<Skeleton>`, add fixture, switch to `isLoading` |
| `packages/HimalayaUI/frontend/src/components/NavModal.tsx` | Wrap experiments + samples lists in `<Skeleton>`, add fixtures |

## Deferred

- Page tab transitions (Index → Inspect → Compare) — out of scope for this iteration
- Optimistic updates on mutations — the stagger reveal on refetch data arrival is sufficient for now
