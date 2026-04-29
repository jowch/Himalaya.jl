# Chat @-mention System — Design Spec

## Context

The Himalaya chat (ChatCard) is a per-sample plain-text notebook. Researchers writing notes about SAXS analyses need to reference specific peaks, phase assignments, exposures, and samples in a way that stays connected to the underlying data — so a note like "shoulder at q=1.223 matches Pn3m on JC001-120.dat" becomes navigable and hover-inspectable rather than inert prose. The centralized DB model (all experiments in one DB) makes cross-experiment references feasible and scientifically useful, since researchers compare replicates across experiments.

## Storage

No schema changes. `sample_messages.body TEXT` stores `[[type:id]]` tokens inline:

```
"shoulder at [[peak:42]] matches [[index:17]] on [[exposure:33]], same as [[sample:8]]"
```

`type` is one of: `peak`, `index`, `exposure`, `sample`, `experiment`. `id` is the integer primary key from the corresponding table. The pair `(type, id)` is globally unambiguous. Tokens survive renames (resolved at render time) and handle deletions gracefully (gray-out). No separate mentions table — query-by-mention is explicitly deferred.

## Compose UX

The existing `<textarea>` in `ChatCard.tsx` gains `@`-trigger behavior via a new `MentionCompose` wrapper component. When the user types `@`:

- A `MentionPicker` dropdown appears **above** the textarea (anchored to the compose box, not the caret)
- Pretext (`@chenglou/pretext`) handles textarea auto-resize without triggering DOM reflow
- Typing continues to filter results; `Esc` dismisses; `↑↓` navigate; `Enter` or click selects
- On selection, `[[type:id]]` is inserted at the `@` position in the raw body string
- Submit path is unchanged — the raw token string goes to the existing `usePostSampleMessage` mutation

### Picker — search scope and ranking

Default `@<query>` searches client-side against already-loaded TanStack Query data, grouped in this order:

1. **Indices** — active exposure
2. **Peaks** — active exposure
3. **Exposures** — active sample
4. **Samples** — active experiment

Type-prefix escape hatch (optional, narrows scope explicitly):

| Prefix | Scope |
|---|---|
| `@index:` | indices, active exposure |
| `@peak:` | peaks, active exposure |
| `@exposure:` | exposures, active sample |
| `@sample:` | samples, active experiment |
| `@experiment:` | all experiments in DB |

Picker rows show phase color matching existing `phases.ts` palette for index results. No ID shown — the picker is the only place the user interacts with the raw token.

## Resolution Strategy — Hybrid (client-first, lazy fallback)

`useMentionResolution(tokens: string[])` accepts an array of `[[type:id]]` tokens extracted from a message body and returns `Map<token, ResolvedEntity | "loading" | "dead">`.

1. **Cache hit** — if TanStack Query already has the entity loaded (common case: current exposure's peaks/indices), resolve instantly with no fetch
2. **Cache miss** — fire `GET /api/{type}s/:id` for the unknown entity; show a neutral loading state in the chip until resolved
3. **404** — mark as `"dead"`; chip renders grayed-out with "no longer exists" tooltip

## Rendered Chips

`renderMentions(body: string, resolutionMap: Map)` is a pure function that splits body text on `[[type:id]]` tokens and returns an array of plain string nodes and `MentionChip` React elements. No side effects; straightforward to unit-test.

### Chip visual design

**Default state** — underline + ambient tinted background per entity type:

| Type | Color | Background |
|---|---|---|
| peak | `#7cb8e8` | `#15222e` |
| index | `#b5a0d8` | `#1e1828` |
| exposure | `#88c0a8` | `#162018` |
| sample | `#c0b878` | `#201c10` |
| dead | `#484848` | `#181818` |

**Hover state** — brightened color, 2px underline, slightly deeper background.

**Chip display text by entity type:**

| Type | Chip text | Example |
|---|---|---|
| peak | `q = {q}` | `q = 1.223` |
| index | `{phase} · {score}` | `Pn3m · 0.91` |
| exposure | `{filename}` | `JC001-120.dat` |
| sample | `{name}` | `copper-B` |
| experiment | `{name}` | `Run 2024-03` |
| dead | original text, grayed | `Pn3m · 0.91` |

### Tooltip content (floating above chip, `position: absolute; bottom: 100%`)

Tooltip shows only what is **not already visible in the chip**. Single compact line, `·` separated.

| Type | Tooltip fields |
|---|---|
| peak | `source · prominence` |
| index (cubic) | `q₁ · d · R² · ngc` (q₁ = `predicted_q[1]`) |
| index (non-cubic) | `q₁ · d · R²` (q₁ = `predicted_q[1]`) |
| exposure | `status` |
| sample | `experiment` (only when cross-experiment) |
| dead | `no longer exists` |

Cubic phases for the `ngc` conditional: Pn3m, Im3m, Ia3d, Fm3m, Fd3m. Same set as the existing phase categorization in `phases.ts`.

### Hover behavior — Zustand wiring

- **Index chip hover** → sets `hoveredIndexId` (already exists in `state.ts`) — trace viewer and Miller plot highlight fire automatically
- **Peak chip hover** → sets `hoveredPeakId` (new, mirrors `hoveredIndexId` pattern) — trace viewer highlights the corresponding tick
- **Exposure/sample chip hover** → tooltip only, no highlight state

### Click behavior — navigation

- **Same-experiment chips** → calls `setActiveSampleId` + `setActiveExposureId` (existing Zustand actions). Navigation feels like a selection change, not a page transition.
- **Cross-experiment chips** → tooltip context shown; click is a no-op. Cross-experiment navigation is deferred (the current UI has no experiment-switcher).

## New Backend Routes

Direct entity-lookup routes needed for cache-miss resolution. Existing routes are exposure-scoped (`GET /api/exposures/:id/peaks`); these are entity-scoped:

| Route | Notes |
|---|---|
| `GET /api/peaks/:id` | New |
| `GET /api/index_groups/:id` | New |
| `GET /api/exposures/:id` | Confirm whether exists |
| `GET /api/samples/:id` | Confirm whether exists |
| `GET /api/experiments/:id` | Confirm whether exists |

Response shapes reuse existing JSON serialization from `json.jl`. Only called on cache misses — same-exposure entities resolve instantly from already-loaded queries.

## Out of Scope (deferred)

- Query-by-mention ("which notes reference this peak?") — needs a `message_mentions` join table; adds complexity not worth it yet
- Cross-experiment navigation on chip click — needs an experiment-switcher in the UI
- Mention editing (changing a tag after insertion) — backspace removes the whole token; re-type to re-insert
- Mention notifications or highlights in the message list
- Rich text beyond mentions (bold, code, links)

## Files to Create / Modify

**New:**
- `frontend/src/components/MentionCompose.tsx` — compose wrapper with `@` trigger + pretext sizing
- `frontend/src/components/MentionPicker.tsx` — autocomplete dropdown
- `frontend/src/components/MentionChip.tsx` — inline chip with hover/click wiring
- `frontend/src/lib/renderMentions.tsx` — pure parser/renderer function
- `frontend/src/hooks/useMentionResolution.ts` — hybrid resolution hook

**Modified:**
- `frontend/src/components/ChatCard.tsx` — swap raw textarea for `MentionCompose`; wrap `MessageRow` body with `renderMentions`
- `frontend/src/state.ts` — add `hoveredPeakId: number | null` + `setHoveredPeakId` action
- `frontend/src/queries.ts` — add per-entity query hooks (`usePeak(id)`, `useIndexGroup(id)`, etc.)
- `packages/HimalayaUI/src/routes_*.jl` — add direct-lookup routes
- `frontend/package.json` — add `@chenglou/pretext` dependency

## Verification

1. `npm test` — unit tests for `renderMentions` (token parsing, dead-token handling), `useMentionResolution` (cache hit, cache miss, 404), `MentionChip` (per-type render, hover state)
2. Playwright E2E — type `@`, verify picker appears; select an index; verify chip renders in sent message; hover chip and verify tooltip; click chip and verify active exposure changes
3. Manual — cross-sample reference: write note on sample A mentioning an exposure from sample B; navigate to sample B; verify the exposure chip resolves correctly
4. Manual — dead reference: write note referencing an index; reanalyze to clear it; verify chip grays out with "no longer exists" tooltip
