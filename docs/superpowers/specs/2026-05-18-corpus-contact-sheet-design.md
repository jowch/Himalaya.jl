# Corpus contact-sheet sample table (I1.4) — design

**Issue:** #160 (I1.4) · **Date:** 2026-05-18 · **Milestone:** HimalayaUI workflow-driven redesign

## Summary

Build the corpus contact-sheet sample table at `/samples` — the redesign's
primary triage surface, replacing the Inspect page. The view renders the whole
corpus as a table; each row shows a sample's identity, an exposure-thumbnail
strip, a Kept count, a free-form Tags column, and a Status. A `?beamtime=<id>`
URL query filters the table to one experiment, and a boneyard skeleton covers
the corpus query's loading state.

This issue builds the **table view only**. Culling interactions (#162) and tag
mutations (#159) wire in separately; this view renders their affordances inert.

## Rationale

The redesign's first surface is a contact sheet over the whole corpus: the user
scans samples, triages exposures, and tags as they go. It replaces the Inspect
page, whose single-sample framing does not scale to corpus-wide triage. The
contact sheet is the entry point of the three-surface workflow model
(redesign-notes §4 — Triage → Index → Series).

## Verified context

The dependencies are merged and the route slot exists:

- **#155 (I1.1) — app shell.** `CorpusShell` renders `CorpusTopbar` plus an
  `<Outlet/>`; `AppRoutes` mounts `<Route path="/samples" element={<SamplesPage/>}/>`
  under the `CorpusShell` layout route. `SamplesPage.tsx` is a placeholder body
  this issue rewrites.
- **#156 (I1.2) — corpus query layer.** `useCorpusSamples()` fetches
  `GET /api/samples` and returns `CorpusSample[]`. `CorpusSample` carries
  `id, experiment_id, name, display_name, notes, tags, q_units` — **identity and
  tags only.** It does **not** carry exposures, exposure status, or any indexed
  phase. The hook deliberately exposes no `?experiment_id=` filter (#156 spec,
  "Verified context").
- **`CorpusTopbar.tsx`** renders an inert "Beamtime" chip. Its own comment
  states `?beamtime=` URL query state is "owned by the `/samples` route (#160)."

Consequences that shape the design:

1. **The exposure strip and Kept count need a second data source.** `CorpusSample`
   has no exposures; `ThumbnailGallery` / `DetectorImage size="thumb"` require
   full `Exposure` objects (`image_path`, `image_version`, `status`,
   `selected`). Each row therefore fetches its own exposures.
2. **The Status column has no data source in #160.** A real phase call needs
   per-sample index queries — a dependency #160 does not declare. Status renders
   a fixed "Not indexed" placeholder behind a typed seam.
3. **`?beamtime=` filters client-side.** `useCorpusSamples()` fetches the whole
   corpus and exposes no filter; the view filters the returned rows.

## Design

Three files: one rewrite, one new component, one small modification.

| File | Action |
|---|---|
| `packages/HimalayaUI/frontend/src/pages/SamplesPage.tsx` | Rewrite — the contact-sheet container |
| `packages/HimalayaUI/frontend/src/components/ContactSheetRow.tsx` | New — one sample row |
| `packages/HimalayaUI/frontend/src/components/CorpusTopbar.tsx` | Modify — make the Beamtime chip a functional experiment picker |
| `packages/HimalayaUI/frontend/test/contact-sheet.test.tsx` | New — Vitest |

### 1. `SamplesPage` — the contact-sheet container

`SamplesPage` owns the corpus query, the `?beamtime=` URL state, the
client-side filter, and the boneyard skeleton.

- Calls `useCorpusSamples()` once.
- Reads `?beamtime=` via `useSearchParams()`. The value is an `experiment_id`
  (experiment is the corpus's provenance facet — redesign-notes §5 architecture
  decision 1). Absent / unparseable → no filter.
- Filters: `samples.filter(s => beamtime === undefined || s.experiment_id === beamtime)`.
- Renders, top to bottom:
  - **Header** — a "Contact sheet" kicker and a sub-line naming the active
    scope: the selected experiment's name when `?beamtime=` is set, otherwise
    "the whole corpus". (No "N / M screened" progress bar — see Out of scope.)
  - **Column header** — Sample · Exposures · Kept · Tags · Status.
  - **Rows** — one `<ContactSheetRow>` per filtered sample, keyed by `sample.id`.
  - **Keyboard legend** — the static mockup legend strip.
- The rows region is wrapped in a boneyard `<Skeleton loading={query.isLoading}>`
  (acceptance: "a boneyard skeleton shows while the corpus query is loading").
  On `query.isError`, render an inline error state instead of an empty table.

### 2. `ContactSheetRow` — one sample row

Props: `{ sample: CorpusSample }`. The row owns its own exposure query so the
table renders incrementally as rows resolve.

- Calls `useExposures(sample.id)` — the existing hook, keyed
  `queryKeys.exposures(sampleId)`. This is the **same cache entry** the culling
  work (#162) and the loupe (#161) read, so no migration is needed when those
  land. Trade-off: one request per rendered row. Accepted — TanStack Query
  dedupes and caches; row virtualization is out of scope (a future concern if a
  corpus grows large).

The five cells:

| Cell | Source | Content |
|---|---|---|
| **Sample** | `CorpusSample` | `display_name ?? name ?? "#" + id`, with the numeric `id` as a secondary line. Identity only — no screened mark (see Out of scope). |
| **Exposures** | `useExposures` | A horizontal strip of `DetectorImage size="thumb"`, one per exposure, in `id` order. Rejected exposures (`status === "rejected"`) render greyed with an ✕; the representative (`selected === true`) renders marked. While `exposuresQuery.isLoading`, the cell shows a lightweight inline placeholder (not a boneyard skeleton — the page-level skeleton already covers initial load). |
| **Kept** | `useExposures` | `kept / total`, where `kept = exposures.filter(e => e.status !== "rejected").length` (every frame is kept until dropped — redesign-notes §6); a "N dropped" sub-label when `total - kept > 0`. Placeholder dash while loading. |
| **Tags** | `CorpusSample.tags` | The sample's `tags` as read-only free-form chips, plus an inert `+` button. The mutation round-trip is #159. |
| **Status** | — (seam) | A fixed "Not indexed" chip. Rendered through a typed `status` value with a `TODO` comment marking where the real phase call wires in — #160 declares no index-data dependency, and no issue is yet scoped for the phase wiring. |

The exposure strip is built directly from `DetectorImage size="thumb"` rather
than reusing `ThumbnailGallery` wholesale: `ThumbnailGallery`'s
`selectedId` / `onSelect` single-select viewer semantics belong to the loupe,
not the multi-select cull strip (#162). The thumbnail-select affordance and the
Tags `+` button render as **inert placeholders** — full visual structure, no
wiring — so #162 and #159 add behaviour with minimal rework.

### 3. `CorpusTopbar` — functional Beamtime chip

The spec mandates a working beamtime filter (master plan §4.1, §11; redesign-notes
line 236 — "a beamtime filter chip in the sample table's topbar"). #155 shipped
the chip as a placeholder; #160 makes it functional.

- The chip becomes a dropdown over `useExperiments()` plus an "All experiments"
  option.
- Selecting an experiment writes `?beamtime=<experiment_id>` to the URL via
  `setSearchParams`; "All experiments" clears the param.
- The chip label reflects the current selection (the experiment name, or
  "Beamtime" when unfiltered).
- The URL is the only channel between the topbar and `SamplesPage` — no prop or
  Zustand coupling. This resolves the issue's open question: `?beamtime=` URL
  state is owned by `/samples`; the chip is purely a writer of it.

### 4. Tests — `test/contact-sheet.test.tsx`

Vitest, following the existing query-hook test harness (`makeClient`, a
`mockOnce`-style fetch mock, `QueryClientProvider` wrapper). The fan-out means
the mock must answer both `GET /api/samples` and the per-sample
`GET /api/samples/{id}/exposures` calls.

- **Table renders the corpus** — given a mocked corpus, each sample produces a
  row with all five cells.
- **Kept count** — a sample with mixed exposure statuses renders the correct
  `kept / total` and "N dropped".
- **Exposure strip** — a row renders one thumbnail per exposure; a rejected
  exposure renders in the rejected state.
- **`?beamtime=` filter** — mounting at `/samples?beamtime=<id>` shows only that
  experiment's samples; the param round-trips (selecting in the chip updates the
  URL, and a URL with the param filters the table).
- **Skeleton** — while the corpus query is loading, the boneyard skeleton is
  shown and the table body is not.

## Out of scope

- **The screened mark and the "N / M screened" progress bar.** `screened` is a
  defined concept (redesign-notes line 150 — a per-row progress mark) but no
  issue is scoped to *derive* it, and #160's acceptance criteria do not list it.
  Screened state is the *output of triage* — and triage (culling, representative
  pick) is #162. The whole screened surface (row dot + header progress bar)
  belongs with #162, which produces the state it reflects. #160's Sample cell is
  identity only, matching the issue's wording.
- **Culling / batch reject / representative pick (#162).** The strip renders
  exposure state; selection and drop are inert here.
- **Tag mutation (#159).** The Tags column renders existing chips and an inert
  `+` button; the add/remove round-trip is #159.
- **The loupe view (#161)** — `/samples/loupe/:sampleId`.
- **The real phase call in the Status column** — needs per-sample index data;
  a fixed "Not indexed" placeholder ships behind a typed seam.
- **The Phase 1 cutover and Inspect deletion (#163 / I1.7).**
- **Row virtualization** — the per-row fan-out renders every filtered row; a
  future concern if a corpus grows large.

## Acceptance criteria

- [ ] `/samples` renders the whole corpus as a contact-sheet table.
- [ ] Each row shows identity, an exposure-thumbnail strip, a Kept count, Tags
      chips, and a Status.
- [ ] The Kept count is `kept / total` with non-rejected exposures counted as
      kept, and a "N dropped" sub-label when any are dropped.
- [ ] `?beamtime=<experiment_id>` filters the table to one experiment and
      round-trips through the URL; the topbar chip reads and writes it.
- [ ] A boneyard skeleton shows while the corpus query is loading.
- [ ] Thumbnail-select and the Tags `+` button render as inert affordances
      (no wiring).
- [ ] The Status column renders "Not indexed" behind a typed seam.
- [ ] Vitest covers the table, the Kept count, the exposure strip, the
      `?beamtime=` filter, and the skeleton.
- [ ] `npm run build` (tsc + vite) passes.

## References

- Issue #160 (I1.4) and issue map: `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign-issue-map.md`
- Master plan §4 / §4.2: `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign.md`
- Corpus query layer (#156): `docs/superpowers/specs/2026-05-17-corpus-sample-query-layer-design.md`
- Mockup: `docs/redesign-mockups/sample-table.html`
- Design record: `docs/redesign-notes.md` (§4 workflow, §5 architecture, line 150 screened mark)
- Reused components: `ThumbnailGallery.tsx`, `DetectorImage.tsx`, the boneyard `<Skeleton>` (`packages/HimalayaUI/docs/boneyard.md`)
- Existing hooks: `useCorpusSamples` / `useExposures` / `useExperiments` (`queries.ts`)
