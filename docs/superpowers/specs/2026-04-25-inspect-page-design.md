# Inspect Page Design

**Date:** 2026-04-25
**Status:** Ready for implementation

---

## Context

SAXS experiments typically produce multiple exposures per sample — different beam positions on the same physical specimen. Before indexing, the user needs to review all exposures and:

1. Reject ones that are obviously bad (flare, missed sample, beam artefacts)
2. Accept ones that look clean
3. Mark one for indexing (with automatic fallback to the first accepted)

Currently, HimalayaUI has no UI for this workflow. The `exposures.selected` boolean handles the "use for indexing" designation but there is no accept/reject status, no image display, and no sample metadata editor. Tag management UI was dropped from the Index page redesign and has no home. This page provides all of it.

---

## Tab Navigation

A third tab — **Inspect** — is added to the `TabRocker` before the existing **Index** tab. The `activePage` Zustand state gains `"inspect"` as a valid value. The page is always scoped to the current `activeSampleId`.

---

## Layout and Responsive Reflow

Three breakpoints using CSS grid. The large layout uses `grid-cols-[28fr_22fr_50fr]` — wider metadata and image columns relative to a compact gallery strip, rather than mirroring the Index page's symmetric `22fr 56fr 22fr`.

### Large (≥ 1400px) — three columns

```
┌───────────────┬──────────────┬───────────────────────────┐
│   Metadata    │   Gallery    │       Image View           │
│    (28fr)     │   (22fr,     │        (50fr,              │
│               │  2-col grid) │      full-height)          │
└───────────────┴──────────────┴───────────────────────────┘
```

Metadata is wider than the Index chat card to accommodate name/notes/tags comfortably. The gallery is a compact 2-column portrait grid. The image view takes the majority of the width to display the portrait detector image well.

### Medium (1100px – 1400px) — two columns, 60/40

```
┌──────────────────────────┬─────────────────┐
│ Metadata (top)           │                 │
├──────────────────────────┤   Image View    │
│ Thumbnail Gallery        │   (40%, full-   │
│ (auto-fill grid,         │    height)      │
│  scrollable)             │                 │
└──────────────────────────┴─────────────────┘
```

### Small (< 1100px) — single column, stacked

```
┌──────────────────────────────┐
│ Sample Metadata              │
├──────────────────────────────┤
│ Horizontal scrollable strip  │  ← thumbnails scroll left/right
├──────────────────────────────┤
│ Image View + controls        │
└──────────────────────────────┘
```

---

## Components

### `InspectPage.tsx`

Composition root. Reads `activeSampleId` from Zustand, calls `useExposures(sampleId)`. On mount, if `activeSampleId` is undefined, opens the `NavModal` (same behaviour as `IndexPage` — reuse the existing `openNavModal` Zustand action). Manages local `selectedExposureId` state (the exposure currently shown in the image view — distinct from the "use for indexing" flag). Default selection on mount: the indexing-marked exposure (`selected = true`) if one exists, otherwise the first non-rejected exposure, otherwise the first exposure. This state does not persist across tab switches. Renders the three cards in the appropriate grid layout. Mirrors the structure of `IndexPage.tsx`.

### `ThumbnailGallery.tsx`

Displays all exposures for the sample as a portrait-thumbnail grid. Each cell:
- Shows the 2D TIFF rendered via `<DetectorImage>` at thumbnail size
- Visual states:
  - **Blue ring** — currently selected in the image view
  - **Green ring** — marked for indexing (`selected = true`)
  - **Dimmed + dashed border** — rejected
  - Default — unreviewed
- Clicking a cell sets `selectedExposureId` locally (does not change the indexing selection)
- Label below each thumbnail: short filename + status icon

In the three-column layout the grid is two columns wide. In medium it is auto-fill. In small it becomes a horizontal scrolling row.

### `DetectorImageCard.tsx`

Right-column card (or bottom card on small). Displays the selected exposure at full card height.

- Header: filename
- `<DetectorImage>` component filling available height at portrait aspect ratio
- Controls at the bottom:
  - **Accept** — sets `status = "accepted"`; no-op if already accepted (active state shown)
  - **Reject** — sets `status = "rejected"`; opens an inline text field for a rejection note (stored as `exposure_tag` with `key = "rejection_reason"`, `source = "manual"`); already-rejected shows the note with an edit affordance
  - **Use for indexing** — calls `PATCH /api/exposures/:id/select`; active state shown when `selected = true`; disabled when exposure is rejected

Auto-fallback rule: if no exposure has `selected = true` for a sample, the CLI `analyze` command selects the first accepted exposure before running the analysis pipeline. This logic lives in `cli.jl`, not in `pipeline.jl`. No UI change needed.

### `DetectorImage.tsx`

Reusable component used by both the card and the gallery thumbnails.

- Fetches `GET /api/exposures/:id/image` (see Backend) which returns a normalized grayscale PNG
- Renders to a `<canvas>` element
- Applies a theme-aware colormap: reads `--color-bg` (maps to intensity 0) and `--color-fg` (maps to intensity 255) from computed CSS, builds a linear LUT, and paints each pixel
- Reapplies the LUT when the theme changes (listen for `data-theme` attribute mutation on `<html>`)
- Accepts a `size` prop (`"thumb" | "full"`); constructs the URL directly as `/api/exposures/${id}/image` (full) or `/api/exposures/${id}/image?thumb=1` (thumb) — no separate query hook needed
- Only renders when `image_path` is non-null (caller passes this as a prop from the exposures list); shows a placeholder otherwise

### `SampleMetadataCard.tsx`

Left-column card (top-left on medium, top on small).

Fields:
- **name** — editable inline text; `PATCH /api/samples/:id` on blur
- **notes** — editable multiline textarea; `PATCH /api/samples/:id` on blur
- **label** — read-only chip (the beamline identifier from the manifest)
- **Tags** — key/value tag list; manifest-sourced tags shown as read-only chips; manual tags have a delete button; `+ add tag` affordance creates a new `sample_tag` via `POST /api/samples/:id/tags`; deletion via `DELETE /api/samples/:id/tags/:tag_id`
- **Summary line** — `N exposures · N accepted · N rejected` computed from the exposures list; not a separate fetch

---

## Backend Changes

### DB migration — `exposures` table

```sql
ALTER TABLE exposures ADD COLUMN status TEXT
  CHECK (status IN ('accepted', 'rejected'));

ALTER TABLE exposures ADD COLUMN image_path TEXT;
```

`status` is nullable (null = unreviewed). `image_path` is nullable (null = no image found).

### Image discovery at ingestion (`pipeline.jl`)

During `init` and `analyze`, after resolving the `.dat` file for each exposure, search for a matching TIFF:

1. Strip the `.dat` extension from the filename, try appending `.tiff` and `.tif` in the same directory.
2. If found, store the absolute path in `image_path`.
3. If not found, `image_path` remains null and the image card shows a "no image" placeholder.

The search pattern is beamline-specific by convention; no config knob yet — this is a known limitation noted for future extension.

### `GET /api/exposures/:id/image`

- Reads `image_path` from DB; returns 404 if null
- Loads the TIFF using `FileIO` / `ImageIO`
- Applies log normalization: `I_norm = log1p.(I) / maximum(log1p.(I))`
- Encodes as 8-bit grayscale PNG
- If `?thumb=1` is present, resizes to fit within 128×128 before encoding (preserving aspect ratio)
- Sets `Cache-Control: max-age=3600` — TIFF files are immutable once ingested
- Content-Type: `image/png`

### `PATCH /api/exposures/:id`

New or extended route accepting `{ status: "accepted" | "rejected" | null }`. Validates the status value. Separate from the existing `PATCH /api/exposures/:id/select`.

### Index page filtering

`GET /api/samples/:id/exposures` gains an optional `?exclude_rejected=true` query param. When set, exposures with `status = "rejected"` are omitted. The Index page passes this param; the Inspect page does not.

---

## Frontend State and Queries

### Zustand (`state.ts`)

- `activePage` type extended: `"index" | "compare" | "inspect"`
- No other state additions — `selectedExposureId` within Inspect is local component state (it does not need to persist across tab switches)

### TanStack Query additions (`queries.ts`)

- `useExposures(sampleId, { excludeRejected?: boolean })` — passes `?exclude_rejected=true` when flag is set; existing Index page call gains the flag
- No new image query hook — `DetectorImage` constructs the URL directly from its props

### Mutations

- `useSetExposureStatus()` — `PATCH /api/exposures/:id` with `{ status }`; invalidates `queryKeys.exposures(sampleId)`
- Rejection note: written as an exposure tag via the existing `useAddExposureTag()` mutation pattern

---

## Index Page Change

`useExposures` call in `IndexPage` / `App.tsx` gains `excludeRejected: true`. No other Index page changes.

---

## Verification

1. **Unit tests (Vitest):** `SampleMetadataCard` renders name/notes/label/tags; inline edit calls mutation on blur; tag add/delete calls correct endpoints. `ThumbnailGallery` applies correct ring classes per status. `DetectorImage` canvas colorization applies LUT using mock CSS variables.
2. **Backend tests (Julia):** Image route returns 404 for null `image_path`; returns PNG for valid path; `?thumb=1` produces smaller image. `PATCH /api/exposures/:id` rejects invalid status values. Ingestion locates `.tiff`/`.tif` sibling files and stores absolute path.
3. **E2E (Playwright):** Navigate to Inspect tab; click a thumbnail → image card updates; click Reject → note field appears → note saved as tag; click Accept → status chip updates; Inspect tab shows rejected exposure dimmed; switch to Index tab → rejected exposure absent from selector.
4. **Manual:** Toggle theme; confirm detector image colormap updates immediately without reload.
