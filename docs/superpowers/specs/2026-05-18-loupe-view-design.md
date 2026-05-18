# Loupe view — design (#161 / I1.5)

**Status:** v1, 2026-05-18 — brainstormed and approved.
**Issue:** #161 — Build the sample loupe view (I1.5).
**Parent plan:** [`2026-05-17-himalaya-ui-redesign.md`](../plans/2026-05-17-himalaya-ui-redesign.md) §4.2 · issue map [`…-issue-map.md`](../plans/2026-05-17-himalaya-ui-redesign-issue-map.md) I1.5.
**Mockup:** `docs/redesign-mockups/sample-table.html` — the `loupe-shell` markup.

---

## 1. What this is

A focused single-sample inspection surface at `/samples/loupe/:sampleId` — the close-up to
the contact sheet's corpus overview. A centered plate: a full detector image, a centered
exposure-thumbnail strip, and a metadata sidebar (per-exposure facts, a Verdict box, a
Representative box, a Sample-tags section, a keyboard legend).

It is the new home for what the Inspect page rendered through `DetectorImageCard` and
`SampleMetadataCard`. The loupe **rebuilds** that content into the mockup's layout reusing
the leaf primitives `DetectorImage` and `ThumbnailGallery` — it does **not** re-compose the
two Inspect cards as-is (their internal layout fights the mockup). The Inspect page and its
cards stay untouched and functional until #163 retires Inspect wholesale.

## 2. Dependencies and scope

**Depends on (both merged):** the corpus app shell (#155 / I1.1) and the corpus query
layer (#156 / I1.2). **Blocks:** the Phase 1 cutover (#163).

**In scope:**

- New loupe view at `/samples/loupe/:sampleId`, mounted under the `CorpusShell` route group.
- Full detector image (`DetectorImage size="full"`) + centered exposure-thumbnail strip.
- Metadata sidebar: "This exposure" meta-list, Verdict box, Representative box, a
  read-only Sample-tags section, a keyboard legend.
- Frame flipping between exposures (strip click + arrow keys).
- Per-exposure drop/restore and set-representative, round-tripping through the queue via
  the existing exposure hooks.
- A boneyard skeleton while loading.
- Vitest coverage for the loupe.

**Out of scope:**

- The contact-sheet table (#160) — nothing links *into* the loupe yet; direct-nav /
  deep-link only until #160 lands.
- The corpus sample-tag add/remove round-trip (#159) — the loupe's Sample-tags section is
  **read-only display** in this issue; #159 wires editing.
- Inspect deletion (#163).
- Sample `display_name` / `notes` editing — moves to the focus workspace (#179 / I4.2, a
  Phase 4 issue, not a Phase 1 sibling); the loupe is read-only for sample-level identity.
- "Index this exposure" navigation — a focus-workspace concern (#179 / I4.2, Phase 4).
- A reject-reason picker — the loupe uses the mockup's plain drop/restore toggle (§6).

## 3. Why "Full" interaction scope

The mockup's loupe carries a Verdict box (drop/restore an exposure) and a Representative
box (pick the indexing frame). The I1.5 acceptance criteria explicitly require only
frame-flipping; I1.6 (#162) owns culling but its card and acceptance criteria are entirely
**contact-sheet** scoped ("reject-only culling on contact-sheet rows", "the multi-select
lives in the table"). The loupe's *own* per-exposure controls are therefore claimed by no
other issue. Building them here completes the mockup and is the natural reading of
"absorbs `DetectorImageCard` internals" (that card already owns status mutation).

Representative-pick (`useSelectExposure`) does appear in master plan §4.2's culling
bullet, which I1.6 implements — but I1.6 wires it into *contact-sheet rows* and the loupe
wires it into the *loupe sidebar*. These are per-surface controls over the same hook by
design, not a contradiction: each surface owns its own affordance, neither owns the
other's.

The exposure hooks (`useSetExposureStatus`, `useSelectExposure`) key their cache writes on
`sampleId` — the same key `useExposures(sampleId)` reads — so drop/restore and
set-representative round-trip correctly in the loupe with **no dependency on the #159
corpus-tag cache work**.

## 4. Architecture — files

Three new files; one existing file modified.

| File | Responsibility |
|---|---|
| `src/pages/LoupePage.tsx` (new) | Route entry. Reads `:sampleId` via `useParams`; fetches data; owns active-exposure state, keyboard handling, and back-nav; renders loading / not-found; lays out the plate. |
| `src/components/LoupeFrame.tsx` (new) | The big `DetectorImage size="full"` with a "Dropped" overlay when the exposure is rejected, plus the centered exposure-thumbnail strip (reusing `ThumbnailGallery`). |
| `src/components/LoupeSidebar.tsx` (new) | The sidebar: "This exposure" meta-list, Verdict box, Representative box, read-only Sample-tags section, keyboard legend. |
| `src/components/AppRoutes.tsx` (modify) | Add `<Route path="/samples/loupe/:sampleId" element={<LoupePage />} />` inside the existing `<CorpusShell>` route group. This is the hoisted route table I1.1 (#155) created for later issues to register slots into — an append-only one-line add, not flagged in issue-map §3 shared-file contention. |
| `src/bones/loupe.bones.json` (new) | Boneyard skeleton fixture for the loupe plate (captured per the boneyard workflow). |

The loupe is **URL-owned**: `sampleId` comes from `useParams`, never from the Zustand
`activeSampleId` (master plan §2.3 — new surfaces own their URL via `useParams`/`useNavigate`).

## 5. Data flow

- `:sampleId` from `useParams<{ sampleId: string }>()`, parsed to a number.
- `useCorpusSamples()` → find the sample by `id`. A `CorpusSample` supplies `name`,
  `display_name`, `experiment_id`, `tags`, and `q_units`.
- `useExposures(sampleId)` → the exposure list for the sample.
- `useExperiment(sample.experiment_id)` → the experiment name for the head's sub-line.
- **Active exposure** — local `useState`, default-selected by the InspectPage rule:
  (1) the exposure with `selected: true`, else (2) the first `status === "accepted"`, else
  (3) the first exposure. The selection **resets when `sampleId` changes**.
- **Loading** — a boneyard skeleton (`loupe.bones.json`, gated on `query.isLoading`) shows
  while either the corpus-samples or the exposures query is loading.
- **Not found** — if the parsed `sampleId` is absent from the corpus list, render a small
  "Sample not found · back to the sheet" panel with a link to `/samples` (not a skeleton,
  not a crash).

## 6. Surface and interactions

Layout follows the `loupe-shell` markup in `sample-table.html`: a centered plate (max
~1100px) — back button → head (sample name + a mono sub-line with experiment name) →
two-column body.

**Left column** — `LoupeFrame`:

- The big square detector frame: `DetectorImage size="full"` for the active exposure. A
  "Dropped" badge overlays the frame when the active exposure's `status === "rejected"`.
- The centered exposure-thumbnail strip — one thumbnail per exposure (`ThumbnailGallery`),
  the active one highlighted.

**Right column** — `LoupeSidebar`:

- **"This exposure" meta-list** — key/value rows for the active exposure: Filename, Kind,
  Frame (*i of N*), Status. (`q_units` rides on the corpus sample for Phase 3
  normalization; it is not surfaced in the I1.5 meta-list.)
- **Verdict box** — current keep/drop state + a toggle. Drop/restore is a plain toggle
  (`status` `'rejected'` ⇄ `null`) via `useSetExposureStatus` — **no reject-reason
  picker** (the mockup's simplification; `DetectorImageCard`'s reason picker is not
  carried over). If a rejected exposure already carries `rejection_reason` tags, they are
  displayed; the loupe does not author new ones.
- **Representative box** — marks the active exposure as the indexing representative via
  `useSelectExposure`; reflects whether the active exposure is already representative.
- **Sample-tags section** — the sample's `tags` rendered as **read-only** chips. The
  add/remove round-trip is #159's deliverable; this section is the placeholder it wires.
- **Keyboard legend** — documents the shortcuts below.

**Interactions:**

| Input | Action |
|---|---|
| Strip thumbnail click | Set the active exposure. |
| **← / →** | Flip to the previous / next exposure. |
| **X** | Drop / restore the active exposure (`useSetExposureStatus`). |
| **R** | Set the active exposure as representative (`useSelectExposure`). |
| **Esc**, back button | `navigate('/samples')`. |

The keyboard handler ignores key events while a text input or textarea is focused (no
text inputs exist in the loupe in I1.5, but #159 will add tag inputs — the guard is in
from the start). Mutations round-trip through the mutation queue via the existing hooks;
optimistic updates land on the `["exposures", sampleId]` cache the loupe reads.

## 7. Not file-per-exposure

The deferred derived-exposure feature must remain possible, so the loupe must not assume
one file per exposure. Concretely: the loupe iterates `exposures` by `id` and renders
every exposure regardless of `kind` (`"file"` / `"averaged"` / `"background_subtracted"`);
nothing keys off `filename` being non-null (an `averaged` exposure may have
`filename: null`). The strip length, frame-flipping, and the meta-list all derive from the
exposure list, not from files.

## 8. Skeleton loading

A new `src/bones/loupe.bones.json` fixture captures the loupe plate skeleton, registered
through the boneyard registry and rendered via `<Skeleton>` gated on the load state
(`corpusSamplesQ.isLoading || exposuresQ.isLoading`) — the same pattern InspectPage uses
for its cards. The exact capture follows the boneyard workflow
(`packages/HimalayaUI/docs/boneyard.md`); the detailed plan pins the capture step.

## 9. Testing (Vitest)

`LoupePage` tests:

- Renders the sample identified by the route param.
- Shows the boneyard skeleton while loading.
- Renders the "Sample not found" fallback for an id absent from the corpus.
- Frame-flip works via both strip click and the arrow keys.
- **X** calls `useSetExposureStatus` (drop and restore); **R** calls `useSelectExposure`.
- The Sample-tags section renders the sample's tags as chips.
- A mixed-`kind` exposure list (including an exposure with `filename: null`) renders every
  exposure in the strip and flips through all of them — the file-per-exposure regression
  guard.

Boneyard skeletons gate on `query.isLoading`. The Phase 1 Playwright mocked spec
(loupe-flip and tag) lands with the Phase 1 cutover (#163 / master plan §4.3), not here.

## 10. Open questions

None. Two design points settled during brainstorming: (a) the reject-reason picker is
dropped in favor of the mockup's plain drop/restore toggle; (b) sample `display_name` /
`notes` editing stays out of the loupe — it lands in the focus workspace (#179).
