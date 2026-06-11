# Editable sample tags — design spec (Layer 1)

**Date:** 2026-06-11
**Status:** approved direction (brainstormed with Jonathan), pending spec review → implementation plan.
**Origin:** Loupe re-score #3 minted **LO-TAGDUP** (P1) — a sample carrying both a manual and a series-scoping tag with the same key+value (`dose=10` manual id 3 + `dose=10` scoping id 7) renders two byte-identical removable pills, and the loupe's `findSampleTagId` matches key+value only, so `×` always deletes the first id and the scoping tag is unremovable.

## The reframing

The duplicate/conflict is a *symptom*. The disease: **sample tags are insert/delete-only — never editable.** Because a human can't fix a value, scoping has no choice but to insert (creating redundancy), and a wrong value can only be deleted-and-retyped. Every tag value is human-derived (typed directly, or assigned through the scoping worksheet), so the fix is to **put editing back in the human's hands**: make tags editable, let scoping update an existing tag instead of inserting a twin, and the duplicate/conflict tangle dissolves.

This spec is **Layer 1**: per-sample tag management. A **Layer 2** corpus-wide vocabulary manager (rename-is-merge across all samples, global duplicate hunting) is explicitly **deferred** to its own spec — Layer 1's corpus-aware suggestion counts already give the user duplicate-spotting at the point of edit, which is the bulk of the value.

## Goals

- A human can **edit** a sample tag's key and value in place (not delete-then-retype).
- Editing surfaces the **existing corpus vocabulary with usage counts**, so a stray value (`dose=10.0`, 2 samples) reads as obviously wrong next to the canonical one (`dose=10`, 44 samples). This is the anti-fragmentation lever.
- The duplicate manual/scoping `dose=10` becomes two **visible, individually-removable** rows — no provenance gymnastics, because everything is now human-editable.
- Scoping **upserts** rather than inserting duplicates.
- Inline quick-add and inline delete (`×`) on the pills **stay** (polished); a dedicated **modal** is the editing home.

## Non-goals (deferred / out of scope)

- **Layer 2** corpus vocabulary manager (global rename-is-merge, cross-sample duplicate hunting as a dedicated page).
- **Exposure** tags (same insert/delete-only model, not the pain point; left untouched).
- Re-homing scoping tags out of `sample_tags` into series-member attributes (a larger re-architecture; noted, not done).
- Undo/redo and versioning of tag edits beyond the immediate optimistic-queue rollback (Layer-4 territory).

## Research grounding (full briefs in session)

Four parallel research tracks (productivity tools, scientific/media tools, anti-patterns + a11y, internal references) converged on:

1. **Show the vocabulary as data with usage counts** (Lightroom keyword list, Linear label table, AEM taxonomy). Counts make a stray value legible without any algorithm.
2. **Autocomplete key and value independently, from the corpus** (Obsidian, Apple Photos, Mendeley). The single biggest lever against fragmentation. The value list is **scoped to the active key**.
3. **Rename-is-merge** (Zotero). Not built in Layer 1 (that's Layer 2's global form), but the per-sample edit + create-as-typed is its local seed.
4. **Anti-patterns to avoid:** the cramped `key|value|🗑` spreadsheet, hard-delete with no undo, modal-in-modal. **A11y bar:** WAI-ARIA combobox + dialog patterns, live-region announcements, 24px targets.

## Architecture

### Backend (Julia / Oxygen.jl / SQLite)

Three cohesive changes (all tag-vocabulary), each behind the existing `with_idempotency` + `apply_event!` + SSE-broadcast machinery.

1. **New route `PATCH /api/samples/{id}/tags/{tag_id}`** — updates `key` and/or `value` of an existing `sample_tags` row (scoped by `sample_id` for safety). Validates shape before the transaction opens (clean 400, never an unguarded 500), mirrors the existing `POST .../tags` / `DELETE .../tags/{tag_id}` handlers in `routes_samples.jl`. Returns the updated `SampleTag`. Editing keeps the row's `id` stable (so optimistic reconciliation and provenance survive).

2. **New event kind `edit_tag`** (non-view-producing, like `add_tag`/`remove_tag` → `apply_event!` returns `nothing`; it logs to `user_actions` and broadcasts over SSE). Scaffold via the `new-event-kind` skill (`/new-event-kind edit_tag --no-view`). Payload carries `{tag_id, key, value, prev_key, prev_value}` for the audit trail. **Round-trip test** required (rebuild-from-log parity) per the event-log contract.

3. **Scoping batch becomes an upsert.** `POST /api/samples/tags/batch` today inserts N rows. Change it to **upsert on `(sample_id, key)`**: if the sample already has a tag with that key, UPDATE its value (and stamp `source='scoping'`); else INSERT. This means re-scoping never duplicates, and scoping a sample that already has a manual `dose=10` updates that one tag instead of minting a twin. Emits `edit_tag` for updates, `add_tag` for inserts. **This adopts a "key is single-valued per sample" model** for the scoping path — see Open Decisions.

4. **Extend `GET /api/sample-tags`** to return a per-`(key,value)` **sample count**: `{key, value, count}` instead of bare `{key, value}`. The frontend derives per-key counts by summing. This feeds the suggestion dropdowns' counts. (Additive field; existing consumers ignore it.)

### Frontend (React / TS / TanStack Query / Zustand / the mutation queue)

- **`TagSuggest`** (new `src/print/ui` primitive) — the shared corpus-aware suggestion combobox, used for **both** the key field and the value field. WAI-ARIA combobox + listbox (`role=combobox`, `aria-expanded`, `aria-activedescendant`; `role=listbox`/`option`). Props: a mode (`key` | `value`), the suggestion list (`{label, count}[]`), the current text, and `onCommit`. Renders suggestions with counts, de-emphasizes low-count entries, shows a **create-as-typed** option, and a **did-you-mean** nudge when the typed text is a near-match (edit-distance ≤ 2 or same digits/different unit) to an existing higher-count entry. Key mode lists corpus keys ranked by sample count; value mode lists values **scoped to the active key**.

- **`ManageTagsModal`** (new `src/print/components`) — built on the **custom-index-modal shell** (`ModalShell size="md"` + `ModalHead` Kicker `"Sample tags"` / serif `h2` = sample display name + `ModalFooter`). Body = a list of editable tag rows, each `[key TagSuggest] [value TagSuggest] [delete IconButton]`, plus an inline **add-a-tag** row at the bottom. Footer = a quiet note ("Tags apply to the whole sample. Saved as you go.") + a `Done` button. Flat modal — **no nested dialogs**; delete confirms/nudges are inline expansions. Initial focus to the first row's value (or the add row when empty); focus-trap + Escape via `ModalShell`'s `useFocusTrap`.

- **Data + writes:**
  - `useCorpusSampleTags()` (existing) extended to the counted shape; drives both `TagSuggest` lists.
  - **`useEditSampleTag(sampleId)`** — new `useQueueMutation` over `api.editSampleTag(id, tagId, {key, value})`, `client_op_id`-keyed, optimistic cache patch, SSE own-op confirmation, foreign `edit_tag` replay-as-rerun. Wire `edit_tag` into `applyRemoteToCache` (mirror the `add_tag`/`remove_tag` branches).
  - Existing `useAddCorpusSampleTag` / `useRemoveCorpusSampleTag` reused for inline add/delete.

- **Boundary adapters carry the id.** `loupeAdapters.toLoupeTags` currently drops `id` + `source`; the modal (and the inline delete) need the exact `id`. Carry `id` (and `source` for the faint provenance marker) onto the view-model the modal consumes so removal/edit target the exact record. The inline pill delete is also fixed to be **id-exact** as part of this (the literal LO-TAGDUP correctness fix), independent of the modal.

- **Loupe integration.** The loupe's "Sample tags" section keeps its pills with inline quick-add (`+ tag`) and inline delete (`×`); a **"Manage"** affordance (next to the section kicker, mirroring the scoping/custom-index header pattern) opens `ManageTagsModal`. The pills are the at-a-glance view; the modal is the editing home.

## Data flow

1. Loupe loads `sample.tags: SampleTag[]` (each with `id`, `key`, `value`, `source`).
2. Pills render from the tags (inline add/delete via the existing queue hooks, now id-exact).
3. "Manage" opens `ManageTagsModal` with the same tags + `useCorpusSampleTags()` counts.
4. Editing a value → `TagSuggest` (value mode, scoped to the row's key) → `onCommit` → `useEditSampleTag.mutate({tagId, key, value})` → optimistic patch → `PATCH` → `edit_tag` event → SSE confirms own-op / broadcasts to peers.
5. Deleting a row → `useRemoveCorpusSampleTag.mutate(tagId)`.
6. Adding a row → `useAddCorpusSampleTag.mutate({key, value})`.
7. Scoping commit (elsewhere) now upserts, so it never creates a duplicate of an existing key.

## Error handling

- Writes go through `useQueueMutation`: 4xx → assertive validation toast (`buildValidationMessage`); 5xx/network → auto-retry + `InfrastructureBanner`; optimistic rollback on terminal failure (per the established Honest-Surface-State machinery).
- Exact `(key,value)` duplicate add/edit → inline non-destructive rejection in the modal (don't mint a second identical row), with an `aria-invalid` + alert, cleared on edit (mirrors `TagEditor`'s empty-key rejection).
- Edit that collides with another existing tag on the same sample (same resulting `key,value`) → inline "this sample already has dose=10" nudge; no silent merge in Layer 1.

## Accessibility (the explicit bar)

- `ManageTagsModal`: `role=dialog`, `aria-modal`, `aria-labelledby` → the `h2`; background `inert`; focus trap; Escape closes; focus returns to the "Manage" trigger.
- `TagSuggest`: full WAI-ARIA combobox (input `role=combobox` + `aria-expanded`/`aria-controls`/`aria-activedescendant`; list `role=listbox`, items `role=option`/`aria-selected`); Down/Up navigate, Enter commits, Escape collapses without clearing; visible label, not placeholder-only; "create" option visually distinct.
- Live region: announce add/edit/remove ("Tag dose=10 updated to dose=12") via a polite `role=status`; imperatively set text on removal (don't rely on `aria-relevant`).
- Delete buttons: real `<button>`, `aria-label="Remove dose=10"` (names the tag), ≥24×24 CSS px target (WCAG 2.5.8).
- Contrast/focus per the house token system + base focus-visible rule (F-STATE).

## Testing (six-layer contract rule)

- **Backend:** route handler test (`PATCH` happy + 400 + scoped-by-sample); `edit_tag` `apply_event!` + **rebuild-from-log round-trip** parity; scoping-batch upsert test (existing key updates, new key inserts, no duplicate); `/api/sample-tags` counted-shape test.
- **Frontend unit (Vitest):** `TagSuggest` (key/value modes, counts, create-as-typed, did-you-mean, ARIA semantics); `ManageTagsModal` (rows render with ids, edit/add/delete wiring, exact-dup rejection, flat-modal/focus); adapters carry id/source; `useEditSampleTag` queue contract (optimistic/own-op/foreign-replay) paired with the backend per the contract-testing doc.
- **E2E (Playwright, mocked):** open Manage from the loupe, edit a value, delete the redundant duplicate row; assert rendered semantics (data-*/text/aria), never class strings.
- **Live render-verify** at :5182 against the real corpus (the duplicate `dose=10` sample) before the work is marked done.

## Open decisions (ratify at spec review)

1. **Single-valued-key model for scoping upsert.** The scoping upsert keys on `(sample_id, key)`, which treats a key as single-valued per sample (correct for ordering variables like `dose`/`temperature`). The **modal** itself does *not* enforce uniqueness — it allows repeated keys generally (multi-value), but rejects an exact `(key,value)` duplicate. Confirm this split is right, or decide keys must be globally unique per sample (which would make the modal enforce it too).
2. **Does scoping overwrite a *manual* same-key tag?** The upsert as specified updates *any* same-key tag (manual or scoping) to the scoped value and stamps `source='scoping'`. Alternative: scoping only upserts within `source='scoping'` and leaves manual tags alone (so a first-time scope of a manually-tagged sample still creates one scoping twin, reconciled in the modal). The spec adopts the former (no twin ever); confirm.

## Components & files (anticipated)

- Backend: `packages/HimalayaUI/src/routes_samples.jl` (+ `PATCH`, upsert), `events.jl` (`edit_tag`), a round-trip test, `/api/sample-tags` count.
- Frontend new: `src/print/ui/TagSuggest.tsx` (+ stories/test), `src/print/components/ManageTagsModal.tsx` (+ test).
- Frontend touched: `src/api.ts` (`editSampleTag`, counted corpus tags), `src/queries.ts` (`useEditSampleTag`, counted `useCorpusSampleTags`), `src/lib/queue/applyRemoteToCache.ts` (`edit_tag` branch), `src/print/pages/loupeAdapters.ts` (carry id/source), `src/print/pages/LoupePage.tsx` + `src/print/components/LoupeSidePanel.tsx` (Manage affordance + modal, id-exact inline delete).
- Closes: **LO-TAGDUP** (P1). The inline id-exact delete is the literal correctness fix; the modal + edit path is the full "editable tags" resolution.
