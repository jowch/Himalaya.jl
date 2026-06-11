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

2. **New event kind `edit_tag`** (non-view-producing, exactly like `add_tag`/`remove_tag` → a `kind == "edit_tag" && return nothing` no-op guard in the `apply_event!` dispatcher; the route `UPDATE`s the `sample_tags` base table directly, then logs to `user_actions` and broadcasts over SSE). Scaffold via `/new-event-kind edit_tag --no-view`. Payload `{tag_id, key, value, experiment_id}` (the `experiment_id`, as `add_tag`/`remove_tag` carry, lets `applyRemoteToCache` invalidate the right samples key). **No rebuild-from-log round-trip test** — `sample_tags` is not a materialized/rebuilt view, so there is nothing to assert parity against (the round-trip contract is for view-producing kinds only). Required test instead: an `apply_event!` test asserting the `user_actions` row + the broadcast fire and **no** view write — the coverage `add_tag` already has.

3. **Scoping batch becomes an upsert.** `POST /api/samples/tags/batch` today inserts N rows **unconditionally** — which is the root of the duplicate: the warm-path worksheet already *reads* each sample's existing value for the ordering key (`proposeOrdering.ts:48-52` does `s.sample.tags.find(t => t.key === orderingKey)` and uses `tag?.value`), so it pre-fills `dose=10` from the existing tag, then on commit inserts a *second* `source='scoping'` row with that same `10`. Change the batch to **upsert on `(sample_id, key)`**: if the sample already has a tag with that key, leave it untouched when the value is unchanged (the common "just use it" case → no write, no `edit_tag` event), or UPDATE the existing row in place (keeping its `id`, leaving its `source`) when the worksheet value was deliberately edited; INSERT `source='scoping'` only for samples that lacked the key. Re-scoping is then idempotent and never mints a twin.

4. **Extend `GET /api/sample-tags`** (lives in **`routes_picker.jl`**, not routes_samples.jl) to return a per-`(key,value)` **sample count**: rewrite its `SELECT DISTINCT key, value` to `SELECT key, value, COUNT(DISTINCT sample_id) AS count … GROUP BY key, value`. Give the experiment-scoped sibling `GET /api/experiments/{eid}/sample-tags` the same treatment for symmetry. The `count` field is **purely additive and consumed only by the edit modal's `TagSuggest`** — it feeds the suggestion ranking (key list ranked by summed sample count = samples carrying the key under the single-valued rule; value list ranked by per-value sample count). **`proposeOrdering` is deliberately left UNCHANGED** (it ignores `count`). The two rankings differ *on purpose, because the objectives differ*: scoping ranks keys by **number of distinct values** — a proxy for "a scientific parameter that varies across the experiment," the right signal for choosing an *ordering variable*; the edit modal ranks by **sample frequency** — the right signal for consistency and duplicate-prevention. Record this so a future reader doesn't "fix" `proposeOrdering` to sample-count and break the ordering proposal.

### Frontend (React / TS / TanStack Query / Zustand / the mutation queue)

- **`TagSuggest`** (new `src/print/ui` primitive — genuinely new; the existing `TagEditor` only does a shallow `knownKeys` chip row, no listbox / counts / value autocomplete). The shared corpus-aware suggestion combobox, used for **both** the key field and the value field. WAI-ARIA combobox + listbox (`role=combobox`, `aria-expanded`, `aria-activedescendant`; `role=listbox`/`option`). Build its dropdown panel + option styling on the existing `Menu`/`Popover` primitives' visual vocabulary (closed-look — don't re-author panel appearance). Props: a mode (`key` | `value`), the suggestion list (`{label, count}[]`), the current text, `onCommit`. Renders suggestions with counts, de-emphasizes low-count entries, and a **create-as-typed** option. Key mode lists corpus keys ranked by sample count; value mode lists values **scoped to the active key**. Counts come from the counted `/api/sample-tags` list only — tolerate an absent/0 count when fed pairs synthesized elsewhere. **Deferred:** the fancy "did-you-mean" near-match nudge (edit-distance / unit-aware matching) is its own sub-feature and is **not** in Layer 1 — the usage counts already make a stray value legible (`10.0 (2)` next to `10 (44)`), which is the anti-fragmentation value; the heuristic nudge is a clean follow-up once the counts are in.

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
- **Duplicate-key rejection (single-valued-key rule).** Adding a key the sample already has, or editing a key so it collides with another row on the same sample → inline non-destructive rejection ("This sample already has a `dose` tag — edit that one instead"), `aria-invalid` + alert, cleared on edit (mirrors `TagEditor`'s empty-key rejection). No silent merge in Layer 1. The one pre-existing live duplicate (sample 3) is shown as two rows so the user can delete the redundant one — the rejection only prevents *new* collisions.

## Accessibility (the explicit bar)

- `ManageTagsModal`: built on `ModalShell` → it already provides `role=dialog`, `aria-modal`, focus trap (`useFocusTrap`), and Escape-closes; `aria-labelledby` → the `h2`. The house pattern is `aria-modal` + focus-trap (not `inert`) — and `suppressGlobalKeys` keys off the open `aria-modal` dialog, so the loupe's window X/R/arrow handlers are suppressed for free while the modal is open. ModalShell does **not** restore focus to the opener, so the modal captures the "Manage" trigger ref and restores focus to it on close itself.
- `TagSuggest`: full WAI-ARIA combobox (input `role=combobox` + `aria-expanded`/`aria-controls`/`aria-activedescendant`; list `role=listbox`, items `role=option`/`aria-selected`); Down/Up navigate, Enter commits, Escape collapses without clearing; visible label, not placeholder-only; "create" option visually distinct.
- Live region: announce add/edit/remove ("Tag dose=10 updated to dose=12") via a polite `role=status`; imperatively set text on removal (don't rely on `aria-relevant`).
- Delete buttons: real `<button>`, `aria-label="Remove dose=10"` (names the tag), ≥24×24 CSS px target (WCAG 2.5.8).
- Contrast/focus per the house token system + base focus-visible rule (F-STATE).

## Testing (six-layer contract rule)

- **Backend:** route handler test (`PATCH` happy + 400 + scoped-by-sample); `edit_tag` `apply_event!` + **rebuild-from-log round-trip** parity; scoping-batch upsert test (existing key updates, new key inserts, no duplicate); `/api/sample-tags` counted-shape test.
- **Frontend unit (Vitest):** `TagSuggest` (key/value modes, counts, create-as-typed, did-you-mean, ARIA semantics); `ManageTagsModal` (rows render with ids, edit/add/delete wiring, exact-dup rejection, flat-modal/focus); adapters carry id/source; `useEditSampleTag` queue contract (optimistic/own-op/foreign-replay) paired with the backend per the contract-testing doc.
- **E2E (Playwright, mocked):** open Manage from the loupe, edit a value, delete the redundant duplicate row; assert rendered semantics (data-*/text/aria), never class strings.
- **Live render-verify** at :5182 against the real corpus (the duplicate `dose=10` sample) before the work is marked done.

## Resolved decisions (ratified with Jonathan 2026-06-11)

1. **A key is single-valued per sample — enforced, no multi-value.** A sample never has two tags with the same key. Enforcement lives in the write paths: the modal's add-row and the inline `+ tag` quick-add both **reject a key already present on the sample** (and offer to edit the existing one instead); a key-edit that would collide with another row on the same sample is rejected inline (no silent merge in Layer 1). The backend `POST .../tags` and `PATCH .../tags/:id` return a 409/400 as a backstop. No DB `UNIQUE(sample_id, key)` constraint is added (the one pre-existing duplicate in live data — sample 3 — would violate it; that row is reconciled by the user in the modal, not by a migration).
2. **Scoping "just uses" the existing tag.** Resolved by the verification above: scoping already reads the existing value (`proposeOrdering.ts:51`); the bug is the unconditional insert on commit. The upsert (change #3) leaves an existing tag untouched when the value is unchanged and only updates it (in place, same `id`, same `source`) on a deliberate worksheet edit. There is no "overwrite a manual tag" surprise, because the worksheet is pre-filled from that very tag.

## Build order & risk (from spec review 2026-06-11)

Five slices, smallest-correct-first; the upsert (the risk locus) goes last behind the most tests:

1. **id-exact inline delete + adapter carries the id** — the literal LO-TAGDUP correctness fix, and it *stands alone*: it ships the P1 close on its own. Mechanically it IS the adapter change — `toLoupeTags` must emit a loupe-local view-model that carries `id` (+`source`), and `onRemoveTag` receives that `id` instead of re-deriving via `findSampleTagId` (which becomes dead and is deleted). Once this lands, sample 3's two `dose=10` pills are individually removable.
2. **Backend `PATCH .../tags/:id` + `edit_tag` event** (+ the apply_event! test; in-txn validation; 409 on a duplicate-key collision via `WHERE sample_id=? AND key=<new> AND id<>?`; `changes()`-guarded 404 on a non-matching `(tag_id, sample_id)`).
3. **Counted `GET /api/sample-tags`** (both the corpus + experiment routes in `routes_picker.jl`). `proposeOrdering` is **left unchanged** — the `count` field is consumed only by the modal's ranking; the two rankings differ by objective (see backend change #4).
4. **`TagSuggest` + `ManageTagsModal`** + the loupe "Manage" affordance + the inline `+ tag` duplicate-key rejection.
5. **Scoping batch upsert** — the **highest-risk** step: it changes a path with a subtle existing idempotency invariant (one `client_op_id` over N rows, the `idx_events_unique_op` uniqueness, the duplicate-`sample_id`-in-batch guard which must stay) and now emits a *heterogeneous* mix per batch (`add_tag` for a new key, `edit_tag` for a changed value, **zero events** for an unchanged value — a behavior change worth its own assertion). The upsert must `SELECT id, value FROM sample_tags WHERE sample_id=? AND key=?` in-txn to know the current value (the batch payload carries none), with a deterministic pick if >1 row exists; and it must use the upsert path, **not** the duplicate-key-rejecting insert path (so it never 409s against the tag it's updating).

## Implementation-plan notes (queue layer — six-layer trap)

`edit_tag` is not two files — the plan must touch the full queue chain (per `lib/queue/AGENTS.md`):
- `lib/queue/types.ts` — add `edit_tag` to the `OpKind` union.
- `lib/queue/mutators/trivial.ts` — new `editSampleTagMutator`: optimistic patch is an **in-place, id-stable map** over the `queryKeys.corpusSamples` cache (the cache the loupe reads), keyed on `tag_id` (NOT `(key,value)` — the add mutator's `(key,value)` placeholder-match is the very thing this feature kills); plus its `synthesizeFromSse` reading `{tag_id, key, value}`.
- `lib/queue/mutatorRegistry.ts` — `edit_tag` arms in **both** `resolveMutator` and `resolveMutatorForEvent` (else own-op confirmation and foreign replay-as-rerun don't dispatch).
- `lib/queue/applyRemoteToCache.ts` — `case "edit_tag"` in the same **invalidation** arm as `add_tag`/`remove_tag` (the tag foreign-event path is pure `invalidateQueries`, not a row-merge): invalidate `corpusSamples` + `corpusSampleTags` + `corpusPickerSamples` (an edited value changes the corpus vocabulary the suggestion dropdown reads).
- `useEditSampleTag` scope carries `sampleId` only (no `experimentId`) so it routes through the corpus mutator.

## Components & files (anticipated)

- Backend: `routes_samples.jl` (`PATCH .../tags/:id`, scoping-batch upsert), `events.jl` (`edit_tag` no-op guard), `routes_picker.jl` (counted `/api/sample-tags` + experiment sibling), backend tests.
- Frontend new: `src/print/ui/TagSuggest.tsx` (+ stories/test), `src/print/components/ManageTagsModal.tsx` (+ test).
- Frontend touched: `src/api.ts` (`editSampleTag`, counted corpus tags), `src/queries.ts` (`useEditSampleTag`, counted `useCorpusSampleTags`), the four `lib/queue/*` files above, `src/print/ui/TagEditor.tsx` + `TagList.tsx` (thread the sample's existing keys for inline duplicate-key rejection — a new prop), `src/print/pages/loupeAdapters.ts` (loupe-local `LoupeTag` view-model carrying id/source — do **not** widen `ui/tag.ts`; delete `findSampleTagId`), `src/print/pages/LoupePage.tsx` + `src/print/components/LoupeSidePanel.tsx` (Manage affordance + modal, id-exact inline delete, reuse `lib/announce`).
- Closes: **LO-TAGDUP** (P1). Slice 1 is the literal correctness fix; slices 2-5 are the full "editable tags" resolution.
