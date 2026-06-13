# Holistic interaction + visual audit — 2026-06-13

A full "what a user sees / touches" pass across all six greenfield surfaces, for
correctness and polish. Method: live walk of every surface against the writable
dev harness (:5182) via Playwright + two code-level deep-dive subagents (tag
model, keyboard/nav). Jonathan seeded two examples (tag-editor model; Focus-vs-
Loupe sample-switch + `[`/`]`) and asked to be holistic, not fixated on those.

Status of this doc: **findings only** (no fixes yet). Prioritized for a fix loop.
Each finding is code-verified; over-claims caught during verification are noted.

---

## Theme A — the tag-editing model is incoherent across THREE surfaces (Jonathan's seed)

There are three tag-entry behaviors — loupe inline (`TagEditor` via `TagList`),
the `ManageTagsModal`, and the samples sheet (`SampleTableRow` → `TagList`) — and
they diverge. Confirmed bugs + model gaps:

| # | Sev | Finding | Where |
|---|-----|---------|-------|
| A1 | P1 | Two divergent add paths: inline commits on **Enter** + rejects empty key loudly; the modal add-row **ignores Enter** (key/value `onCommit` only set state, never `handleAdd`) and **silently** no-ops an empty key. | `TagEditor.tsx:85-88` vs `ManageTagsModal.tsx:148,259-261` |
| A2 | P1 | Editing an existing tag is **modal-only**. Inline you can only add/remove, so changing `dose=10`→`12` means delete+re-add, which churns the tag id and breaks provenance. | `LoupeSidePanel.tsx:122-129` |
| A3 | P1 | **Samples sheet renders `editable` TagList with NO `onAdd`/`onRemove` wired** — a lying affordance: no `+ tag`, pills have no ×, yet the roving-grid still lets you "Enter the cell" to edit nothing. | `SampleTableRow.tsx:379`, `SamplesPage.tsx` (no tag hooks) |
| A4 | P2 | `knownKeys` (inline key suggestions) is a **dead prop** — no consumer passes it; `TagList` doesn't even forward it. So inline add offers zero key assist while the modal has a ranked combobox. | `TagEditor.tsx:14-17`, `TagList.tsx` |
| A5 | P2 | Single-valued-key check validates against the **server `tags` snapshot, not the live drafts**; two in-modal rows targeting the same new key don't see each other → backend 409 → generic toast, draft left dirty (spec wants inline rejection). | `ManageTagsModal.tsx:107,126,149` |
| A6 | P2 | `commitEdit` fires a mutation on **every field blur** incl. no-op (rewrites tag to itself → event-log + SSE churn) and empty-value (stores `value:""`). | `ManageTagsModal.tsx:136` |
| A7 | P2 | **key-only vs key=value has two encodings** (`undefined` in `Tag` vs `""` in `SampleTag`/DB). `TagPill` checks `value !== undefined`; `toLoupeTags` checks truthiness — they disagree, so a `value:""` tag reaching TagPill renders a broken empty-value pill. | `TagPill.tsx:33`, `loupeAdapters.ts:96-102` |
| A8 | P2 | Tag **provenance (`source` manual/scoping) is invisible** (pill never shows it) and editing a scoping tag silently keeps `source:"scoping"` while the value no longer reflects the rule. Optimistic vs SSE-synth source also disagree. | `TagPill.tsx`, `trivial.ts:408-409,433` |
| A9 | P3 | Modal drafts shadow concurrent peer edits while open (LWW, acceptable but undocumented). | `ManageTagsModal.tsx:71-80` |

**The model should be:** one shared `TagRowEditor` (key + value `TagSuggest`
comboboxes, corpus-frequency-ranked), used identically in (a) inline add, (b)
inline click-pill-to-edit, (c) the modal (N rows + an add row). It owns: trim,
loud-inline empty-key rejection, one canonical "no value = omit" representation
enforced at the API boundary, Enter-to-commit, duplicate validation against the
**live draft set** with 409s mapped back to the originating row's inline alert.
Provenance rendered as a faint marker + protected/converted on edit. The samples
sheet either adopts the editor (full parity) or drops `editable` to be honestly
read-only.

---

## Theme B — keyboard + inter-sample navigation is inconsistent (Jonathan's seed)

| # | Sev | Finding | Where |
|---|-----|---------|-------|
| B1 | P1 | **Focus has no `[`/`]` sample nav** (and no window keydown handler at all). Loupe binds `←→ x k r [ ] Esc`. Focus only gets sample nav via the *global* `,`/`.`. So `[`/`]` learned on Loupe do nothing on Focus (exactly Jonathan's report). | `FocusPage.tsx` (none) vs `LoupePage.tsx:331-360` |
| B2 | P1 | **Focus has no `Escape`-to-back.** Every other navigable surface closes/returns on Esc. | `FocusPage.tsx` |
| B3 | P1 | The Focus topbar **stepper is click-only** + has no `aria-keyshortcuts`; the global `,`/`.` that drive it are unhinted. | `CorpusTopbar.tsx:207-236` |
| B4 | P1 | **Two parallel inter-sample navigators with different affordances:** Loupe's own `‹ Prev [ / Next ]` (hinted, `[`/`]`) shows only on Loupe; the topbar stepper (`,`/`.`, unhinted) shows only on Focus. No single consistent sample-switch. | `LoupePage.tsx` + `CorpusTopbar.tsx` |
| B5 | P2 | Samples `X`/`K`/`Esc` work only with a selection and have **no legend hint** when nothing is selected (CullBar only appears post-selection). | `SamplesPage.tsx:318-332` |
| B6 | P2 | Scoping binds `⌘Z` (undo) but there is **no visible Undo button/legend** — only a "Discard". | `SeriesScopingPage.tsx:508` |
| B7 | P3 | **Reorder affordance differs:** Builder reorders via explicit `▲▼` buttons (keyboard-accessible, correct disabled-at-ends); Scoping reorders via drag handle only (no `▲▼`). Same action, two affordances; scoping lacks keyboard reorder. | `SeriesBuilderPage` (▲▼) vs `SeriesScopingPage` (drag) |

**Corrected over-claim:** the keyboard subagent reported the **Builder** as "drag
reorder only, WCAG 2.1.1 violation." **Verified false** — the Builder's TRACES
rows have real `Move up`/`Move down` buttons (first row's Move-up correctly
disabled). The keyboard-reorder gap is **Scoping's**, not the Builder's (B7).

**Direction:** give Focus the same `[`/`]` + `Esc` handler Loupe has; add
`aria-keyshortcuts` to the stepper; decide one sample-switch affordance and show
it on both Focus and Loupe; surface a visible Undo on Scoping; add a samples
X/K legend hint; give Scoping `▲▼` reorder to match the Builder.

---

## Theme C — cross-cutting + per-surface polish (my live pass)

| # | Sev | Finding | Where |
|---|-----|---------|-------|
| C1 | P2 | **No sample door is a real `<a href>`** — name, frames, and "Index" are all `<button onClick>`. No middle-click / cmd-click / open-in-new-tab anywhere, on a tool whose whole point is comparing samples. | `SampleTableRow.tsx`, app-wide doors |
| C2 | P3 | Loupe: a **kept** frame still shows both "Restore" and "Drop" — "Restore" is a no-op on an already-kept frame. The action set should reflect the current call (kept → offer Drop; dropped → offer Restore). | `LoupeSidePanel`/loupe status panel |
| C3 | P3 | Focus right rail is very sparse in the no-peaks empty state (large empty third). Honest but could carry more onboarding/guidance density. | `FocusPage.tsx` rail |

**Verified clean:** Focus log axes are correct 1-2-5 decades (0.01/0.02/0.05/0.1/0.2
· 20/50/100/200/500/1000/2000) — not irregular. Folio is well-composed. Scoping
a11y is solid (pencil = "Edit value for HEPES Only", Skip = "Skip this read: N",
reorder handles labeled). TOP-CHIPMOBILE already fixed this pass.

---

## Suggested fix order

1. **Theme A — tag model unification** (the biggest, Jonathan's seed). Start with
   the confirmed bugs A1/A3 (modal Enter-to-add + sheet lying affordance), then the
   shared `TagRowEditor` refactor (A2/A4), then the model-coherence items
   (A5–A8). A9 doc-only.
2. **Theme B — keyboard/nav parity.** B1/B2 (Focus `[`/`]` + `Esc`) are ~10 min and
   close the reported bug. Then B3/B4 (stepper hints + unify), B6 (Scoping Undo),
   B5 (samples legend), B7 (Scoping `▲▼`).
3. **Theme C** — C1 (anchor doors) is the highest-value polish; C2/C3 are P3.
