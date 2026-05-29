# HimalayaUI — Functional Interaction Sweep (2026-05-29)

> **Method.** This is a *functional* sweep, not a visual audit: every interaction was traced to what it actually does. Two halves — (1) a **9-agent code inventory + reachability synthesis** across all surfaces (workflow `functional-interaction-sweep`), and (2) a **live click-through** driving the running app against the sanitized copy DB, to confirm the high-stakes cases empirically. Lens: the IA is ~70% right and stays (no new theme, no rewrite); the goal is to **minimize friction between intent and result** and to make the app actually functional.

> **Scope note.** The contact-sheet→loupe *triage* layer is wired and works. The breakage is concentrated where triage hands off to **indexing** and **series**, and in a set of controls that *look* interactive but do nothing.

## The one-sentence finding

The corpus is a viewing gallery with **no door into the work**: the primary indexing workspace (`/sample/:id`) is unreachable by clicking, the loupe and the ⌘K palette both dead-end short of it, and several of the loudest controls (Confirm series, Cross-experiment filter, "+ Peak" mode, mention chips) do nothing — so the app *reads* as broken even though the rendering and backend are sound.

## Live-verified (clicked in the running app)

| Interaction | Verified result |
|---|---|
| Corpus sample name / #id / screened-mark | **Inert** — no click target into the sample (confirmed: 0 `/sample/:id` links on the corpus) |
| Corpus thumbnail double-click | → `/samples/loupe/:id` (**loupe, never focus**) |
| Loupe → indexing | **No link to `/sample/:id` anywhere on the loupe**; only Back-to-sheet / Series |
| Top-nav "Index" tab | **`[disabled]`** on every route |
| Beamtime filter | Works (corpus 139 → 70) |
| Thumbnail select → cull bar → Drop | Works; Drop persists (grease-pencil ✕). **No un-drop on the sheet.** |
| Tag add → chip → remove | Works; inline key/value form; chip shows value (key in `title`) |
| Focus candidate-row select | Works → updates the PHASE CALL rail |
| Focus `+ Peak` arm | Button turns terracotta, **but no crosshair over the plot** (armed mode is cosmetic — see B below) |
| Focus log/lin, Auto-fit + zoom-indicator | Work (auto-fit shows "zoomed · reset") |
| Folio series card | Works → opens builder |
| `/compare/*` | Redirects to `/series` (standalone Compare page is gone) |
| @-mention chat | **Absent** on focus and in notes — the subsystem is orphaned (mounted by no page) |

Everything below is the cross-surface synthesis from the code inventory (9 agents, all 9 surfaces), which independently reached the same reachability conclusion by grepping every `navigate()` call site, and adds the builder/scoping/chat/global-chrome detail the live pass didn't reach.

---

# HimalayaUI — Cross-Surface UX Synthesis ("The Print")

## 1. Reachability map

**Start: `/samples` (the corpus, ~139 sample rows).** What a click can actually reach:

```
/samples (contact sheet)
 ├─ frame thumbnail · single-click ........ batch-select (row-local) → CullBar [Drop|Clear]
 ├─ frame thumbnail · DOUBLE-click ........ /samples/loupe/:id     ← ONLY mouse path off a row
 ├─ frame thumbnail · Enter/Space ......... /samples/loupe/:id     (keyboard: open only, no select)
 ├─ +tag / tag × .......................... inline tag CRUD (real, stays on page)
 ├─ beamtime <select> ..................... ?beamtime= filter (this page only)
 ├─ Samples tab ........................... /samples (self)
 ├─ Series tab ............................ /series
 ├─ Index tab ............................. DISABLED (no `to`) ✗
 └─ Loupe segment ......................... DISABLED ("Open a sample to use the loupe") ✗

/samples/loupe/:id (loupe)
 ├─ switch exposure (click/←→), Drop/Restore (X), Set rep (R), tag CRUD
 ├─ Back / Esc ............................ /samples (DROPS ?beamtime=)
 └─ → /sample/:id (Index) ................. DOES NOT EXIST ✗  (sidebar says rep "carries forward
                                            to the Index stage" but offers no on-ramp)

/sample/:id (Focus / Index workspace) — REACHABLE ONLY BY TYPING THE URL or legacy redirect
 ├─ topbar sample stepper ‹ › ............. /sample/:sibling  (works — the one real inter-sample nav)
 ├─ candidate row toggle, Auto-fit, Re-analyze, SpeculativeBuilder, notes
 └─ ...

/series (folio)
 ├─ + New series ......................... /series/new (scoping)
 ├─ series card (whole-card button) ....... /series/:id (builder)
 └─ Cross-experiment chip ................. ALWAYS returns [] ✗

/series/new (scoping) — sole entry is folio "+ New series"; carries NO selection
/series/:id (builder) — Confirm series CTA = no-op ✗
```

**The confirmed gap (corpus → sample), and it is worse than reported.** From `/samples` there is **no discoverable path into a sample at all**:
- Sample name + `#id` are inert `<span>`s — no onClick/href/role/focus (`ContactSheetRow.tsx:347-354`). The single most-targeted element on the row is dead.
- "Not indexed" status chip is an inert `<span>` (`ContactSheetRow.tsx:510-515`) — it literally invites indexing and offers no way to index.
- Screened-mark dot looks like a checkbox but is inert/derived (`ContactSheetRow.tsx:323-346`).
- The only mouse path is **double-click a thumbnail** (`ContactSheetRow.tsx:92,236-238`) — undiscoverable, documented only in the sticky footer legend — and it lands in the **loupe**, never the Index workspace.

**The Index workspace `/sample/:id` has no front door anywhere in the corpus→loupe path.** Grep-confirmed: the only `navigate('/sample/:id')` call sites are the topbar stepper (`CorpusTopbar.tsx:251,270`), which renders *only once you are already on `/sample/:id`*. The loupe — the place a scientist lands after the hidden double-click — also dead-ends: its only exit is "Back to the sheet" (`LoupePage.tsx:137`); the sidebar promises the rep "carries forward to the Index stage" (`LoupeSidebar.tsx:311`) but provides no link. **The Samples→Index→Series workflow stepper is not traversable from step one.** The Index tab is permanently disabled on every route except `/sample/:id`, where it's an inert label (`CorpusTopbar.tsx:139-161`).

**Other dead ends / unreachable features:**
- **NavModal (`/` or ⌘K)** — the app's "jump to any sample" picker — is a silent no-op off `/sample/:id`. `commitSample()` sets store state and closes, but never navigates (`NavModal.tsx:137-144`); the route sync is one-way URL→store (`useSyncActiveSampleFromRoute.ts:14-17`). From `/samples` or `/series` you pick a sample and *nothing visible happens*. Its only visible opener (TitleButton) is **dead code, mounted nowhere** — so even the picker's entry is invisible.
- **Entire chat + @-mention subsystem is orphaned** — `ChatCard` is mounted by no page (only referenced in its own file, its test, and a retirement comment). A fully backend-wired composer (`POST /api/samples/:id/messages`) with zero entry points. The focus workspace uses a plain `FocusNotesMargin` textarea instead.
- **Standalone Compare page** — deleted; `/compare/*` redirect to `/series`. The legacy plot (`MultiTracePlot`) is live as the builder's render core, but its peak-edit cycle (`onPeakClick`) is never wired in the builder, so clicking a peak is inert.
- **Rename** — `useUpdateSample({display_name})` exists and is wired (used only for notes), but there is no rename affordance anywhere in the UI.
- **Cross-experiment filter chip** (`folioFilter.ts:42-43`) always returns `[]` — a live-looking control that can never succeed.
- **Scoping "Ordered by" field** (`ScopingOrderField.tsx`) — styled as a dropdown with a ▾ caret, but a static div; the "one real decision on this surface" cannot be exercised.

## 2. Cross-surface friction patterns (ranked by friction added)

**P0 — Primary entities aren't clickable / primary nav silently goes nowhere.** The single dominant theme. The corpus has no clickable path into a sample (name, #id, status chip all inert). The loupe has no path onward to indexing. NavModal sample-commit is a no-op off the focus route. The Index tab is never navigable. Rename doesn't exist despite a wired mutator. Across surfaces, the most intent-aligned target a scientist reaches for is dead. *This is the friction that defines the product.* Cites: `ContactSheetRow.tsx:347-354,510-515`; `LoupePage.tsx:137` + `LoupeSidebar.tsx:311`; `NavModal.tsx:137-144`; `CorpusTopbar.tsx:139-161`.

**P0 — Affordance lies: `cursor-pointer` / dropdown styling / drag grips with no behavior.** Controls that *look* interactive and aren't, actively teaching users the app is broken:
- MentionChip `cursor-pointer` + no onClick (`MentionChip.tsx:196`).
- Reflection rows `cursor-pointer` + no onClick, no keyboard (`FocusReflectionsTable.tsx:142-166`).
- Scoping "Ordered by" ▾ dropdown that's a static div (`ScopingOrderField.tsx:10-26`).
- Scoping drag grip ⠿ `cursor-grab` with no DnD handlers (`ScopingRow.tsx:28-34`).
- Series-builder row overflow button + drag-cue, both inert (`RowActionZone.tsx:8`).
- "+ Peak" toggle whose armed state has no behavioral effect (`PlotCard.tsx:498-513`).

**P1 — Dead/no-op primary CTAs.** Accent-colored buttons (the single terracotta interaction accent in this design system) that do nothing: "Confirm series" autogroup CTA (`SeriesBuilderPage.tsx:430`) and the Cross-experiment chip that always empties the grid. A dead primary CTA erodes trust in *every* control.

**P1 — Overloaded, unlabeled gestures on canvases (esp. the trace).** Peak editing stacks add/remove/exclude on one click, plus wheel-zoom, dblclick-reset, brush-zoom — none labeled, none with cursor affordance. Same single click does three different things depending on what's under the cursor and the peak's source (`TraceViewer.tsx:264-291`). Contact-sheet thumbnails overload click=select / dblclick=open / shift=range. The sticky footer legend carries the entire discoverability burden.

**P1 — Destructive/primary affordances hidden behind hover.** "+tag" reveal is opacity-0 until hover (`ContactSheetRow.tsx:483-503`); frame selection is a hidden affordance; the "+ Add speculative" CTA is buried in a collapsed `<details>` that only opens once speculatives already exist (`PhasePanel.tsx:359-390`) — inverted from when the user needs it.

**P2 — Keyboard paths missing or invisible.** Candidate rows reachable but no focus-visible ring (`PhasePanel.tsx:93-121`); reflection rows not tab-reachable at all; filmstrip cells are clickable divs with no tabIndex/role; contact-sheet thumbs let keyboard *open* but not *select* (cull flow is mouse-only); the `,`/`.` sample shortcut is dead on `/sample/:id` though onboarding teaches it (`useGlobalShortcuts.ts:48-61`); the documented ←→ tab nav is bound to nothing.

**P2 — Inconsistent dialog/state behavior.** Notes drawer is `role=dialog` with no focus trap and no Esc (`FocusWorkspaceLayout.tsx:135-166`) while SpeculativeBuilder and modals elsewhere have both. Queue mutations (loupe Drop/Set rep) give no in-flight feedback. Notes save silently on blur with no confirmation.

**P2 — Modal where inline would do.** SpeculativeBuilder occludes the very plot it's fitting against (`SpeculativeBuilder.tsx:113-282`). Scoping Confirm modal restates what the footer note already says (two clicks for one non-destructive action).

**P2 — Persistent controls that no-op on the current surface.** Beamtime chip rendered globally but honored only on `/samples`; changeable but inert on `/series` and `/sample/:id`, and dropped on stage-switch. Offset slider stays enabled in heatmap mode where it does nothing. Cross-trace tracking checkbox toggles with no visible effect when no member has a confirmed phase.

**P3 — Stale copy / naming debt.** Onboarding teaches retired controls (←→ tabs, Inspect/Compare, title-button) (`OnboardingFlow.tsx:20-47`); `MultiTracePlot` still calls itself the Compare-page host; CLAUDE.md describes a three-card Index that no longer exists.

## 3. Broken / dead / disabled — fix-or-kill list

### P0 — Blocks the core workflow or silently misleads on every use

| Item | Location | Disposition |
|---|---|---|
| Corpus row → sample: name/#id inert | `ContactSheetRow.tsx:347-354` | **Quick wire-up** — wrap identity in `<Link to="/samples/loupe/:id">` (or `/sample/:id`). |
| "Not indexed" status chip inert | `ContactSheetRow.tsx:510-515` | **Quick wire-up** — link → `/sample/:id`. Highest-leverage single change: the most intent-aligned corpus→Index door. |
| Loupe → Index has no link | `LoupePage.tsx` / `LoupeSidebar.tsx:311-321` | **Quick wire-up** — add "Open in Focus / Index this sample" in the Representative section, enabled once a rep exists; give it a key. Closes the inspect→index loop. |
| NavModal sample-commit no-op | `NavModal.tsx:137-144` | **Quick wire-up** — add `useNavigate`, `navigate('/sample/:id')` in `commitSample`. Small fix converts the app's fastest intent→result path from dead-end to live. |
| Cross-experiment filter chip always `[]` | `folioFilter.ts:41-44` | **Rethink or hide** — needs a `spans_multiple_experiments` flag on `SeriesSummary` (backend). Hide until then; a control that universally empties is worse than absent. |
| "Confirm series" CTA = no-op | `SeriesBuilderPage.tsx:430` | **Rethink** — wire to persist/confirm grouping, or delete. Dead terracotta primary CTA. |
| "+ Peak" armed mode no behavioral effect | `PlotCard.tsx:498-513` / `TraceViewer.tsx:264-291` | **Rethink** — make armed mode real (interior click adds, suppress remove/exclude) or delete the button. A mode that isn't a mode is a trap. |

### P1 — Dead affordances that look interactive

| Item | Location | Disposition |
|---|---|---|
| MentionChip `cursor-pointer`, no onClick | `MentionChip.tsx:196` | **Quick wire-up** (if chat revived) — onClick → navigate by entity type; else remove cursor-pointer. |
| Scoping "Ordered by" fake dropdown | `ScopingOrderField.tsx:10-26` | **Quick win or strip** — wire to corpus tag-key picker (freq map already computed) or remove the ▾/dropdown styling. |
| Scoping drag grip ⠿, no DnD | `ScopingRow.tsx:28-34` | **Strip** (order is machine-derived) or implement DnD. |
| Builder row overflow + drag-cue inert | `RowActionZone.tsx:8` | **Build menu or delete**; drop false drag-cue. |
| "Variable" sort = A-Z title sort | `folioFilter.ts:48-50` | **Quick fix** — relabel "Title", or surface `ordering_variable` on the listing and sort for real. A label that contradicts behavior erodes trust. |
| `+N new match` pill never rendered | `SeriesFolioCard.tsx:133-150` | **Defer/strip** — never receives a non-zero prop; if lit later, wire to an add-matching action, not a passive label. |

### P2 — Disabled-that-shouldn't-be / wiring gaps

| Item | Location | Disposition |
|---|---|---|
| `,`/`.` shortcut dead on `/sample/:id` | `useGlobalShortcuts.ts:48-61` + `AppRoutes.tsx:79` | **Quick wire-up** — bind to `navigate('/sample/:sibling)` using the stepper's sibling derivation. |
| Reflection row `cursor-pointer`, no onClick/keyboard | `FocusReflectionsTable.tsx:142-166` | **Quick wire-up** — give the row a job (zoom trace to q / toggle exclude) + tabIndex/onFocus, or drop cursor-pointer. |
| Notes drawer no focus-trap / no Esc | `FocusWorkspaceLayout.tsx:135-166` | **Quick wire-up** — add `useFocusTrap` + Esc, matching SpeculativeBuilder. |
| MentionChip stuck "loading" on non-404 | `useMentionResolution.ts:43-45` | **Quick fix** (if chat revived) — distinct "unavailable"/retry state. |
| Index tab never navigable | `CorpusTopbar.tsx:139-161` | **Rethink** — navigate to last/first sample's focus, or drop the tab and keep it only as the active marker. |
| Loupe Back drops `?beamtime=` | `LoupePage.tsx:137` | **Quick fix** — preserve query, mirroring `CorpusTopbar.tsx:196`. |
| Beamtime chip no-op off `/samples` | `CorpusTopbar.tsx:165-179` | **Quick fix** — hide where it does nothing, or make `/series` honor it. |
| Empty states with no action (folio no-series / no-match; builder no-members; scoping cold-corpus) | `SeriesFolioPage.tsx:164-171`; `SeriesBuilderPage.tsx:339-349`; `SeriesScopingPage.tsx:269-277` | **Quick wire-up** — embed the CTA (New series / Clear filters / Edit recipe / link to contact sheet). |

### P3 — Cleanup / honesty

- Onboarding stale copy (←→ tabs, Inspect/Compare, title-button) — `OnboardingFlow.tsx:20-47`. **Rewrite.**
- TitleButton dead code — **delete or re-home** as the NavModal's visible opener.
- `MultiTracePlot` stale Compare naming/docstrings — **rename** to series-builder vocabulary.
- Speculative delete 🗑 emoji — `PhasePanel.tsx:159-170`. **Replace** with ink-stroke SVG (only emoji in the rail, off-brand).
- Wordmark not a home link — `CorpusTopbar.tsx:95`. **Wrap** in `<Link to="/samples">`.
- Loupe "integration"/"collected" permanent em-dashes — `LoupeSidebar.tsx:241-242`. **Plumb or hide.**
- **Potential data bug (verify):** scoping `toggleFlag` clears the flag without checking `value!==''` (`SeriesScopingPage.tsx:179-183`), and `canBuild` only checks `!flagged` — an accepted-but-empty ordering value can pass the gate and write an empty-string tag. Confirm backend rejects empty values; if not, it's a wiring bug, not friction.

## 4. Redesign requirements (the four named targets)

### A. Tag component (corpus "+tag" chips + add form)
**What's wrong:** Mechanically sound and fully wired (real queue mutators) — the friction is discoverability and form symmetry, not breakage. The "+" reveal is opacity-0 until hover/focus-within when chips exist (`ContactSheetRow.tsx:483-503`) — invisible on touch, easy to miss. The key input has no Enter/Esc handler though the value input does (`ContactSheetRow.tsx:435-481`; same defect in `LoupeSidebar.tsx:136-159`) — submitting from the key field silently fails. There's no tag/value autocomplete over the corpus, so tagging is free-text-prone to drift (relevant downstream: scoping groups on tag *keys*).
**Low-friction shape:** Low-contrast **always-visible "+"** (quiet, not hover-gated). Add the same `onKeyDown` (Enter/Esc) to the key input. Add **key/value autocomplete** drawn from existing corpus tag keys so tags stay consistent enough to drive scoping. Keep the chip-body-inert / ×-only-removes model (correct).
**Verdict: modernize-in-place.** It's the most finished component in the audit; tighten discoverability + form parity + autocomplete.

### B. Plot / trace components (Focus TraceViewer / PlotCard — peak add/move/delete)
**What's wrong:** Peak editing is entirely gesture-based and unlabeled — the heart of "plot components need redesign." Empty-click=add, click-manual=delete, click-auto=toggle-exclude (`TraceViewer.tsx:264-291`), wheel=zoom, dblclick=reset; none discoverable, none with cursor feedback. The "+ Peak" toggle is a *false mode* (P0 above). Two "reset" verbs disagree: dblclick resets to full q-range while "zoomed · reset" / Auto-fit re-FIT to features, and the ZoomIndicator labeled "reset" actually calls `onFitFeatures` (`PlotCard.tsx:480,567-582`). There is **no peak-move** at all — only add/delete/exclude — and the owner explicitly named "move" as friction.
**Low-friction shape:**
- **Make "+ Peak" a real mode** with a crosshair cursor + inline "click to place a peak" hint while armed; suppress the remove/exclude branch in that mode; disarm on place.
- **Disambiguate triangle clicks** with hover tooltips ("click to exclude" on auto, "click to delete" on manual) and visually distinguish auto vs manual vs in-flight (the outlined optimistic triangle is a good guard — make it dashed for legibility, `TraceViewer.tsx:111-130`).
- **Add peak-move**: drag a triangle along q (the missing verb), with a snap-to-local-max assist.
- **Reconcile the two resets**: relabel the ZoomIndicator "fit"; reserve "reset" for true full-range. Add a one-time zoom/brush hint.
- Add focus-visible rings on candidate rows / scale toggle for keyboard parity.
**Verdict: modernize-in-place** — the rendering and mutation plumbing are solid; the redesign is the *interaction grammar* (modes, labels, move, reset coherence).

### C. Old compare plot component
**Reachability:** The standalone Compare page is **gone** — `ComparePage/Edit/Compare.tsx` deleted; `/compare/*` redirect to `/series` (`AppRoutes.tsx:123-124`). The component itself, `MultiTracePlot.tsx`, is **live** — the render core of `SeriesBuilderPage` (`:8,371`), `SeriesMiniWaterfall`, and figure export. The fold already happened.
**What's broken in its new home:** Wheel/dblclick/brush zoom are undiscoverable, brush draws no live selection rectangle (`MultiTracePlot.tsx:435-577`); the per-member peak-edit cycle (`onPeakClick`, `:477-502`) is shipped-but-inert because the builder never passes `onPeakClick` — clicking a peak does nothing; view state (grouping/representation/offset/scale/tracking) is local-only and lost on navigate-away; `Adjust` and `Edit` are duplicate entry points to the same `onStartEdit`; offset slider stays live in heatmap mode where it's a no-op; stale "Compare page" naming throughout.
**Verdict: MODERNIZE-IN-PLACE, do NOT delete.** Concrete requirements: (1) own/wire-or-remove the "Confirm series" no-op CTA; (2) decide explicitly to **wire** `onPeakClick` into the builder's edit mode (so per-trace peak curation works) **or strip** the hit-test plumbing as dead; (3) surface zoom/brush/reset (visible reset chip, live brush rectangle, crosshair); (4) persist view-composition (`view_grouping_mode` exists); (5) rename Compare-era identifiers.

### D. @-mention chips (never designed — full interaction + chip design)
**Reality is worse than "never designed":** the *entire chat + mention subsystem is orphaned* — `ChatCard` is mounted nowhere; backend message endpoints are live with no reachable UI; the focus workspace shipped the weaker `FocusNotesMargin` (single textarea, no threading/author/timestamp/citation) instead. **Decision is binary and must come first** (this is a `/shape` question, not a `/craft` one):
- **(A) Revive** — re-home a threaded notes surface on `/sample/:id` (replacing/augmenting FocusNotesMargin), because the scientist genuinely wants citable notes ("clean Pn3m here, see @peak q=0.064"). Then the chip requirements below apply.
- **(B) Retire** — delete ChatCard/MentionCompose/MentionPicker/MentionChip/useMentionResolution/renderMentions + the message hooks; keep FocusNotesMargin.

**If revived, chip + interaction requirements:**
1. **Click-to-navigate is the whole payoff** — chips carry `cursor-pointer` but no onClick (`MentionChip.tsx:196`). Wire: `@sample`→`/sample/:id`, `@exposure`→loupe, `@index`/`@peak`→focus + highlight. A non-navigable citation is a dead reference.
2. **Re-token to paper+ink+terracotta** — chip colors and tooltip are hardcoded dark-theme hex (`#252525`/`#3a3a3a`, `MentionChip.tsx:22-42,211-216`) — leftover "Darkroom" styling that bypasses design tokens and will clash on warm paper.
3. **Surface the grammar** — the `@` trigger and especially the powerful `:type` filter (`@peak:`, `@index:`, `MentionPicker.tsx:35-39`) are invisible; add a visible "@ mention" affordance and render kind headers/tabs in the picker.
4. **Cite-from-context** — let a user cite the peak/index they're looking at on the trace rather than re-searching in the picker.
5. **Fix degraded states** — non-404 errors should be a distinct "unavailable"/retry state, not permanent "loading" (`useMentionResolution.ts:43-45`); preserve the good 404 "no longer exists" handling.
6. **Retarget or drop the comparison arm** — the comparison-mention + hash-drift machinery targets a retired concept (`MentionChip.tsx:124-132`); retarget at Series or delete.
7. Co-mount with the trace plot so the hover-halo cross-highlight (which *does* work, `MentionChip.tsx:144-156`) is visible.
8. Make "Sign in to post…" an actual link to OnboardingFlow, not a dead placeholder.

**Verdict: shape-then-decide.** Lean revive (the citable-notes value is real and the backend is built), but the chips MUST be navigable + on-palette or the system is decorative. Either way, don't leave a backend-wired composer with zero entry points.

## 5. Recommended sequencing (mapped to /impeccable commands)

**Wire immediately — cheap, high-friction-relief (no IA change, mostly `/live` + `/clarify`):**
1. **Corpus→sample door** (`/live`): make "Not indexed" chip + sample name links into `/sample/:id` (or loupe). *Single highest-leverage change in the audit* — it gives the workflow a front door. `ContactSheetRow.tsx:347-354,510-515`.
2. **Loupe→Index "Open in Focus"** (`/live`): closes the inspect→index loop. `LoupeSidebar.tsx:314`.
3. **NavModal `useNavigate` in commitSample** (`/live`): converts the app's fastest jump-to-sample from dead-end to live. `NavModal.tsx:137-144`.
4. **Kill the lying controls** (`/live` + `/clarify`): hide Cross-experiment chip; relabel "Variable"→"Title"; remove `cursor-pointer` from inert chips/rows; strip the scoping drag grip + builder drag-cue; either wire or remove the "+ Peak" mode and "Confirm series" CTA. *Lying controls erode trust in every other control — do these before any polish.*
5. **Embed empty-state CTAs** (`/clarify`): folio, builder, scoping cold-corpus. `SeriesFolioPage.tsx:164-171`, etc.
6. **Quick keyboard/dialog fixes** (`/harden`): `,`/`.` on focus route; notes-drawer focus-trap + Esc; key-input Enter/Esc; candidate-row focus-visible ring.
7. **Stale-copy sweep** (`/clarify`): rewrite onboarding to current IA; delete/re-home TitleButton; rename MultiTracePlot's Compare vocabulary. (Update CLAUDE.md to match.)

**Critique first, then shape before building (need a design decision, not just a wire-up):**
8. **`/critique` the trace peak-editing grammar** → **`/shape`** the mode/move/disambiguation/reset model → **`/craft`**. This is the owner's named "plot redesign"; the gesture overload and missing peak-move need a designed model, not a patch. (Target B.)
9. **`/shape` the @-mention revive-or-retire decision** (Target D) — binary, blocks all chip work. If revive: `/shape` the threaded-notes surface on `/sample/:id`, then `/craft` chips (navigate + tokens + grammar) and `/colorize`/re-token off the dark palette.
10. **`/shape` the series-builder peak-edit decision** (Target C): wire `onPeakClick` into edit mode or strip it; then `/live` the zoom/brush discoverability + view-state persistence.
11. **`/shape` selection hand-off from contact sheet → `/series/new`** — the root cause behind scoping's "You selected N samples" lie, whole-corpus loose-matches, and the missing contact-sheet entry point. One plumbing change fixes four findings.
12. **`/harden`** pass last: queue-mutation in-flight feedback (loupe), focus-visible across remaining controls, keyboard reach for filmstrip/reflection rows, scoping empty-tag data-bug guard.

**Sequencing rationale:** steps 1-3 restore the broken core loop (corpus→sample→index) at near-zero cost; step 4 stops the trust erosion from fake affordances; 5-7 are cheap clarity wins. Only then spend design effort (8-11) on the genuinely under-specified surfaces — the trace grammar and the mention system — where building before shaping would just re-pour concrete. Tag (A) is `/craft`-in-place and can slot in anywhere after step 4.

**The through-line for the owner's one goal (minimize intent→result friction):** the corpus is a viewing gallery with no door into the work; the loupe and NavModal both dead-end short of indexing; and several of the loudest controls do nothing. Fix the doors (1-3), kill the lies (4), then design the two surfaces that were never finished (trace grammar, mentions). None of this requires a new theme or IA rewrite — it's wiring the doors the existing IA already implies.