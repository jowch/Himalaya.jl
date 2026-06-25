# Interaction Architecture — Design Spec

**Date:** 2026-06-25
**Status:** Design approved; implementation not started.
**Scope:** The *interaction* layers of The Print (HimalayaUI frontend) — **not** visual design.

## 1. Problem

The app's interaction model is the culmination of partially-applied iterations. Three layers are each fragmented:

- **Dock** is a dumb `children` container; **6 pages each hand-build their own** dock contents (`ExperimentCorpusPage`, `FocusPage`, `LoupePage`, `SeriesScopingPage`, `SeriesBuilderPage`, `GroupingReviewPage`).
- **Keyboard** has a partial registry (`shortcuts.ts` — a global vocabulary feeding legend + aria) but pages still bind handlers per-page via `useShortcuts(bindings)`, and **~25 components carry ad-hoc `onKeyDown`**. Registry ids already strain (`nextDetail` means "frame" on Corpus/Loupe but "candidate" on Focus).
- **Cursor / current target** is reinvented in **5 pages** (`useState(index)` each). The dock stepper and keyboard both drive it, re-implemented each time.

The concrete symptom that motivated this: **on the Corpus you can only move the cursor with the dock stepper — clicking a row doesn't park the cursor there.** The interaction logic is stranded in the dock's local wiring; the table is non-functional without the stepper. This is not patchable; it needs a designed, consistent interaction architecture.

## 2. Goals & non-goals

**Goals**
- One consistent interaction layer the whole app shares; **each page declares what it needs** (its target + its actions) and the chrome derives.
- **Mouse and keyboard parity**: the cursor and every action are reachable by both, and they never diverge.
- Retire the half-applied remnants (page-interpreted shortcut ids, ad-hoc keydown, per-page dock construction, per-page cursor state).
- Accessible (WCAG 2.1.x) from the ground up, not bolted on.

**Non-goals**
- Visual redesign (colors, type, layout). The Dock's *appearance* and the existing `ui/` primitives stay.
- The parallel UX-consistency cleanup (error/toast channels, number/unit formatting, empty-states, modals, phase-color) — tracked separately from the UX-fragmentation inventory.

## 3. Principles

1. **Single source of truth, drift killed by elimination.** A page declares its target + actions *once*; the Dock, the keyboard layer, the on-screen hints/legend, and a future command palette all *derive* from that declaration. They cannot drift because there is nothing to synchronize.
2. **Input parity.** The cursor is input-agnostic (click, arrow, stepper all move the one cursor). Every target has a *primary action* invoked identically by `Enter` / double-click / the dock's primary button. Every action is both a clickable dock button and a bound shortcut.
3. **A + C target model.** (A) Each page is always parked on exactly **one** current item — a cursor; multi-select is a **mode** layered on top, not a co-equal. (C) Pages have **different target shapes** (Corpus = sample; Loupe = (sample, frame); Scoping = ordered members) but expose **one shared contract** (`current / moveBy / activate / actions`).
4. **Hybrid action vocabulary.** A small **fixed core** of cross-page gestures (nav, undo/redo, open Focus/Loupe, help, find) is identical everywhere; **page-local verbs** (Drop, Merge, Add peak) are declared per-page. Both flow through one `Action` shape.
5. **Build the thin layer; adopt no heavy dep.** Prior-art research (kbar, cmdk, VS Code, react-hotkeys-hook, XState, TanStack Table, React Aria) returned **study-don't-adopt / avoid** across the board — the libraries are palette-centric or would fight our Zustand/d3/design-system. We borrow the *patterns*; the registry is ~60 lines.

## 4. Core model

### 4.1 The `Action`

The atom every surface derives from:

```ts
interface Action {
  id: ActionId;
  label: string;                       // dock button · legend · aria-keyshortcuts · palette
  keys?: string[];                     // normalized combos ("x", "Mod+z", "ArrowDown")
  group: "Navigate" | "Act" | "Screen" | "Edit";
  enabled?: () => boolean;             // pure closure read via stateRef (never stale)
  run: (e?: KeyboardEvent | MouseEvent) => void;
  dock?: boolean | "primary";          // show in dock; "primary" = Enter / double-click target
  mode?: ModeKind;                     // live only in this mode; omit = always live
  glyph?: ReactNode;
}
```

`enabled()` is a **pure synchronous closure read through a `stateRef`** (`stateRef.current = latestPageState` on every render) so it is never stale. The Dock pulls it at render (button greys); the keyboard pulls it before firing (shortcut inert). One reference, no drift. This is the same ref discipline already in `useShortcuts.ts`.

### 4.2 Hybrid vocabulary — `core()` / `page()`

Two constructors make "fixed core" vs "page-local" explicit at the call site:

- `core(id, {...})` — `id` drawn from a fixed `CORE` registry (today's `shortcuts.ts`, evolved): `openFocus`, `openLoupe`, `undo`, `redo`, `help`, `find`, `back`, nav. **Keys are fixed app-wide** (`openFocus` is always `Enter`; `undo` always `Mod+Z`). The page supplies only the *handler* (and optionally `enabled`/`dock`).
- `page(id, {...})` — page-local verb (`cull`, `merge`, `addPeak`) with its own `label`/`keys`/`group`. Same `Action` shape. A build-time check rejects a page key colliding with a core key.

### 4.3 The page declaration — `usePageActions`

```ts
usePageActions({ cursor, actions: Action[] })
```

Registers the page's cursor + actions into the shell registry on mount; clears on unmount. **Registration is per-page-mount (lifecycle); within-page gating is via `enabled()`/`mode`** — the sweet spot between kbar's pure-lifecycle model and VS Code's always-registered-`when` model.

Worked example — Corpus, ~10 lines replacing the hand-built dock + per-page keydown + stepper wiring:

```ts
const cursor = useListCursor({ ids: sampleIds, onActivate: openFocus });
usePageActions({
  cursor,
  actions: [
    core("openFocus", { run: openFocus, dock: "primary" }),       // Enter · dbl-click · dock "Focus ↵"
    core("openLoupe", { run: openLoupe, dock: true }),
    page("cull", { label: "Cull", keys: ["x"], group: "Act",
                   enabled: () => cursor.selected.size > 0, run: enterCull, dock: true, mode: "selection" }),
  ],
});
// rows:  <SampleRow {...cursor.rowProps(id)} />
```

`prevSample`/`nextSample` are **not** declared — they come free from `cursor` (stepper + arrows + click all drive it).

## 5. Cursor & roving-tabindex

`useListCursor` is the **C** contract and the engine of input parity:

```ts
interface ListCursor {
  cursorId: number | null;          // stored as ID, never an index (SSE-safe)
  selected: Set<number>;            // selection-as-mode (A); empty by default
  setCursor(id): void;              // a row click lands here
  moveBy(delta): void;              // arrows AND the dock stepper both call this
  activate(): void;                 // Enter / double-click → onActivate(cursorId)
  toggleSelect(id?): void;          // selection mode
  rowProps(id): {...};              // tabIndex, ref, onClick, role, aria-current, data-cursored
  stepperProps(): {...};            // feeds DockStepper (count, prev/next, disabled)
}
useListCursor({ ids, onActivate }): ListCursor
```

**Roving tabindex — the cursor _is_ DOM focus.** Exactly one row has `tabIndex=0` (the cursor); the rest `-1`. `rowProps` wires that plus the ref. Moving the cursor calls `.focus()` on the new row — yielding **free scroll-into-view, `:focus-visible`, and "row N of M" for screen readers**, and deleting the manual `scrollIntoView` effects. The container is `role="grid" aria-multiselectable` (selection is a mode, so multiselectable is always declared).

**Every input path converges on the same two functions:**

| Gesture | Calls | Result |
|---|---|---|
| Click a row | `setCursor(id)` | cursor + DOM focus move there |
| `↑` / `↓` | `moveBy(±1)` | same |
| Dock stepper ◂ ▸ | `moveBy(±1)` | **same function** — stepper is a thin consumer of `stepperProps()` |
| `Enter` / double-click / dock primary | `activate()` | runs the page's primary action |
| `Space` (selection mode) | `toggleSelect()` | cursor keeps moving; nav is the base layer |

**Three load-bearing decisions:**
- **Cursor-as-ID** (not index) — an SSE row insert/remove before the cursor keeps the *same item* cursored. (Corpus has the index bug latently today.)
- **Cursor = DOM focus** — one substrate serves mouse + keyboard + a11y + scroll + focus-visible.
- **Selection is a `Set` orthogonal to the cursor** — never collapsed (the bug every grid library warns about).

## 6. Keyboard layer — one listener, three rungs

A precedence ladder with `e.defaultPrevented` as the cross-rung hand-off:

- **Rung 1 — widget `onKeyDown`** (element-local, React synthetic). Legit composite-widget navigation that stays local: `Menu`/`SegmentedControl` arrows, `TagEditor`, `useInlineEdit` Enter/Escape, all text inputs. They `preventDefault()` keys they own; never `stopPropagation`. **~15 of the 25 sites stay unchanged.**
- **Rung 2 — document listener** (modal layer). `ModalShell` owns Escape + focus-trap and `preventDefault()`s the Escape it consumes. **~3 sites consolidate.**
- **Rung 3 — window listener** (the registry). **The single global keyboard entry point.** Matches the event against registered actions and dispatches `run()`. **The ~6 page-level keydown blocks + `useShortcuts`/`useReorderShortcuts`/`useGlobalShortcuts` collapse here.**

**Rung 3 guards, in order:**
1. `if (e.defaultPrevented) return` — a widget/modal already handled it. **This is the microtask-race fix** — guard on the *signal*, never on `querySelector('[aria-modal]')`. (The repo already learned this via the `jsdom_dispatch_false_green` lesson.)
2. `if (isTyping(e.target) && isBareKey(e)) return` — don't fire `x`/`k` while typing; `Mod`-chords and `Escape` still pass.
3. Match → check `enabled()` → `run()`.

**WCAG 2.1.4 (Character Key Shortcuts):** bare single-key actions (`x`,`k`,`r`,`l`,`p`) fire only when focus is inside the page's interactive region — and since **cursor = DOM focus**, "a row is focused" *is* that scope. `Mod`-chords and `Escape` stay global. The legend (derived from the registry) satisfies the discoverability requirement.

## 7. The Dock as a pure renderer

The Dock moves **into the app shell** as persistent chrome and renders **from the active page's declaration**. One shell-level **`InteractionProvider`** wraps the routed pages and owns the registry store (Zustand), the rung-3 window listener, and the Dock + legend/help (+ future palette).

A page mounts → `usePageActions` populates the store → the shell's Dock and listener react. Page unmounts → store clears.

**Dock renders left→right, fully derived:**
- **Up-link** ← `core("back")` → `DockUpLink`
- **Stepper** ← `cursor.stepperProps()` (only if a cursor was declared) → `DockStepper`
- **Action buttons** ← actions with `dock:true`, grouped by `group`, each showing its key hint, `disabled = !enabled()`, `onClick = run`
- **Primary button** ← the `dock:"primary"` action, prominent + `↵`

The recently-extracted `DockStepper` / `DockUpLink` primitives become the shell's rendering vocabulary, now fed by the contract instead of hand-wired per page.

**Ownership line:**
- **Shell owns:** Dock chrome, keyboard listener, registry, legend, palette.
- **Page owns:** cursor state, handlers, data — and *declares* them in. A page's JSX stops containing `<Dock>`.

## 8. Interaction modes

Replace per-page boolean soup (`addArmed`, `previewIndexId`, `editing`, `selectionMode`…) with **one discriminated union per page**:

```ts
type CorpusMode =
  | { kind: "browse" }
  | { kind: "selection" }
  | { kind: "edit"; field: "name" };
```

Two consistency rules (no statechart — a union + a tiny escape stack; XState was study-don't-adopt):
1. **Nav is always the base layer.** Entering selection/edit never freezes `↑/↓`/click. Modes *widen* or *scope* actions; they don't replace navigation.
2. **One Escape ladder.** `Escape` unwinds the **innermost mode only** and stops (`edit → browse`, `selection → browse`); only at `browse` does it fall through to `core("back")`. The shell coordinates this (same `defaultPrevented` hand-off), so every page's Escape feels identical.

**Selection is the one mode standardized across pages** (Corpus cull, Loupe frame-cull, Scoping multi-select) — same enter/toggle/clear/Escape grammar, driven by `cursor.selected`. An action declares `mode?` when mode-specific; the page derives its action list from its current mode.

## 9. Migration & remnant removal

**Build once, prove on one page, then roll.** The `InteractionProvider` coexists with un-migrated pages: its Dock renders only for a page that called `usePageActions`; un-migrated pages keep their own `<Dock>` until their turn. Ships **page-by-page, each behind a green build + a live walk-through**.

**Order (value × risk):**
1. **Corpus** — the proof (click-to-cursor, roving tabindex, cull-as-mode, dock-from-declaration).
2. **Loupe** — compound `(sample, frame)` target (stress-tests **C**).
3. **Focus** — candidate cursor + peak-**edit** modes (stress-tests the mode stack).
4. **Scoping → Builder** — ordered members + reorder + multi-select.
5. **Grouping** — folds + rename (already on `useInlineEdit`); last.

**Remnant kill-list (deleted, not patched):**
- `shortcuts.ts` **page-interpreted ids** → explicit `core()`/`page()` declarations (the shared CORE vocabulary survives, evolved).
- `useShortcuts` · `useReorderShortcuts` · `useGlobalShortcuts` → absorbed into the one window listener; deleted at the last migration.
- The **~6 page-level `onKeyDown` blocks** → declared actions.
- **Per-page hand-built `<Dock>`** → shell renders from declaration.
- **Per-page `useState(index)` cursors** → `useListCursor` (ID-based).
- **Manual `scrollIntoView`** in cursored lists → free from roving tabindex (`safeScrollIntoView` stays for non-cursor uses).
- Retired "Locked decisions" comments + dead axis remnants in `shortcuts.ts`.

Each page migration is one reviewable PR; the branch stays unmerged for review per project convention.

## 10. Testing

**Centerpiece — a shared cursor-contract test harness**, parameterized over *any* page's `ListCursor`, making input parity mechanical:
- click-sets-cursor **==** arrow-moves-cursor **==** stepper-moves-cursor (same end state)
- `Enter` **==** double-click **==** dock primary (same `onActivate` call + arg)
- **SSE-safety:** insert/remove an id before the cursor → same item stays cursored
- **roving invariant:** exactly one row `tabIndex=0` at all times
- **orthogonality:** moving the cursor never mutates `selected`; `toggleSelect` never moves the cursor

**Keyboard-ladder tests** pin precedence + the `defaultPrevented` hand-off, including the **microtask-race regression** (a `preventDefault`'d document listener with no dialog must not let the page's `dismiss` fire). ⚠️ Because **jsdom false-greens dispatch tests** (repo's `jsdom_dispatch_false_green` lesson), the ladder + parity invariants also get a **real-browser Playwright check**, not jsdom alone.

**Per-page declaration tests** assert each migrated page registers the expected actions, the Dock renders them, and `enabled()` tracks state.

**Per-migration live walk-through** (the julia-backend + Vite + Playwright harness): real-browser verification that click→cursor, `Enter`→activate, stepper==arrows, and focus/scroll/a11y hold before the page's PR merges.

Slots into the existing six-layer contract-testing philosophy (`docs/contract-testing.md`) — the cursor harness *is* a contract test for the interaction layer.

## 11. Prior-art grounding (what we borrowed, adopted nothing)

- **Registry shape** ← kbar's `useRegisterActions` (`{ id, keys, label, section, handler, enabled, priority }`) — *pattern only*; kbar itself is palette-centric and brings its own store.
- **`enabled()`/`when` enablement** ← VS Code (always-evaluable predicate, not lifecycle) — gives no-stale-handler gating.
- **Roving tabindex + role=grid + aria-multiselectable** ← ARIA APG grid/listbox patterns — the a11y substrate that doubles as the mouse-parity + scroll substrate.
- **Cursor + selection as two orthogonal vars** ← Gmail/Linear list models.
- **Three-rung `defaultPrevented` precedence ladder** ← the synthesis of VS Code `when`-context + the repo's own modal-Escape lesson.
- **Mode = discriminated union, Escape pops innermost** ← Figma/tldraw tool-state, minus the statechart weight.

## 12. Open questions / risks

- **Command palette** is left as a near-free future surface (the registry feeds it) — out of scope for the first build; revisit after Corpus + Loupe land.
- **Compound cursors (Loupe `(sample, frame)`)** need the `ListCursor` contract to compose two axes; validate the contract shape on Loupe before committing all pages.
- **Reorder (Scoping/Builder)** interacts with the cursor + selection model; confirm the drag-reorder hook coexists with roving tabindex during that migration.
- **Risk:** roving-tabindex wiring is the main added effort vs. the window-keydown status quo; mitigated by extracting it once in `useListCursor` and proving it on Corpus first.
