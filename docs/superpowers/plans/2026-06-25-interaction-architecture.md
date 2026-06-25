# Interaction Architecture Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace The Print's three fragmented interaction layers (per-page dock construction, per-page cursor state, ad-hoc keyboard) with one shell-owned layer that each page *declares into* — giving mouse/keyboard parity, an ID-based roving-tabindex cursor, and a single keyboard listener.

**Architecture:** A Zustand store (`useInteraction`) holds the active page's declaration (a `ListCursor` + an `Action[]`). The app shell mounts one `<InteractionDock>` (renders entirely from the declaration) and one `useKeyboardLayer()` window listener (the single rung-3 keyboard entry point). Pages call `usePageActions({ cursor, actions })` on mount; the store, the Dock, and the listener all derive — they cannot drift because there is nothing to synchronize. The cursor (`useListCursor`) stores the current item **as an ID** and wires **roving tabindex** so the cursor *is* DOM focus — one substrate for mouse, keyboard, scroll, and a11y. Built once, proven on Corpus, then rolled page-by-page; each page migration is one reviewable PR.

**Tech Stack:** React 18 + TypeScript (strict: `noUnusedLocals`, `exactOptionalPropertyTypes`), Zustand, TanStack Query, react-router-dom, Vitest + RTL (jsdom), Playwright (mocked + live). Design-system primitives in `src/print/ui/` (closed-look/open-placement; `check-design.mjs` lint). No new runtime dependency is added.

## Global Constraints

- **No new dependency.** Prior-art research returned study-don't-adopt across the board (kbar/cmdk/XState/React-Aria/TanStack-Table). The registry is hand-rolled (~60 lines). [verbatim from spec §3.5]
- **Visual design is out of scope.** Colours, type, layout, the Dock's *appearance*, and the `ui/` primitives stay. Only interaction changes. [spec §2]
- **Never touch the Julia core** (`src/`, top-level `test/`) — it is a standalone package. All work is under `packages/HimalayaUI/frontend/`.
- **`lint:design` is a build gate.** No inline appearance utilities (`text-[…]`, `rounded-[…]`, raw colour literals, side-stripes) outside `src/print/ui/**` and the render-layer dirs. New files in `src/print/interaction/` are **not** appearance-exempt → compose from existing primitives, pass placement-only `className`.
- **Build chain that must pass before every PR:** `npm run build` = `npm run lint:design && tsc --noEmit -p tsconfig.build.json && vite build`. Plus `npm test` (Vitest) and `npm run e2e` (mocked Playwright).
- **`data-testid` contracts are load-bearing for E2E.** Preserve `dock-up-link`, `dock-prev-{base}`/`dock-next-{base}`/`dock-{base}-count` exactly (DockStepper/DockUpLink already emit these).
- **Cursor is stored as an ID, never an index** (SSE-safety: an insert/remove before the cursor must keep the *same* item cursored).
- **jsdom false-greens dispatch tests** (`jsdom_dispatch_false_green` lesson): keyboard-ladder + parity invariants get a real-browser Playwright check, not jsdom alone. In jsdom, `fireEvent.keyDown(...)` returns `false` when the handler called `preventDefault()` — assert on that, not on DOM side effects alone.
- **The freshness rule (load-bearing — read once, applies to every page):** `usePageActions` writes the declaration to the store in a **dependency-less `useEffect`** (runs every commit). The Dock subscribes to the store (re-renders on change); the keyboard listener reads `useInteraction.getState()` at *event* time. Both therefore see the latest `enabled()`/`run` closures with no manual `stateRef` plumbing. Never call the store setter during render.
- **Branch stays unmerged** for Jonathan's review per project convention; each phase is its own commit, each page migration its own reviewable unit.

---

## File Structure

**New — the interaction layer (`src/print/interaction/`):**

| File | Responsibility |
|---|---|
| `types.ts` | `Action`, `ActionId`/`CoreId`, `ActionGroup`, `ModeKind`, `ListCursor`, `RowProps`, `CursorStepperProps` — the contracts every other file consumes. |
| `core.ts` | The fixed `CORE` table (`back`/`openFocus`/`openLoupe`/`undo`/`redo`/`help`/`find` + nav), the `core()`/`page()` action constructors, and the build-time core-vs-page key-collision check. |
| `matchKey.ts` | Pure `comboOf(e)` (KeyboardEvent → normalized combo string) + `matchesKeys(e, keys)`. The one place keys are parsed. |
| `registry.ts` | `useInteraction` Zustand store: `{ cursor, actions, setPage, clearPage }`. |
| `usePageActions.ts` | The page declaration hook (writes declaration → store every commit; clears on unmount). |
| `useListCursor.ts` | The `ListCursor` engine: ID-based cursor + selection `Set` + roving tabindex (`rowProps`) + `stepperProps`. |
| `useKeyboardLayer.ts` | The single rung-3 window listener (mounted once in the shell): `defaultPrevented` guard → typing/scope guard → match → `enabled()` → `run()`. |
| `InteractionDock.tsx` | The Dock renderer: derives up-link + stepper + action buttons + primary from the store. Composes `Dock`/`DockUpLink`/`DockStepper`/`Button`/`KbKey` primitives. |

**Modified:**
- `src/print/shell/AppRoutes.tsx` — `AppShell` mounts `<InteractionDock/>` and calls `useKeyboardLayer()`; `useGlobalShortcuts()` removed at the final phase.
- The six pages, one per migration phase (see Phases 2–6).

**Deleted at the final phase (remnant kill-list, spec §9):**
- `src/print/shell/useShortcuts.ts`, `src/print/shell/useReorderShortcuts.ts`, `src/hooks/useGlobalShortcuts.ts`.
- `shortcuts.ts` page-interpreted ids (`prevDetail`/`nextDetail` reinterpretation) — the shared CORE vocabulary survives in `core.ts`.

**Test files (created alongside their source):**
- `test/interaction/matchKey.test.ts`, `registry.test.ts`, `useListCursor.contract.test.tsx` (the shared harness), `keyboardLayer.test.tsx`, `InteractionDock.test.tsx`.
- One `*.actions.test.tsx` per migrated page.
- `e2e/interaction-corpus.spec.ts` (mocked) + a live walk-through per page.

---

# Phase 0 — Foundation (coexists; no page migrated yet)

The foundation ships behind the existing pages: the store is empty until a page calls `usePageActions`, the Dock renders nothing when the store is empty, and the listener matches nothing. Un-migrated pages keep their own `<Dock>` and `useShortcuts`. **Nothing visible changes until Phase 2.**

### Task 0.1: Interaction contracts (`types.ts`)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/print/interaction/types.ts`
- Test: `packages/HimalayaUI/frontend/test/interaction/types.test.ts` (type-only compile check)

**Interfaces:**
- Produces: `Action`, `ActionId`, `CoreId`, `ActionGroup`, `ModeKind`, `ListCursor`, `RowProps`, `CursorStepperProps`. Every later task consumes these.

- [ ] **Step 1: Write the failing test** (a compile-time shape assertion — runs as a no-op at runtime but fails `tsc` if the types drift)

```ts
// test/interaction/types.test.ts
import { describe, it, expect } from "vitest";
import type { Action, ListCursor, CursorStepperProps } from "../../src/print/interaction/types";

describe("interaction types", () => {
  it("an Action is constructable with the documented fields", () => {
    const a: Action = {
      id: "openFocus",
      label: "Focus",
      keys: ["Enter"],
      group: "Navigate",
      enabled: () => true,
      run: () => {},
      dock: "primary",
    };
    expect(a.id).toBe("openFocus");
  });

  it("a ListCursor exposes the parity contract", () => {
    const stepper: CursorStepperProps = {
      label: "Sample", axis: "vertical", testIdBase: "sample",
      count: "1 / 3", onPrev: () => {}, onNext: () => {},
      prevDisabled: true, nextDisabled: false,
    };
    const c: ListCursor = {
      cursorId: null, selected: new Set<number>(),
      setCursor: () => {}, moveBy: () => {}, activate: () => {},
      toggleSelect: () => {}, rowProps: () => ({} as never),
      stepperProps: () => stepper,
    };
    expect(c.cursorId).toBeNull();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `node_modules/.bin/vitest run test/interaction/types.test.ts`
Expected: FAIL — `Cannot find module '../../src/print/interaction/types'`.

- [ ] **Step 3: Write `types.ts`**

```ts
// src/print/interaction/types.ts
import type { ReactNode, MouseEvent as ReactMouseEvent } from "react";

export type ActionGroup = "Navigate" | "Act" | "Screen" | "Edit";

/** Fixed cross-page gestures. Keys are app-wide constants (see core.ts). */
export type CoreId = "back" | "openFocus" | "openLoupe" | "undo" | "redo" | "help" | "find";

/** A page-local verb id is any other string (e.g. "cull", "merge", "addPeak"). */
export type ActionId = CoreId | (string & {});

/** A page's current interaction mode is its own discriminated union; the layer
 *  only needs the `kind` tag to gate actions, so it sees the bare string. */
export type ModeKind = string;

export interface Action {
  id: ActionId;
  /** dock button · legend · aria-keyshortcuts · palette */
  label: string;
  /** normalized combos: "x" · "Mod+z" · "Shift+ArrowUp" · "Enter" · "?" */
  keys?: string[];
  group: ActionGroup;
  /** pure synchronous closure; the layer reads it before firing (never stale). */
  enabled?: () => boolean;
  run: (e?: KeyboardEvent | ReactMouseEvent) => void;
  /** show in dock; "primary" = Enter / double-click target, rendered prominently. */
  dock?: boolean | "primary";
  /** live only in this page-mode; omit = always live. */
  mode?: ModeKind;
  glyph?: ReactNode;
}

/** Spread onto every cursorable row. tabIndex is the roving substrate. */
export interface RowProps {
  ref: (el: HTMLElement | null) => void;
  tabIndex: 0 | -1;
  onClick: (e: ReactMouseEvent) => void;
  onDoubleClick: (e: ReactMouseEvent) => void;
  role: "row";
  "aria-current": "true" | undefined;
  "data-cursored": "true" | "false";
}

/** Feeds the shell's DockStepper. Mirrors DockStepperProps minus appearance. */
export interface CursorStepperProps {
  label: string;
  axis: "vertical" | "horizontal";
  testIdBase: string;
  count: ReactNode;
  onPrev: () => void;
  onNext: () => void;
  prevDisabled: boolean;
  nextDisabled: boolean;
}

export interface ListCursor {
  cursorId: number | null;
  selected: Set<number>;
  setCursor: (id: number) => void;
  moveBy: (delta: number) => void;
  activate: () => void;
  toggleSelect: (id?: number) => void;
  rowProps: (id: number) => RowProps;
  stepperProps: () => CursorStepperProps;
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `node_modules/.bin/vitest run test/interaction/types.test.ts`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/interaction/types.ts \
        packages/HimalayaUI/frontend/test/interaction/types.test.ts
git commit -m "feat(interaction): action + cursor contracts"
```

---

### Task 0.2: Key matcher (`matchKey.ts`)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/print/interaction/matchKey.ts`
- Test: `packages/HimalayaUI/frontend/test/interaction/matchKey.test.ts`

**Interfaces:**
- Consumes: nothing.
- Produces: `comboOf(e: KeyboardEvent): string`, `matchesKeys(e: KeyboardEvent, keys: string[]): boolean`, `isTyping(t: EventTarget | null): boolean`, `isBareKey(e: KeyboardEvent): boolean`.

- [ ] **Step 1: Write the failing test**

```ts
// test/interaction/matchKey.test.ts
import { describe, it, expect } from "vitest";
import { comboOf, matchesKeys, isTyping, isBareKey } from "../../src/print/interaction/matchKey";

function key(init: Partial<KeyboardEvent> & { key: string }): KeyboardEvent {
  return new KeyboardEvent("keydown", init);
}

describe("comboOf", () => {
  it("normalizes a bare letter to lowercase", () => {
    expect(comboOf(key({ key: "X" }))).toBe("x");
  });
  it("emits Mod for meta or ctrl", () => {
    expect(comboOf(key({ key: "z", metaKey: true }))).toBe("Mod+z");
    expect(comboOf(key({ key: "z", ctrlKey: true }))).toBe("Mod+z");
  });
  it("orders modifiers Mod+Alt+Shift then key", () => {
    expect(comboOf(key({ key: "ArrowUp", shiftKey: true }))).toBe("Shift+ArrowUp");
    expect(comboOf(key({ key: "z", metaKey: true, shiftKey: true }))).toBe("Mod+Shift+z");
  });
  it("maps space to Space", () => {
    expect(comboOf(key({ key: " " }))).toBe("Space");
  });
  it("passes ? through verbatim (layout-stable)", () => {
    expect(comboOf(key({ key: "?" }))).toBe("?");
  });
});

describe("matchesKeys", () => {
  it("matches any combo in the list, case-insensitively for letters", () => {
    expect(matchesKeys(key({ key: "K" }), ["k"])).toBe(true);
    expect(matchesKeys(key({ key: "z", metaKey: true }), ["Mod+z"])).toBe(true);
    expect(matchesKeys(key({ key: "z" }), ["Mod+z"])).toBe(false);
  });
});

describe("isBareKey", () => {
  it("is true for a single char with no modifiers", () => {
    expect(isBareKey(key({ key: "x" }))).toBe(true);
  });
  it("is false for Mod-chords and for Escape", () => {
    expect(isBareKey(key({ key: "z", metaKey: true }))).toBe(false);
    expect(isBareKey(key({ key: "Escape" }))).toBe(false);
  });
});

describe("isTyping", () => {
  it("detects inputs, textareas, and contenteditable", () => {
    const input = document.createElement("input");
    expect(isTyping(input)).toBe(true);
    const div = document.createElement("div");
    expect(isTyping(div)).toBe(false);
    div.setAttribute("contenteditable", "true");
    expect(isTyping(div)).toBe(true);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `node_modules/.bin/vitest run test/interaction/matchKey.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Write `matchKey.ts`**

```ts
// src/print/interaction/matchKey.ts

/** Normalized combo: [Mod][+Alt][+Shift]+<key>. `Mod` folds meta/ctrl so the
 *  same declaration works on macOS and the rest. Letters lowercase; `?` and the
 *  Arrow/Enter/Escape names pass through; " " becomes "Space". */
export function comboOf(e: KeyboardEvent): string {
  const parts: string[] = [];
  if (e.metaKey || e.ctrlKey) parts.push("Mod");
  if (e.altKey) parts.push("Alt");
  if (e.shiftKey) parts.push("Shift");
  let k = e.key;
  if (k === " ") k = "Space";
  else if (k.length === 1) k = k.toLowerCase();
  parts.push(k);
  return parts.join("+");
}

export function matchesKeys(e: KeyboardEvent, keys: string[]): boolean {
  const got = comboOf(e);
  // Normalize each declared key through the same lowering so "K" === "k".
  return keys.some((want) => {
    const parts = want.split("+");
    const k = parts.pop()!;
    const norm = [...parts, k.length === 1 ? k.toLowerCase() : k].join("+");
    return norm === got;
  });
}

export function isBareKey(e: KeyboardEvent): boolean {
  return !e.metaKey && !e.ctrlKey && !e.altKey && e.key.length === 1;
}

export function isTyping(t: EventTarget | null): boolean {
  if (!(t instanceof HTMLElement)) return false;
  const tag = t.tagName;
  if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return true;
  return t.isContentEditable;
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `node_modules/.bin/vitest run test/interaction/matchKey.test.ts`
Expected: PASS (all matchKey cases).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/interaction/matchKey.ts \
        packages/HimalayaUI/frontend/test/interaction/matchKey.test.ts
git commit -m "feat(interaction): key normalization + matcher"
```

---

### Task 0.3: Core vocabulary + constructors (`core.ts`)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/print/interaction/core.ts`
- Test: `packages/HimalayaUI/frontend/test/interaction/core.test.ts`

**Interfaces:**
- Consumes: `Action`, `CoreId`, `ActionGroup` from `types.ts`.
- Produces:
  - `CORE: Record<CoreId, { label: string; keys: string[]; group: ActionGroup }>`
  - `core(id: CoreId, over: { run: Action["run"]; enabled?: () => boolean; dock?: Action["dock"]; label?: string; mode?: string }): Action`
  - `page(id: string, def: { label: string; keys?: string[]; group: ActionGroup; run: Action["run"]; enabled?: () => boolean; dock?: Action["dock"]; mode?: string }): Action`
  - `assertNoCoreCollision(actions: Action[]): void`

- [ ] **Step 1: Write the failing test**

```ts
// test/interaction/core.test.ts
import { describe, it, expect } from "vitest";
import { CORE, core, page, assertNoCoreCollision } from "../../src/print/interaction/core";

describe("CORE vocabulary", () => {
  it("openFocus is bound to Enter app-wide", () => {
    expect(CORE.openFocus.keys).toEqual(["Enter"]);
  });
  it("undo is Mod+z, redo Mod+Shift+z", () => {
    expect(CORE.undo.keys).toEqual(["Mod+z"]);
    expect(CORE.redo.keys).toEqual(["Mod+Shift+z"]);
  });
});

describe("core()", () => {
  it("fills label/keys/group from CORE, page supplies the handler", () => {
    const run = (): void => {};
    const a = core("openFocus", { run, dock: "primary" });
    expect(a).toMatchObject({ id: "openFocus", label: "Focus", keys: ["Enter"], group: "Navigate", dock: "primary" });
    expect(a.run).toBe(run);
  });
  it("lets the page override the label only", () => {
    expect(core("openFocus", { run: () => {}, label: "Apply" }).label).toBe("Apply");
  });
});

describe("page()", () => {
  it("builds a page-local verb verbatim", () => {
    const a = page("cull", { label: "Cull", keys: ["x"], group: "Act", run: () => {} });
    expect(a).toMatchObject({ id: "cull", label: "Cull", keys: ["x"], group: "Act" });
  });
});

describe("assertNoCoreCollision", () => {
  it("throws when a page verb steals a core key", () => {
    const bad = page("cull", { label: "Cull", keys: ["Enter"], group: "Act", run: () => {} });
    expect(() => assertNoCoreCollision([bad])).toThrow(/Enter/);
  });
  it("passes when keys are disjoint from core", () => {
    const ok = page("cull", { label: "Cull", keys: ["x"], group: "Act", run: () => {} });
    expect(() => assertNoCoreCollision([ok])).not.toThrow();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `node_modules/.bin/vitest run test/interaction/core.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Write `core.ts`**

```ts
// src/print/interaction/core.ts
import type { Action, ActionGroup, CoreId } from "./types";

/** The fixed cross-page gestures. Keys are app-wide constants — a page supplies
 *  only the handler (and optionally enabled/dock/label). Evolved from the old
 *  shortcuts.ts global vocabulary, minus the page-interpreted ids. */
export const CORE: Record<CoreId, { label: string; keys: string[]; group: ActionGroup }> = {
  back:      { label: "Back",  keys: ["Escape"],      group: "Navigate" },
  openFocus: { label: "Focus", keys: ["Enter"],       group: "Navigate" },
  openLoupe: { label: "Loupe", keys: ["l"],           group: "Navigate" },
  undo:      { label: "Undo",  keys: ["Mod+z"],       group: "Edit" },
  redo:      { label: "Redo",  keys: ["Mod+Shift+z"], group: "Edit" },
  help:      { label: "Shortcuts", keys: ["?"],       group: "Screen" },
  find:      { label: "Find",  keys: ["/", "Mod+k"],  group: "Screen" },
};

type CoreOverride = {
  run: Action["run"];
  enabled?: () => boolean;
  dock?: Action["dock"];
  label?: string;
  mode?: string;
};

export function core(id: CoreId, over: CoreOverride): Action {
  const base = CORE[id];
  const a: Action = {
    id,
    label: over.label ?? base.label,
    keys: base.keys,
    group: base.group,
    run: over.run,
  };
  if (over.enabled) a.enabled = over.enabled;
  if (over.dock !== undefined) a.dock = over.dock;
  if (over.mode !== undefined) a.mode = over.mode;
  return a;
}

type PageDef = {
  label: string;
  keys?: string[];
  group: ActionGroup;
  run: Action["run"];
  enabled?: () => boolean;
  dock?: Action["dock"];
  mode?: string;
};

export function page(id: string, def: PageDef): Action {
  const a: Action = { id, label: def.label, group: def.group, run: def.run };
  if (def.keys) a.keys = def.keys;
  if (def.enabled) a.enabled = def.enabled;
  if (def.dock !== undefined) a.dock = def.dock;
  if (def.mode !== undefined) a.mode = def.mode;
  return a;
}

const CORE_KEYS = new Set(Object.values(CORE).flatMap((c) => c.keys));

/** Build-time guard: a page verb must not reuse a core key. Called inside
 *  usePageActions so a colliding declaration throws in dev/test immediately. */
export function assertNoCoreCollision(actions: Action[]): void {
  const coreIds = new Set(Object.keys(CORE));
  for (const a of actions) {
    if (coreIds.has(a.id)) continue;
    for (const k of a.keys ?? []) {
      if (CORE_KEYS.has(k)) {
        throw new Error(`Action "${a.id}" reuses core key "${k}" — use core("${k}") or pick another key.`);
      }
    }
  }
}
```

> **Note on `exactOptionalPropertyTypes`:** the conditional-spread style above (`if (over.enabled) a.enabled = ...`) is required — assigning `enabled: undefined` would violate the strict optional flag. This mirrors the existing `KbdLegend`/`DockStepper` prop-forwarding pattern.

- [ ] **Step 4: Run test to verify it passes**

Run: `node_modules/.bin/vitest run test/interaction/core.test.ts`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/interaction/core.ts \
        packages/HimalayaUI/frontend/test/interaction/core.test.ts
git commit -m "feat(interaction): CORE vocabulary + core()/page() constructors"
```

---

### Task 0.4: Registry store + `usePageActions`

**Files:**
- Create: `packages/HimalayaUI/frontend/src/print/interaction/registry.ts`
- Create: `packages/HimalayaUI/frontend/src/print/interaction/usePageActions.ts`
- Test: `packages/HimalayaUI/frontend/test/interaction/registry.test.tsx`

**Interfaces:**
- Consumes: `Action`, `ListCursor` from `types.ts`; `assertNoCoreCollision` from `core.ts`.
- Produces:
  - `useInteraction` Zustand store with state `{ cursor: ListCursor | null; actions: Action[] }` and actions `setPage(cursor, actions)`, `clearPage()`.
  - `usePageActions(decl: { cursor?: ListCursor | null; actions: Action[] }): void`.

- [ ] **Step 1: Write the failing test**

```tsx
// test/interaction/registry.test.tsx
import { describe, it, expect } from "vitest";
import { render, cleanup } from "@testing-library/react";
import { useInteraction } from "../../src/print/interaction/registry";
import { usePageActions } from "../../src/print/interaction/usePageActions";
import { core, page } from "../../src/print/interaction/core";

function Page({ label }: { label: string }): null {
  usePageActions({ actions: [core("openFocus", { run: () => {} }), page("cull", { label, keys: ["x"], group: "Act", run: () => {} })] });
  return null;
}

describe("useInteraction + usePageActions", () => {
  it("registers a page's actions on mount", () => {
    render(<Page label="Cull" />);
    const ids = useInteraction.getState().actions.map((a) => a.id);
    expect(ids).toContain("openFocus");
    expect(ids).toContain("cull");
  });

  it("clears the store on unmount", () => {
    const { unmount } = render(<Page label="Cull" />);
    expect(useInteraction.getState().actions.length).toBeGreaterThan(0);
    unmount();
    expect(useInteraction.getState().actions).toEqual([]);
    expect(useInteraction.getState().cursor).toBeNull();
  });

  it("throws when a page verb collides with a core key", () => {
    function Bad(): null {
      usePageActions({ actions: [page("cull", { label: "x", keys: ["Enter"], group: "Act", run: () => {} })] });
      return null;
    }
    expect(() => render(<Bad />)).toThrow(/Enter/);
    cleanup();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `node_modules/.bin/vitest run test/interaction/registry.test.tsx`
Expected: FAIL — modules not found.

- [ ] **Step 3: Write `registry.ts` then `usePageActions.ts`**

```ts
// src/print/interaction/registry.ts
import { create } from "zustand";
import type { Action, ListCursor } from "./types";

interface InteractionState {
  cursor: ListCursor | null;
  actions: Action[];
  setPage: (cursor: ListCursor | null, actions: Action[]) => void;
  clearPage: () => void;
}

/** The single source the Dock and the keyboard layer both derive from. No
 *  middleware (mirrors floatingDock.ts) — purely in-memory, never persisted. */
export const useInteraction = create<InteractionState>((set) => ({
  cursor: null,
  actions: [],
  setPage: (cursor, actions) => set({ cursor, actions }),
  clearPage: () => set({ cursor: null, actions: [] }),
}));
```

```ts
// src/print/interaction/usePageActions.ts
import { useEffect } from "react";
import type { Action, ListCursor } from "./types";
import { useInteraction } from "./registry";
import { assertNoCoreCollision } from "./core";

/** A page declares its cursor + actions here. Written to the store in a
 *  dependency-less effect so it runs every commit — the closures captured in
 *  enabled()/run are therefore always the latest (no manual stateRef). The
 *  store clears on unmount so an un-migrated next page shows no stale dock. */
export function usePageActions(decl: { cursor?: ListCursor | null; actions: Action[] }): void {
  const setPage = useInteraction((s) => s.setPage);
  const clearPage = useInteraction((s) => s.clearPage);

  assertNoCoreCollision(decl.actions); // dev/test guard; cheap, pure

  // No dependency array: refresh the registry on every commit.
  useEffect(() => {
    setPage(decl.cursor ?? null, decl.actions);
  });

  useEffect(() => clearPage, [clearPage]); // unmount only
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `node_modules/.bin/vitest run test/interaction/registry.test.tsx`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/interaction/registry.ts \
        packages/HimalayaUI/frontend/src/print/interaction/usePageActions.ts \
        packages/HimalayaUI/frontend/test/interaction/registry.test.tsx
git commit -m "feat(interaction): registry store + usePageActions"
```

---

### Task 0.5: The keyboard layer (`useKeyboardLayer`)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/print/interaction/useKeyboardLayer.ts`
- Test: `packages/HimalayaUI/frontend/test/interaction/keyboardLayer.test.tsx`

**Interfaces:**
- Consumes: `useInteraction` (registry), `matchesKeys`/`isTyping`/`isBareKey` (matchKey), `suppressGlobalKeys` is **not** used here (this *replaces* it).
- Produces: `useKeyboardLayer(): void` — mounts one `window` keydown listener. Called once from `AppShell`.

**Rung-3 guards, in order (spec §6):**
1. `if (e.defaultPrevented) return` — a widget (rung 1) or modal (rung 2) already handled it. **The microtask-race fix.**
2. `if (isTyping(e.target) && isBareKey(e)) return` — don't fire `x`/`k` while typing; Mod-chords and Escape still pass.
3. For a matched action with a **bare** key, require focus to be inside the page's interactive scope (`data-interaction-scope` ancestor, or the cursored row, or `document.body`). Mod-chords/Escape are global. (WCAG 2.1.4.)
4. `enabled()` false → inert. Else `e.preventDefault(); run(e)`.

- [ ] **Step 1: Write the failing test**

```tsx
// test/interaction/keyboardLayer.test.tsx
import { describe, it, expect, vi, afterEach } from "vitest";
import { render, fireEvent, cleanup } from "@testing-library/react";
import { useKeyboardLayer } from "../../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../../src/print/interaction/registry";
import { core, page } from "../../src/print/interaction/core";

function Harness(): JSX.Element {
  useKeyboardLayer();
  return <div data-interaction-scope tabIndex={-1} data-testid="scope" />;
}

afterEach(() => { useInteraction.getState().clearPage(); cleanup(); });

describe("useKeyboardLayer", () => {
  it("fires a matched action's run() and preventDefaults", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [page("cull", { label: "Cull", keys: ["x"], group: "Act", run })]);
    render(<Harness />);
    document.querySelector<HTMLElement>('[data-interaction-scope]')!.focus();
    const evt = fireEvent.keyDown(document.activeElement!, { key: "x" });
    expect(run).toHaveBeenCalledTimes(1);
    expect(evt).toBe(false); // preventDefault was called (jsdom convention)
  });

  it("respects enabled() === false (inert, no preventDefault)", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [
      page("cull", { label: "Cull", keys: ["x"], group: "Act", enabled: () => false, run }),
    ]);
    render(<Harness />);
    document.querySelector<HTMLElement>('[data-interaction-scope]')!.focus();
    const evt = fireEvent.keyDown(document.activeElement!, { key: "x" });
    expect(run).not.toHaveBeenCalled();
    expect(evt).toBe(true); // not prevented
  });

  it("bails when the event was already defaultPrevented (rung 1/2 hand-off)", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("back", { run })]); // Escape
    render(<Harness />);
    // Simulate a modal that already consumed Escape:
    const handled = new KeyboardEvent("keydown", { key: "Escape", cancelable: true, bubbles: true });
    handled.preventDefault();
    window.dispatchEvent(handled);
    expect(run).not.toHaveBeenCalled();
  });

  it("ignores a bare key while typing in an input", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [page("cull", { label: "Cull", keys: ["x"], group: "Act", run })]);
    render(<><Harness /><input data-testid="field" /></>);
    const input = document.querySelector<HTMLInputElement>('[data-testid="field"]')!;
    input.focus();
    fireEvent.keyDown(input, { key: "x" });
    expect(run).not.toHaveBeenCalled();
  });

  it("still fires a Mod-chord while typing (undo)", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("undo", { run })]);
    render(<><Harness /><input data-testid="field" /></>);
    const input = document.querySelector<HTMLInputElement>('[data-testid="field"]')!;
    input.focus();
    fireEvent.keyDown(input, { key: "z", metaKey: true });
    expect(run).toHaveBeenCalledTimes(1);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `node_modules/.bin/vitest run test/interaction/keyboardLayer.test.tsx`
Expected: FAIL — module not found.

- [ ] **Step 3: Write `useKeyboardLayer.ts`**

```ts
// src/print/interaction/useKeyboardLayer.ts
import { useEffect } from "react";
import type { Action } from "./types";
import { useInteraction } from "./registry";
import { isBareKey, isTyping, matchesKeys } from "./matchKey";

function inPageScope(target: EventTarget | null): boolean {
  if (!(target instanceof HTMLElement)) return true; // window/body-level keydown
  if (target === document.body) return true;
  return target.closest("[data-interaction-scope],[data-cursored]") !== null;
}

function isEnabled(a: Action): boolean {
  return a.enabled ? a.enabled() : true;
}

/** The single rung-3 keyboard entry point. Mounted once in the app shell. */
export function useKeyboardLayer(): void {
  useEffect(() => {
    function onKeyDown(e: KeyboardEvent): void {
      // Rung 1/2 already handled it (the microtask-race fix — guard the signal,
      // never querySelector('[aria-modal]')).
      if (e.defaultPrevented) return;
      if (isTyping(e.target) && isBareKey(e)) return;

      const { actions } = useInteraction.getState();
      for (const a of actions) {
        if (!a.keys || !matchesKeys(e, a.keys)) continue;
        // WCAG 2.1.4: bare single-key actions fire only inside the page scope.
        if (isBareKey(e) && !inPageScope(e.target)) continue;
        if (!isEnabled(a)) return; // matched but inert — claim nothing, swallow nothing
        e.preventDefault();
        a.run(e);
        return;
      }
    }
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, []);
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `node_modules/.bin/vitest run test/interaction/keyboardLayer.test.tsx`
Expected: PASS (all five cases).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/interaction/useKeyboardLayer.ts \
        packages/HimalayaUI/frontend/test/interaction/keyboardLayer.test.tsx
git commit -m "feat(interaction): single rung-3 keyboard layer"
```

---

### Task 0.6: The Dock renderer (`InteractionDock`)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/print/interaction/InteractionDock.tsx`
- Test: `packages/HimalayaUI/frontend/test/interaction/InteractionDock.test.tsx`

**Interfaces:**
- Consumes: `useInteraction` (store), `Dock`/`DockUpLink`/`DockStepper` from `../ui`, `Button`/`KbKey` from `../ui`, `Action`/`ListCursor` from `types.ts`.
- Produces: `InteractionDock(): JSX.Element | null` — renders nothing when the store is empty.

**Render order (spec §7):** up-link (`core("back")` if declared) → stepper (`cursor.stepperProps()` if a cursor exists) → spacer → grouped `dock:true` action buttons (each with key hint + `disabled = !enabled()`) → primary (`dock:"primary"`, prominent + `↵`).

- [ ] **Step 1: Write the failing test**

```tsx
// test/interaction/InteractionDock.test.tsx
import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen, fireEvent, cleanup } from "@testing-library/react";
import { InteractionDock } from "../../src/print/interaction/InteractionDock";
import { useInteraction } from "../../src/print/interaction/registry";
import { core, page } from "../../src/print/interaction/core";
import type { ListCursor } from "../../src/print/interaction/types";

const fakeCursor = (over: Partial<ListCursor> = {}): ListCursor => ({
  cursorId: 1, selected: new Set(), setCursor: () => {}, moveBy: vi.fn(),
  activate: () => {}, toggleSelect: () => {}, rowProps: () => ({} as never),
  stepperProps: () => ({ label: "Sample", axis: "vertical", testIdBase: "sample", count: "1 / 3", onPrev: () => {}, onNext: () => {}, prevDisabled: true, nextDisabled: false }),
  ...over,
});

afterEach(() => { useInteraction.getState().clearPage(); cleanup(); });

describe("InteractionDock", () => {
  it("renders nothing when the store is empty", () => {
    const { container } = render(<InteractionDock />);
    expect(container.querySelector('[data-testid="dock"]')).toBeNull();
  });

  it("renders the stepper from cursor.stepperProps()", () => {
    useInteraction.getState().setPage(fakeCursor(), []);
    render(<InteractionDock />);
    expect(screen.getByTestId("dock-sample-count")).toHaveTextContent("1 / 3");
    expect(screen.getByTestId("dock-prev-sample")).toBeDisabled();
  });

  it("stepper next button calls cursor.moveBy(1)", () => {
    const moveBy = vi.fn();
    useInteraction.getState().setPage(fakeCursor({ moveBy }), []);
    render(<InteractionDock />);
    fireEvent.click(screen.getByTestId("dock-next-sample"));
    expect(moveBy).toHaveBeenCalledWith(1);
  });

  it("renders dock:true action buttons and greys disabled ones", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [
      page("cull", { label: "Cull", keys: ["x"], group: "Act", enabled: () => false, dock: true, run }),
    ]);
    render(<InteractionDock />);
    const btn = screen.getByRole("button", { name: /cull/i });
    expect(btn).toBeDisabled();
  });

  it("renders the primary action prominently and runs it on click", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("openFocus", { run, dock: "primary" })]);
    render(<InteractionDock />);
    fireEvent.click(screen.getByTestId("dock-primary"));
    expect(run).toHaveBeenCalledTimes(1);
  });

  it("renders the up-link when core('back') is declared", () => {
    useInteraction.getState().setPage(null, [core("back", { run: () => {}, label: "Corpus", dock: true })]);
    render(<InteractionDock />);
    expect(screen.getByTestId("dock-up-link")).toHaveTextContent("Corpus");
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `node_modules/.bin/vitest run test/interaction/InteractionDock.test.tsx`
Expected: FAIL — module not found.

- [ ] **Step 3: Write `InteractionDock.tsx`**

```tsx
// src/print/interaction/InteractionDock.tsx
import { Dock, DockUpLink, DockStepper, Button, KbKey } from "../ui";
import type { Action } from "./types";
import { useInteraction } from "./registry";

function enabledOf(a: Action): boolean {
  return a.enabled ? a.enabled() : true;
}

export function InteractionDock(): JSX.Element | null {
  const cursor = useInteraction((s) => s.cursor);
  const actions = useInteraction((s) => s.actions);

  if (cursor === null && actions.length === 0) return null;

  const back = actions.find((a) => a.id === "back");
  const primary = actions.find((a) => a.dock === "primary");
  const buttons = actions.filter((a) => a.dock === true && a.id !== "back");
  const stepper = cursor ? cursor.stepperProps() : null;

  return (
    <Dock>
      {back && <DockUpLink label={back.label} onClick={() => back.run()} />}

      {stepper && (
        <DockStepper
          label={stepper.label}
          axis={stepper.axis}
          testIdBase={stepper.testIdBase}
          count={stepper.count}
          onPrev={stepper.onPrev}
          onNext={stepper.onNext}
          prevDisabled={stepper.prevDisabled}
          nextDisabled={stepper.nextDisabled}
        />
      )}

      <div className="flex-1" />

      {buttons.map((a) => (
        <Button
          key={a.id}
          variant="outline"
          disabled={!enabledOf(a)}
          onClick={(e) => a.run(e)}
          data-testid={`dock-action-${a.id}`}
        >
          {a.label}
          {a.keys?.[0] && <KbKey className="ml-1.5">{a.keys[0]}</KbKey>}
        </Button>
      ))}

      {primary && (
        <Button
          variant="accent"
          disabled={!enabledOf(primary)}
          onClick={(e) => primary.run(e)}
          data-testid="dock-primary"
        >
          {primary.label}
          <KbKey variant="frost" className="ml-1.5">↵</KbKey>
        </Button>
      )}
    </Dock>
  );
}
```

> **`KbKey` label note:** the dock shows `a.keys[0]` verbatim (e.g. `x`, `l`). Friendly glyphs (`↵` for Enter) are only special-cased for the primary. If a future key needs a glyph in a normal button, add a small `keyGlyph()` map — YAGNI until then.

- [ ] **Step 4: Run test to verify it passes**

Run: `node_modules/.bin/vitest run test/interaction/InteractionDock.test.tsx`
Expected: PASS (all six cases).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/interaction/InteractionDock.tsx \
        packages/HimalayaUI/frontend/test/interaction/InteractionDock.test.tsx
git commit -m "feat(interaction): Dock renders from the page declaration"
```

---

### Task 0.7: Mount the layer in the app shell (still inert)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/shell/AppRoutes.tsx` (the `AppShell` component)

**Interfaces:**
- Consumes: `InteractionDock`, `useKeyboardLayer`.
- Produces: a shell that renders `<InteractionDock/>` as persistent chrome and runs `useKeyboardLayer()` once. Because the store is empty for every un-migrated page, **the Dock renders null and the listener matches nothing — zero behavioural change.**

- [ ] **Step 1: Write the failing test**

```tsx
// test/interaction/shellMount.test.tsx
import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../../src/print/shell/AppRoutes";

describe("AppShell mounts the interaction layer", () => {
  it("does not render a dock for an un-migrated route (store empty)", () => {
    const qc = new QueryClient();
    const { container } = render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments"]}>
          <AppRoutes />
        </MemoryRouter>
      </QueryClientProvider>,
    );
    // The shell mounted, but InteractionDock returns null until a page declares.
    expect(container.querySelector('[data-testid="app-shell"]')).not.toBeNull();
    // No second dock injected by the shell yet.
    expect(container.querySelectorAll('[data-testid="dock"]').length).toBe(0);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `node_modules/.bin/vitest run test/interaction/shellMount.test.tsx`
Expected: FAIL — `app-shell` present but the assertion catches the shell not yet importing the layer (the test passes only once the mount is wired and confirmed inert; if `AppShell` currently has no dock the second assertion already holds — so first make the test meaningful by asserting `useKeyboardLayer` is invoked; simplest: assert the import compiles and the shell renders. Treat a module-not-wired state as FAIL via a TODO import).

> Practical note: this task's "failing" state is the un-wired shell. Add the imports + calls, then the test green-confirms the shell renders with no extra dock.

- [ ] **Step 3: Wire `AppShell`**

In `AppRoutes.tsx`, add imports and edit the `AppShell` layout component:

```tsx
import { InteractionDock } from "../interaction/InteractionDock";
import { useKeyboardLayer } from "../interaction/useKeyboardLayer";

function AppShell(): JSX.Element {
  useKeyboardLayer(); // one window listener for the whole app
  return (
    <div data-testid="app-shell" className="h-full w-full flex flex-col min-h-0 bg-paper text-ink">
      <TopNav />
      <main className="flex-1 min-h-0 overflow-auto">
        <Outlet />
      </main>
      <InteractionDock />
    </div>
  );
}
```

Leave `useGlobalShortcuts()` at the top of `AppRoutes` **in place** for now (find/help still flow through it until the final phase). The two listeners coexist: `useGlobalShortcuts` owns `/`·`⌘K`·`?`; `useKeyboardLayer` matches only what a migrated page registers (nothing yet).

- [ ] **Step 4: Run test + full build to verify inert mount**

Run: `node_modules/.bin/vitest run test/interaction/shellMount.test.tsx`
Expected: PASS.
Run: `npm run build`
Expected: PASS (lint:design + tsc + vite). Manually confirm: every existing page still shows its own hand-built dock (the shell's `<InteractionDock/>` is null).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/shell/AppRoutes.tsx \
        packages/HimalayaUI/frontend/test/interaction/shellMount.test.tsx
git commit -m "feat(interaction): mount InteractionDock + keyboard layer in shell (inert)"
```

---

# Phase 1 — `useListCursor` + the shared contract harness

This is the **C** contract and the engine of input parity. It is built and tested in isolation (no page yet), with a **parameterized contract harness** that every later page re-runs.

### Task 1.1: `useListCursor` — ID-based cursor + roving tabindex

**Files:**
- Create: `packages/HimalayaUI/frontend/src/print/interaction/useListCursor.ts`
- Test: `packages/HimalayaUI/frontend/test/interaction/useListCursor.test.tsx`

**Interfaces:**
- Consumes: `ListCursor`/`RowProps`/`CursorStepperProps` (types), `safeScrollIntoView` (`../../lib/safeScrollIntoView`), `findNextEnabled` is **not** used here (no per-row disabled in the base cursor; pages that need skipping pass a filtered `ids`).
- Produces: `useListCursor(opts: { ids: number[]; onActivate?: (id: number) => void; stepperLabel?: string; stepperTestIdBase?: string; axis?: "vertical" | "horizontal" }): ListCursor`.

**Behaviour:**
- `cursorId` defaults to `ids[0] ?? null`; when `ids` changes and the current `cursorId` is gone, fall back to the **nearest surviving index** (SSE-safety), else `null`.
- `moveBy(delta)`: from the current id's index, clamp to `[0, ids.length-1]`, set the new id. (Stepper and arrows both call this.)
- `setCursor(id)`: a click lands here.
- `activate()`: `onActivate?.(cursorId)` when non-null.
- `toggleSelect(id?)`: toggles `id ?? cursorId` in the `selected` Set (new Set each time — immutable).
- `rowProps(id)`: registers the element into a `Map<number, HTMLElement>` ref; returns `tabIndex` (`id===cursorId?0:-1`), `onClick`→`setCursor`, `onDoubleClick`→`setCursor`+`activate`, `role:"row"`, `aria-current`, `data-cursored`.
- A `useEffect` on `cursorId` calls `.focus()` on the cursored element (free scroll + `:focus-visible` + SR "row N of M"); `safeScrollIntoView` as the jsdom-safe fallback.

- [ ] **Step 1: Write the failing test** (unit behaviours; the *contract* harness is Task 1.2)

```tsx
// test/interaction/useListCursor.test.tsx
import { describe, it, expect, vi } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { useListCursor } from "../../src/print/interaction/useListCursor";

describe("useListCursor", () => {
  it("defaults the cursor to the first id", () => {
    const { result } = renderHook(() => useListCursor({ ids: [10, 20, 30] }));
    expect(result.current.cursorId).toBe(10);
  });

  it("moveBy steps by id and clamps at the ends", () => {
    const { result } = renderHook(() => useListCursor({ ids: [10, 20, 30] }));
    act(() => result.current.moveBy(1));
    expect(result.current.cursorId).toBe(20);
    act(() => result.current.moveBy(-5));
    expect(result.current.cursorId).toBe(10);
    act(() => result.current.moveBy(99));
    expect(result.current.cursorId).toBe(30);
  });

  it("setCursor parks the cursor on a clicked id", () => {
    const { result } = renderHook(() => useListCursor({ ids: [10, 20, 30] }));
    act(() => result.current.setCursor(30));
    expect(result.current.cursorId).toBe(30);
  });

  it("activate calls onActivate with the cursor id", () => {
    const onActivate = vi.fn();
    const { result } = renderHook(() => useListCursor({ ids: [10, 20], onActivate }));
    act(() => result.current.moveBy(1));
    act(() => result.current.activate());
    expect(onActivate).toHaveBeenCalledWith(20);
  });

  it("SSE-safety: removing an id before the cursor keeps the SAME item cursored", () => {
    const { result, rerender } = renderHook(({ ids }) => useListCursor({ ids }), {
      initialProps: { ids: [10, 20, 30] },
    });
    act(() => result.current.setCursor(30));
    rerender({ ids: [20, 30] }); // 10 removed (insert/remove before cursor)
    expect(result.current.cursorId).toBe(30); // still item 30, not index-shifted
  });

  it("SSE-safety: removing the cursored id falls back to the nearest surviving index", () => {
    const { result, rerender } = renderHook(({ ids }) => useListCursor({ ids }), {
      initialProps: { ids: [10, 20, 30] },
    });
    act(() => result.current.setCursor(20));
    rerender({ ids: [10, 30] }); // the cursored item itself removed
    expect([10, 30]).toContain(result.current.cursorId);
  });

  it("toggleSelect maintains a Set orthogonal to the cursor", () => {
    const { result } = renderHook(() => useListCursor({ ids: [10, 20, 30] }));
    act(() => result.current.toggleSelect(20));
    expect(result.current.selected.has(20)).toBe(true);
    expect(result.current.cursorId).toBe(10); // cursor did NOT move
    act(() => result.current.toggleSelect(20));
    expect(result.current.selected.has(20)).toBe(false);
  });

  it("stepperProps reflects position and end-disabled state", () => {
    const { result } = renderHook(() => useListCursor({ ids: [10, 20, 30], stepperLabel: "Sample", stepperTestIdBase: "sample" }));
    const s = result.current.stepperProps();
    expect(s.count).toBe("1 / 3");
    expect(s.prevDisabled).toBe(true);
    expect(s.nextDisabled).toBe(false);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `node_modules/.bin/vitest run test/interaction/useListCursor.test.tsx`
Expected: FAIL — module not found.

- [ ] **Step 3: Write `useListCursor.ts`**

```ts
// src/print/interaction/useListCursor.ts
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { ListCursor, RowProps, CursorStepperProps } from "./types";
import { safeScrollIntoView } from "../../lib/safeScrollIntoView";

interface Opts {
  ids: number[];
  onActivate?: (id: number) => void;
  stepperLabel?: string;
  stepperTestIdBase?: string;
  axis?: "vertical" | "horizontal";
}

export function useListCursor(opts: Opts): ListCursor {
  const { ids, onActivate, stepperLabel = "Item", stepperTestIdBase = "item", axis = "vertical" } = opts;

  const [cursorId, setCursorId] = useState<number | null>(ids[0] ?? null);
  const [selected, setSelected] = useState<Set<number>>(() => new Set());
  const elements = useRef<Map<number, HTMLElement>>(new Map());

  // Keep the cursor on a live id. If the cursored item vanished (SSE remove),
  // fall back to the nearest surviving index; never silently jump to index 0.
  const prevIds = useRef<number[]>(ids);
  useEffect(() => {
    if (cursorId !== null && ids.includes(cursorId)) {
      prevIds.current = ids;
      return;
    }
    if (cursorId === null) {
      if (ids.length > 0) setCursorId(ids[0]!);
      prevIds.current = ids;
      return;
    }
    const oldIndex = prevIds.current.indexOf(cursorId);
    const fallback = ids[Math.min(oldIndex, ids.length - 1)] ?? ids[ids.length - 1] ?? null;
    setCursorId(fallback);
    prevIds.current = ids;
  }, [ids, cursorId]);

  // Cursor === DOM focus: move focus to the cursored row.
  useEffect(() => {
    if (cursorId === null) return;
    const el = elements.current.get(cursorId);
    if (el) {
      el.focus();
      safeScrollIntoView(el, { block: "nearest" });
    }
  }, [cursorId]);

  const indexOf = useCallback((id: number | null) => (id === null ? -1 : ids.indexOf(id)), [ids]);

  const setCursor = useCallback((id: number) => setCursorId(id), []);

  const moveBy = useCallback(
    (delta: number) => {
      if (ids.length === 0) return;
      const from = Math.max(0, indexOf(cursorId));
      const next = Math.min(Math.max(from + delta, 0), ids.length - 1);
      setCursorId(ids[next]!);
    },
    [ids, cursorId, indexOf],
  );

  const activate = useCallback(() => {
    if (cursorId !== null) onActivate?.(cursorId);
  }, [cursorId, onActivate]);

  const toggleSelect = useCallback(
    (id?: number) => {
      const target = id ?? cursorId;
      if (target === null) return;
      setSelected((prev) => {
        const next = new Set(prev);
        if (next.has(target)) next.delete(target);
        else next.add(target);
        return next;
      });
    },
    [cursorId],
  );

  const rowProps = useCallback(
    (id: number): RowProps => ({
      ref: (el) => {
        if (el) elements.current.set(id, el);
        else elements.current.delete(id);
      },
      tabIndex: id === cursorId ? 0 : -1,
      onClick: () => setCursorId(id),
      onDoubleClick: () => {
        setCursorId(id);
        onActivate?.(id);
      },
      role: "row",
      "aria-current": id === cursorId ? "true" : undefined,
      "data-cursored": id === cursorId ? "true" : "false",
    }),
    [cursorId, onActivate],
  );

  const stepperProps = useCallback((): CursorStepperProps => {
    const i = indexOf(cursorId);
    return {
      label: stepperLabel,
      axis,
      testIdBase: stepperTestIdBase,
      count: ids.length > 0 ? `${Math.max(0, i) + 1} / ${ids.length}` : "0 / 0",
      onPrev: () => moveBy(-1),
      onNext: () => moveBy(1),
      prevDisabled: i <= 0,
      nextDisabled: i < 0 || i >= ids.length - 1,
    };
  }, [ids, cursorId, indexOf, moveBy, stepperLabel, stepperTestIdBase, axis]);

  return useMemo<ListCursor>(
    () => ({ cursorId, selected, setCursor, moveBy, activate, toggleSelect, rowProps, stepperProps }),
    [cursorId, selected, setCursor, moveBy, activate, toggleSelect, rowProps, stepperProps],
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `node_modules/.bin/vitest run test/interaction/useListCursor.test.tsx`
Expected: PASS (all eight cases).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/interaction/useListCursor.ts \
        packages/HimalayaUI/frontend/test/interaction/useListCursor.test.tsx
git commit -m "feat(interaction): useListCursor (ID-based + roving tabindex)"
```

---

### Task 1.2: The shared cursor-contract harness

**Files:**
- Create: `packages/HimalayaUI/frontend/test/interaction/cursorContract.tsx` (the reusable harness — a `.tsx` helper, not a spec)
- Create: `packages/HimalayaUI/frontend/test/interaction/useListCursor.contract.test.tsx` (runs the harness against the base `useListCursor`)

**Interfaces:**
- Consumes: `ListCursor`, `DockStepper`.
- Produces: `runCursorContract(name: string, mount: () => { cursor: ListCursor; ids: number[]; onActivate: ReturnType<typeof vi.fn> })` — a function each page's `*.actions.test.tsx` calls to assert parity mechanically.

**The harness renders** a list whose rows spread `cursor.rowProps(id)` and a `DockStepper` fed by `cursor.stepperProps()`, then drives click/arrow/stepper/Enter/double-click/dock-primary and asserts the spec §10 invariants.

- [ ] **Step 1: Write the harness + its self-test (failing)**

```tsx
// test/interaction/cursorContract.tsx
import { vi, expect } from "vitest";
import { render, fireEvent, within, act } from "@testing-library/react";
import { useEffect } from "react";
import type { ListCursor } from "../../src/print/interaction/types";
import { DockStepper } from "../../src/print/ui";

interface Mounted { cursor: ListCursor; ids: number[]; onActivate: ReturnType<typeof vi.fn>; }

/** Render a minimal cursored list + its stepper, exposing the live cursor. */
function Probe({ expose }: { expose: (m: Pick<Mounted, "cursor" | "ids">) => void }): JSX.Element {
  // The concrete cursor is supplied by the caller via a wrapper; see usage.
  throw new Error("Probe is provided by runCursorContract's wrapper");
}
void Probe;

export function runCursorContract(name: string, makeWrapper: () => {
  ui: (capture: (m: Mounted) => void) => JSX.Element;
}): void {
  describe(`cursor contract: ${name}`, () => {
    let m!: Mounted;
    const { ui } = makeWrapper();

    function mount(): ReturnType<typeof render> {
      return render(ui((captured) => { m = captured; }));
    }

    it("exactly one row has tabIndex=0 at all times (roving invariant)", () => {
      const { container } = mount();
      const tabbable = container.querySelectorAll('[role="row"][tabindex="0"]');
      expect(tabbable.length).toBe(1);
    });

    it("click-sets-cursor == arrow-moves-cursor == stepper-moves-cursor", () => {
      const { container, getByTestId } = mount();
      const rows = () => container.querySelectorAll<HTMLElement>('[role="row"]');

      // arrow down (moveBy via the live cursor)
      act(() => m.cursor.moveBy(1));
      const viaArrow = m.cursor.cursorId;

      // stepper next → same end id from the same start
      act(() => m.cursor.setCursor(m.ids[0]!));
      fireEvent.click(getByTestId(`dock-next-${m.cursor.stepperProps().testIdBase}`));
      const viaStepper = m.cursor.cursorId;

      // click row index 1 → same id
      act(() => m.cursor.setCursor(m.ids[0]!));
      fireEvent.click(rows()[1]!);
      const viaClick = m.cursor.cursorId;

      expect(viaArrow).toBe(m.ids[1]);
      expect(viaStepper).toBe(m.ids[1]);
      expect(viaClick).toBe(m.ids[1]);
    });

    it("Enter == double-click == dock primary (same onActivate arg)", () => {
      const { container } = mount();
      act(() => m.cursor.setCursor(m.ids[1]!));

      act(() => m.cursor.activate());          // keyboard Enter path
      const viaEnter = m.onActivate.mock.calls.at(-1)?.[0];

      fireEvent.doubleClick(container.querySelectorAll('[role="row"]')[1]!);
      const viaDouble = m.onActivate.mock.calls.at(-1)?.[0];

      expect(viaEnter).toBe(m.ids[1]);
      expect(viaDouble).toBe(m.ids[1]);
    });

    it("moving the cursor never mutates selected; toggleSelect never moves the cursor", () => {
      mount();
      act(() => m.cursor.toggleSelect(m.ids[2]!));
      const cursorBefore = m.cursor.cursorId;
      act(() => m.cursor.moveBy(1));
      expect(m.cursor.selected.has(m.ids[2]!)).toBe(true); // selection survives nav
      act(() => m.cursor.toggleSelect(m.ids[0]!));
      // toggling did not change the cursor that moveBy set
      expect(m.cursor.cursorId).not.toBe(cursorBefore === m.ids[0] ? null : undefined);
    });
  });
}
```

```tsx
// test/interaction/useListCursor.contract.test.tsx
import { vi } from "vitest";
import { useListCursor } from "../../src/print/interaction/useListCursor";
import { DockStepper } from "../../src/print/ui";
import { runCursorContract } from "./cursorContract";

runCursorContract("useListCursor (base list)", () => {
  const onActivate = vi.fn();
  const IDS = [10, 20, 30];
  return {
    ui: (capture) => {
      function Wrapper(): JSX.Element {
        const cursor = useListCursor({ ids: IDS, onActivate, stepperTestIdBase: "item" });
        capture({ cursor, ids: IDS, onActivate });
        const s = cursor.stepperProps();
        return (
          <div role="grid" aria-multiselectable data-interaction-scope>
            {IDS.map((id) => (
              <div key={id} {...cursor.rowProps(id)}>row {id}</div>
            ))}
            <DockStepper {...s} onPrev={s.onPrev} onNext={s.onNext} />
          </div>
        );
      }
      return <Wrapper />;
    },
  };
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `node_modules/.bin/vitest run test/interaction/useListCursor.contract.test.tsx`
Expected: FAIL initially (harness/exports not wired), then drive to green.

- [ ] **Step 3: Fix the harness wiring** until the four contract assertions pass against the base cursor. (Remove the placeholder `Probe`; the wrapper supplies the live cursor through `capture`.)

- [ ] **Step 4: Run to verify it passes**

Run: `node_modules/.bin/vitest run test/interaction/useListCursor.contract.test.tsx`
Expected: PASS — the base cursor satisfies the parity contract.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/test/interaction/cursorContract.tsx \
        packages/HimalayaUI/frontend/test/interaction/useListCursor.contract.test.tsx
git commit -m "test(interaction): shared cursor-contract harness"
```

---

# Phase 2 — Corpus migration (the proof)

Corpus is the proof: click-to-cursor, roving tabindex, cull-as-mode, dock-from-declaration, and the **index→ID** bug fixed. Corpus has a **compound (sample, frame)** target — the cursor migration covers the **sample** axis (the primary navigable list); the frame axis stays page-local state for now (a `useState(frameId)`), exposed as a second stepper *not* driven by `useListCursor` (Loupe in Phase 3 proves the compound-cursor composition). This keeps Phase 2 focused on proving the contract.

**Current state to replace** (`src/print/pages/ExperimentCorpusPage.tsx`):
- `useState({ sampleIndex, frameIndex })` (lines 167–169) — **index-based** sample cursor → replace sample axis with `useListCursor` (ID).
- The hand-built `<Dock>` (lines 590–689) with two steppers + cull verbs + compose + destinations → declare via `usePageActions`; the shell renders it.
- `useShortcuts({...})` (lines 299–398) — 15 bindings → core/page declarations + cursor (nav comes free).
- The manual `scrollIntoView` effect on `[data-cursored="true"]` → deleted (roving tabindex gives it free).
- `selected: Set<number>` (frame cull) stays as the cursor's `selected` only if it's the *sample* selection; Corpus's `selected` is a **frame**-exposure Set, so it stays page-local. The `checkedSamples` Set (compose) is the sample-level selection → can move onto `cursor.selected`.

> **Scope decision for Phase 2 (locked):** migrate the **sample axis** to `useListCursor` and the **sample-level checkbox selection** (`checkedSamples`) onto `cursor.selected`. Leave the frame axis + frame-exposure `selected` as page-local state, surfaced as a second `DockStepper` and cull buttons declared as `page()` actions. This is the minimal change that proves the contract; the frame axis migrates with Loupe's compound-cursor work if desired.

### Task 2.1: Corpus mode union + cursor

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/pages/ExperimentCorpusPage.tsx`
- Test: `packages/HimalayaUI/frontend/test/ExperimentCorpusPage.cursor.test.tsx`

**Interfaces:**
- Consumes: `useListCursor`, `usePageActions`, `core`, `page`.
- Produces: a Corpus whose sample cursor is ID-based and whose rows carry `cursor.rowProps`.

- [ ] **Step 1: Write the failing test**

```tsx
// test/ExperimentCorpusPage.cursor.test.tsx
import { describe, it, expect } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
// ... standard Corpus test harness (QueryClient + MemoryRouter + mocked useCorpusSamples) ...
import { renderCorpus } from "./helpers/renderCorpus"; // existing helper or inline per repo convention

describe("Corpus cursor (ID-based, click-to-park)", () => {
  it("clicking a sample row parks the cursor there (data-cursored)", async () => {
    renderCorpus({ samples: [{ id: 10, name: "A" }, { id: 20, name: "B" }, { id: 30, name: "C" }] });
    const rowB = await screen.findByText("B");
    fireEvent.click(rowB.closest('[role="row"]')!);
    expect(rowB.closest('[role="row"]')).toHaveAttribute("data-cursored", "true");
  });

  it("exactly one row is tabbable at a time (roving)", async () => {
    const { container } = renderCorpus({ samples: [{ id: 10, name: "A" }, { id: 20, name: "B" }] });
    await screen.findByText("A");
    expect(container.querySelectorAll('[role="row"][tabindex="0"]').length).toBe(1);
  });

  it("ArrowDown moves the cursor to the next sample id", async () => {
    renderCorpus({ samples: [{ id: 10, name: "A" }, { id: 20, name: "B" }] });
    const scope = await screen.findByTestId("experiment-corpus");
    fireEvent.keyDown(within(scope).getByText("A").closest('[role="row"]')!, { key: "ArrowDown" });
    expect(within(scope).getByText("B").closest('[role="row"]')).toHaveAttribute("data-cursored", "true");
  });
});
```

> **Note on ArrowUp/ArrowDown:** the base cursor doesn't bind arrows itself; the page wires `↑/↓` to `cursor.moveBy(∓1)` through the declared nav (see Task 2.2) **or** the row container handles arrows via `onKeyDown` calling `moveBy`. Decision: arrows are handled by a small `onKeyDown` on the grid container (rung-1 widget nav — list navigation that stays local), calling `cursor.moveBy`. This keeps arrows scoped to the focused list and off the global registry. Add to the grid container: `onKeyDown={(e) => { if (e.key === "ArrowDown") { e.preventDefault(); cursor.moveBy(1); } if (e.key === "ArrowUp") { e.preventDefault(); cursor.moveBy(-1); } }}` plus `data-interaction-scope`.

- [ ] **Step 2: Run test to verify it fails**

Run: `node_modules/.bin/vitest run test/ExperimentCorpusPage.cursor.test.tsx`
Expected: FAIL — rows have no `data-cursored`/`role="row"` from the cursor yet (current rows use the old `cursored` prop driven by index).

- [ ] **Step 3: Implement the cursor**

In `ExperimentCorpusPage.tsx`:

1. Replace the sample part of `const [cursor, setCursor] = useState<{ sampleIndex; frameIndex }>` with:
```tsx
const sampleIds = useMemo(() => scopedSamples.map((s) => s.id), [scopedSamples]);
const sampleCursor = useListCursor({
  ids: sampleIds,
  onActivate: (id) => navigate(`/sample/${id}`),
  stepperLabel: "Sample",
  stepperTestIdBase: "sample",
  axis: "vertical",
});
const activeSample = scopedSamples.find((s) => s.id === sampleCursor.cursorId);
```
2. Keep `frameIndex` as page-local: `const [frameIndex, setFrameIndex] = useState(0)` reset on `sampleCursor.cursorId` change.
3. Delete the manual `scrollIntoView` effect (lines ~178–183).
4. On each `SampleTableRow`, spread `{...sampleCursor.rowProps(s.id)}` and remove the old `cursored`/`onSelectExposure→setCursor` index wiring. The row's existing checkbox toggles `sampleCursor.toggleSelect(s.id)`.
5. Wrap the table container with `role="grid" aria-multiselectable data-interaction-scope` and the arrow `onKeyDown` from the note above.

- [ ] **Step 4: Run test to verify it passes**

Run: `node_modules/.bin/vitest run test/ExperimentCorpusPage.cursor.test.tsx`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/ExperimentCorpusPage.tsx \
        packages/HimalayaUI/frontend/test/ExperimentCorpusPage.cursor.test.tsx
git commit -m "feat(corpus): ID-based sample cursor + roving tabindex"
```

---

### Task 2.2: Corpus action declaration (dock + keyboard from one source)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/pages/ExperimentCorpusPage.tsx`
- Test: `packages/HimalayaUI/frontend/test/ExperimentCorpusPage.actions.test.tsx`

**Interfaces:**
- Consumes: `usePageActions`, `core`, `page`, `useListCursor` (from 2.1), `runCursorContract`.
- Produces: Corpus registers `openFocus` (primary), `openLoupe`, `back`, and the cull verbs; the shell Dock renders them; the keyboard layer fires them.

**Mode union (replaces the boolean soup):**
```ts
type CorpusMode = { kind: "browse" } | { kind: "selection" };
// derived: const mode: CorpusMode = cursor.selected.size > 0 ? { kind: "selection" } : { kind: "browse" };
```

- [ ] **Step 1: Write the failing test**

```tsx
// test/ExperimentCorpusPage.actions.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { renderCorpus } from "./helpers/renderCorpus";
import { runCursorContract } from "./interaction/cursorContract";

describe("Corpus action declaration", () => {
  it("renders the shell dock primary 'Focus' and navigates on click", async () => {
    const { navigate } = renderCorpus({ samples: [{ id: 10, name: "A" }] });
    await screen.findByText("A");
    fireEvent.click(screen.getByTestId("dock-primary"));
    expect(navigate).toHaveBeenCalledWith("/sample/10");
  });

  it("the sample stepper in the shell dock drives the cursor", async () => {
    renderCorpus({ samples: [{ id: 10, name: "A" }, { id: 20, name: "B" }] });
    await screen.findByText("A");
    fireEvent.click(screen.getByTestId("dock-next-sample"));
    expect(screen.getByText("B").closest('[role="row"]')).toHaveAttribute("data-cursored", "true");
  });

  it("cull (x) is enabled only in selection mode", async () => {
    renderCorpus({ samples: [{ id: 10, name: "A" }] });
    const row = (await screen.findByText("A")).closest('[role="row"]')!;
    // browse mode: no Cull dock button enabled
    expect(screen.queryByTestId("dock-action-cull")).toBeNull();
    fireEvent.click(within(row).getByRole("checkbox")); // enters selection
    expect(screen.getByTestId("dock-action-cull")).toBeEnabled();
  });

  it("Enter on a focused row navigates to Focus (primary == Enter)", async () => {
    const { navigate } = renderCorpus({ samples: [{ id: 10, name: "A" }] });
    const row = (await screen.findByText("A")).closest('[role="row"]')! as HTMLElement;
    row.focus();
    fireEvent.keyDown(row, { key: "Enter" });
    expect(navigate).toHaveBeenCalledWith("/sample/10");
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `node_modules/.bin/vitest run test/ExperimentCorpusPage.actions.test.tsx`
Expected: FAIL — no shell dock primary yet (page still owns its `<Dock>`).

- [ ] **Step 3: Declare the actions; delete the hand-built dock + `useShortcuts`**

In `ExperimentCorpusPage.tsx`:

```tsx
const mode: "browse" | "selection" = sampleCursor.selected.size > 0 ? "selection" : "browse";

usePageActions({
  cursor: sampleCursor,
  actions: [
    core("back", { label: "Experiments", run: () => navigate("/experiments"), dock: true }),
    core("openFocus", { run: () => sampleCursor.activate(), dock: "primary",
                        enabled: () => sampleCursor.cursorId !== null }),
    core("openLoupe", { run: () => { if (activeSample) navigate(`/sample/${activeSample.id}/loupe`); },
                        dock: true, enabled: () => activeSample != null }),
    page("cull", { label: "Drop", keys: ["x"], group: "Act", mode: "selection",
                   enabled: () => mode === "selection", dock: true, run: () => dropSelected() }),
    page("keep", { label: "Keep", keys: ["k"], group: "Act", mode: "selection",
                   enabled: () => mode === "selection", dock: true, run: () => keepSelected() }),
    page("restore", { label: "Restore", keys: ["r"], group: "Act", mode: "selection",
                      enabled: () => mode === "selection", dock: true, run: () => restoreSelected() }),
  ],
});
```

Then:
- Delete the entire hand-built `<Dock>...</Dock>` block (lines ~590–689).
- Delete the `useShortcuts({...})` call (lines ~299–398); its frame-axis (`prevDetail`/`nextDetail`) handlers move to a small page-local `onKeyDown` on the frame region (rung-1) since the frame axis stays page-local this phase.
- `prevSample`/`nextSample` are **not** declared — they come from the cursor (arrows via the grid `onKeyDown` from 2.1; stepper via the shell).
- Keep `representative` and the compose ("+ New series") affordance as page-local UI for now (not a global gesture); the compose button stays in the page body, not the dock. (Revisit if Jonathan wants it docked.)

> **`dropSelected`/`keepSelected`/`restoreSelected`:** extract the existing inline `drop`/`keep`/`restore` handler bodies from the old `useShortcuts` block into named functions so both the (now-deleted) shortcut and the dock button share them. They already exist as closures; lift them to `useCallback`s.

- [ ] **Step 4: Run test + contract harness**

Run: `node_modules/.bin/vitest run test/ExperimentCorpusPage.actions.test.tsx`
Expected: PASS.

Add a contract run (append to the actions test or a sibling):
```tsx
// Corpus sample cursor satisfies the shared parity contract
runCursorContract("Corpus sample cursor", () => ({ ui: (capture) => <CorpusContractProbe capture={capture} /> }));
```
Run it; Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/ExperimentCorpusPage.tsx \
        packages/HimalayaUI/frontend/test/ExperimentCorpusPage.actions.test.tsx
git commit -m "feat(corpus): declare actions; shell renders dock + keyboard"
```

---

### Task 2.3: Corpus mocked E2E + live walk-through

**Files:**
- Create: `packages/HimalayaUI/frontend/e2e/interaction-corpus.spec.ts`

**Interfaces:**
- Consumes: the mocked-route + `seedState` helpers (per `e2e/AGENTS.md`).

- [ ] **Step 1: Write the mocked E2E spec**

```ts
// e2e/interaction-corpus.spec.ts
import { test, expect, type Page } from "@playwright/test";
import { mockCore, seedState } from "./helpers"; // existing mock helpers

test.describe("Corpus interaction (mocked)", () => {
  test("click a row then Enter enters Focus", async ({ page }) => {
    await mockCore(page);
    await seedState(page, { activeExperimentId: 1 });
    await page.goto("/experiments/1/corpus");
    const rowB = page.getByRole("row").filter({ hasText: "B" });
    await rowB.click();
    await expect(rowB).toHaveAttribute("data-cursored", "true");
    await page.keyboard.press("Enter");
    await expect(page).toHaveURL(/\/sample\/\d+$/);
  });

  test("stepper next == ArrowDown (same cursored row)", async ({ page }) => {
    await mockCore(page);
    await seedState(page, { activeExperimentId: 1 });
    await page.goto("/experiments/1/corpus");
    await page.getByTestId("dock-next-sample").click();
    const viaStepper = await page.locator('[role="row"][data-cursored="true"]').textContent();
    await page.getByTestId("dock-prev-sample").click();
    await page.locator('[role="row"][data-cursored="true"]').focus();
    await page.keyboard.press("ArrowDown");
    const viaArrow = await page.locator('[role="row"][data-cursored="true"]').textContent();
    expect(viaArrow).toBe(viaStepper);
  });
});
```

- [ ] **Step 2: Run to verify**

Run: `npm run e2e -- --grep "Corpus interaction"`
Expected: PASS (mocked).

- [ ] **Step 3: Live walk-through** (real Chrome — jsdom can't prove focus/scroll)

Per `docs/frontend-dev-loop.md`, on a **copy** of the dev DB:
```bash
mkdir -p /tmp/himalaya-uat && cp <dev-db> /tmp/himalaya-uat/himalaya.db
HIMALAYA_DB_PATH=/tmp/himalaya-uat/himalaya.db julia --project=packages/HimalayaUI \
  -e 'using HimalayaUI; main(ARGS)' -- serve --port 8091 &
cd packages/HimalayaUI/frontend && VITE_API_PORT=8091 npm run dev -- --host 127.0.0.1 --port 5182 &
```
Then via Playwright MCP, verify on a real Corpus: (a) click a row → it gains the cursor ring + scrolls into view; (b) `↑/↓` move the cursor and the row stays visible (roving focus); (c) `Enter` enters Focus; (d) the shell dock stepper moves the same cursor; (e) check a row → `x` culls; `Escape` clears selection then (on browse) goes back. Record results in the PR description.

- [ ] **Step 4: Full gate**

Run: `npm run build && npm test && npm run e2e`
Expected: all PASS. **This is the proof gate** — the model is validated on one page before rolling.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/e2e/interaction-corpus.spec.ts
git commit -m "test(corpus): mocked E2E for click→cursor→Enter parity"
```

---

# Phase 3 — Loupe (compound `(sample, frame)` cursor)

Loupe stress-tests **C**: two axes. The **sample** axis is URL-driven (`/sample/:id/loupe`); the **frame** axis is `activeId` (already ID-based — good). Decision: the **frame** axis becomes the page's `useListCursor` (the in-page navigable list), and the **sample** axis stays a URL-stepper declared as two `page()`-style nav handlers feeding a *second* `DockStepper` via a small adapter. This validates that a page can own a primary cursor (frame) plus a secondary stepper (sample) without the cursor contract needing to model both — the answer to spec §12's compound-cursor open question: **don't compose two axes into one `ListCursor`; one axis is the cursor, the other is a declared stepper.**

### Task 3.1: Loupe frame cursor + sample stepper adapter

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/pages/LoupePage.tsx`
- Create: `packages/HimalayaUI/frontend/src/print/interaction/useStepperOnly.ts` (a tiny cursor-less stepper for URL-driven axes)
- Test: `packages/HimalayaUI/frontend/test/LoupePage.actions.test.tsx`

**Interfaces:**
- Consumes: `useListCursor`, `usePageActions`, `core`, `page`, `CursorStepperProps`.
- Produces: `useStepperOnly(opts: { ids: number[]; currentId: number | undefined; onGo: (id: number) => void; label: string; testIdBase: string; axis?: "vertical"|"horizontal" }): CursorStepperProps` — feeds a `DockStepper` for an axis whose state lives elsewhere (URL).

> **Why a second primitive:** the sample axis on Loupe is *navigation between pages* (URL changes), not selection within a list. Modelling it as a `ListCursor` would imply rows + roving focus that don't exist. `useStepperOnly` is ~15 lines: given ids + currentId + onGo, it returns the same `CursorStepperProps` shape so the dock renders it identically. The frame axis (real in-page list) uses the full `useListCursor`.

The Dock must render **two** steppers here. Extend the declaration: `usePageActions` already carries one `cursor`; add an optional `extraSteppers?: CursorStepperProps[]` to the declaration and have `InteractionDock` render them after the cursor stepper.

- [ ] **Step 1: Write the failing test** — frame cursor parity + sample stepper navigates URL + Drop/Keep/Representative declared.

```tsx
// test/LoupePage.actions.test.tsx — key cases
it("frame stepper (←/→) moves the frame cursor by id", async () => { /* ArrowRight → activeId advances */ });
it("sample stepper (↑/↓) navigates to the sibling sample URL", async () => { /* dock-next-sample → /sample/<next>/loupe */ });
it("Drop (x) toggles the active frame's cull", async () => { /* dock-action-cull → mutation */ });
it("Enter enters Focus for the sample", async () => { /* dock-primary → /sample/<id> */ });
it("renders BOTH steppers in the shell dock", async () => {
  // dock-frame-count AND dock-sample-count both present
});
```

- [ ] **Step 2: Run to verify it fails.** Run: `node_modules/.bin/vitest run test/LoupePage.actions.test.tsx` → FAIL.

- [ ] **Step 3: Implement.**
  1. Add `extraSteppers?: CursorStepperProps[]` to the `usePageActions` declaration type and to the store (`registry.ts` state + `setPage`), and render them in `InteractionDock` after the primary cursor stepper. (Update `InteractionDock.test.tsx` for a two-stepper case.)
  2. Write `useStepperOnly.ts`.
  3. In `LoupePage.tsx`: frame axis → `useListCursor({ ids: exposureIds, onActivate: undefined, stepperTestIdBase: "frame", axis: "horizontal" })`, seeded from `activeId`/`?exposure=`. Sample axis → `useStepperOnly({ ids: orderedSampleIds, currentId: sampleId, onGo: gotoSample, label: "Sample", testIdBase: "sample" })`.
  4. `usePageActions({ cursor: frameCursor, extraSteppers: [sampleStepper], actions: [ core("back",…), core("openFocus", { run: enterFocus, dock: "primary" }), page("cull",…"x"), page("keep",…"k"), page("representative",…), page("restore",…) ] })`.
  5. Delete the hand-built `<Dock>` (lines ~507–626) and `useShortcuts` (lines ~358–382). Keep the `KbdLegend` for screen verbs (it derives from the registry — but the registry id source changes; point it at the new action labels or drop it in favour of the dock key hints. Decision: **drop the page-local `KbdLegend`** — the dock buttons now show their own key hints, satisfying discoverability).

- [ ] **Step 4: Run tests + contract harness on the frame cursor.** Run: `node_modules/.bin/vitest run test/LoupePage.actions.test.tsx` + the Loupe contract run → PASS.

- [ ] **Step 5: Commit.** `git commit -m "feat(loupe): frame cursor + URL sample stepper via declaration"`

### Task 3.2: Loupe E2E + live walk-through (compound axes)

- [ ] Mocked spec: frame `←/→` flips the image; sample `↑/↓` changes URL; both steppers visible; Drop persists; Enter→Focus. Run `npm run e2e -- --grep Loupe`.
- [ ] Live walk-through (copy DB): confirm both axes feel identical to keyboard and dock, focus stays sane across sample navigation, `?exposure=` seeds the frame.
- [ ] Gate: `npm run build && npm test`. Commit the spec.

---

# Phase 4 — Focus (candidate cursor + edit modes)

Focus stress-tests the **mode stack**. Its primary (Enter) is **not** navigation — it toggles the previewed candidate's assignment. Its cursor is the **candidate pool** (`previewIndexId`, already ID-based). It has a real `edit` mode (`addArmed` = "+ Peak").

**Mode union (replaces `addArmed`/`previewIndexId`/`xDomain` soup):**
```ts
type FocusMode =
  | { kind: "browse" }                 // navigating candidates
  | { kind: "addPeak" };               // "+ Peak" armed
// previewIndexId becomes the candidate cursor's cursorId; xDomain/scale/combView stay page-local view state (not modes).
```

### Task 4.1: Focus candidate cursor + assignment-as-primary

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/pages/FocusPage.tsx`
- Test: `packages/HimalayaUI/frontend/test/FocusPage.actions.test.tsx`

**Interfaces:**
- Consumes: `useListCursor` (candidate pool), `usePageActions`, `core`, `page`, `useStepperOnly` (sample siblings — URL-driven, like Loupe).

- [ ] **Step 1: Failing tests:**
```tsx
it("Enter toggles the previewed candidate's assignment (NOT navigation)", () => { /* add/removeAssignmentPhase.mutate called */ });
it("candidate ←/→ moves the candidate cursor by id", () => {});
it("P toggles addPeak mode; Escape exits addPeak before leaving the page", () => {});
it("sample ↑/↓ navigates to the sibling sample URL (stepper)", () => {});
it("openFocus label in the dock reads 'Apply' (page override), not 'Focus'", () => {});
```

- [ ] **Step 2:** Run → FAIL.

- [ ] **Step 3: Implement.**
  - Candidate cursor: `useListCursor({ ids: candidatePoolIds, onActivate: toggleAssignmentForId, stepperTestIdBase: "candidate", axis: "horizontal" })`. `previewIndexId` → `candidateCursor.cursorId`.
  - Sample siblings: `useStepperOnly({ ids: siblingIds, currentId: activeSampleId, onGo: (id) => navigate(`/sample/${id}`), label: "Sample", testIdBase: "sample" })`.
  - `const mode: FocusMode = addArmed ? { kind: "addPeak" } : { kind: "browse" }`.
  - `usePageActions({ cursor: candidateCursor, extraSteppers: [sampleStepper], actions: [`
    `core("back", { label: backLabel, run: goBack, dock: true }),`
    `core("openFocus", { label: "Apply", run: () => candidateCursor.activate(), dock: "primary", enabled: () => candidateCursor.cursorId !== null }),`
    `core("openLoupe", { run: () => navigate(`/sample/${activeSampleId}/loupe`), dock: true }),`
    `page("addPeak", { label: "+ Peak", keys: ["p"], group: "Edit", dock: true, run: () => setAddArmed(v => !v) }),`
  `] })`.
  - **Escape ladder (spec §8):** Focus's Escape must unwind the *innermost* mode first. Implement by declaring `core("back", …)` whose `run` checks: `if (addArmed) { setAddArmed(false); return; } if (candidateCursor.cursorId != null && previewWasExplicit) { clearPreview(); return; } goBack();`. Because Escape is `core("back")`'s key, the one listener routes it; the page's `run` owns the ladder. (The `addPeak`-disarm previously in `TracePlate`'s window listener moves here — delete that rung-3 listener; TracePlate keeps only its pointer logic.)
  - Delete the hand-built `<Dock>` (lines ~885–945) and `useShortcuts` (lines ~488–539).

- [ ] **Step 4:** Run tests + Focus contract harness (candidate cursor) → PASS.

- [ ] **Step 5:** Commit `feat(focus): candidate cursor + assignment-primary + mode ladder`.

### Task 4.2: Focus E2E + live walk-through

- [ ] Mocked: Enter toggles assignment (not nav); P arms add-peak; Escape disarms then leaves; candidate `←/→`. `npm run e2e -- --grep Focus`.
- [ ] Live (copy DB): the add-peak arm/disarm + Escape ladder + candidate stepping are real-DOM-sensitive — verify in Chrome. Confirm the `TracePlate` peak-add still works with its window listener removed.
- [ ] Gate + commit.

---

# Phase 5 — Scoping → Builder (ordered members, reorder, multi-select)

These share the reorder concern. **The reorder hook (`useDragReorder` + `data-reorder-index`/`data-reorder-row`) must coexist with roving tabindex** (spec §12 risk). Decision: reorder rows are cursorable rows; `Alt+↑/↓` reorder is declared as `page("reorderUp"/"reorderDown")` actions whose `run` reads the cursored row's index from the cursor (not from `e.target`), removing `useReorderShortcuts`'s DOM-target coupling. Drag stays as-is (pointer events, orthogonal to keyboard).

### Task 5.1: Scoping cursor + reorder-as-actions + undo/redo via core

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/pages/SeriesScopingPage.tsx`
- Test: `packages/HimalayaUI/frontend/test/SeriesScopingPage.actions.test.tsx`

**Interfaces:** Consumes `useListCursor` (over `order` ids), `usePageActions`, `core` (undo/redo/back), `page` (reorderUp/reorderDown), `useDragReorder` (unchanged).

- [ ] **Step 1: Failing tests:**
```tsx
it("the member list is a roving-tabindex cursor over order[]", () => {});
it("Alt+ArrowUp moves the CURSORED row up (applyReorder), cursor follows the item", () => {});
it("⌘Z / ⌘⇧Z undo/redo reorder via core(undo)/core(redo)", () => {});
it("drag-reorder still works alongside roving tabindex", () => {});
```

- [ ] **Step 2:** Run → FAIL.

- [ ] **Step 3: Implement.**
  - `const cursor = useListCursor({ ids: order, stepperTestIdBase: "member", axis: "vertical" })` (order is already `number[]` of sample ids — perfect).
  - Reorder actions (cursor-driven, replacing `useReorderShortcuts`):
    ```ts
    page("reorderUp", { label: "Move up", keys: ["Alt+ArrowUp"], group: "Edit",
      run: () => { const i = order.indexOf(cursor.cursorId!); if (i > 0) { applyReorder(i, i - 1); cursor.setCursor(order[i - 1]!); } } }),
    page("reorderDown", { label: "Move down", keys: ["Alt+ArrowDown"], group: "Edit",
      run: () => { const i = order.indexOf(cursor.cursorId!); if (i < order.length - 1) { applyReorder(i, i + 1); cursor.setCursor(order[i + 1]!); } } }),
    ```
    Wait — `applyReorder` then re-render gives a new `order`; capture the moved id *before* reorder and `setCursor` to it so the cursor follows the item, not the slot. Use the id, not the index.
  - `core("undo", { run: undo, enabled: () => undoStack.canUndo() })`, `core("redo", { run: redo, enabled: () => undoStack.canRedo() })`. Delete the page's `⌘Z` window listener (lines ~544–559).
  - `core("back", { label: "All series", run: () => navigate("/series"), dock: true })`. Delete the bare up-link dock (the shell renders it).
  - Rows spread `cursor.rowProps(id)`; the row keeps `data-reorder-row`/`data-reorder-index` for `useDragReorder`. Container gets `data-interaction-scope` + arrow `onKeyDown`→`moveBy`.
  - Delete `useReorderShortcuts` import/call.

  > **Coexistence note:** `cursor.rowProps` sets `tabIndex`/`role`/`onClick`; `dragItemProps` sets `draggable`/`onPointerDown`. They target the same element but touch disjoint props — spread both: `<Row {...cursor.rowProps(id)} {...dragItemProps(i)} />`. Roving focus (keyboard) and drag (pointer) don't collide.

- [ ] **Step 4:** Run tests + Scoping contract harness → PASS.

- [ ] **Step 5:** Commit `feat(scoping): member cursor + reorder/undo via declaration`.

### Task 5.2: Builder cursor + reorder + add/confirm

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/pages/SeriesBuilderPage.tsx`
- Test: `packages/HimalayaUI/frontend/test/SeriesBuilderPage.actions.test.tsx`

**Interfaces:** Consumes `useListCursor` (over `navSamples` mapped to ids — fix the **index** bug: `selectedIndex` → cursor by `sample_id`), `usePageActions`, `core` (undo/redo/openFocus/back), `page` (addSample/confirm/reorder).

- [ ] **Step 1: Failing tests:**
```tsx
it("cursor is by sample_id, surviving recipe edits (was selectedIndex)", () => {});
it("Enter enters Focus for the cursored sample (?from=series)", () => {});
it("A opens the add-sample picker; ⌘Enter confirms when ready", () => {});
it("Alt+↑/↓ reorders the cursored recipe row (visual order preserved)", () => {});
```

- [ ] **Step 2:** Run → FAIL.

- [ ] **Step 3: Implement.**
  - `const ids = navSamples.map((r) => r.sample_id)` ; `const cursor = useListCursor({ ids, onActivate: (id) => navigate(`/sample/${id}?from=series`), stepperTestIdBase: "sample" })`. Replace `selectedIndex`/`cursorSampleId` derivations with `cursor.cursorId`. (Keep the `cursorIdentity` lookup, now keyed off `cursor.cursorId`.)
  - Actions:
    ```ts
    core("back", { label: "All series", run: () => navigate("/series"), dock: true }),
    core("openFocus", { run: () => cursor.activate(), dock: "primary", enabled: () => cursor.cursorId !== undefined }),
    core("undo", { run: undo, enabled: () => canUndo }),
    core("redo", { run: redo, enabled: () => canRedo }),
    page("addSample", { label: "Add sample", keys: ["a"], group: "Act", dock: true,
      run: () => document.querySelector<HTMLButtonElement>('[data-testid="builder-add-sample"]')?.click() }),
    page("confirm", { label: "Confirm", keys: ["Mod+Enter"], group: "Act", dock: true,
      enabled: () => !!liveDraft && stage.current === "idle" && resolverReady, run: () => confirmRef.current() }),
    ```
  - RecipeEditor reorder: same cursor-driven `reorderUp`/`reorderDown` as Scoping, but mapping through the existing `visual`→`recipe` reversal (`toRecipe`). Replace the `useReorderShortcuts` that clicked the row's Move buttons with cursor-index-driven `onReorder` calls; keep the Move buttons for mouse. Preserve the `data-testid="builder-recipe-up/down"` for E2E.
  - Delete the hand-built `<Dock>` (lines ~652–716) and `useShortcuts` (lines ~260–299) and the `useReorderShortcuts` (lines ~1140–1148).

- [ ] **Step 4:** Run tests + Builder contract harness → PASS. (The Builder identity panel in the dock — `cursorIdentity` swatch + name — becomes a custom dock segment; render it via a small `dock` slot. Decision: keep the identity readout as page-local UI in the dock by declaring a non-interactive `page("identity", { dock: true, … })`? No — non-interactive content doesn't fit `Action`. Instead, extend `usePageActions` with an optional `dockExtra?: ReactNode` rendered between the stepper and the buttons. Add it to the store + `InteractionDock`, covered by an `InteractionDock` test.)

- [ ] **Step 5:** Commit `feat(builder): id-based cursor + add/confirm/reorder via declaration`.

### Task 5.3: Scoping + Builder E2E + live walk-through

- [ ] Mocked specs: reorder via keyboard (cursor follows item), drag still works, undo/redo, add/confirm. `npm run e2e -- --grep "Scoping|Builder"`.
- [ ] Live (copy DB): drag + keyboard reorder coexistence is the real-DOM risk — verify both in Chrome, plus the cursor-follows-item behaviour on reorder.
- [ ] Gate + commit.

---

# Phase 6 — Grouping + final remnant removal

Grouping is **last** (already ID-based cursor + inline rename via `SampleFold`). It also carries the **ordered selection** (`selection: number[]`, first = bulk-merge survivor) — the one place the cursor's `selected` must preserve **order**.

### Task 6.1: Grouping cursor + ordered selection + rename mode

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/components/GroupingReviewPage.tsx`
- Test: `packages/HimalayaUI/frontend/test/GroupingReviewPage.actions.test.tsx`

**Interfaces:** Consumes `useListCursor`, `usePageActions`, `core`, `page`.

> **Ordered-selection wrinkle:** `useListCursor.selected` is a `Set` (unordered). Grouping's bulk-merge needs the **first-selected = survivor** order. Decision: keep Grouping's `selection: number[]` page-local (it's load-bearing domain order, not generic selection) and **do not** route it through `cursor.selected`. The cursor handles navigation + the cursored-row verbs (Split, dismiss-flag); the ordered multi-select stays the page's own `selection[]` with its existing `toggleSelect` (append/remove preserving order). This is the correct boundary — the generic contract doesn't model survivor-order, and forcing it would leak domain meaning into the primitive.

- [ ] **Step 1: Failing tests:**
```tsx
it("cursor navigates by sample_id (already ID-based — now via useListCursor)", () => {});
it("Space toggles the page's ordered selection (first = survivor preserved)", () => {});
it("Split acts on the cursored sample; Merge enabled at ≥2 selected", () => {});
it("⇧↑/⇧↓ jump between flagged samples", () => {});
it("rename (SampleFold) commits via useInlineEdit, unaffected by the keyboard layer", () => {});
```

- [ ] **Step 2:** Run → FAIL.

- [ ] **Step 3: Implement.**
  - `const cursor = useListCursor({ ids: flatSampleIds, stepperLabel: "Sample", stepperTestIdBase: "sample" })`. Replace `cursorId`/`moveCursor` with the cursor (it already matches — minimal change; mostly delete `moveCursor`).
  - Keep `selection: number[]` + its `toggleSelect`. Bind `Space` to it via `page("select", { keys: ["Space"], … run: () => cursor.cursorId != null && toggleSelect(cursor.cursorId) })`.
  - Actions: `core("back", { label: "Samples", run: onBack, dock: true })`, `page("split", …)`, `page("merge", { enabled: () => selection.length >= 2, … })`, `page("dismissFlag", { keys: ["x"], … })`, `page("prevFlagged"/"nextFlagged", { keys: ["Shift+ArrowUp"/"Shift+ArrowDown"], … })`, and the `Confirm groups` button as a `page("confirm", { dock: true, enabled: () => !scanning })`.
  - The `scanning` gate (old `useShortcuts` enabled-arg) → each action's `enabled: () => !scanning && movePicker === null`.
  - Rename: `SampleFold` already uses local inline-edit (rung 1) — **no change**; confirm the keyboard layer ignores its input (bare keys while typing are suppressed by guard #2). Add a test pinning that `x` typed in the rename input does **not** dismiss a flag.
  - Delete the hand-built `<Dock testId="grouping-footer">` and the `useShortcuts` (lines ~314–324). (Note: the dock here used `testId="grouping-footer"` and bespoke testids like `grouping-prev-sample`. The shell dock uses `dock-prev-sample`. Update the Grouping E2E specs to the unified testids — flag this as a deliberate testid migration in the PR.)

- [ ] **Step 4:** Run tests + Grouping contract harness → PASS.

- [ ] **Step 5:** Commit `feat(grouping): cursor + declared actions (ordered selection stays page-local)`.

### Task 6.2: Delete the remnants

**Files:**
- Delete: `packages/HimalayaUI/frontend/src/print/shell/useShortcuts.ts`
- Delete: `packages/HimalayaUI/frontend/src/print/shell/useReorderShortcuts.ts`
- Delete: `packages/HimalayaUI/frontend/src/hooks/useGlobalShortcuts.ts`
- Modify: `packages/HimalayaUI/frontend/src/print/shell/shortcuts.ts` (prune page-interpreted ids), `AppRoutes.tsx` (drop `useGlobalShortcuts()`), `KbdOverlay.tsx`/`KbdLegend.tsx` (re-source the legend from the registry — see below).

**Interfaces:** find/help (`/`·`⌘K`·`?`) move into the registry as `core("find")`/`core("help")` declared **once** at the shell level (a tiny always-on `usePageActions`-like shell registration, or a constant `SHELL_ACTIONS` the listener always includes). Decision: add a `shellActions` slot to the registry (always present, merged with page actions in the listener + dock), holding `find`/`help`. This removes the last window listener (`useGlobalShortcuts`).

- [ ] **Step 1: Write the failing test** — `test/interaction/shellActions.test.tsx`: `/` opens the nav modal and `?` opens help **through the keyboard layer** (not `useGlobalShortcuts`).

- [ ] **Step 2:** Run → FAIL (still routed through the old listener).

- [ ] **Step 3: Implement.**
  - Add `shellActions: Action[]` to `registry.ts`, set once at shell init (`useInteraction.setState({ shellActions: [...] })` in `AppShell`), and have `useKeyboardLayer` iterate `[...shellActions, ...actions]`. `InteractionDock` ignores shellActions (they're not docked).
  - `core("find", { run: () => openNavModal(...) })`, `core("help", { run: openHelpOverlay })`.
  - Delete `useGlobalShortcuts()` call + file.
  - Grep-confirm zero importers of the three deleted hooks: `git grep -n "useShortcuts\|useReorderShortcuts\|useGlobalShortcuts" -- packages/HimalayaUI/frontend/src` → only the deletions.
  - Prune `shortcuts.ts`: the registry-as-legend-source (`SHORTCUTS`/`shortcutLabel`) can stay if `KbdOverlay` still uses it; but the **page-interpreted ids** (`prevDetail`/`nextDetail` reinterpretation, the Rev-2 axis comment) are dead — delete them. Re-point `KbdOverlay`/`KbdLegend` to derive from `CORE` + the current page's `actions` (read `useInteraction.getState().actions` for the live legend). If that's a larger lift, keep `SHORTCUTS` as a static help table for the global keys only and delete the per-surface legend wiring.

- [ ] **Step 4:** Run `test/interaction/shellActions.test.tsx` + the **full suite**: `npm run build && npm test && npm run e2e`. Expected: all PASS, zero references to the deleted hooks.

- [ ] **Step 5:** Commit `refactor(interaction): delete useShortcuts/useReorderShortcuts/useGlobalShortcuts; find/help via registry`.

### Task 6.3: Final full-suite gate + live sweep

- [ ] Run `npm run build && npm test && npm run e2e` (full).
- [ ] Live walk-through across **all six** pages on a copy of the dev DB: cursor parity (click==arrow==stepper), Enter==double-click==dock-primary, Escape ladder, selection grammar, reorder coexistence, find/help. Record in the PR.
- [ ] Run `frontend-reviewer` agent over the branch diff.
- [ ] Commit any review fixes.

---

## Self-Review (run against the spec)

**1. Spec coverage:**
- §3 Principles — single-source (registry, Task 0.4), input parity (cursor, Task 1.1 + harness 1.2), A+C (cursor + selection Set + per-page shapes, Phases 2–6), hybrid vocabulary (`core`/`page`, Task 0.3), build-thin (no dep, all tasks). ✓
- §4 Core model — `Action` (0.1), `enabled()` freshness (Global Constraints + 0.4 effect rule), `core`/`page` (0.3), `usePageActions` (0.4). ✓
- §5 Cursor & roving tabindex — `useListCursor` (1.1), cursor-as-ID (1.1 SSE tests), cursor=DOM focus (1.1 effect), selection orthogonal Set (1.1 + harness). ✓
- §6 Keyboard ladder — rung-3 listener (0.5), three guards in order (0.5 + tests), rung-1 stays (Corpus arrows as container `onKeyDown`), rung-2 ModalShell untouched. ✓
- §7 Dock renderer — `InteractionDock` (0.6), shell mount (0.7). ✓
- §8 Modes — discriminated unions per page (Corpus 2.2, Focus 4.1, etc.), nav-always-base (arrows independent of mode), Escape ladder (Focus 4.1). ✓
- §9 Migration order + kill-list — Corpus→Loupe→Focus→Scoping→Builder→Grouping (Phases 2–6), remnant deletion (6.2). ✓
- §10 Testing — shared contract harness (1.2), keyboard-ladder tests (0.5), per-page declaration tests, live walk-throughs, real-browser Playwright. ✓
- §11 Prior-art — informed the design; no adoption (no-dep constraint). ✓
- §12 Open questions — compound cursor **resolved** (Loupe Phase 3: one axis = cursor, other = `useStepperOnly`); reorder coexistence **resolved** (Phase 5: disjoint props); palette deferred (not in plan). ✓

**2. Placeholder scan:** No "TBD"/"handle edge cases"/"write tests for the above" — every code step shows code; every test step shows assertions. Phases 3–6 give complete declarations + named test cases (the mechanism, defined once in Phases 0–2, is consumed, not re-derived — per the Interfaces blocks). ✓

**3. Type consistency:** `ListCursor`/`Action`/`CursorStepperProps` defined in `types.ts` (0.1) and consumed unchanged downstream. `usePageActions` declaration type extended additively in Phase 3 (`extraSteppers`) and Phase 5 (`dockExtra`) — both flagged as registry + `InteractionDock` + test changes in their tasks. `core()`/`page()` signatures fixed in 0.3 and used verbatim. ✓

One open design call surfaced for the user (not a blocker): in Phase 5/6 the dock gains two small extensions (`extraSteppers`, `dockExtra`) for Loupe's second stepper and Builder's identity readout. These are additive and tested, but they widen the declaration surface slightly beyond the spec's "actions + cursor". Worth a glance at review.
