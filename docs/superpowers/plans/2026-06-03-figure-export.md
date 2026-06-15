# Figure Export Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the greenfield "The Print" figure-export control — a `useFigureExport` logic hook + a presentational `ExportButton` split button (Copy → clipboard; chevron → Download PNG/SVG).

**Architecture:** Logic/presentation split per the Print contract. `useFigureExport` is the sole consumer of the existing `src/lib/figure-export/*` engine; `ExportButton` composes the `Button`/`IconButton`/`Menu` primitives and owns only its menu-open UI state (the `Field` precedent). The clean scientific render is the engine adapter's job (wired at Layer 4) — these units are style-agnostic.

**Tech Stack:** React 18 + TypeScript (strict, `exactOptionalPropertyTypes: true`), Vitest + Testing Library (JSDOM), Storybook (CSF3), Tailwind v4. Spec: `docs/superpowers/specs/2026-06-03-figure-export-design.md`.

---

## ⚠ STANDING CONSTRAINTS (every task, non-negotiable)

- **Branch `worktree-greenfield-ui-rebuild` stays UNMERGED and UNPUSHED.** Do **NOT** run `superpowers:finishing-a-development-branch`, do **NOT** offer merge/PR/push, even after the final task. The rebuild lands as one piece later.
- **Commit ONLY the explicitly-named files** in each task's commit step. **NEVER** `git add -A` or `git add .`.
- **NEVER stage** `src/bones/registry.ts`, `src/bones/contact-sheet.bones.json`, or anything under `docs/superpowers/plans/` (this plan file included — it stays untracked).
- **Commit message's exact last line:** `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
- **Typecheck command:** `npx tsc --noEmit -p tsconfig.build.json` (NOT plain `tsc`).
- **Work from** `packages/HimalayaUI/frontend/`. All paths below are relative to that dir.
- **Dispatch implementers SEQUENTIALLY** — concurrent `git commit` collides on `.git/index.lock`.

## Design-guard & test rules (load-bearing)

- `src/print/components/**` is **NOT** design-guard-exempt. `ExportButton` may use only **token color classes** (`border-hair-strong`, `bg-hair-strong`, `bg-paper-sunk`, `text-ink`…), **named type roles** (`text-meta`…), **layout/positioning utilities** (`relative`, `inline-flex`, `items-stretch`, `w-px`, `right-0`, `top-full`, `overflow-hidden`), and `rounded`/`rounded-md`. **NO** inline appearance literals (`text-[…]`, `rounded-[…]`, `bg-[#…]`, raw `oklch(`/hex, side-stripe borders). `npm run lint:design` enforces this and must stay green. (`useFigureExport.ts` has no JSX → guard-neutral, like the sibling `useDragReorder.ts`.)
- **Class-agnostic tests:** assert via `data-testid` / ARIA role+name / text only. Never assert class strings.
- House `cx` helper is copied locally into each `.tsx` file (no shared export):
  ```ts
  function cx(...parts: Array<string | false | null | undefined>): string {
    return parts.filter(Boolean).join(" ");
  }
  ```

## File Structure

| File | Responsibility |
|---|---|
| `src/print/components/useFigureExport.ts` (create) | The logic hook. Sole consumer of `lib/figure-export/*` + `lib/toast`. Returns handlers + disabled/pending flags. |
| `src/print/components/ExportButton.tsx` (create) | Presentational split button. Composes `Button`+`IconButton`+`Menu`; owns only menu-open state. |
| `src/print/components/ExportButton.stories.tsx` (create) | Storybook stories with spy/no-op handlers. |
| `test/print-components/useFigureExport.test.tsx` (create) | Hook test — `renderHook` + `vi.spyOn` the engine modules + `setToastImpl`. |
| `test/print-components/ExportButton.test.tsx` (create) | Component test — RTL, class-agnostic. |
| `docs/greenfield-component-ledger.md` (modify) | Flip `ExportButton` ⬜ → ✅ (final task). |

### Verified API reference (do not re-derive — these are confirmed against live source)

- **Engine** (`src/lib/figure-export/`): `buildExportPng(spec: ExportSpec, scale=2): Promise<Blob>` and `canExportPng(): boolean` (both in `renderer`); `buildExportSvg(spec): SVGSVGElement` (`renderer`); `canCopyPngToClipboard(): boolean` + `copyPngToClipboard(blob): Promise<void>` (`clipboard`); `downloadBlob(blob, filename): void` (`download`); `buildFilename(stem, "png"|"svg"): string` (`filename`); `ExportSpec` (`types`).
- **Toast** (`src/lib/toast`): `showToast(msg: string, kind?: "info"|"success"|"warning"|"error"): void`; `setToastImpl(fn|null)` for tests.
- **`Button`** (`../ui`): `variant?: "solid"|"accent"|"ghost"|"danger"|"outline"|"ghostInverse"`; extends `ButtonHTMLAttributes` (so `data-testid`, `aria-label`, `disabled`, `onClick` pass through). Base already includes `disabled:opacity-45 disabled:cursor-not-allowed`.
- **`IconButton`** (`../ui`): `{ label: string /* = aria-label, REQUIRED */, tone?, dismiss?, children? }` + `ButtonHTMLAttributes` minus `aria-label` (so `data-testid`, `aria-haspopup`, `aria-expanded`, `disabled`, `onClick` pass through).
- **`Menu<T extends string>`** (`../ui`): `{ open, options: ReadonlyArray<MenuOption<T>>, onSelect: (v:T)=>void, onClose: ()=>void, "aria-label": string, activeValue?, className? }`. `MenuOption<T> = { value: T; label: ReactNode; disabled?: boolean }`. Root is `<div role="menu" data-testid="menu" class="card absolute z-20 mt-1 min-w-[180px] py-1 {className}">`; items are `role="menuitem"`. Menu owns Escape (→ `onClose`) + arrow-nav + click-select-then-close. It does **NOT** own outside-click — the consumer does (the `Field` pattern).

---

## Task 1: `useFigureExport` hook

**Files:**
- Create: `src/print/components/useFigureExport.ts`
- Test: `test/print-components/useFigureExport.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/print-components/useFigureExport.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { useFigureExport } from "../../src/print/components/useFigureExport";
import * as renderer from "../../src/lib/figure-export/renderer";
import * as clipboard from "../../src/lib/figure-export/clipboard";
import * as download from "../../src/lib/figure-export/download";
import * as filename from "../../src/lib/figure-export/filename";
import * as toastModule from "../../src/lib/toast";
import type { ExportSpec } from "../../src/lib/figure-export/types";

const fakeSpec: ExportSpec = {
  title: { primary: "P" }, width: 100, height: 80,
  plot: { marks: [], width: 80, height: 60 },
};
const specThunk = (): ExportSpec => fakeSpec;

let toastSpy: ReturnType<typeof vi.fn>;

beforeEach(() => {
  toastSpy = vi.fn();
  toastModule.setToastImpl(toastSpy);
  // JSDOM has no OffscreenCanvas → canExportPng() is false by default. Default
  // it true for the bulk of tests; capability-off cases set it explicitly.
  vi.spyOn(renderer, "canExportPng").mockReturnValue(true);
  vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(true);
});
afterEach(() => {
  toastModule.setToastImpl(null);
  vi.restoreAllMocks();
});

describe("useFigureExport", () => {
  it("onCopy renders PNG, writes clipboard, emits success toast", async () => {
    const blob = new Blob([new Uint8Array([1])], { type: "image/png" });
    vi.spyOn(renderer, "buildExportPng").mockResolvedValue(blob);
    const write = vi.spyOn(clipboard, "copyPngToClipboard").mockResolvedValue();

    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "the trace"));
    act(() => { result.current.onCopy(); });

    await waitFor(() => expect(write).toHaveBeenCalledWith(blob));
    expect(toastSpy).toHaveBeenCalledWith(
      expect.stringContaining("Copied the trace to clipboard"), "success",
    );
  });

  it("onCopy failure emits an error toast", async () => {
    vi.spyOn(renderer, "buildExportPng").mockRejectedValue(new Error("boom"));
    vi.spyOn(clipboard, "copyPngToClipboard").mockResolvedValue();

    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "the trace"));
    act(() => { result.current.onCopy(); });

    await waitFor(() =>
      expect(toastSpy).toHaveBeenCalledWith(
        expect.stringContaining("Couldn't copy figure"), "error",
      ),
    );
  });

  it("onDownloadPng renders PNG and downloads with the built filename", async () => {
    const blob = new Blob([new Uint8Array([2])], { type: "image/png" });
    vi.spyOn(renderer, "buildExportPng").mockResolvedValue(blob);
    vi.spyOn(filename, "buildFilename").mockReturnValue("stem-2099-01-01.png");
    const dl = vi.spyOn(download, "downloadBlob").mockImplementation(() => {});

    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "the trace"));
    act(() => { result.current.onDownloadPng(); });

    await waitFor(() => expect(dl).toHaveBeenCalledWith(blob, "stem-2099-01-01.png"));
  });

  it("onDownloadSvg serializes the SVG and downloads it", async () => {
    const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    vi.spyOn(renderer, "buildExportSvg").mockReturnValue(svg);
    vi.spyOn(filename, "buildFilename").mockReturnValue("stem-2099-01-01.svg");
    const dl = vi.spyOn(download, "downloadBlob").mockImplementation(() => {});

    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "the trace"));
    act(() => { result.current.onDownloadSvg(); });

    await waitFor(() => {
      expect(dl).toHaveBeenCalledTimes(1);
      const [blobArg, nameArg] = dl.mock.calls[0];
      expect(blobArg).toBeInstanceOf(Blob);
      expect((blobArg as Blob).type).toBe("image/svg+xml");
      expect(nameArg).toBe("stem-2099-01-01.svg");
    });
  });

  it("copyDisabled is true when clipboard is unsupported", () => {
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(false);
    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "x"));
    expect(result.current.copyDisabled).toBe(true);
  });

  it("copyDisabled and pngDisabled are true when PNG export is unsupported", () => {
    vi.spyOn(renderer, "canExportPng").mockReturnValue(false);
    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "x"));
    expect(result.current.copyDisabled).toBe(true);
    expect(result.current.pngDisabled).toBe(true);
  });

  it("guards against a double-invocation while a render is in flight", async () => {
    // buildExportPng never resolves → the first call holds the in-flight guard.
    const build = vi.spyOn(renderer, "buildExportPng").mockReturnValue(new Promise<Blob>(() => {}));
    vi.spyOn(clipboard, "copyPngToClipboard").mockResolvedValue();

    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "x"));
    act(() => { result.current.onCopy(); result.current.onCopy(); });

    expect(build).toHaveBeenCalledTimes(1);
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npm test -- useFigureExport`
Expected: FAIL — `useFigureExport` is not exported / file does not exist.

- [ ] **Step 3: Write the implementation**

Create `src/print/components/useFigureExport.ts`:

```ts
import { useCallback, useRef, useState } from "react";
import type { ExportSpec } from "../../lib/figure-export/types";
import { buildExportPng, buildExportSvg, canExportPng } from "../../lib/figure-export/renderer";
import { canCopyPngToClipboard, copyPngToClipboard } from "../../lib/figure-export/clipboard";
import { downloadBlob } from "../../lib/figure-export/download";
import { buildFilename } from "../../lib/figure-export/filename";
import { showToast } from "../../lib/toast";

export interface UseFigureExport {
  onCopy: () => void;
  onDownloadPng: () => void;
  onDownloadSvg: () => void;
  copyDisabled: boolean;
  pngDisabled: boolean;
  pending: boolean;
}

/**
 * The figure-export logic hook — the ONLY consumer of `lib/figure-export/*`.
 * `spec` is a thunk evaluated at click time so it captures fresh state. The
 * clean scientific styling lives in the engine adapter that builds the
 * `ExportSpec`; this hook is style-agnostic.
 */
export function useFigureExport(
  spec: () => ExportSpec,
  filenameStem: string,
  ariaContext: string,
): UseFigureExport {
  const [pending, setPending] = useState(false);
  // A synchronous in-flight ref: `pending` state flips async, so two clicks in
  // one tick would both pass a state-only guard (the legacy FigureExportControls
  // double-fire window). The ref closes it; the state still drives the UI.
  const inFlight = useRef(false);

  const run = useCallback(
    async (work: () => Promise<void> | void, failMsg: string): Promise<void> => {
      if (inFlight.current) return;
      inFlight.current = true;
      setPending(true);
      try {
        await work();
      } catch (err) {
        showToast(failMsg, "error");
        // eslint-disable-next-line no-console
        console.warn(err);
      } finally {
        inFlight.current = false;
        setPending(false);
      }
    },
    [],
  );

  const onCopy = useCallback(() => {
    void run(async () => {
      const blob = await buildExportPng(spec());
      await copyPngToClipboard(blob);
      showToast(`Copied ${ariaContext} to clipboard`, "success");
    }, "Couldn't copy figure. Try Download instead.");
  }, [run, spec, ariaContext]);

  const onDownloadPng = useCallback(() => {
    void run(async () => {
      const blob = await buildExportPng(spec());
      downloadBlob(blob, buildFilename(filenameStem, "png"));
    }, "Couldn't render figure for download.");
  }, [run, spec, filenameStem]);

  const onDownloadSvg = useCallback(() => {
    void run(() => {
      const svg = buildExportSvg(spec());
      const xml = new XMLSerializer().serializeToString(svg);
      const blob = new Blob([xml], { type: "image/svg+xml" });
      downloadBlob(blob, buildFilename(filenameStem, "svg"));
    }, "Couldn't render figure for download.");
  }, [run, spec, filenameStem]);

  const canCopy = canCopyPngToClipboard();
  const canPng = canExportPng();
  return {
    onCopy,
    onDownloadPng,
    onDownloadSvg,
    copyDisabled: !canCopy || !canPng || pending,
    pngDisabled: !canPng || pending,
    pending,
  };
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npm test -- useFigureExport`
Expected: PASS (7/7).

- [ ] **Step 5: Lint + typecheck**

Run: `npm run lint:design && npx tsc --noEmit -p tsconfig.build.json`
Expected: both exit 0.

- [ ] **Step 6: Commit**

```bash
git add src/print/components/useFigureExport.ts test/print-components/useFigureExport.test.tsx
git commit -m "$(cat <<'EOF'
feat(print): useFigureExport hook — figure-export logic over lib/figure-export

Sole consumer of the export engine. Copy → PNG → clipboard; download PNG/SVG;
capability gating (canExportPng/canCopyPngToClipboard); in-flight ref guard;
error toasts. Style-agnostic — the clean render is the adapter's job.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: `ExportButton` presentational split button

**Files:**
- Create: `src/print/components/ExportButton.tsx`
- Test: `test/print-components/ExportButton.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/print-components/ExportButton.test.tsx`:

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ExportButton } from "../../src/print/components/ExportButton";

function setup(overrides: Partial<React.ComponentProps<typeof ExportButton>> = {}) {
  const props = {
    onCopy: vi.fn(),
    onDownloadPng: vi.fn(),
    onDownloadSvg: vi.fn(),
    ariaContext: "the trace",
    ...overrides,
  };
  render(<ExportButton {...props} />);
  return props;
}

describe("ExportButton", () => {
  it("renders the Copy action with an aria-label built from ariaContext", () => {
    setup();
    expect(screen.getByRole("button", { name: /copy the trace to clipboard/i }))
      .toBeInTheDocument();
  });

  it("does not show the download menu until the chevron is clicked", () => {
    setup();
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    const menu = screen.getByRole("menu", { name: /download formats/i });
    expect(menu).toBeInTheDocument();
    expect(screen.getByRole("menuitem", { name: /download as png/i })).toBeInTheDocument();
    expect(screen.getByRole("menuitem", { name: /download as svg/i })).toBeInTheDocument();
  });

  it("primary click calls onCopy", () => {
    const props = setup();
    fireEvent.click(screen.getByTestId("export-copy"));
    expect(props.onCopy).toHaveBeenCalledTimes(1);
  });

  it("selecting PNG calls onDownloadPng and closes the menu", () => {
    const props = setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    fireEvent.click(screen.getByRole("menuitem", { name: /download as png/i }));
    expect(props.onDownloadPng).toHaveBeenCalledTimes(1);
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
  });

  it("selecting SVG calls onDownloadSvg", () => {
    const props = setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    fireEvent.click(screen.getByRole("menuitem", { name: /download as svg/i }));
    expect(props.onDownloadSvg).toHaveBeenCalledTimes(1);
  });

  it("disables Copy when copyDisabled", () => {
    setup({ copyDisabled: true });
    expect(screen.getByTestId("export-copy")).toBeDisabled();
  });

  it("disables the PNG menuitem when pngDisabled (SVG stays enabled)", () => {
    setup({ pngDisabled: true });
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    expect(screen.getByRole("menuitem", { name: /download as png/i })).toBeDisabled();
    expect(screen.getByRole("menuitem", { name: /download as svg/i })).not.toBeDisabled();
  });

  it("the page-level disabled prop ORs into Copy and the trigger", () => {
    setup({ disabled: true });
    expect(screen.getByTestId("export-copy")).toBeDisabled();
    expect(screen.getByRole("button", { name: /download formats/i })).toBeDisabled();
  });

  it("Escape closes the open menu", () => {
    setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    const menu = screen.getByRole("menu", { name: /download formats/i });
    fireEvent.keyDown(menu, { key: "Escape" });
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
  });

  it("an outside pointerdown closes the open menu", () => {
    setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    expect(screen.getByRole("menu", { name: /download formats/i })).toBeInTheDocument();
    fireEvent.pointerDown(document.body);
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npm test -- ExportButton`
Expected: FAIL — `ExportButton` does not exist.

- [ ] **Step 3: Write the implementation**

Create `src/print/components/ExportButton.tsx`:

```tsx
import { useEffect, useRef, useState } from "react";
import { Button, IconButton, Menu } from "../ui";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ExportButtonProps {
  onCopy: () => void;
  onDownloadPng: () => void;
  onDownloadSvg: () => void;
  /** Copy unavailable (no clipboard / no PNG renderer / a render in flight). */
  copyDisabled?: boolean;
  /** PNG download unavailable (no PNG renderer / a render in flight). */
  pngDisabled?: boolean;
  /** A render is in flight (blocks SVG too). */
  pending?: boolean;
  /** Page-level gate (data not ready). ORs into every action. */
  disabled?: boolean;
  /** Fills aria-labels: "Copy {ariaContext} to clipboard". */
  ariaContext: string;
  /** PLACEMENT ONLY. */
  className?: string;
}

/**
 * The figure-export split button (mockup: series-builder rail-foot "Copy as
 * PNG"): a bordered group with a primary **Copy** action and a `▾` chevron that
 * opens a two-item download menu (PNG / SVG). Presentational — every side
 * effect is a prop (wire it with `useFigureExport`). Owns ONLY its menu-open
 * state + outside-pointerdown dismissal (the `Field` precedent; `Menu` owns
 * Escape + arrow-nav).
 */
export function ExportButton({
  onCopy,
  onDownloadPng,
  onDownloadSvg,
  copyDisabled = false,
  pngDisabled = false,
  pending = false,
  disabled = false,
  ariaContext,
  className,
}: ExportButtonProps): JSX.Element {
  const [open, setOpen] = useState(false);
  const wrapRef = useRef<HTMLSpanElement | null>(null);

  // Outside-pointerdown closes the menu (mirrors Field/Popover). Bound only
  // while open.
  useEffect(() => {
    if (!open) return;
    const onPointerDown = (e: PointerEvent): void => {
      if (wrapRef.current && !wrapRef.current.contains(e.target as Node)) {
        setOpen(false);
      }
    };
    document.addEventListener("pointerdown", onPointerDown);
    return () => document.removeEventListener("pointerdown", onPointerDown);
  }, [open]);

  const copyOff = disabled || copyDisabled;
  const pngOff = disabled || pngDisabled;
  const svgOff = disabled || pending;

  return (
    <span
      ref={wrapRef}
      data-testid="export-button"
      className={cx(
        "relative inline-flex items-stretch border border-hair-strong rounded overflow-hidden",
        className,
      )}
    >
      <Button
        variant="ghost"
        data-testid="export-copy"
        aria-label={`Copy ${ariaContext} to clipboard`}
        disabled={copyOff}
        onClick={onCopy}
      >
        Copy
      </Button>
      <span className="w-px bg-hair-strong" aria-hidden="true" />
      <IconButton
        label="Download formats"
        data-testid="export-menu-trigger"
        aria-haspopup="menu"
        aria-expanded={open}
        disabled={disabled || pending}
        onClick={() => setOpen((o) => !o)}
      >
        ▾
      </IconButton>
      <Menu<"png" | "svg">
        open={open}
        options={[
          { value: "png", label: "Download as PNG", disabled: pngOff },
          { value: "svg", label: "Download as SVG", disabled: svgOff },
        ]}
        onSelect={(v) => {
          if (v === "png") onDownloadPng();
          else onDownloadSvg();
        }}
        onClose={() => setOpen(false)}
        aria-label="Download formats"
        className="right-0 top-full"
      />
    </span>
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npm test -- ExportButton`
Expected: PASS (10/10).

- [ ] **Step 5: Lint + typecheck**

Run: `npm run lint:design && npx tsc --noEmit -p tsconfig.build.json`
Expected: both exit 0. (If `lint:design` flags `border-hair-strong`/`bg-hair-strong`, that would contradict their use in `RailBack`/`Dock`/`ScopePlate` — it will not; do not swap to a literal.)

- [ ] **Step 6: Commit**

```bash
git add src/print/components/ExportButton.tsx test/print-components/ExportButton.test.tsx
git commit -m "$(cat <<'EOF'
feat(print): ExportButton — presentational figure-export split button

Copy primary + chevron menu (Download PNG/SVG) over Button/IconButton/Menu.
Owns only menu-open state + outside-pointerdown dismissal (Field precedent);
all side effects are props. Placement-only; lint:design clean.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: `ExportButton` Storybook stories

**Files:**
- Create: `src/print/components/ExportButton.stories.tsx`

- [ ] **Step 1: Write the stories**

Create `src/print/components/ExportButton.stories.tsx`. (No test precedes this — stories are verified by `build-storybook` + visual check. Match the `Meta`/`StoryObj` import style of a sibling story such as `MemberRow.stories.tsx` if it differs from below.)

```tsx
import type { Meta, StoryObj } from "@storybook/react";
import { ExportButton } from "./ExportButton";

const noop = (): void => {};

const meta: Meta<typeof ExportButton> = {
  title: "Components/ExportButton",
  component: ExportButton,
  args: {
    onCopy: noop,
    onDownloadPng: noop,
    onDownloadSvg: noop,
    ariaContext: "the LL37 series",
  },
};
export default meta;

type Story = StoryObj<typeof ExportButton>;

/** Default — Copy enabled; click the chevron to reveal Download PNG / SVG. */
export const Default: Story = {};

/** Clipboard unavailable (e.g. insecure origin) — Copy is disabled, downloads remain. */
export const CopyDisabled: Story = { args: { copyDisabled: true } };

/** A render is in flight — every action is blocked. */
export const Pending: Story = { args: { pending: true, copyDisabled: true, pngDisabled: true } };

/** Page-level gate (data not ready) — fully disabled. */
export const FullyDisabled: Story = { args: { disabled: true, copyDisabled: true, pngDisabled: true } };
```

- [ ] **Step 2: Build Storybook to verify the stories compile**

Run: `npm run build-storybook`
Expected: exit 0 (Storybook globs `src/print/**/*.stories.tsx` automatically).

- [ ] **Step 3: Typecheck**

Run: `npx tsc --noEmit -p tsconfig.build.json`
Expected: exit 0.

- [ ] **Step 4: Commit**

```bash
git add src/print/components/ExportButton.stories.tsx
git commit -m "$(cat <<'EOF'
test(print): ExportButton Storybook stories (default / disabled / pending)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Ledger flip + full gate

**Files:**
- Modify: `docs/greenfield-component-ledger.md`

- [ ] **Step 1: Flip the ledger entry ⬜ → ✅**

In `docs/greenfield-component-ledger.md`, the "Other overlays" table row for `ExportButton` currently ends with:
`| series-plot | ⬜ (design approved 2026-06-03) |`
Change the status cell to:
`| series-plot | ✅ | \`components/ExportButton.tsx\` + \`components/useFigureExport.ts\` |`

And in the "Modals/overlays" summary line, change
`\`ExportButton\` (+\`useFigureExport\`) ⬜ **design approved 2026-06-03**`
to
`\`ExportButton\` (+\`useFigureExport\`) ✅ **built 2026-06-03**`.

- [ ] **Step 2: Run the full gate**

Run each; all must pass:
```bash
npm test -- print-components        # the two new suites + the rest of print-components green
npm run lint:design                  # placement-only proven
npx tsc --noEmit -p tsconfig.build.json
npm run build-storybook              # exit 0
```
Expected: all green / exit 0.

- [ ] **Step 3: Visual verification in Storybook**

Build is already done. Serve and inspect (Playwright MCP), then clean up artifacts:
```bash
python3 -m http.server 6017 --directory storybook-static
```
Navigate to `http://localhost:6017/iframe.html?id=components-exportbutton--default&viewMode=story`. Confirm: the split button reads as one bordered control with a hairline divider; **Copy** on the left, `▾` on the right; clicking `▾` opens a right-aligned menu with "Download as PNG" / "Download as SVG". Check `--copy-disabled` and `--fully-disabled` stories show the disabled treatment. (Benign favicon.ico 404 in console is expected. The production Storybook build needs a full pointer-event sequence to open the menu programmatically — open it by clicking in the live page instead.)

Then kill the server and remove build artifacts:
```bash
# stop the http.server job; then:
rm -rf storybook-static .playwright-mcp
```

- [ ] **Step 4: Commit**

```bash
git add docs/greenfield-component-ledger.md
git commit -m "$(cat <<'EOF'
docs(print): ledger — ExportButton + useFigureExport built (Layer 3)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 5: STOP.** Do **not** run `finishing-a-development-branch`, do **not** merge/push/PR. Report completion; the branch stays unmerged per the standing constraints.

---

## Self-Review (completed by plan author)

- **Spec coverage:** Export model (content-WYSIWYG / style-transformed, no preview) → realized by the style-agnostic hook + `CLEAN_SCIENTIFIC` note (Task 1 docstring). Split button (Copy + PNG/SVG menu) → Task 2. `useFigureExport` contract (handlers + `copyDisabled`/`pngDisabled`/`pending`) → Task 1. `ExportButton` props/ownership → Task 2. Error & capability handling → Task 1 tests (copy/png gating, error toasts) + Task 2 (disabled propagation). Testing (hook + component + stories, class-agnostic) → Tasks 1–3. Ledger bookkeeping → Task 4. Conflict modal is a non-goal — correctly **not** a task here.
- **Placeholder scan:** none — every code/test/command step is concrete.
- **Type consistency:** `UseFigureExport` fields (`onCopy`/`onDownloadPng`/`onDownloadSvg`/`copyDisabled`/`pngDisabled`/`pending`) match `ExportButtonProps` and the spread wiring; `Menu<"png"|"svg">` generic matches the `onSelect` switch; engine signatures match the verified API reference.
- **Spec deviations (deliberate, within the approved spec):** (1) hook colocated in `components/` (not a new `hooks/` dir) to match the sole `useDragReorder` precedent — the spec permitted "or colocated"; (2) the dropdown is asserted by `role="menu"` + `aria-label="Download formats"` because the `Menu` primitive's root testid is fixed at `"menu"` (no `testId` prop) — the spec's `export-menu` testid is superseded; `export-button`/`export-copy`/`export-menu-trigger` testids are preserved.
