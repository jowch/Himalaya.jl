# Greenfield UI Primitive Sweep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring the seeded `src/print/ui/` design-system primitives to The Print's mockup fidelity — harden the four primitives with mockup gaps, build the ~20 atomic/layout/form primitives the mockups need, and verify every one against the **impeccable** quality bar in Storybook — so Phase 2 component rebuilds compose from a complete, impeccable primitive set.

**Architecture:** Each primitive is a closed-look component in `src/print/ui/` (appearance internal via tokenized Tailwind utilities + a local `cx()` join; consumer `className` is placement-only). Each is TDD'd against `data-*` attributes (never class strings), gets a Storybook story rendering all of its states, and passes an explicit **impeccable pass** (a per-primitive checklist distilled from the `impeccable` skill + `DESIGN.md` named rules) before commit. All work happens in the greenfield tree only — the retiring `src/components/ui/` is untouched and deleted at Phase-4 cutover.

**Tech Stack:** React 18, TypeScript strict (`exactOptionalPropertyTypes`), TailwindCSS v4 (`@theme` tokens in `src/styles.css`), Vitest + @testing-library/react, Storybook v10 (`@storybook/react-vite`). Design enforced by `scripts/check-design.mjs` (`lint:design`, fails build on inline appearance utilities / legacy imports outside `print/ui/`).

**Spec:** `docs/superpowers/specs/2026-05-31-greenfield-ui-rebuild-design.md` (this is the back half of Phase 1 — primitive hone + expand).
**Inputs:** `docs/greenfield-rebuild/design-system-catalog.md` (the gap list), `docs/greenfield-rebuild/surface-map.md`, `DESIGN.md`, `PRODUCT.md`, and the `impeccable` skill (`~/.claude/skills/impeccable`).

**Working directory for all commands:** `packages/HimalayaUI/frontend`. Paths below are relative to it unless noted (docs live at repo root `docs/`).

---

## Conventions every task follows

These hold for **every** primitive task (2–25). Stated once here; each task gives only its concrete code (interface, implementation, story, test cases) — not re-pasted boilerplate.

**File set per new primitive `X`:**
- Create `src/print/ui/X.tsx` — the primitive.
- Create `src/print/ui/X.stories.tsx` — Storybook story (mirrors `Button.stories.tsx`: `satisfies Meta<typeof X>`, one named export per state/variant).
- Create `test/print-ui/X.test.tsx` — Vitest unit test importing from `../../src/print/ui/X`.
- Modify `src/print/ui/index.ts` — add `export { X } from "./X";` (+ `export type { ... }` for any public unions), keeping the existing alphabter-agnostic grouping.

**The closed-look idiom (house style, verified against the seeded primitives):**
- Appearance lives in the primitive via Tailwind utilities that resolve to `@theme` tokens (`bg-paper`, `text-ink-soft`, `border-hair-strong`, `rounded-sm`, `rounded-full`, `text-sm`, …). Never inline `text-[Npx]` / `rounded-[Npx]` / raw `oklch(...)` / `#hex` in a `print/ui` consumer-facing string — the guard *allows* it in `print/ui/**` but the impeccable pass forbids off-scale values; use named roles.
- Phase/dynamic colors come from `phaseColor(phase)` (import from `../../phases`) applied via inline `style`, exactly as the seeded `PhaseChip`/`ScoreBar` do. Never duplicate phase hexes.
- Variant → class maps use a local `const fooClass: Record<Variant,string> = {...}` plus a local `cx(...parts)` join helper:
  ```ts
  function cx(...parts: Array<string | false | null | undefined>): string {
    return parts.filter(Boolean).join(" ");
  }
  ```
- The consumer `className` prop is appended LAST so placement utilities win position; it never carries appearance.
- Set `data-*` attributes for every variant/state/tone so tests assert on them (e.g. `data-active`, `data-tone`, `data-armed`). Tests must NOT assert on Tailwind class strings.

**The per-primitive step template (each task's steps follow this shape):**
1. **Write the failing test** — `test/print-ui/X.test.tsx` asserting render + each state's `data-*` + behavior (onClick/onChange) + a11y (`role`/`aria-*`). Imports: `import { render, screen, fireEvent } from "@testing-library/react"; import { describe, it, expect, vi } from "vitest"; import { X } from "../../src/print/ui/X";`.
2. **Run it; verify it fails** — `npm test -- X` → FAIL (cannot resolve `../../src/print/ui/X`).
3. **Implement** `src/print/ui/X.tsx` to the interface + spec given in the task; add the barrel export.
4. **Run it; verify it passes** — `npm test -- X` → PASS.
5. **Write the Storybook story** `src/print/ui/X.stories.tsx` covering every state the impeccable pass requires (default + hover/focus/active where interactive + disabled + each variant).
6. **Impeccable pass** — (a) run the deterministic scan `npx impeccable detect src/print/ui/X.tsx` → expect exit 0 / no anti-pattern hits (fix any); (b) verify against `docs/greenfield-rebuild/impeccable-primitive-checklist.md` (Task 1) and fix any miss; (c) run the two guards: `npm run lint:design` (exit 0) and `npx tsc --noEmit -p tsconfig.build.json` (exit 0).
7. **Commit** — `git add src/print/ui/X.tsx src/print/ui/X.stories.tsx src/print/ui/index.ts test/print-ui/X.test.tsx && git commit -m "feat(print-ui): <X> primitive at mockup fidelity + impeccable pass"`.

**Tokens:** This plan adds **no** `@theme` tokens by default. Phase colors stay single-sourced in `phases.ts`. `--color-unindexed`, `--color-peak-manual`, and the large serif display step are needed by Phase-2 renderers / Phase-3 pages, not these primitives, and are deferred there. The one exception (a serif size for `EmptyState`) is handled inside Task 19 with an explicit decision.

**DEFERRED to Phase 2 (do NOT build here):** `Sparkline`, `MiniWaterfall` (figure renderers needing real trace data via `lib/plot/sparkline.ts` + the intensity model), and all composite surface components (`DataTable`, `SampleRow`, `ThumbnailStrip`, `Thumbnail`, `FolioCard`, `PhasecallBlock`, `CandidateRow`, `CombPanel`, `CustomIndexSheet`, `ScopingPlate`, `ScopingSampleRow`, `NotesMargin`, `ExposureSwitcher`, `RegionBracket`, `PhaseReadingPanel`, `MemberList`, `VerdictCard`, `RepresentativeBox`, `RailbackTab`, `FloatingDock`, `FloatingCullBar`, `ExportButton`, `Stepper`). These assemble the primitives below.

---

### Task 1: Impeccable primitive checklist (criteria baked in from the start)

Author the per-primitive impeccable bar so every later task verifies against a fixed, written standard rather than ad-hoc judgement.

**Files:**
- Create: `docs/greenfield-rebuild/impeccable-primitive-checklist.md`

- [ ] **Step 1: Confirm the impeccable context loads** (so the checklist reflects the real register/laws)

Run (from `packages/HimalayaUI/frontend`):
```bash
IMPECCABLE_CONTEXT_DIR="$(git rev-parse --show-toplevel)" node ~/.claude/skills/impeccable/scripts/load-context.mjs | head -5
```
Expected: JSON with `"hasProduct": true` and `"hasDesign": true`. (Register is **product**.) If either is false, STOP and report — do not invent criteria.

- [ ] **Step 2: Write the checklist**

Create `docs/greenfield-rebuild/impeccable-primitive-checklist.md` with this exact content:

```markdown
# Impeccable primitive checklist (greenfield src/print/ui)

Every `src/print/ui/` primitive must satisfy this before commit. Distilled from the
`impeccable` skill (register: **product**) + `DESIGN.md` named rules + `PRODUCT.md`.
Source of truth for the per-task "Impeccable pass" step.

## A. Second-Channel rule (HARD — color-blind requirement)
- No meaning rests on hue alone. Any color-coded element pairs hue with a second
  channel: shape, label, position, or pattern. Phase/status discs use the Semantic
  Dot pattern (color + `role="img"` + `aria-label`). Verify the encoding still reads
  in grayscale.

## B. Semantic color + Phase-carries-surface
- Every color is a label (interaction / provenance / phase / status). No decorative color.
- Phase color via `phaseColor()` only (never a literal). Terracotta `accent` is the one
  rationed interaction mark — used for the single primary/grease-pencil action, focus
  rings, the live-edit/reject mark. Not for decoration or inactive states.

## C. Interaction states (product register — ship them all)
- Every interactive primitive defines: default, hover, focus-visible, active, and
  disabled where applicable. Focus is a visible `focus-visible` outline in `accent`,
  never removed without replacement. No heavy/full-saturation accent on inactive states.

## D. Touch target
- Interactive targets are ≥44×44px effective, OR the primitive documents (in a code
  comment) that it is a known-dense control whose hit area is enlarged by its container.
  (The compact toggles were a flagged audit defect — do not reproduce silently.)

## E. Voice discipline (DESIGN.md type rules)
- Serif (Newsreader, via `.text-display`/`.text-headline`/`.text-headline-lg`) = titles
  on a plate only. Sans = chrome/prose. Monospace = measured values only (q, score,
  counts, ids, lattice). No serif in labels/buttons/data; no off-scale `text-[Npx]`.

## F. Tokenized appearance (Fixed-Scale + radius)
- Sizes/colors/radii come from named scale roles / `@theme` tokens. Radius is
  `rounded-sm` (5px) / `rounded-full`. No arbitrary `rounded-[Npx]`, no raw color literal
  in a consumer string. (`lint:design` enforces outside `print/ui/`; we hold `print/ui/`
  to the same bar by choice.)

## G. Absolute bans (impeccable shared + product)
- No side-stripe accent border (`border-l/r` >1px as a colored accent) — full border +
  leading icon/word instead. No gradient text. No decorative glassmorphism/`backdrop-blur`.
  No display font in UI labels. No bounce/elastic motion; transitions are 120ms ease-out
  on color/background/border/opacity only — never on layout properties. Respect
  `prefers-reduced-motion`. No em dashes in any rendered copy.

## H. Flat-except-the-plate
- Only `Card` (the plate) and plate-like surfaces carry the Plate-Lift shadow. Other
  primitives are flat: tonal steps (paper → paper-sunk → plate) + hairlines do separation.

## I. Closed-look API hygiene
- Appearance internal; consumer `className` is placement-only and appended last. Variants
  are a `Record<Variant,string>` map + local `cx()`. Public union types exported from the
  barrel. `data-*` attributes set for every variant/state so tests avoid class strings.

## J. Deterministic scan (impeccable detect)
- `npx impeccable detect src/print/ui/<Name>.tsx` returns no anti-pattern hits. This is
  the deterministic floor (side-stripes, gradient text, glass, off-scale values, provider
  tells); a clean scan is necessary but NOT sufficient — A–I still apply by inspection.

## Per-primitive sign-off
For each primitive, confirm: ☐ A ☐ B ☐ C ☐ D ☐ E ☐ F ☐ G ☐ H ☐ I ☐ J, and that the
Storybook story renders every state listed in C. A primitive that legitimately has no
interactive state (e.g. a static badge) marks C/D N/A with a one-line reason.
```

- [ ] **Step 3: Commit**
```bash
git add docs/greenfield-rebuild/impeccable-primitive-checklist.md
git commit -m "docs(rebuild): impeccable primitive checklist — the per-primitive quality bar"
```

---

## Hardening the four gapped primitives (Tasks 2–5)

These modify the **seeded `src/print/ui/`** copies only (the old `src/components/ui/` is left to rot until cutover). Exact mockup specimens verified 2026-05-31. New tests live in `test/print-ui/` and import from `src/print/ui/`.

### Task 2: Button — `armed` toggle state

The focus-plot "+ Peak" tool button has an armed/active state: terracotta fill, paper text. Add an `armed` boolean that renders the accent-filled active look, distinct from the `accent` variant (which is for primary actions, not toggles) by setting `data-armed` and `aria-pressed`.

**Files:**
- Modify: `src/print/ui/Button.tsx`
- Test: `test/print-ui/Button.test.tsx`

- [ ] **Step 1: Write the failing test** — `test/print-ui/Button.test.tsx`:
```tsx
import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Button } from "../../src/print/ui/Button";

describe("<Button> armed state", () => {
  it("defaults to not-armed and sets no aria-pressed", () => {
    render(<Button>+ Peak</Button>);
    const b = screen.getByRole("button", { name: "+ Peak" });
    expect(b.getAttribute("data-armed")).toBe(null);
    expect(b.getAttribute("aria-pressed")).toBe(null);
  });
  it("reflects armed on data-armed and aria-pressed", () => {
    render(<Button armed>+ Peak</Button>);
    const b = screen.getByRole("button", { name: "+ Peak" });
    expect(b.getAttribute("data-armed")).toBe("true");
    expect(b.getAttribute("aria-pressed")).toBe("true");
  });
  it("keeps the variant data-attr when armed", () => {
    render(<Button variant="ghost" armed>x</Button>);
    expect(screen.getByRole("button").getAttribute("data-variant")).toBe("ghost");
  });
  it("still fires onClick when armed", () => {
    const onClick = vi.fn();
    render(<Button armed onClick={onClick}>x</Button>);
    fireEvent.click(screen.getByRole("button"));
    expect(onClick).toHaveBeenCalledTimes(1);
  });
});
```

- [ ] **Step 2: Run; verify fail** — `npm test -- Button` → FAIL (`armed` not a prop; `data-armed` null when armed).

- [ ] **Step 3: Implement.** Edit `src/print/ui/Button.tsx`. Add `armed` to the interface and an armed-class branch. Final file:
```tsx
import type { ButtonHTMLAttributes } from "react";

export type ButtonVariant = "solid" | "accent" | "ghost" | "danger";

interface ButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: ButtonVariant;
  /** Toggled/active tool-button state (focus-plot "+ Peak" armed): terracotta
   *  fill, paper text, `aria-pressed`. Distinct from `variant="accent"` (a
   *  primary action, not a toggle). */
  armed?: boolean;
}

const variantClass: Record<ButtonVariant, string> = {
  solid:
    "bg-ink border border-ink text-paper hover:brightness-110 " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  accent:
    "bg-accent border border-accent text-paper hover:brightness-110 " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  ghost:
    "text-ink-soft hover:text-ink hover:bg-paper-sunk border border-transparent " +
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent",
  danger:
    "text-ink-soft hover:text-error border border-transparent " +
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent",
};

// Armed overrides the resting look with the terracotta active fill (mockup
// .tool-btn.armed). Kept as an override layer so any variant can be armed.
const armedClass =
  "bg-accent border border-accent text-paper hover:brightness-110";

export function Button({
  variant = "ghost",
  armed = false,
  className = "",
  children,
  ...props
}: ButtonProps): JSX.Element {
  return (
    <button
      data-variant={variant}
      data-armed={armed ? "true" : undefined}
      aria-pressed={armed ? true : undefined}
      className={`rounded-md px-2.5 py-1 transition-colors ${armed ? armedClass : variantClass[variant]} ${className}`}
      {...props}
    >
      {children}
    </button>
  );
}
```
(Matches the seed exactly — explicit `children`, NO `type` default — plus the `armed`/`armedClass` addition. `rounded-md` == `rounded-sm` == 5px per the collapsed scale; keep `rounded-md` to match the seed. Do NOT introduce a `type` default; the seed has none and Task scope is only `armed`.)

- [ ] **Step 4: Run; verify pass** — `npm test -- Button` → PASS.

- [ ] **Step 5: Story** — append to `src/print/ui/Button.stories.tsx`:
```tsx
export const Armed: Story = { args: { variant: "ghost", armed: true, children: "+ Peak" } };
```

- [ ] **Step 6: Impeccable pass** — checklist C (armed = the "active" state; default/hover/focus already present), B (accent rationed: armed is a genuine active toggle, allowed), D (if "+ Peak" is dense, note container enlarges hit area). Run `npx impeccable detect src/print/ui/Button.tsx` (clean), then `npm run lint:design` + `npx tsc --noEmit -p tsconfig.build.json` → exit 0.

- [ ] **Step 7: Commit**
```bash
git add src/print/ui/Button.tsx src/print/ui/Button.stories.tsx test/print-ui/Button.test.tsx
git commit -m "feat(print-ui): Button armed toggle state (focus-plot + Peak) + impeccable pass"
```

---

### Task 3: Card — `draft` variant

The folio draft card (`.card.is-draft`) uses a dashed `hair-strong` border and removes the Plate-Lift shadow. Add a `draft` boolean.

**Files:**
- Modify: `src/print/ui/Card.tsx`
- Test: `test/print-ui/Card.test.tsx`

- [ ] **Step 1: Write the failing test** — `test/print-ui/Card.test.tsx`:
```tsx
import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Card } from "../../src/print/ui/Card";

describe("<Card> draft variant", () => {
  it("is not draft by default", () => {
    render(<Card data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-draft")).toBe(null);
  });
  it("marks data-draft when draft", () => {
    render(<Card draft data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-draft")).toBe("true");
  });
  it("still honors elevated alongside the existing API", () => {
    render(<Card elevated data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-elevated")).toBe("true");
  });
});
```

- [ ] **Step 2: Run; verify fail** — `npm test -- Card` → FAIL (`draft` not a prop).

- [ ] **Step 3: Implement.** Edit `src/print/ui/Card.tsx`: add `draft?: boolean` to `CardOwnProps`, set `data-draft`, and when draft, swap to a dashed `hair-strong` border and drop the shadow. Add to the own-props type:
```tsx
type CardOwnProps<T extends CardElement> = {
  as?: T;
  elevated?: boolean;
  border?: "hair" | "strong";
  /** Folio draft/recipe card: dashed hair-strong border, no Plate-Lift shadow,
   *  figure dimmed by the consumer. (mockup .card.is-draft) */
  draft?: boolean;
};
```
The live `Card` computes a single `appearance` string: `elevated ? "card" : `rounded-md bg-plate ${borderClass[border]}`` (where `.card` is the Plate-Lift CSS class in `styles.css`, and the module-level `borderClass` record maps `hair`/`strong`). A draft card is a FLAT card (already shadowless) with a dashed `hair-strong` border, and draft overrides elevation (folio draft cards never lift). So change ONLY the `appearance` ternary and the data-attrs — do not invent a `shadow-*` utility (none exists; elevation is the `.card` class):
```tsx
const appearance = draft
  ? "rounded-md bg-plate border border-dashed border-hair-strong"
  : elevated
    ? "card"
    : `rounded-md bg-plate ${borderClass[border]}`;
```
and in the assembled `props`, set `data-draft` when draft and keep `data-elevated` only for a non-draft elevated card:
```tsx
    ...(draft ? { "data-draft": "true" } : {}),
    ...(elevated && !draft ? { "data-elevated": "true" } : {}),
```
Preserve everything else verbatim (the `as` polymorphism, the `type="button"` injection for `as="button"`, the `...rest` spread, the `cx`-free string join).

- [ ] **Step 4: Run; verify pass** — `npm test -- Card` → PASS.

- [ ] **Step 5: Story** — append a `Draft` story to `src/print/ui/Card.stories.tsx` (create the stories file if absent, mirroring `Button.stories.tsx`): export `Default`, `Elevated`, `Draft`.

- [ ] **Step 6: Impeccable pass** — checklist H (draft removes shadow → still flat-except-plate compliant), F (dashed border is structural, not a side-stripe ban violation — it is a full dashed border, allowed). Run `npx impeccable detect` on this task's modified `src/print/ui/*.tsx` (clean), then guards (`lint:design` + `tsc -p tsconfig.build.json`) exit 0.

- [ ] **Step 7: Commit**
```bash
git add src/print/ui/Card.tsx src/print/ui/Card.stories.tsx test/print-ui/Card.test.tsx
git commit -m "feat(print-ui): Card draft variant (folio recipe cards) + impeccable pass"
```

---

### Task 4: SegmentedControl — `xs` size

The focus-plot mini comb toggle (`.seg.mini`) is a step smaller than `sm`: `4px 9px` padding, `10px` text. Add `"xs"` to `SegmentedSize`.

**Files:**
- Modify: `src/print/ui/SegmentedControl.tsx`
- Test: `test/print-ui/SegmentedControl.test.tsx`

- [ ] **Step 1: Write the failing test** — `test/print-ui/SegmentedControl.test.tsx`:
```tsx
import { render, screen } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SegmentedControl } from "../../src/print/ui/SegmentedControl";

const opts = [
  { value: "comb", label: "Comb" },
  { value: "resid", label: "Resid" },
] as const;

describe("<SegmentedControl> xs size", () => {
  it("reflects size=xs on data-size", () => {
    render(
      <SegmentedControl aria-label="view" options={opts} value="comb" onChange={() => {}} size="xs" />,
    );
    expect(screen.getByRole("group").getAttribute("data-size")).toBe("xs");
  });
  it("still toggles via onChange at xs", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl aria-label="view" options={opts} value="comb" onChange={onChange} size="xs" />,
    );
    // second segment
    const btns = screen.getAllByRole("button");
    btns[1].click();
    expect(onChange).toHaveBeenCalledWith("resid");
  });
});
```

- [ ] **Step 2: Run; verify fail** — `npm test -- SegmentedControl` → FAIL (`"xs"` not assignable to `size`).

- [ ] **Step 3: Implement.** Edit `src/print/ui/SegmentedControl.tsx`: extend the size union and the per-size button-padding/text map.
```tsx
export type SegmentedSize = "xs" | "sm" | "md";
```
Find the existing size→class map (the seed maps `sm`/`md` to button padding + text-size utilities). Add an `xs` entry using the smallest on-scale step: `"px-2 py-1 text-xs"` (text-xs = 10.5px ≈ mockup 10px; the 4px/9px padding maps to `py-1 px-2`). Keep `sm`/`md` unchanged. Ensure `data-size={size}` already flows (it does).

- [ ] **Step 4: Run; verify pass** — `npm test -- SegmentedControl` → PASS.

- [ ] **Step 5: Story** — append `MiniXs` to `src/print/ui/SegmentedControl.stories.tsx` (create if absent), showing `size="xs"` alongside `sm`/`md`.

- [ ] **Step 6: Impeccable pass** — checklist D is the sharp one: `xs` is below 44px. Add a code comment that `xs` is a known-dense in-panel toggle whose row provides the effective hit area, mirroring the audit's accepted density exception; confirm focus-visible still renders at `xs`. Run `npx impeccable detect` on this task's modified `src/print/ui/*.tsx` (clean), then guards (`lint:design` + `tsc -p tsconfig.build.json`) exit 0.

- [ ] **Step 7: Commit**
```bash
git add src/print/ui/SegmentedControl.tsx src/print/ui/SegmentedControl.stories.tsx test/print-ui/SegmentedControl.test.tsx
git commit -m "feat(print-ui): SegmentedControl xs size (mini comb toggle) + impeccable pass"
```

---

### Task 5: PhaseChip — coexistence (two-phase) label

The builder renders a coexistence chip as `"Pn3m + Lam"` using the **dominant** phase's color (tint fill + border + text), not a split. Add an optional `coexistWith` second phase; when present, render `"<short(phase)> + <short(coexistWith)>"` with the dominant (`phase`) styling. Short-name mapping mirrors the builder.

**Files:**
- Modify: `src/print/ui/PhaseChip.tsx`
- Test: `test/print-ui/PhaseChip.test.tsx`

- [ ] **Step 1: Write the failing test** — `test/print-ui/PhaseChip.test.tsx`:
```tsx
import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { PhaseChip } from "../../src/print/ui/PhaseChip";

describe("<PhaseChip> coexistence", () => {
  it("renders single phase by default", () => {
    render(<PhaseChip phase="Pn3m" />);
    expect(screen.getByTestId("phase-chip").textContent).toBe("Pn3m");
  });
  it("renders two-phase label with dominant first", () => {
    render(<PhaseChip phase="Pn3m" coexistWith="Lamellar" />);
    const chip = screen.getByTestId("phase-chip");
    expect(chip.textContent).toBe("Pn3m + Lam");
    expect(chip.getAttribute("data-coexist")).toBe("true");
  });
});
```

- [ ] **Step 2: Run; verify fail** — `npm test -- PhaseChip` → FAIL (`coexistWith` not a prop).

- [ ] **Step 3: Implement.** Edit `src/print/ui/PhaseChip.tsx`. Add a short-name map + `coexistWith` prop; when set, render `"<phase> + <short>"` keeping the dominant-phase `variantStyle`. Add:
```tsx
const PHASE_SHORT: Record<string, string> = {
  Pn3m: "Pn3m", Im3m: "Im3m", Ia3d: "Ia3d", Fm3m: "Fm3m",
  Fd3m: "Fd3m", Hexagonal: "Hex", Lamellar: "Lam", Square: "Sq",
};
function short(p: string): string { return PHASE_SHORT[p] ?? p; }
```
Add `coexistWith?: string;` to the interface. In render, compute `const text = coexistWith ? `${short(phase)} + ${short(coexistWith)}` : phase;` and set `data-coexist={coexistWith ? "true" : undefined}`. Color/fill/border stay driven by `phaseColor(phase)` (dominant). Render `{text}` as children. Keep `data-testid="phase-chip"`, `data-variant`, `data-size`.

- [ ] **Step 4: Run; verify pass** — `npm test -- PhaseChip` → PASS.

- [ ] **Step 5: Story** — append `Coexistence: Story = { args: { phase: "Pn3m", coexistWith: "Lamellar" } }` to the existing/new `PhaseChip.stories.tsx`.

- [ ] **Step 6: Impeccable pass** — checklist A is central: the chip is a monospace **label** (second channel = the phase name text), so phase identity never rests on hue alone — confirm the text channel carries both phases. E (mono = data voice, correct). Run `npx impeccable detect` on this task's modified `src/print/ui/*.tsx` (clean), then guards (`lint:design` + `tsc -p tsconfig.build.json`) exit 0.

- [ ] **Step 7: Commit**
```bash
git add src/print/ui/PhaseChip.tsx src/print/ui/PhaseChip.stories.tsx test/print-ui/PhaseChip.test.tsx
git commit -m "feat(print-ui): PhaseChip coexistence two-phase label + impeccable pass"
```

---

## Building the new primitives (Tasks 6–25)

Each follows the per-primitive step template (Conventions section). Specimens are exact mockup values (verified 2026-05-31), mapped to Print token utilities. Where a value is genuinely off the named scale, the task says so and uses the nearest role.

### Task 6: Badge — inline mono count

The notes-button count (`.btn .badge`) is a plain mono span, `text-xs`, `ink-faint`, left margin — no fill/border. Reusable wherever a button carries a count. (This replaces the catalog's mis-filed "IconButton badge" — IconButton itself already matches the mockups and is not changed.)

- **Interface:** `{ children: ReactNode; className?: string }` (extends `Omit<HTMLAttributes<HTMLSpanElement>, "children">`).
- **Render:** `<span data-testid="badge" className={cx("font-mono text-xs text-ink-faint ml-1.5", className)}>{children}</span>`.
- **Test cases:** renders children; has `data-testid="badge"`; forwards a placement `className`.
- **Story:** `Default` (`children: "1"`), and a `WithinButton` story showing `<Button>Notes <Badge>3</Badge></Button>`.
- **Impeccable pass:** C/D N/A (non-interactive). E: mono = a measured count, correct. No fill (flat).

### Task 7: TopBar — fixed 56px header shell

`.topbar`: 56px tall, flex, `gap-4`, `px-6`, `border-b border-hair`, paper background, `flex-shrink-0`.

- **Interface:** `{ wordmark?: ReactNode; children?: ReactNode; rightSlot?: ReactNode; className?: string }`. Layout: wordmark, then `children` (stage tabs / facet / stepper), then a `flex-1` spacer, then `rightSlot`.
- **Render:** a `<header data-testid="topbar" className={cx("h-14 flex-shrink-0 flex items-center gap-4 px-6 border-b border-hair bg-paper", className)}>` with `{wordmark}{children}<div className="flex-1" />{rightSlot}`. (`h-14` = 56px.)
- **Test cases:** renders wordmark + children + rightSlot in DOM; `data-testid="topbar"`; is a `banner` landmark (`<header>` → assert `screen.getByRole("banner")`).
- **Story:** `Default` with a wordmark + a placeholder tabs node + a rightSlot button.
- **Impeccable pass:** standard nav pattern (product permission), flat + hairline (H), banner landmark (a11y).

### Task 8: StageTabs — Samples/Index/Series tabs with dot

`.stage-tab`: uppercase 11px, `letter-spacing 0.04em`, `font-weight 600`, `ink-faint`, `px-[11px] py-1.5 rounded-sm`; active = `ink` text on `paper-sunk`. `.dot`: 5px circle, `hair-strong` inactive / `accent` active, `mr-1.5`.

- **Types:** `export type StageKey = "samples" | "index" | "series";`
- **Interface:** `{ active: StageKey; onChange: (s: StageKey) => void; className?: string }`. Render the three tabs from a fixed `[{key,label}]` array.
- **Render:** container `role="tablist"`, each tab a `<button role="tab" aria-selected={active===key} data-active={...} data-stage={key} className="text-sm font-semibold uppercase tracking-wide px-3 py-1.5 rounded-sm ...">` with a leading `Dot` (reuse the seeded `Dot` primitive: `<Dot tone={active===key?"accent":"neutral"} size="xs" />`). Use `tracking-wide` or the nearest tracking utility; text via `text-sm` (11.5px ≈ 11px) — note the −0.5px deviation in a comment.
- **Test cases:** three tabs render with labels SAMPLES/INDEX/SERIES; active tab has `data-active="true"` + `aria-selected="true"`; clicking a tab calls `onChange` with its key; the active tab's Dot is accent-toned (assert via the Dot's `data-tone="accent"` within the active tab).
- **Story:** one per active stage (`Samples`, `Index`, `Series`).
- **Impeccable pass:** A (active state carries position + dot + text weight, not hue alone), C (hover/focus/active), D (44px: `py-1.5 px-3` ~ check; enlarge if needed), reuse `Dot` (DRY).

### Task 9: FacetChip — pill dropdown trigger

`.facet-chip`: `flex items-center gap-[7px]`, `border border-hair-strong`, `bg-plate`, `rounded-full`, `px-[11px] py-[5px]`; hover `bg-paper-sunk`. Label `.fc-k` (ink, 11px, 600) + chevron `.fc-chev` (ink-faint, "▾").

- **Interface:** `{ label: string; onClick?: () => void; className?: string }`.
- **Render:** `<button data-testid="facet-chip" className="inline-flex items-center gap-2 border border-hair-strong bg-plate rounded-full px-3 py-1 text-sm font-semibold text-ink hover:bg-paper-sunk ...">{label}<span aria-hidden className="text-ink-faint">▾</span></button>`.
- **Test cases:** renders label; `data-testid`; click fires onClick; chevron is `aria-hidden`.
- **Story:** `Default` (`label: "Beamtime"`).
- **Impeccable pass:** rounded-full is intentional pill (not "rounded for charm" — it is a standard filter affordance), C (hover/focus), distinct from FilterChip (Task 21) which is a toggle.

### Task 10: KbKey — keyboard key badge

`.k`: mono, `text-xs` (10.5px), `bg-plate`, `border border-hair-strong`, **2px bottom border** (3D), `rounded-sm`, `px-1.5 py-px`, `ink-soft`.

- **Interface:** `{ children: ReactNode; className?: string }` → renders `<kbd>`.
- **Render:** `<kbd data-testid="kbkey" className={cx("font-mono text-xs bg-plate border border-hair-strong border-b-2 rounded-sm px-1.5 py-px text-ink-soft", className)}>{children}</kbd>`. (The `border-b-2` is a uniform-color thicker bottom border, NOT a colored side-stripe — allowed; it is the physical-key affordance.)
- **Test cases:** renders as `<kbd>` with children; `data-testid`.
- **Story:** `Default` (`children: "⌘K"`), `Letter` (`children: "X"`).
- **Impeccable pass:** G — confirm the `border-b-2` is hairline-colored (not an accent stripe), so it is not a side-stripe ban violation; it is a 3D key, a recognized affordance.

### Task 11: KbLegend — shortcut legend row

`.kb-legend`: `flex flex-wrap gap-5`, `text-sm` (11.5px), `ink-faint`. Each item = a `KbKey` + description text.

- **Types:** `export interface Shortcut { keyLabel: string; description: string; }`
- **Interface:** `{ shortcuts: Shortcut[]; className?: string }`.
- **Render:** `<div data-testid="kb-legend" className={cx("flex flex-wrap gap-5 text-sm text-ink-faint", className)}>` mapping each to `<span className="inline-flex items-center"><KbKey>{s.keyLabel}</KbKey>{s.description}</span>`. Reuse Task 10 `KbKey`.
- **Test cases:** renders N items; each shows key + description; `data-testid`.
- **Story:** `Default` with 3 shortcuts (e.g. `←/→ navigate`, `X reject`, `Esc close`).
- **Impeccable pass:** DRY (composes KbKey), flat, no interactive state (C N/A).

### Task 12: TagPill — single tag pill

`.tag`: `text-xs` (10.5px), `font-semibold`, `ink-soft`, `bg-plate`, `border border-hair`, `rounded-full`, `px-[9px] py-px`, `whitespace-nowrap`. Optional remove affordance.

- **Interface:** `{ children: ReactNode; onRemove?: () => void; className?: string }`.
- **Render:** a `<span data-testid="tag-pill" className="inline-flex items-center gap-1 text-xs font-semibold text-ink-soft bg-plate border border-hair rounded-full px-2 py-px whitespace-nowrap ...">{children}{onRemove && <button aria-label="Remove tag" onClick={onRemove} className="text-ink-faint hover:text-accent">×</button>}</span>`.
- **Test cases:** renders children; `data-testid`; no remove button unless `onRemove`; remove button has `aria-label` and fires `onRemove`.
- **Story:** `Default`, `Removable`.
- **Impeccable pass:** C (the remove ×: hover→accent, focus), B (accent only on the remove hover — rationed grease-pencil).

### Task 13: TagList — wrapping tag list + add invite

`.tags`: `flex items-center flex-wrap gap-[5px]`. Holds `TagPill`s. `.tag-add`: dashed `hair-strong` pill, `ink-faint`, hover→accent; hidden until row hover when tags exist, always visible (`invite`) when empty.

- **Interface:** `{ tags: string[]; onAdd?: () => void; onRemove?: (tag: string) => void; editable?: boolean; className?: string }`.
- **Render:** `<div data-testid="tag-list" className={cx("flex items-center flex-wrap gap-1.5 group/tags", className)}>` mapping `tags` → `<TagPill onRemove={editable && onRemove ? () => onRemove(t) : undefined}>{t}</TagPill>`, then if `onAdd` an add button: `<button data-testid="tag-add" aria-label="Add tag" onClick={onAdd} className={cx("text-xs font-semibold rounded-full border border-dashed border-hair-strong px-2 py-px text-ink-faint hover:text-accent hover:border-accent transition-colors", tags.length>0 && "opacity-0 group-hover/tags:opacity-100")}>+ tag</button>`. (Reveal-on-hover via Tailwind `group` — appearance is a tokenized opacity transition, no layout animation.)
- **Test cases:** renders all tags; add button present only with `onAdd`, has `aria-label`, fires `onAdd`; when `editable` + `onRemove`, each pill removable; empty list still shows the add invite.
- **Story:** `WithTags`, `Empty` (invite), `Editable`.
- **Impeccable pass:** C (add hover/focus), G (opacity transition only, no layout anim), composes TagPill.

### Task 14: EmptyState — centered empty message

`.empty`: `text-center py-[70px] px-5`, `ink-faint`. `.empty h2`: Newsreader 500, 21px, `ink-soft`. Body: `text-base`, `ink-faint`.

- **Interface:** `{ title: string; body?: ReactNode; className?: string }`.
- **Render:** `<div data-testid="empty-state" className={cx("text-center py-16 px-5", className)}><h2 className="text-headline text-ink-soft mb-1.5">{title}</h2>{body && <div className="text-base text-ink-faint">{body}</div>}</div>`. **Decision (token):** the mockup h2 is 21px serif; the nearest existing serif role is `.text-headline` (19px). Use `.text-headline` and note the −2px; do NOT add a new token for a 2px delta (YAGNI; the audit warns against one-off sizes). If review insists on exact 21px, add `--text-empty: 21px` + a `.text-empty` serif role class in `styles.css` and use it — but default to reuse.
- **Test cases:** renders title as a heading (`getByRole("heading", { name })`); renders body when provided; no body node when omitted; `data-testid`.
- **Story:** `FilterNoMatch` (`title: "No series match"`, body copy), `TitleOnly`.
- **Impeccable pass:** E (serif title = correct, it is a plate title moment), product "empty states that teach" — body should guide, not just say "nothing here" (story copy reflects this), no em dashes.

### Task 15: ScreenedMark — completion badge

`.screened-mark`: 13px circle, `flex-shrink-0`; empty = `border-[1.5px] border-hair-strong` transparent; `.done` = `bg-ink border-ink` + white check (inline SVG path, stroke `#faf6ee`≈paper). Second channel: shape (filled + check) not just color.

- **Interface:** `{ done: boolean; className?: string }`.
- **Render:** a `<span data-testid="screened-mark" role="img" aria-label={done ? "Screened" : "Not screened"} data-done={done?"true":undefined} className={cx("inline-block w-[13px] h-[13px] rounded-full flex-shrink-0 border", done ? "bg-ink border-ink" : "border-hair-strong", className)}>` — when `done`, render an inline check `<svg viewBox="0 0 13 13" className="w-full h-full"><path d="M3.4 6.8l2.1 2.1 4.2-4.6" fill="none" stroke="var(--color-paper)" strokeWidth={1.7} strokeLinecap="round" strokeLinejoin="round"/></svg>`. (Use `var(--color-paper)` for the stroke — a token, not a literal; the 13px size is a fixed semantic badge dimension, acceptable as an arbitrary size on a non-text element — note it; if the guard flags `w-[13px]`, it won't, since `print/ui/` is exempt, but keep the comment.)
- **Test cases:** `role="img"` with the right `aria-label` for each state; `data-done` only when done; the check SVG present only when done.
- **Story:** `Screened` (done), `Unscreened`.
- **Impeccable pass:** A is the point — completion carries shape (fill + checkmark) + `aria-label`, not hue alone. F: stroke via token.

### Task 16: RejectOverlay — grease-pencil ✕

The hand-skewed terracotta ✕ over rejected detector frames. `viewBox="0 0 100 100"`, two `<line>`s with rotate(±2), `stroke=accent`, `stroke-width=7`, round caps. Scales with container via `inset-0`.

- **Interface:** `{ className?: string }` (purely decorative overlay; absolutely positioned by consumer, or fills its relative parent).
- **Render (verbatim geometry):**
```tsx
export function RejectOverlay({ className = "" }: { className?: string }): JSX.Element {
  return (
    <svg
      data-testid="reject-overlay"
      aria-hidden="true"
      viewBox="0 0 100 100"
      preserveAspectRatio="none"
      className={cx("absolute inset-0 w-full h-full pointer-events-none", className)}
    >
      <line x1="16" y1="20" x2="86" y2="82" stroke="var(--color-accent)" strokeWidth={7} strokeLinecap="round" transform="rotate(-2 50 50)" />
      <line x1="84" y1="18" x2="14" y2="84" stroke="var(--color-accent)" strokeWidth={7} strokeLinecap="round" transform="rotate(2 50 50)" />
    </svg>
  );
}
```
- **Test cases:** renders an `aria-hidden` svg; `data-testid`; two `<line>` elements.
- **Story:** `OverThumbnail` — render inside a 62px `relative` dark box; `OverFrame` — inside a larger box (proves it scales).
- **Impeccable pass:** A (rejection carries the ✕ shape + frame dimming, not hue alone — the dimming is applied by the consumer; note that in the story), B (accent = the reject grease-pencil mark, a rationed semantic use), `aria-hidden` correct (the state is announced elsewhere).

### Task 17: ProgressBar — accent capsule track

`.pbar`: `h-1` (4px), `rounded-full`, `bg-hair`, `overflow-hidden`; fill `.pbar i`: `h-full bg-accent`, width = value/total %.

- **Interface:** `{ value: number; total: number; label?: string; className?: string }`.
- **Render:** `<div data-testid="progressbar" role="progressbar" aria-valuenow={value} aria-valuemin={0} aria-valuemax={total} aria-label={label} className={cx("h-1 rounded-full bg-hair overflow-hidden", className)}><span className="block h-full bg-accent" style={{ width: `${total>0 ? Math.min(100,(value/total)*100) : 0}%` }} /></div>`.
- **Test cases:** `role="progressbar"` with `aria-valuenow/min/max`; fill width style reflects ratio (assert inline `style.width`); clamps at 100% when value>total; 0% when total=0.
- **Story:** `Partial` (value 8/12), `Full`, `Empty`.
- **Impeccable pass:** A (progress carries position/width, not hue), B (accent = the one progress indicator, semantic), C N/A (non-interactive but has ARIA value).

### Task 18: BonnetBadge — Gauss–Bonnet flag pill

`.bonnet`: `inline-flex items-center gap-[3px]`, `text-[9px]` (below scale — use nearest), `font-bold uppercase tracking`, `text-accent`, border = `color-mix(accent 40%, hair)`, bg = `color-mix(accent 8%, transparent)`, `rounded-full`, `px-1.5 py-px`, leading `⭙`.

- **Interface:** `{ className?: string }` (decorative status pill; no props).
- **Render:** `<span data-testid="bonnet-badge" className={cx("inline-flex items-center gap-0.5 text-xs font-bold uppercase tracking-wide rounded-full px-1.5 py-px", className)} style={{ color: "var(--color-accent)", borderWidth: 1, borderStyle: "solid", borderColor: "color-mix(in oklab, var(--color-accent) 40%, var(--color-hair))", background: "color-mix(in oklab, var(--color-accent) 8%, transparent)" }}><span aria-hidden>⭙</span>Bonnet</span>`. (The 9px mockup size is below the `text-xs` 10.5px floor; use `text-xs` and note the +1.5px — do not add a sub-xs token. The `color-mix` lives in inline `style` with token vars, mirroring how `PhaseChip` does dynamic color; not a raw literal.)
- **Test cases:** renders "Bonnet" text; `data-testid`; the `⭙` glyph is `aria-hidden`.
- **Story:** `Default`, and `InCandidateRow` (inline next to placeholder text).
- **Impeccable pass:** B (accent-tinted, but this is a genuine semantic flag — the Gauss–Bonnet ratio signal; the text "Bonnet" is the second channel so it survives grayscale), A satisfied by the word.

### Task 19: GripHandle — drag handle

`.grip`: `⋮⋮` (two `⋮` U+22EE), `text-base`, `leading-none`, `tracking-[-2px]`, `hair-strong` at rest → `ink-faint` on row hover, `cursor-grab`, `flex-shrink-0`.

- **Interface:** `{ className?: string }` (decorative; drag wiring is the consumer's). Hover reveal is driven by the consumer row's `group` — the primitive exposes the hover class via `group-hover`.
- **Render:** `<span data-testid="grip-handle" aria-hidden="true" className={cx("select-none leading-none text-base text-hair-strong group-hover:text-ink-faint cursor-grab flex-shrink-0 tracking-tighter", className)}>⋮⋮</span>`. (`tracking-tighter` approximates the −2px column packing; note the approximation.)
- **Test cases:** renders `⋮⋮`; `aria-hidden`; `data-testid`.
- **Story:** `Default`, and `InRow` (inside a `group` flex row with hover note).
- **Impeccable pass:** `aria-hidden` correct (drag is not keyboard-primary here; note that reorder must have a keyboard path at the *consumer* level — flag for Phase 2), G (color transition only).

### Task 20: SearchInput — search field

`.search`: `flex items-center gap-[7px]`, `bg-plate`, `border border-hair-strong`, `rounded-sm`, `px-[11px] py-1.5`, `min-w-[230px]`; input transparent, `text-base`, `ink`, placeholder `ink-faint`; leading magnifier SVG.

- **Interface:** `{ value: string; onChange: (v: string) => void; placeholder?: string; className?: string }`.
- **Render:** a `<div data-testid="search-input" className="flex items-center gap-2 bg-plate border border-hair-strong rounded-sm px-3 py-1.5 focus-within:border-accent ...">` with the magnifier `<svg width={14} height={14} viewBox="0 0 14 14" fill="none" aria-hidden><circle cx="6" cy="6" r="4.3" stroke="var(--color-ink-faint)" strokeWidth={1.5}/><line x1="9.2" y1="9.2" x2="12.5" y2="12.5" stroke="var(--color-ink-faint)" strokeWidth={1.5} strokeLinecap="round"/></svg>` and `<input value={value} onChange={(e)=>onChange(e.target.value)} placeholder={placeholder} className="flex-1 bg-transparent border-none outline-none text-base text-ink placeholder:text-ink-faint" />`. Add `focus-within:border-accent` for the focus ring.
- **Test cases:** renders an input with the placeholder; typing fires `onChange` with the new value (`fireEvent.change`); icon is `aria-hidden`; reflects `value`.
- **Story:** `Empty`, `WithValue`.
- **Impeccable pass:** C (focus-within ring in accent), D (44px: `py-1.5` + min-width row — fine), product standard search affordance.

### Task 21: FilterChip — toggle pill

`.chip`: `text-sm` (11.5px), `font-semibold`, `border border-hair-strong`, `bg-plate`, `rounded-full`, `px-3 py-1`, `ink-soft`; hover border→ink-faint; `.on` = `bg-ink text-paper border-ink` (ink-on-paper inversion). Distinct from FacetChip (a dropdown trigger).

- **Interface:** `{ label: string; active: boolean; onClick: () => void; className?: string }`.
- **Render:** `<button data-testid="filter-chip" aria-pressed={active} data-active={active?"true":"false"} onClick={onClick} className={cx("text-sm font-semibold rounded-full px-3 py-1 border transition-colors", active ? "bg-ink text-paper border-ink" : "bg-plate text-ink-soft border-hair-strong hover:border-ink-faint", className)}>{label}</button>`.
- **Test cases:** `aria-pressed` reflects `active`; `data-active` string reflects state; click fires `onClick`; label renders.
- **Story:** `Off`, `On`.
- **Impeccable pass:** A (active carries inversion = fill + text-weight contrast, not hue alone; it is ink/paper, not a color anyway), C (hover/focus/active), D.

### Task 22: Slider — range input

`input[type=range]`: track `h-[3px]` `bg-hair-strong` `rounded-full`; thumb 16px accent disc, 2.5px plate ring + 1px hair-strong outer ring. Includes label row + value.

- **Interface:** `{ value: number; min: number; max: number; step?: number; onChange: (v: number) => void; label?: string; valueDisplay?: ReactNode; className?: string }`.
- **Render:** a wrapper with an optional label/value row (`flex justify-between items-baseline`) and the `<input type="range" data-testid="slider" ... />`. The custom thumb/track needs real CSS (pseudo-elements `::-webkit-slider-thumb` / `::-moz-range-thumb` can't be done with Tailwind utilities). **Author these in `src/styles.css` `@layer components`** as a `.print-range` class (styles.css is not guard-scanned, and this is the legitimate home for pseudo-element styling), then apply `className="print-range"` on the input. Add to `styles.css`:
```css
@layer components {
  .print-range { -webkit-appearance: none; appearance: none; width: 100%; height: 3px;
    background: var(--color-hair-strong); border-radius: 9999px; outline: none; }
  .print-range::-webkit-slider-thumb { -webkit-appearance: none; appearance: none;
    width: 16px; height: 16px; border-radius: 50%; background: var(--color-accent);
    cursor: pointer; border: 2.5px solid var(--color-plate);
    box-shadow: 0 0 0 1px var(--color-hair-strong); }
  .print-range::-moz-range-thumb { width: 16px; height: 16px; border: 2.5px solid var(--color-plate);
    border-radius: 50%; background: var(--color-accent); box-shadow: 0 0 0 1px var(--color-hair-strong); }
  .print-range:focus-visible::-webkit-slider-thumb { outline: 2px solid var(--color-accent); outline-offset: 2px; }
}
```
- **Test cases:** renders a `slider` role input with `min`/`max`/`step`/`value`; change fires `onChange` with the numeric value; label + valueDisplay render when given.
- **Story:** `Default` (offset 0.4–1.4), `WithValue`.
- **Impeccable pass:** C (focus-visible thumb outline added), B (accent thumb = the live control, semantic), motion none (G).

### Task 23: ToggleSwitch — pill switch

32×18 track, `rounded-full`; off `bg-hair-strong`, on `bg-accent`; 14px plate disc, 2px inset, translateX(14) on. Label beside it. Pseudo-free: use a styled checkbox or a button-driven span.

- **Interface:** `{ checked: boolean; onChange: (v: boolean) => void; label: string; className?: string }`.
- **Render (button + spans, no pseudo-element so it's testable & tokenized):**
```tsx
<label data-testid="toggle-switch" className={cx("inline-flex items-center gap-2.5 cursor-pointer", className)}>
  <button type="button" role="switch" aria-checked={checked} aria-label={label}
    onClick={() => onChange(!checked)}
    className={cx("relative w-8 h-[18px] rounded-full transition-colors flex-shrink-0",
      checked ? "bg-accent" : "bg-hair-strong")}>
    <span className={cx("absolute top-0.5 left-0.5 w-3.5 h-3.5 rounded-full bg-plate shadow-sm transition-transform",
      checked && "translate-x-3.5")} />
  </button>
  <span className="text-base font-semibold text-ink">{label}</span>
</label>
```
(`w-8`=32, `h-[18px]` arbitrary height OK in print/ui; `w-3.5`=14, `translate-x-3.5`=14. The transition is `transform`/`color` only — allowed.)
- **Test cases:** `role="switch"` with `aria-checked` reflecting `checked`; click fires `onChange(!checked)`; label renders and is the `aria-label`.
- **Story:** `Off`, `On`.
- **Impeccable pass:** A (switch carries position + role, not hue), C (focus on the switch button — add `focus-visible:outline` accent), D (the 32×18 control is small — the full `<label>` row is the hit target; note it), G (transform transition only).

### Task 24: MetaList — mono key/value list

`.meta-list`: `flex flex-col gap-1.5`; row `flex justify-between font-mono text-sm` (11.5px); key `ink-faint`, value `ink`.

- **Types:** `export interface MetaEntry { key: string; value: ReactNode; }`
- **Interface:** `{ entries: MetaEntry[]; className?: string }`.
- **Render:** `<dl data-testid="meta-list" className={cx("flex flex-col gap-1.5 font-mono text-sm", className)}>` mapping each to `<div className="flex justify-between"><dt className="text-ink-faint">{e.key}</dt><dd className="text-ink">{e.value}</dd></div>`. (Use `<dl>/<dt>/<dd>` for semantic key/value.)
- **Test cases:** renders each key + value; uses `dt`/`dd` (assert `getByRole("term")`/`definition` or query `dt`); `data-testid`.
- **Story:** `Default` with frame/integration/collected/signal rows (the last value a placeholder for SignalBars).
- **Impeccable pass:** E (mono = measured values, correct), semantic `<dl>` (a11y), flat.

### Task 25: SignalBars — 5-bar strength indicator

`.signal-bars`: `inline-flex gap-0.5`; 5 bars `w-[5px] h-[11px] rounded-[1px]`; active `bg-ink-soft`, inactive `bg-hair-strong`; active count = round(value/max*5). Second channel: the *count* of filled bars (position), not hue.

- **Interface:** `{ value: number; max?: number; className?: string }` (`max` default 5; `value` is 0–max, or a 0–1 fraction × max — define: `value` is on the same scale as `max`).
- **Render:** compute `const on = Math.max(0, Math.min(5, Math.round((value / (max ?? 5)) * 5)));` then `<span data-testid="signal-bars" role="img" aria-label={`Signal ${on} of 5`} className={cx("inline-flex gap-0.5 items-end", className)}>` with 5 `<i className={cx("w-[5px] h-[11px] rounded-[1px]", i < on ? "bg-ink-soft" : "bg-hair-strong")} />`.
- **Test cases:** renders exactly 5 bars; `on` count reflects value/max (e.g. value 4 max 5 → 4 active — assert by counting elements whose class includes the active token via `data-on` attr instead: set `data-on={i<on?"true":"false"}` and count); `role="img"` with the `aria-label`; clamps to 0..5. (Set `data-on` so the test avoids class assertions.)
- **Story:** `Strong` (4/5), `Weak` (1/5), `Full`.
- **Impeccable pass:** A (strength carries the count/position + `aria-label`, never hue alone — ink-soft vs hair-strong is tonal, not chromatic, and the count is the real channel), C N/A.

---

### Task 26: Storybook review assembly + whole-set impeccable gate

Confirm the full primitive set is browsable in Storybook and clears the deterministic scan as a set, then launch it for the user's review. This is the plan's review gate — by its end, every primitive built here has a story and renders with live Print tokens.

**Files:**
- (No new source.) Verification + a short index doc.
- Create: `docs/greenfield-rebuild/primitive-storybook-index.md`

- [ ] **Step 1: Confirm every primitive has a story**

Run (from `packages/HimalayaUI/frontend`):
```bash
echo "primitives:"; ls src/print/ui/*.tsx | grep -v '\.stories\.tsx$' | xargs -n1 basename | sed 's/\.tsx$//' | sort
echo "stories:"; ls src/print/ui/*.stories.tsx | xargs -n1 basename | sed 's/\.stories\.tsx$//' | sort
```
Expected: every primitive (except the pure-barrel `index.ts`, `peakMark.ts`, and any non-visual helper) has a matching `*.stories.tsx`. The set must include the 4 hardened (Button, Card, SegmentedControl, PhaseChip) + 20 new (Badge, TopBar, StageTabs, FacetChip, KbKey, KbLegend, TagPill, TagList, EmptyState, ScreenedMark, RejectOverlay, ProgressBar, BonnetBadge, GripHandle, SearchInput, FilterChip, Slider, ToggleSwitch, MetaList, SignalBars) + the pre-existing seeded stories. List any primitive missing a story and add it before proceeding.

- [ ] **Step 2: Whole-set deterministic scan**

Run: `npx impeccable detect src/print/ui/` → expect exit 0 / no hits across the entire primitive tree. Fix any straggler (the per-task passes should already keep this clean).

- [ ] **Step 3: Build Storybook headlessly**

Run: `npm run build-storybook`
Expected: PASS — emits `storybook-static/`. Confirm: `test -d storybook-static && echo OK`.

- [ ] **Step 4: Write the review index**

Create `docs/greenfield-rebuild/primitive-storybook-index.md` listing each primitive, its Storybook title (`ui/<Name>`), the states/variants its story shows, and a one-line note on the mockup it came from. This is the reviewer's map.

- [ ] **Step 5: Launch Storybook for review**

Run (foreground, the user opens it): `npm run storybook` → http://localhost:6006. Tell the user it's up; the sidebar `ui/*` entries are the full Phase-1 primitive set rendering on the warm-paper field with real Print tokens. (Leave running until the user is done reviewing; stop on request.)

- [ ] **Step 6: Commit the index**
```bash
git add docs/greenfield-rebuild/primitive-storybook-index.md
git commit -m "docs(rebuild): primitive Storybook index — Phase-1 review map"
```

---

## Self-Review

**1. Spec coverage** (against `2026-05-31-greenfield-ui-rebuild-design.md` Phase 1 "hone + expand … each primitive verified in Storybook" + `design-system-catalog.md`):
- Harden the 4 gapped existing primitives (Button armed, Card draft, SegmentedControl xs, PhaseChip coexist) → Tasks 2–5. ✓ (The catalog's 5th gap, "IconButton notes badge," is correctly re-resolved to a new `Badge` atom (Task 6) since the specimen is a Button+span, not an IconButton change.)
- Build the true atomic/layout/form primitives → Tasks 6–25 (Badge, TopBar, StageTabs, FacetChip, KbKey, KbLegend, TagPill, TagList, EmptyState, ScreenedMark, RejectOverlay, ProgressBar, BonnetBadge, GripHandle, SearchInput, FilterChip, Slider, ToggleSwitch, MetaList, SignalBars). ✓
- Each verified in Storybook → every task Step 5; full set browsable + launched for review → Task 26. ✓
- Impeccable criteria built in from the start → Task 1 checklist (incl. the `npx impeccable detect` deterministic floor, section J) + per-task Impeccable pass step (detect + checklist + guards) + whole-set detect gate (Task 26) + (execution) impeccable-review stage. ✓
- DEFERRED correctly (not missing): `Sparkline`/`MiniWaterfall` (renderers needing the intensity model / real trace data via `lib/plot/sparkline.ts`) → Phase 2; all composites (DataTable, FolioCard, CombPanel, PhasecallBlock, CandidateRow, ScopingPlate, NotesMargin, ExposureSwitcher, Stepper, …) → Phase 2; token additions (`--color-unindexed`, `--color-peak-manual`, large serif display step) → the Phase-2/3 surfaces that need them; phase colors stay single-sourced in `phases.ts` (no `@theme` duplication). ✓

**2. Placeholder scan:** Every code task gives a concrete interface + implementation + story + test cases. The few "match the seed's existing utility" notes (Card elevation class, SegmentedControl size-map entry) are explicit instructions to read one named line, not vague TODOs — acceptable because the seeded files are the live source and the exact target line is named. No "add error handling"/"similar to Task N"/"TBD".

**3. Type/name consistency:** `StageKey` (Task 8) used only there. `Shortcut` (Task 11) consumed by KbLegend. `MetaEntry` (Task 24) consumed by MetaList. `Badge` (Task 6) composed by TopBar/Button usage examples. `KbKey` (Task 10) composed by KbLegend (Task 11). `TagPill` (Task 12) composed by TagList (Task 13). `Dot` (seeded) reused by StageTabs. All composition references point at primitives built in an earlier task number (KbKey<KbLegend, TagPill<TagList, Dot<StageTabs). Barrel export form (`export { X } from "./X"`) matches the seeded `index.ts`. `data-*` attribute names are unique per primitive and asserted in tests (no class-string assertions).

**4. Token discipline:** Only `styles.css` `@layer components` additions are the `.print-range` block (Task 22) and the optional `.text-empty` (Task 19, default = reuse `.text-headline`). No `@theme` color tokens added; `styles.css` is not guard-scanned so these are permitted; both are real-CSS needs (pseudo-elements / a serif role) that cannot be expressed as Tailwind utilities.

**5. Ordering:** Composition deps go forward-only (KbKey→KbLegend, TagPill→TagList). Tasks 2–5 (hardening) precede the new builds; Task 1 (checklist) precedes everything since every later Impeccable-pass step references it.
