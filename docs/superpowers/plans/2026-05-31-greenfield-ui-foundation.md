# Greenfield UI Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the clean-room `src/print/` greenfield tree — a runnable second Vite entry, a mechanical `no-legacy-import` build guard, the old UI excluded from the new build, the design-system primitives seeded from today's clean extraction, and Storybook rendering them — so all later component/page rebuilds happen in an environment where the old UI is unreachable.

**Architecture:** A new `src/print/` subtree holds the new app (`main.tsx`, `App.tsx`, `ui/`, later `components/`/`pages/`). It imports the unchanged logic core (`src/api.ts`, `src/queries.ts`, `src/state.ts`, `src/phases.ts`, `src/lib/**`, `src/hooks/**`) by relative path but is forbidden — by a hard build-failing rule in `scripts/check-design.mjs` — from importing `src/components/**` or `src/pages/**`. The old app (`index.html` → `src/main.tsx`) is left untouched and is dropped from the production build + typecheck; it is deleted wholesale at a later cutover.

**Tech Stack:** React 18, Vite, TypeScript strict, TailwindCSS v4 (`@tailwindcss/vite`), Vitest, Storybook v10 (`@storybook/react-vite`). Zero-dependency Node guard script (`scripts/check-design.mjs`).

**Spec:** `docs/superpowers/specs/2026-05-31-greenfield-ui-rebuild-design.md`

**Working directory for all commands:** `packages/HimalayaUI/frontend` (the frontend package root). Paths below are relative to it unless noted.

---

## File Structure

Created:
- `docs/greenfield-rebuild/inventory-core-purity.md` — ledger: shared core is import-clean; renderer geometry functions to reuse
- `docs/greenfield-rebuild/design-system-catalog.md` — every reusable primitive/pattern extracted from the 7 mockups + gap list
- `docs/greenfield-rebuild/surface-map.md` — each mockup → new page + component breakdown
- `packages/HimalayaUI/frontend/print.html` — second Vite entry
- `packages/HimalayaUI/frontend/src/print/main.tsx` — new app bootstrap (mounts `#print-root`)
- `packages/HimalayaUI/frontend/src/print/App.tsx` — minimal greenfield shell (full routing/SSE wiring deferred to the pages plan)
- `packages/HimalayaUI/frontend/src/print/ui/**` — primitives seeded (copied) from `src/components/ui/**`
- `packages/HimalayaUI/frontend/src/print/ui/Button.stories.tsx` — exemplar Storybook story
- `packages/HimalayaUI/frontend/.storybook/main.ts`, `.storybook/preview.ts` — Storybook config
- `packages/HimalayaUI/frontend/test/print-shell.test.tsx` — smoke test for the shell

Modified:
- `packages/HimalayaUI/frontend/scripts/check-design.mjs` — add `scanLegacyImports`, extend `isExcluded` to `print/ui/`
- `packages/HimalayaUI/frontend/test/check-design.test.ts` — cases for the new rule + exclusion
- `packages/HimalayaUI/frontend/vite.config.ts` — `build.rollupOptions.input` → print entry (ESM-safe path)
- `packages/HimalayaUI/frontend/package.json` — Storybook devDeps + scripts

(`tsconfig.build.json` is intentionally NOT modified — see Task 4 Step 3.)

Scope note — DEFERRED to follow-up plans (do NOT do here): the full primitive-hardening sweep against the catalog, and the rebuild of components/pages. This plan delivers the *foundation only*: a runnable, isolated, Storybook-backed shell with seeded (not-yet-expanded) primitives.

---

### Task 1: Core-purity ledger (de-risk the import boundary)

Confirm the shared core the new tree will import does not itself reach into old UI; record the renderer geometry functions to reuse. If anything is entangled, this surfaces it before later tasks depend on it.

**Files:**
- Create: `docs/greenfield-rebuild/inventory-core-purity.md`

- [ ] **Step 1: Scan shared core for forbidden imports**

Run (from `packages/HimalayaUI/frontend`):

```bash
grep -rnE 'from "(\.\./)+(components|pages)/' \
  src/api.ts src/queries.ts src/state.ts src/phases.ts \
  src/lib src/hooks
```

Expected (KNOWN, verified 2026-05-31): `api.ts`/`queries.ts`/`state.ts`/`phases.ts`/`hooks/**` are clean. **`src/lib/figure-export/**` is NOT** — it imports old renderer components:
- `lib/figure-export/marks/traceExportMarks.ts` → `components/ui/peakMark`
- `lib/figure-export/marks/multiTraceExportMarks.ts` → `components/MemberHeatmapLayer`, `components/ui/peakMark`, `components/CrossTraceTrackingLayer`, `components/RepresentationToggle`
- `lib/figure-export/adapters/multiTraceAdapter.ts` → `components/RepresentationToggle` (type-only)

This coupling is **expected and not a foundation blocker**: `figure-export` reuses renderer *mark* logic, so it is genuinely coupled to the visual layer (it is not pure "shared core"). The foundation `src/print/` tree does **not** import `figure-export`, so the clean room holds for this plan. The coupling is resolved later — when the renderer components are rebuilt in `src/print/`, `figure-export` is repointed at the `print/` versions (note: `peakMark` is seeded into `src/print/ui/` in Task 6, so a `print/` repoint target already exists for it). Record this in the ledger; do **not** try to extract it now. Any *other* hit (outside `figure-export/`) IS a surprise and must be recorded for extraction before Task 5.

- [ ] **Step 2: Enumerate renderer geometry functions to reuse**

Run:

```bash
grep -rnE '^export (function|const) ' src/lib/plot src/lib/qRing.ts src/lib/detectorOrient.ts \
  | sed 's/{$//'
```

- [ ] **Step 3: Write the ledger**

Create `docs/greenfield-rebuild/inventory-core-purity.md` with these sections, filled from Steps 1–2:

```markdown
# Core-purity ledger (greenfield rebuild)

## Shared core importable by src/print/
- src/api.ts, src/queries.ts, src/state.ts, src/phases.ts
- src/lib/** (list any file that is NOT pure logic here; expected: none)
- src/hooks/** (logic hooks)
- src/ErrorBoundary.tsx (shared)

## Forbidden from src/print/ (enforced by no-legacy-import guard)
- src/components/** (incl. old src/components/ui/**)
- src/pages/**

## Known coupling — lib/figure-export → renderer components (NOT pure core)
figure-export imports old renderer components (peakMark, MemberHeatmapLayer,
CrossTraceTrackingLayer, RepresentationToggle). Expected: it reuses renderer mark
logic. NOT a foundation blocker (print/ does not import figure-export). Resolved
when renderers are rebuilt in print/ and figure-export is repointed.
<paste exact Step-1 figure-export hits here>

## Surprise entanglements to extract (must be empty before Task 5)
<paste any Step-1 hits OUTSIDE lib/figure-export/, or "none">

## Renderer geometry functions to reuse (Option 2: rebuild visuals, reuse math)
<paste Step-2 export inventory, grouped by file>
```

- [ ] **Step 4: Commit**

```bash
git add docs/greenfield-rebuild/inventory-core-purity.md
git commit -m "docs(rebuild): core-purity ledger — shared core is import-clean"
```

---

### Task 2: Design-system catalog + surface map (mockup extraction)

Walk all 7 mockups in `docs/redesign-mockups/` and `DESIGN.md`; extract *every* reusable primitive/pattern and map each mockup to its new page + components. This is the foundation knowledge the primitive-hardening plan consumes. (During execution this reading-heavy task may be delegated to a subagent; the deliverable structure is fixed below.)

**Files:**
- Create: `docs/greenfield-rebuild/design-system-catalog.md`
- Create: `docs/greenfield-rebuild/surface-map.md`
- Read: `docs/redesign-mockups/{sample-table,focus-workspace,series-folio,series-scoping,series-builder}.html`, `docs/redesign-mockups/2026-05-29-{focus,series}-plot.html`, `DESIGN.md`

- [ ] **Step 1: List the existing seeded primitives to cross-reference against**

Run (from `packages/HimalayaUI/frontend`):

```bash
ls src/components/ui/*.tsx | xargs -n1 basename
```

These 13 are the baseline. The catalog's job is to confirm each is in the mockups and to list **gaps** — recurring mockup patterns NOT yet covered by a primitive.

- [ ] **Step 2: Write the catalog**

Create `docs/greenfield-rebuild/design-system-catalog.md`:

```markdown
# Design-system catalog (extracted from The Print mockups)

Source mockups: docs/redesign-mockups/*.html + DESIGN.md
Baseline primitives: src/components/ui/** (seeded into src/print/ui/ in Task 6)

## Existing primitives — confirmed in mockups
For each of the 13 baseline primitives: name, where it appears in the mockups,
any variant/state the mockup shows that the current primitive lacks.

| Primitive | Appears in | Variants/states seen | Gap vs current |
|---|---|---|---|
| Button | … | … | … |
| … | | | |

## Gap primitives — recurring patterns with NO current primitive
For each: proposed name, what it is, which mockups use it, a one-line prop sketch.

| Proposed primitive | Description | Mockups | Prop sketch |
|---|---|---|---|
| … | | | |

## Tokens
Any color/spacing/type role used in the mockups not yet in styles.css @theme.
```

- [ ] **Step 3: Write the surface map**

Create `docs/greenfield-rebuild/surface-map.md`:

```markdown
# Surface map (mockup → new page → components)

| Mockup | New page (src/print/pages/) | Key components (src/print/components/) | Renderer(s) |
|---|---|---|---|
| sample-table.html | ContactSheetPage | … | — |
| sample-table.html (loupe section) | LoupePage (own route) | … | DetectorImage |
| focus-workspace.html + focus-plot.html | FocusWorkspacePage | … | TraceViewer, MillerPlot |
| series-folio.html | SeriesFolioPage | … | — |
| series-scoping.html | SeriesScopingPage | … | — |
| series-builder.html + series-plot.html | SeriesBuilderPage | … | MultiTracePlot |
```

- [ ] **Step 4: Commit**

```bash
git add docs/greenfield-rebuild/design-system-catalog.md docs/greenfield-rebuild/surface-map.md
git commit -m "docs(rebuild): design-system catalog + surface map from mockups"
```

---

### Task 3: no-legacy-import guard (TDD)

Add a hard build-failing rule: any file under `src/print/**` importing (relatively) from `src/components/**` or `src/pages/**` is a violation. Also extend the appearance-rule exclusion so `src/print/ui/**` is design-authoring (like `src/components/ui/**`).

**Files:**
- Modify: `scripts/check-design.mjs`
- Test: `test/check-design.test.ts`

- [ ] **Step 1: Write the failing tests**

Add to `test/check-design.test.ts` (it already imports from `../scripts/check-design.mjs`; add `scanLegacyImports` to that import). Append:

```ts
import { scanLegacyImports } from "../scripts/check-design.mjs";

describe("scanLegacyImports", () => {
  it("flags a print/ file importing old components (one level up to src)", () => {
    const v = scanLegacyImports("print/App.tsx", `import { X } from "../components/X";`);
    expect(v).toEqual([
      { rule: "no-legacy-import", file: "print/App.tsx", line: 1, text: "../components/X" },
    ]);
  });

  it("flags a print/ui file importing old components (two levels up to src)", () => {
    const v = scanLegacyImports("print/ui/Button.tsx", `import { X } from "../../components/X";`);
    expect(v.map((h) => h.text)).toEqual(["../../components/X"]);
  });

  it("flags importing old pages", () => {
    const v = scanLegacyImports("print/pages/P.tsx", `import P from "../../pages/Old";`);
    expect(v.map((h) => h.rule)).toEqual(["no-legacy-import"]);
  });

  it("flags importing the OLD ui primitives (under components/)", () => {
    const v = scanLegacyImports("print/ui/Card.tsx", `import { Dot } from "../../components/ui/Dot";`);
    expect(v).toHaveLength(1);
  });

  it("allows importing the shared core", () => {
    const src = [
      `import { api } from "../api";`,
      `import { phaseColor } from "../phases";`,
      `import { foo } from "../lib/queue/types";`,
      `import { useFocusTrap } from "../hooks/useFocusTrap";`,
    ].join("\n");
    expect(scanLegacyImports("print/App.tsx", src)).toEqual([]);
  });

  it("allows print-internal imports (its own components/pages)", () => {
    const v = scanLegacyImports("print/App.tsx", `import { Page } from "./pages/ContactSheet";\nimport { C } from "./components/C";`);
    expect(v).toEqual([]);
  });

  it("ignores non-print files entirely", () => {
    expect(scanLegacyImports("components/Foo.tsx", `import x from "../pages/Bar";`)).toEqual([]);
  });
});

describe("isExcluded via scanContent — print/ui authoring", () => {
  it("excludes src/print/ui from appearance rules", () => {
    expect(scanContent("print/ui/PhaseChip.tsx", `<span className="text-[11px]" />`)).toEqual([]);
  });
  it("still enforces appearance rules on print/ consumers", () => {
    const v = scanContent("print/components/Foo.tsx", `<span className="text-[11px]" />`);
    expect(v.map((h) => h.rule)).toContain("no-arbitrary-text");
  });
});
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `npm test -- check-design`
Expected: FAIL — `scanLegacyImports is not a function` and the print/ui exclusion cases fail.

- [ ] **Step 3: Implement**

In `scripts/check-design.mjs`:

(a) Extend the path import at the top:

```js
import { join, relative, sep, posix } from "node:path";
```

(b) Extend `isExcluded` so `print/ui/` is also design-authoring:

```js
// src/components/ui/** and src/print/ui/** are excluded entirely (where appearance is authored).
function isExcluded(relPath) {
  return relPath.startsWith("components/ui/") || relPath.startsWith("print/ui/");
}
```

(c) Add the new scanner (after `scanContent`):

```js
// Relative import/export specifiers on one line. Group 1 = `... from "X"`, group 2 = dynamic import("X").
const IMPORT_SPEC_RE =
  /(?:import|export)\b[^'"]*?\bfrom\s*["']([^"']+)["']|import\s*\(\s*["']([^"']+)["']\s*\)/g;

// Import-boundary guard: a file under src/print/** may not import (relatively) from the OLD
// src/components/** or src/pages/**. We resolve each relative specifier against the importer's
// dir (POSIX) and flag it if it lands under components/ or pages/. Print-internal imports
// (./components, ./pages — i.e. src/print/components, src/print/pages) resolve under print/ and pass.
// Exported pure for unit testing.
export function scanLegacyImports(relPath, content) {
  if (!relPath.startsWith("print/")) return [];
  const dir = posix.dirname(relPath);
  const out = [];
  const lines = content.split("\n");
  for (let i = 0; i < lines.length; i++) {
    IMPORT_SPEC_RE.lastIndex = 0;
    let m;
    while ((m = IMPORT_SPEC_RE.exec(lines[i])) != null) {
      const spec = m[1] ?? m[2];
      if (!spec || !spec.startsWith(".")) continue; // bare/external specifiers are fine
      const resolved = posix.normalize(posix.join(dir, spec));
      if (resolved.startsWith("components/") || resolved.startsWith("pages/")) {
        out.push({ rule: "no-legacy-import", file: relPath, line: i + 1, text: spec });
      }
    }
  }
  return out;
}
```

(d) Aggregate it in `runCli` — change the per-file push to include both scanners:

```js
  for (const abs of files) {
    const rel = relToSrc(abs);
    const content = readFileSync(abs, "utf8");
    all.push(...scanContent(rel, content));
    all.push(...scanLegacyImports(rel, content));
  }
```

(e) Update the CLI hint string to mention the new rule:

```js
    process.stderr.write(
      "Move the appearance utility into src/components/ui/** (or src/print/ui/**), or use a named " +
        "scale/role token. Raw color literals belong only in the color-authoring files. " +
        "src/print/** may not import from src/components/** or src/pages/** (no-legacy-import).\n",
    );
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `npm test -- check-design`
Expected: PASS (all new + existing cases).

- [ ] **Step 5: Commit**

```bash
git add scripts/check-design.mjs test/check-design.test.ts
git commit -m "feat(guard): no-legacy-import rule + print/ui design-authoring exclusion"
```

---

### Task 4: Second Vite entry (build only the new app)

Add the `print.html` entry and make the production build emit ONLY the new app. The old app is left for dev reference. NOTE: we do **not** change `tsconfig.build.json` in this plan — see Step 3 for why excluding old from the typecheck is both unnecessary now and ineffective (defeated by `lib/figure-export` → components coupling, verified 2026-05-31).

**Files:**
- Create: `print.html`
- Modify: `vite.config.ts`

- [ ] **Step 1: Create `print.html`**

Mirror `index.html`, but a distinct root id and the print entry:

```html
<!doctype html>
<html lang="en">
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>Himalaya — The Print</title>
  </head>
  <body>
    <div id="print-root"></div>
    <script type="module" src="/src/print/main.tsx"></script>
  </body>
</html>
```

- [ ] **Step 2: Point the production build at the print entry**

`vite.config.ts` is ESM (`package.json` has `"type": "module"`), so **`__dirname` is undefined** — do not use it. Use `fileURLToPath(new URL(...))`. Add the import at the top of `vite.config.ts`:

```ts
import { fileURLToPath } from "node:url";
```

Then replace the existing `build` block with:

```ts
  build: {
    outDir: "dist",
    emptyOutDir: true,
    sourcemap: true,
    // The production build emits ONLY the new greenfield app. The old
    // index.html → src/main.tsx is left for dev reference until cutover; it is
    // not built (the rebuild does not preserve old-app usability). In dev,
    // `vite` still serves index.html by path regardless of this input.
    rollupOptions: {
      input: { print: fileURLToPath(new URL("./print.html", import.meta.url)) },
    },
  },
```

- [ ] **Step 3: Do NOT modify `tsconfig.build.json` (rationale)**

Leave `tsconfig.build.json` unchanged. Two reasons (verified 2026-05-31):
- **Unnecessary now:** old code currently typechecks green and this plan only *adds* `src/print/`. `tsc` over old+core+new stays green; old isn't touched.
- **Ineffective anyway:** TS `exclude` does not remove a file that an *included* file imports. `src/lib/figure-export/**` (included shared core) imports old `components/**`, so excluding `src/components` would not actually drop the old typecheck — it would silently pull old components (and their closures) back in.

The typecheck-time exclusion of old becomes relevant only later, when old genuinely rots — and it must be paired with severing/repointing `lib/figure-export`'s component imports. That is a later-phase task, recorded in the Task 1 ledger, not foundation work.

- [ ] **Step 4: Verify (deferred until Task 5 creates the entry files)**

`npm run build` cannot pass until `src/print/main.tsx`/`App.tsx` exist (Task 5). No commit yet — proceed to Task 5, then verify+commit there.

---

### Task 5: Greenfield shell `src/print/` (TDD smoke)

Create the minimal new-app bootstrap and a placeholder shell, proving the entry builds and renders. Full routing/SSE/queue wiring is deferred to the pages plan.

**Files:**
- Create: `src/print/App.tsx`
- Create: `src/print/main.tsx`
- Test: `test/print-shell.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/print-shell.test.tsx`:

```tsx
import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { PrintApp } from "../src/print/App";

describe("PrintApp shell", () => {
  it("renders the greenfield shell landmark", () => {
    render(<PrintApp />);
    expect(screen.getByTestId("print-shell")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- print-shell`
Expected: FAIL — cannot resolve `../src/print/App`.

- [ ] **Step 3: Implement the shell**

Create `src/print/App.tsx` (minimal, provider-free so it renders standalone in the test; appearance uses only tokenized utilities so `lint:design` passes):

```tsx
/**
 * PrintApp — root of the greenfield ("The Print") UI tree.
 *
 * Foundation placeholder. Routing, SSE wiring, and the mutation-queue
 * persistence/rehydrate effects (mirrored from the old src/App.tsx) are added
 * in the pages-assembly plan. For now this proves the src/print/ entry builds,
 * renders, and is reachable via print.html.
 */
export function PrintApp(): JSX.Element {
  return (
    <main data-testid="print-shell" className="min-h-screen bg-paper text-ink">
      <h1 className="text-display">The Print — greenfield shell</h1>
    </main>
  );
}
```

Create `src/print/main.tsx` (mirrors `src/main.tsx` providers; mounts `#print-root`):

```tsx
import "../styles.css";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { BrowserRouter } from "react-router-dom";
import { PrintApp } from "./App";
import { ErrorBoundary } from "../ErrorBoundary";

const queryClient = new QueryClient({
  defaultOptions: {
    queries: { staleTime: 30_000, retry: 1, refetchOnWindowFocus: false },
  },
});

const root = document.getElementById("print-root");
if (!root) throw new Error("#print-root missing");
createRoot(root).render(
  <StrictMode>
    <ErrorBoundary>
      <QueryClientProvider client={queryClient}>
        <BrowserRouter>
          <PrintApp />
        </BrowserRouter>
      </QueryClientProvider>
    </ErrorBoundary>
  </StrictMode>,
);
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npm test -- print-shell`
Expected: PASS.

- [ ] **Step 5: Verify the full build (Task 4 + 5 together)**

Run: `npm run build`
Expected: PASS — `lint:design` clean (incl. no-legacy-import), `tsc` clean on new tree + core, `vite build` emits `dist/print.html`. Confirm: `test -f dist/print.html && echo OK`.

- [ ] **Step 6: Commit**

```bash
git add print.html vite.config.ts src/print/App.tsx src/print/main.tsx test/print-shell.test.tsx
git commit -m "feat(print): greenfield shell + second Vite entry (old UI dropped from build output)"
```

---

### Task 6: Seed design-system primitives into `src/print/ui/`

Copy today's clean `ui/` primitives into the new tree as the starting point (hardening/expansion is a follow-up plan). A straight copy preserves all relative imports — `src/components/ui/` and `src/print/ui/` are the same depth under `src/`, and the primitives only import siblings + `../../{phases,lib/toast,hooks/useFocusTrap}` (verified).

**Files:**
- Create: `src/print/ui/**` (copied from `src/components/ui/**`)

- [ ] **Step 1: Copy the primitives**

Run (from `packages/HimalayaUI/frontend`):

```bash
mkdir -p src/print/ui
cp -R src/components/ui/. src/print/ui/
ls src/print/ui
```

Expected: all primitives + `index.ts` + `peakMark.ts` present.

- [ ] **Step 2: Verify no seeded primitive trips the guards or the typecheck**

Run: `npm run lint:design`
Expected: exit 0 — `print/ui/**` is design-authoring (Task 3 exclusion) so internal appearance utilities are allowed, and none of the primitives import `components/`/`pages/` (they import siblings + shared core).

Run: `npx tsc --noEmit -p tsconfig.build.json`
Expected: exit 0 — `src/print/ui/**` typechecks as part of the new tree.

- [ ] **Step 3: Commit**

```bash
git add src/print/ui
git commit -m "feat(print): seed design-system primitives into src/print/ui from clean extraction"
```

---

### Task 7: Storybook (component lab) + exemplar story

Install Storybook v10 (`@storybook/react-vite`), point it at `src/print/**`, import the Print stylesheet so primitives render with real tokens, and prove the harness with a Button story.

**Files:**
- Modify: `package.json`
- Create: `.storybook/main.ts`, `.storybook/preview.ts`
- Create: `src/print/ui/Button.stories.tsx`

- [ ] **Step 1: Install Storybook**

Run (from `packages/HimalayaUI/frontend`):

```bash
npm install -D storybook@^10 @storybook/react-vite@^10
```

- [ ] **Step 2: Add scripts to `package.json`**

In the `scripts` block, add:

```json
    "storybook": "storybook dev -p 6006 --no-open",
    "build-storybook": "storybook build"
```

- [ ] **Step 3: Configure `.storybook/main.ts`**

Create `.storybook/main.ts` — stories scoped to `src/print/**` only (keeps Storybook clean-room; it never loads old components):

```ts
import { defineMain } from "@storybook/react-vite/node";

// No `core.builder` — the Vite builder is implied by the react-vite framework.
// builder-vite auto-merges this project's vite.config.ts, so @tailwindcss/vite
// and @vitejs/plugin-react apply in stories with no extra wiring.
export default defineMain({
  framework: "@storybook/react-vite",
  stories: ["../src/print/**/*.stories.@(ts|tsx)"],
});
```

- [ ] **Step 4: Configure `.storybook/preview.ts`**

Create `.storybook/preview.ts` — import the Print stylesheet (the existing `@tailwindcss/vite` plugin is reused via builder-vite, so `@theme` tokens compile):

```ts
import { definePreview } from "@storybook/react-vite";
import "../src/styles.css";

export default definePreview({
  parameters: {
    layout: "centered",
  },
});
```

- [ ] **Step 5: Write the exemplar story**

Create `src/print/ui/Button.stories.tsx` (CSF3):

```tsx
import type { Meta, StoryObj } from "@storybook/react-vite";
import { Button } from "./Button";

const meta = {
  title: "ui/Button",
  component: Button,
  args: { children: "Curate" },
} satisfies Meta<typeof Button>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Ghost: Story = { args: { variant: "ghost" } };
export const Solid: Story = { args: { variant: "solid" } };
export const Accent: Story = { args: { variant: "accent" } };
export const Danger: Story = { args: { variant: "danger", children: "Reject" } };
```

- [ ] **Step 6: Verify Storybook builds headlessly**

Run: `npm run build-storybook`
Expected: PASS — emits `storybook-static/`. Confirm: `test -d storybook-static && echo OK`.

(Optional interactive check during execution: `npm run storybook` and open http://localhost:6006 — the four Button variants render on the warm paper field with terracotta accent, confirming Print tokens are live.)

- [ ] **Step 7: Ignore the Storybook build output**

Append to `.gitignore` (frontend package root — verify the path; the repo root `.gitignore` already ignores `node_modules`):

```bash
printf '\n# Storybook\nstorybook-static/\n' >> .gitignore
```

- [ ] **Step 8: Commit**

```bash
git add package.json package-lock.json .storybook .gitignore src/print/ui/Button.stories.tsx
git commit -m "feat(print): Storybook component lab over src/print + Button exemplar story"
```

---

## Self-Review

**1. Spec coverage** (against `2026-05-31-greenfield-ui-rebuild-design.md`):
- Clean-room `src/print/` tree → Task 5. ✓
- Import-boundary guard (`no-legacy-import`) → Task 3. ✓
- `print/ui` as design-authoring → Task 3. ✓
- Separate Vite entry → Task 4. ✓
- Old excluded from build *output* → Task 4 (typecheck-exclusion deferred: unnecessary now + defeated by `lib/figure-export`→components coupling; see Task 4 Step 3). ✓
- Primitives seeded from clean extraction → Task 6. ✓
- Storybook for component fidelity → Task 7. ✓
- Phase 0 inventory (core-purity + catalog + surface map) → Tasks 1–2. ✓
- DEFERRED (correctly out of this plan): primitive hardening/expansion, renderer rebuilds, component/page rebuilds, Playwright page-fidelity, cutover. These are follow-up plans gated on Tasks 1–2's deliverables.

**2. Placeholder scan:** Tasks 1–2 are research deliverables with fixed document templates and the exact gathering commands; the catalog/surface-map *content* is produced during execution by reading the mockups (the only legitimately data-dependent output). All code tasks (3–7) contain complete code. No "TBD"/"add error handling"/"similar to" placeholders.

**3. Type/name consistency:** `PrintApp` (export) used identically in `src/print/App.tsx`, `src/print/main.tsx`, and `test/print-shell.test.tsx`. `scanLegacyImports(relPath, content)` signature matches between `check-design.mjs` and its tests; violation shape `{ rule, file, line, text }` matches `scanContent`'s existing shape. Root id `print-root` consistent between `print.html` and `src/print/main.tsx`. Storybook imports use `@storybook/react-vite` / `@storybook/react-vite/node` consistently.

**4. Ordering:** Task 4 intentionally defers its build verification + commit to Task 5 (the entry files it points at are created there) — called out explicitly in Task 4 Step 4.
