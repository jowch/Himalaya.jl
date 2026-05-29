# Component-Library Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract HimalayaUI's reusable UI patterns into a closed-look/open-placement primitive library, enforced by a design-token guard, so the unenforced design system stops drifting.

**Architecture:** Primitives in `src/components/ui/` own all appearance via semantic props (`variant`/`size`/`tone`/domain props); consumer `className` is placement-only and a `scripts/check-design.mjs` guard (content-hash baseline, per-rule ratchet) bans appearance utilities outside `ui/`. Four sequential phases: foundation (tokens + guard + catalog + adopt-or-delete + side-stripe fix) → core data primitives → chrome primitives → sweep & enforce.

**Tech Stack:** React 18 + TypeScript, Tailwind v4 (`@theme`), Vitest + Testing Library, zero new runtime deps (hand-rolled variant maps, no cva/clsx).

**Source spec:** `docs/superpowers/specs/2026-05-29-component-library-extraction-design.md` · **Migration checklist:** `docs/2026-05-29-library-extraction-inventory.json` (the exhaustive per-primitive call-site lists this plan's migration tasks work through).

---


## Phase 0 — Foundation (tokens, guard, catalog, adopt-or-delete, side-stripe fix)

Phase 0 is the hard prerequisite for every later wave. All commands run from `packages/HimalayaUI/frontend` unless stated. The guard script + baseline live at `packages/HimalayaUI/frontend/scripts/`.

### Task 1: Token foundation — radius scale, scrim, print-accent single-source

Phase 0 token edits in `src/styles.css` `@theme` (spec §3, §8). This is CSS/config — there is no React
component to render, so the test is **not RTL**: it parses the `@theme` source text (mirroring
`test/phases.test.ts`, which reads from source) AND compiles the real stylesheet with Tailwind v4's
programmatic `compile()` to assert the *resolved utility output* (`rounded-md` computes 5px), per the spec
§3 note that property-presence alone is insufficient (a stray `@layer`/specificity issue could shadow the
utility while the property exists).

**Verified facts (confirmed against live source + a Tailwind v4 compile probe — reproduce exactly):**
- `src/styles.css` `@theme` opens at line **9**. Current anchors: `--color-accent: oklch(0.555 0.150 38);`
  at line **20**; a DUPLICATE terracotta literal `--color-print-accent: oklch(0.555 0.150 38);` at line
  **48**; `--radius: 6px;` at line **65** (DEAD — no `var(--radius)` consumer); `.card` rule hardcodes
  `border-radius: 12px;` at line **175**. (Verify each line before editing — the file may have drifted.)
- `tailwindcss` (v4.2.4) exports a `compile(css, opts)` function. With the namespaced radius tokens in
  `@theme`, `rounded-md` compiles to `.rounded-md { border-radius: var(--radius-md); }` AND the compiled
  output emits the `--radius-md: 5px` declaration (Tailwind only emits *used* theme vars, so it appears
  precisely because `rounded-md` references it) — the full var chain resolves to `5px`. Confirmed by probe.
- `compile()` REQUIRES a `loadStylesheet` callback to resolve `@import "tailwindcss"`; resolve it via
  `createRequire(import.meta.url).resolve("tailwindcss/index.css")` and return its file contents. The
  `@import "@fontsource/*"` lines are irrelevant to token resolution — return `{ content: "" }` for them.
- `rounded-full` in Tailwind v4 is a STOCK utility resolving to `calc(infinity * 1px)`, NOT
  `var(--radius-full)`. Defining `--radius-full: 9999px` is still correct per spec §3 (consumers that
  reference the token directly, e.g. `border-radius: var(--radius-full)`, get 9999px) — but do NOT write a
  test asserting `rounded-full` → 9999px; it will fail. Assert only `--radius-full`'s presence/value in the
  `@theme` source.
- Tests may import from `node:fs` — `test/contentHash.test.ts` already does (`import { readFileSync } from
  "node:fs"`). Vitest env is jsdom with `globals: true` (`vitest.config.ts`).

**Files:**
- Modify: `src/styles.css` (lines ~20, ~48, ~65, ~175 — verify before editing)
- Test: `test/theme.test.ts` (new)

- [ ] **Step 1: Write the failing test** — create `test/theme.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import { resolve } from "node:path";
import { compile } from "tailwindcss";

const require = createRequire(import.meta.url);
const STYLES_PATH = resolve(__dirname, "../src/styles.css");
const SRC = readFileSync(STYLES_PATH, "utf-8");

// Extract the contents of the top-level `@theme { ... }` block from source.
// (Mirrors test/phases.test.ts: assert against source text, not a runtime DOM.)
function themeBlock(src: string): string {
  const open = src.indexOf("@theme");
  if (open < 0) throw new Error("no @theme block");
  const brace = src.indexOf("{", open);
  let depth = 0;
  for (let i = brace; i < src.length; i++) {
    if (src[i] === "{") depth++;
    else if (src[i] === "}") {
      depth--;
      if (depth === 0) return src.slice(brace + 1, i);
    }
  }
  throw new Error("unterminated @theme block");
}

// Read a `--name: value;` declaration's value from a CSS text region.
function declValue(css: string, name: string): string | null {
  const m = new RegExp(`${name.replace(/[-]/g, "\\-")}\\s*:\\s*([^;]+);`).exec(css);
  return m ? m[1]!.trim() : null;
}

// Compile the real stylesheet and return the resolved value of the var chain a
// utility's declared property points at (e.g. `border-radius: var(--radius-md)`
// resolved against the emitted `--radius-md: 5px`).
async function resolvedUtilityValue(utility: string, property: string): Promise<string> {
  const compiler = await compile(SRC, {
    base: resolve(__dirname, ".."),
    loadStylesheet: async (id: string, base: string) => {
      if (id === "tailwindcss") {
        const p = require.resolve("tailwindcss/index.css");
        return { base, content: readFileSync(p, "utf-8"), path: p };
      }
      // fontsource + any other import: irrelevant to token resolution
      return { base, content: "", path: id };
    },
  });
  const out = compiler.build([utility]);
  const ruleRe = new RegExp(`\\.${utility.replace(/[\\^$.*+?()[\]{}|]/g, "\\$&")}\\s*\\{([^}]*)\\}`);
  const rule = ruleRe.exec(out);
  if (!rule) throw new Error(`utility .${utility} not generated`);
  const propRe = new RegExp(`${property}\\s*:\\s*([^;]+);`);
  const prop = propRe.exec(rule[1]!);
  if (!prop) throw new Error(`${property} not set by .${utility}`);
  let value = prop[1]!.trim();
  // Follow a single var(--x) hop to its emitted :root declaration.
  const varRef = /^var\((--[\w-]+)\)$/.exec(value);
  if (varRef) {
    const resolved = declValue(out, varRef[1]!);
    if (!resolved) throw new Error(`${varRef[1]} not emitted in compiled output`);
    value = resolved;
  }
  return value;
}

describe("styles.css @theme tokens (Phase 0 token foundation)", () => {
  const theme = themeBlock(SRC);

  it("defines --color-scrim", () => {
    expect(declValue(theme, "--color-scrim")).toBe("oklch(0.05 0 0 / 0.65)");
  });

  it("sources --color-print-accent from --color-accent (single source, no duplicated literal)", () => {
    expect(declValue(theme, "--color-print-accent")).toBe("var(--color-accent)");
  });

  it("drops the dead singular --radius token", () => {
    // No `--radius:` declaration that is not part of `--radius-sm/md/full`.
    expect(/--radius\s*:/.test(theme)).toBe(false);
  });

  it("defines the namespaced radius scale", () => {
    expect(declValue(theme, "--radius-sm")).toBe("5px");
    expect(declValue(theme, "--radius-md")).toBe("5px");
    expect(declValue(theme, "--radius-full")).toBe("9999px");
  });

  it("resolves rounded-md utility to 5px (not stock Tailwind 6px)", async () => {
    expect(await resolvedUtilityValue("rounded-md", "border-radius")).toBe("5px");
  });

  it("resolves rounded-sm utility to 5px (not stock Tailwind 4px)", async () => {
    expect(await resolvedUtilityValue("rounded-sm", "border-radius")).toBe("5px");
  });

  it("no longer hardcodes a 12px corner on .card", () => {
    const cardRule = /\.card\s*\{([^}]*)\}/.exec(SRC);
    expect(cardRule).not.toBeNull();
    expect(/border-radius\s*:\s*12px/.test(cardRule![1]!)).toBe(false);
  });
});
```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/theme.test.ts` → Expected: FAIL (`--color-scrim` absent; `--color-print-accent` equals the literal `oklch(0.555 0.150 38)` not `var(--color-accent)`; singular `--radius:` still present; `--radius-sm/md/full` undefined; `rounded-md` resolves to stock 6px; `.card` still has `border-radius: 12px`).

- [ ] **Step 3: Implement (a) — replace the dead singular `--radius` with the namespaced scale** — in `src/styles.css`, replace the line `  --radius: 6px;` (≈ line 65, verify) with:

```css
  /* Radius scale — make the radius tokens REAL. Tailwind v4 generates
     rounded-* from the --radius-* NAMESPACE, not singular --radius (which was
     dead: no var(--radius) consumer, and rounded-sm/md resolved to Tailwind
     STOCK 4px/6px). Normalize the primary surfaces to the mockups' 5px corner
     (see spec §3/§8). sm and md deliberately both = 5px so the utility-name
     semantics (controls rounded-sm, cards rounded-md) survive for a future
     re-differentiation by editing one token. */
  --radius-sm:   5px;     /* buttons, inputs, chips, segmented controls */
  --radius-md:   5px;     /* cards, plates, sheets, modals */
  --radius-full: 9999px;  /* pills, dots */
```

- [ ] **Step 4: Implement (b) — single-source `--color-print-accent`** — in `src/styles.css`, replace the duplicated literal line `  --color-print-accent: oklch(0.555 0.150 38);` (≈ line 48, verify) with:

```css
  --color-print-accent: var(--color-accent); /* single source — was a duplicated terracotta literal */
```

- [ ] **Step 5: Implement (c) — add `--color-scrim`** — in `src/styles.css`, add directly after the `--color-print-accent` line:

```css
  --color-scrim: oklch(0.05 0 0 / 0.65); /* modal backdrops → bg-scrim (was hand-inlined literals) */
```

- [ ] **Step 6: Implement (d) — fold `.card` 12px corner onto the radius token** — in `src/styles.css` `.card` rule (≈ line 175, verify), replace `    border-radius: 12px;` with:

```css
    border-radius: var(--radius-md); /* 5px — folds onto the radius scale (spec §3/§8 sign-off) */
```

- [ ] **Step 7: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/theme.test.ts` → Expected: PASS (all 7 cases green).

- [ ] **Step 8: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (the `@theme` edits are CSS-only; `vite build` consumes the new tokens; no `tsc` surface changed). NOTE: if the guard task (`lint:design`) has already been wired into `build` by an earlier-merged Phase 0 task, this step also exercises the guard — it must still pass because this task introduces no new `className` appearance.

- [ ] **Step 9: Commit** — `cd packages/HimalayaUI/frontend && git add src/styles.css test/theme.test.ts && git commit -m "feat(tokens): namespaced radius scale (5px), --color-scrim, single-source --color-print-accent

Phase 0 token foundation (spec §3/§8). Replaces dead singular --radius with
--radius-sm/md (both 5px) + --radius-full (9999px); folds .card 12px corner onto
--radius-md; adds --color-scrim; points --color-print-accent at --color-accent.
test/theme.test.ts pins the @theme source + the resolved rounded-md utility (5px).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"`

> Phase 0. Lands the `scripts/check-design.mjs` design-token guard, its seed `design-baseline.json`,
> a fixture unit test, the `lint:design` build wiring, and a WARN-ONLY settings.json hook.
> Depends on the token foundation unit only for the *eventual* `bg-scrim`/radius-scale migrations
> that drain the baseline (Phase 3) — the guard itself is independent and ships its baseline so that
> nothing existing blocks and nothing new lands. All commands run from
> `packages/HimalayaUI/frontend` unless stated otherwise.

The guard is a **deliberately approximate line/token detector** (spec §4). It is authored as an
ES-module with **pure, exported functions** (`scanContent`, `hashViolation`, `loadBaseline`,
`diffBaseline`) plus a CLI entry guarded by `import.meta.url === ...`, so the Vitest fixture test
imports the functions directly (the repo's `package.json` is `"type": "module"`, Node 24).

### Task 2: Author the design guard script (scripts/check-design.mjs)

**Files:**
- Create: `packages/HimalayaUI/frontend/scripts/check-design.mjs`

The script is written first (no test yet) because the next task's fixture test imports its functions;
this task ends with a smoke-run that proves the file parses and runs. Write it IN FULL:

- [ ] **Step 1: Create the script** — write `packages/HimalayaUI/frontend/scripts/check-design.mjs` with exactly this content:

```js
#!/usr/bin/env node
// Zero-dependency design-token guard. Globs src/**/*.{ts,tsx}, excludes src/components/ui/**,
// flags banned appearance utilities / raw color literals against a content-hash baseline.
// Pure functions are exported for unit testing; the CLI runs only when invoked directly.
import { readdirSync, readFileSync, writeFileSync, existsSync } from "node:fs";
import { join, relative, sep } from "node:path";
import { fileURLToPath } from "node:url";
import { createHash } from "node:crypto";

const HERE = fileURLToPath(new URL(".", import.meta.url));
const SRC_DIR = join(HERE, "..", "src");
const BASELINE_PATH = join(HERE, "design-baseline.json");

// --- Ban rules ---------------------------------------------------------------
// Each rule: { id, test(line) -> matched-substring|null }. Rule #3 also consults `relPath`
// (POSIX-normalized, relative to src/) against its allowlist.
//
// Rule #3 allowlist: arbitrary-value color UTILITIES that legitimately need a runtime color.
// (spec §4: rules #1/#2/#4/#5 take NO allowlist.)
const RULE3_ALLOWLIST = new Set([
  "phases.ts",
  "lib/comparison/coloring.ts",
  "components/MemberHeatmapLayer.tsx",
  "components/DetectorImage.tsx",
  "components/FocusDetectorPanel.tsx",
  "main.tsx",
]);
// figure-export/** is allowlisted by prefix (the whole export palette dir).
const RULE3_ALLOW_PREFIXES = ["lib/figure-export/"];

function rule3Allowed(relPath) {
  if (RULE3_ALLOWLIST.has(relPath)) return true;
  return RULE3_ALLOW_PREFIXES.some((p) => relPath.startsWith(p));
}

const COLOR_LITERAL = "(?:oklch|oklab|hsla?|rgba?)\\(|#[0-9a-fA-F]{3,8}\\b";

const RULES = [
  {
    // #1 inline arbitrary type size or color: text-[10px], text-[var(--color-…)]
    id: "no-arbitrary-text",
    test: (line) => {
      const m = /\btext-\[/.exec(line);
      return m ? line.slice(m.index, m.index + 24) : null;
    },
  },
  {
    // #2 inline arbitrary radius: rounded-[3px], rounded-[1.5px], rounded-[7px], …
    id: "no-arbitrary-radius",
    test: (line) => {
      const m = /\brounded-\[/.exec(line);
      return m ? line.slice(m.index, m.index + 20) : null;
    },
  },
  {
    // #3 raw appearance color inside a color UTILITY (not shadow-[…]).
    // shadow-[…] is excluded — 4 legit Plate-Lift rgba shadows live in tsx.
    id: "no-raw-color-utility",
    allowlisted: rule3Allowed,
    test: (line) => {
      const re = new RegExp(
        "\\b(?:bg|text|border|fill|stroke|ring|outline|from|via|to|decoration)-\\[(?:[^\\]]*?)(?:" +
          COLOR_LITERAL +
          ")",
      );
      const m = re.exec(line);
      return m ? line.slice(m.index, m.index + 40) : null;
    },
  },
  {
    // #4 side-stripe: border-l-4 / border-r-4 / border-{l,r}-[…] / border-{l,r}-2..9
    id: "no-side-stripe",
    test: (line) => {
      const m = /\bborder-(?:l|r)-(?:4|[2-9]|\[)/.exec(line);
      return m ? line.slice(m.index, m.index + 16) : null;
    },
  },
  {
    // #5 raw color literal ANYWHERE on the line (covers multi-line style={{…}}).
    // Strip every var(--color-*) token FIRST, then flag a remaining raw literal.
    // This passes legit computed-color sites (their non-var color is a ${…} expression
    // with no raw literal) and flags the raw oklch(0.05…) scrim hand-inlines.
    id: "no-raw-color-literal",
    test: (line) => {
      const stripped = line.replace(/var\(--color-[^)]*\)/g, "");
      const re = new RegExp("(?:" + COLOR_LITERAL + ")");
      const m = re.exec(stripped);
      return m ? stripped.slice(m.index, m.index + 32).trim() : null;
    },
  },
];

// --- Violation discovery -----------------------------------------------------
function listSourceFiles(dir) {
  const out = [];
  for (const ent of readdirSync(dir, { withFileTypes: true })) {
    const full = join(dir, ent.name);
    if (ent.isDirectory()) {
      out.push(...listSourceFiles(full));
    } else if (/\.(ts|tsx)$/.test(ent.name)) {
      out.push(full);
    }
  }
  return out;
}

// POSIX-normalized path relative to src/ (so allowlist + exclusion are OS-independent).
function relToSrc(absPath) {
  return relative(SRC_DIR, absPath).split(sep).join("/");
}

// src/components/ui/** is excluded entirely (it is where appearance is authored).
function isExcluded(relPath) {
  return relPath.startsWith("components/ui/");
}

// Scan one file's text. Returns [{ rule, file, line, text }].
// Exported so the fixture test can drive it without touching disk.
export function scanContent(relPath, content) {
  if (isExcluded(relPath)) return [];
  const violations = [];
  const lines = content.split("\n");
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    for (const rule of RULES) {
      if (rule.allowlisted && rule.allowlisted(relPath)) continue;
      const matched = rule.test(line);
      if (matched != null) {
        violations.push({ rule: rule.id, file: relPath, line: i + 1, text: matched });
      }
    }
  }
  return violations;
}

// Content-hash key: (rule, normalized-violation-text) — NOT (file,line), NOT (file,rule)-count.
// Normalize: collapse all whitespace runs to one space + trim, so reflow/indent moves are no-ops.
export function hashViolation(rule, text) {
  const norm = text.replace(/\s+/g, " ").trim();
  return createHash("sha1").update(rule + " " + norm).digest("hex").slice(0, 16);
}

export function loadBaseline(path = BASELINE_PATH) {
  if (!existsSync(path)) return { hashes: [] };
  return JSON.parse(readFileSync(path, "utf8"));
}

// baseline_new ⊆ baseline_old check + not-in-baseline detection.
// Returns { notInBaseline: [...violations], grew: bool }.
export function diffBaseline(violations, baseline) {
  const baseSet = new Set(baseline.hashes ?? []);
  const seen = new Set();
  const notInBaseline = [];
  for (const v of violations) {
    const h = hashViolation(v.rule, v.text);
    seen.add(h);
    if (!baseSet.has(h)) notInBaseline.push({ ...v, hash: h });
  }
  // grew = the current set introduced a hash the baseline does not contain.
  return { notInBaseline, grew: notInBaseline.length > 0, seenHashes: seen };
}

// --- CLI ---------------------------------------------------------------------
function runCli(argv) {
  const initMode = argv.includes("--init");
  const files = listSourceFiles(SRC_DIR);
  const all = [];
  for (const abs of files) {
    const rel = relToSrc(abs);
    all.push(...scanContent(rel, readFileSync(abs, "utf8")));
  }

  if (initMode) {
    const hashes = [...new Set(all.map((v) => hashViolation(v.rule, v.text)))].sort();
    const baseline = {
      note: "Content-hash baseline for scripts/check-design.mjs. Keyed on (rule, normalized text). Ratchet: only ever shrinks. Regenerate seed with `node scripts/check-design.mjs --init`.",
      generated: new Date().toISOString().slice(0, 10),
      hashes,
    };
    writeFileSync(BASELINE_PATH, JSON.stringify(baseline, null, 2) + "\n");
    process.stderr.write(`[check-design] wrote ${hashes.length} baseline hashes to ${BASELINE_PATH}\n`);
    return 0;
  }

  const baseline = loadBaseline();
  const { notInBaseline } = diffBaseline(all, baseline);
  if (notInBaseline.length > 0) {
    process.stderr.write(`[check-design] ${notInBaseline.length} NEW design violation(s) not in baseline:\n`);
    for (const v of notInBaseline) {
      process.stderr.write(`  ${v.rule}  src/${v.file}:${v.line}  ${JSON.stringify(v.text)}  [${v.hash}]\n`);
    }
    process.stderr.write("Fix the appearance utility (move it into src/components/ui/**) or, if intentional and pre-existing, the baseline is wrong.\n");
    return 2;
  }
  return 0;
}

if (import.meta.url === `file://${process.argv[1]}`) {
  process.exit(runCli(process.argv.slice(2)));
}
```

- [ ] **Step 2: Smoke-run it parses and runs** — `cd packages/HimalayaUI/frontend && node scripts/check-design.mjs; echo "exit=$?"` → Expected: it runs without a syntax error; with no baseline file yet it prints a list of NEW violations to stderr and `exit=2` (that is correct: every current violation is "not in baseline" until the seed is generated in the next task). The point of this step is only that the module loads and the regexes compile.
- [ ] **Step 3: Commit** — `git add packages/HimalayaUI/frontend/scripts/check-design.mjs && git commit -m "feat(guard): zero-dep design-token guard (scan + content-hash baseline diff)"`

### Task 3: Guard fixture unit test (spec §4 cases)

**Files:**
- Test: `packages/HimalayaUI/frontend/test/check-design.test.ts`

The test imports the guard's pure functions and drives each spec §4 fixture case through
`scanContent` / `diffBaseline`. It asserts on `rule` ids and counts — never on regex internals.

- [ ] **Step 1: Write the failing test** — create `packages/HimalayaUI/frontend/test/check-design.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { scanContent, hashViolation, diffBaseline } from "../scripts/check-design.mjs";

// Helper: collect the set of rule ids a single line trips (outside ui/).
function rulesFor(line: string): string[] {
  return scanContent("components/Probe.tsx", line).map((v) => v.rule);
}

describe("check-design guard — ban rules (spec §4)", () => {
  // Rule #1 — arbitrary text size/color
  it("flags text-[10px]", () => {
    expect(rulesFor('<span className="text-[10px] mt-1" />')).toContain("no-arbitrary-text");
  });
  it("passes the on-scale text-base", () => {
    expect(rulesFor('<span className="text-base font-medium" />')).not.toContain("no-arbitrary-text");
  });

  // Rule #2 — arbitrary radius
  it("flags rounded-[3px]", () => {
    expect(rulesFor('<div className="rounded-[3px]" />')).toContain("no-arbitrary-radius");
  });
  it("passes rounded-md", () => {
    expect(rulesFor('<div className="rounded-md bg-plate" />')).not.toContain("no-arbitrary-radius");
  });

  // Rule #4 — side stripe
  it("flags border-l-4", () => {
    expect(rulesFor('<div className="border-l-4 border-error" />')).toContain("no-side-stripe");
  });

  // Rule #3 — raw color in a color utility, shadow-[…] excluded
  it("flags bg-[oklch(...)] in a non-allowlisted file", () => {
    expect(rulesFor('<div className="bg-[oklch(0.05_0_0/0.65)]" />')).toContain("no-raw-color-utility");
  });
  it("passes a shadow-[…rgba…] Plate-Lift literal (rule #3 excludes shadow-)", () => {
    expect(rulesFor('<div className="shadow-[0_8px_26px_-10px_rgba(60,52,40,.34)]" />')).not.toContain("no-raw-color-utility");
  });
  it("allowlists bg-[oklch] inside an allowlisted file (rule #3 only)", () => {
    const v = scanContent("components/DetectorImage.tsx", '<div className="bg-[oklch(0.15_0.01_55)]" />');
    expect(v.map((x) => x.rule)).not.toContain("no-raw-color-utility");
  });

  // Rule #5 — raw color literal anywhere; var(--color-*) stripped first
  it("passes a var-only color-mix string literal", () => {
    const line = '  background: "color-mix(in oklab, var(--color-accent) 9%, transparent)",';
    expect(rulesFor(line)).not.toContain("no-raw-color-literal");
  });
  it("passes a computed style={{ color }} identifier", () => {
    expect(rulesFor("        style={{ color }}")).not.toContain("no-raw-color-literal");
  });
  it("passes a template-literal interpolated color", () => {
    expect(rulesFor("        style={{ background: `color-mix(in oklab, ${x} 20%, transparent)` }}")).not.toContain("no-raw-color-literal");
  });
  it("flags a raw oklch scrim string literal", () => {
    expect(rulesFor('        background: "oklch(0.05 0 0 / 0.65)",')).toContain("no-raw-color-literal");
  });
  it("flags a bare var() with NO raw literal as clean (FocusReflectionsTable/SeriesFolioCard cases)", () => {
    expect(rulesFor('        style={{ background: "var(--color-success)" }}')).not.toContain("no-raw-color-literal");
  });

  // Scope exclusion — src/components/ui/** is never scanned
  it("excludes src/components/ui/** entirely", () => {
    expect(scanContent("components/ui/Toast.tsx", '<div className="border-l-4 text-[10px]" />')).toHaveLength(0);
  });
});

describe("check-design guard — baseline diff", () => {
  it("returns clean (no new) when every violation hash is in the baseline", () => {
    const violations = scanContent("components/Probe.tsx", '<div className="text-[10px]" />');
    const baseline = { hashes: violations.map((v) => hashViolation(v.rule, v.text)) };
    const { notInBaseline, grew } = diffBaseline(violations, baseline);
    expect(notInBaseline).toHaveLength(0);
    expect(grew).toBe(false);
  });

  it("flags a not-in-baseline violation (CI would exit 2)", () => {
    const violations = scanContent("components/Probe.tsx", '<div className="text-[10px]" />');
    const { notInBaseline, grew } = diffBaseline(violations, { hashes: [] });
    expect(notInBaseline).toHaveLength(1);
    expect(grew).toBe(true);
  });

  it("hash is stable across whitespace/indent reflow (move-between-files is a no-op)", () => {
    expect(hashViolation("no-arbitrary-text", "text-[10px]")).toBe(
      hashViolation("no-arbitrary-text", "   text-[10px]   "),
    );
  });
});
```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/check-design.test.ts` → Expected: FAIL only if the script does not yet export the functions; since Unit task 1 already wrote and committed the script, this should actually PASS on first run. If it fails, the failure is a guard-logic bug — fix `scripts/check-design.mjs` (not the test) until the spec §4 cases pass. (TDD note: the script and test are split into two commits per the brief; the test is the executable spec for the §4 case table.)
- [ ] **Step 3: Run, expect pass** — same command → Expected: PASS (all cases green).
- [ ] **Step 4: Commit** — `git add packages/HimalayaUI/frontend/test/check-design.test.ts && git commit -m "test(guard): fixture cases for the 5 ban rules + baseline diff (spec §4)"`

### Task 4: Generate and commit the seed design-baseline.json

**Files:**
- Create: `packages/HimalayaUI/frontend/scripts/design-baseline.json` (machine-generated by `--init`)

The seed is generated from **today's actual violations** by running the guard in `--init` mode (not
hand-written from the inventory — the script is the source of truth so the hash function and the
baseline can never drift apart). The inventory's violation lists (the ~96 `text-[`, the 9 `rounded-[`,
the raw-color/scrim sites, the one enforced side-stripe at `InfrastructureBanner.tsx:55`) are the
*expectation* the count is sanity-checked against.

- [ ] **Step 1: Generate the seed** — `cd packages/HimalayaUI/frontend && node scripts/check-design.mjs --init` → Expected: writes `scripts/design-baseline.json` and prints `[check-design] wrote N baseline hashes …` to stderr. N should be on the order of ~100–130 unique hashes (≈96 `text-[` minus duplicates that collapse by content, +9 `rounded-[`, + the `rounded-lg` sites, + the handful of raw scrim/color sites, + the single `InfrastructureBanner.tsx` `border-l-4`; `Toast.tsx` is in `ui/` and excluded).
- [ ] **Step 2: Sanity-check it is clean against itself** — `node scripts/check-design.mjs; echo "exit=$?"` → Expected: `exit=0` and **no** stderr violation list (every current violation now hashes into the baseline — this proves "nothing existing blocks").
- [ ] **Step 3: Spot-check the expected hits are present** — `node scripts/check-design.mjs --init 2>&1 >/dev/null | head` then confirm the count printed matches step 1; and confirm `InfrastructureBanner.tsx:55` (the one enforced `border-l-4`) is captured and `Toast.tsx:85` is NOT (it is under `ui/`): `node -e "import('./scripts/check-design.mjs').then(async m=>{const {readFileSync}=await import('node:fs');const inf=m.scanContent('components/InfrastructureBanner.tsx', readFileSync('src/components/InfrastructureBanner.tsx','utf8'));console.log('infra side-stripe hits:', inf.filter(v=>v.rule==='no-side-stripe').length);const t=m.scanContent('components/ui/Toast.tsx','x');console.log('ui Toast scanned:', t.length);})"` → Expected: `infra side-stripe hits: 1`, `ui Toast scanned: 0`.
- [ ] **Step 4: Commit** — `git add packages/HimalayaUI/frontend/scripts/design-baseline.json && git commit -m "chore(guard): seed design-baseline.json from current violations (Phase 0 ratchet floor)"`

### Task 5: Wire lint:design into package.json build

**Files:**
- Modify: `packages/HimalayaUI/frontend/package.json` (`scripts` block)

- [ ] **Step 1: Add the script + prepend to build** — in `packages/HimalayaUI/frontend/package.json`, change the `scripts` block so that:
  - `"build"` becomes `"npm run lint:design && tsc --noEmit -p tsconfig.build.json && vite build"` (was `"tsc --noEmit -p tsconfig.build.json && vite build"`).
  - a new key `"lint:design": "node scripts/check-design.mjs"` is added.

  Resulting `scripts` block:
```json
  "scripts": {
    "dev": "vite",
    "build": "npm run lint:design && tsc --noEmit -p tsconfig.build.json && vite build",
    "lint:design": "node scripts/check-design.mjs",
    "preview": "vite preview",
    "test": "vitest run",
    "test:watch": "vitest",
    "e2e": "playwright test",
    "e2e:live": "playwright test -c playwright.live.config.ts"
  }
```

- [ ] **Step 2: Verify the wiring runs and blocks correctly** — `cd packages/HimalayaUI/frontend && npm run lint:design; echo "exit=$?"` → Expected: `exit=0` (baseline is current, nothing new). Then prove it blocks: `printf 'export const X = () => <div className="text-[99px]" />;\n' > src/components/__guardprobe.tsx && (npm run lint:design; echo "exit=$?") ; rm src/components/__guardprobe.tsx` → Expected: prints the new `no-arbitrary-text` violation for `__guardprobe.tsx` and `exit=2`.
- [ ] **Step 3: Confirm full build still passes** — `npm run build` → Expected: `lint:design` passes (exit 0), then tsc + vite build complete. (If vite/tsc are slow, this is the gate the PR must clear; capture once.)
- [ ] **Step 4: Commit** — `git add packages/HimalayaUI/frontend/package.json && git commit -m "build(guard): run lint:design as the blocking first step of build"`

### Task 6: Add the WARN-ONLY guard hook to settings.json

**Files:**
- Modify: `.claude/settings.json` (root; the `hooks.PostToolUse` array)

The CI `lint:design` step is the hard gate (`exit 2`). This hook is **fast feedback only**: on a
`frontend/src/**/*.tsx` (or `.ts`) edit it runs the guard and prints any NEW (not-in-baseline)
violations to stderr, but **always `exit 0`** so it never blocks the author mid-refactor while they
are lowering the baseline (spec §4 "Wiring"). It mirrors the existing tsc/vitest PostToolUse hook
shape (matcher `Edit|Write|MultiEdit`, `$CLAUDE_PROJECT_DIR` rooting, `jq` to read `file_path`).

- [ ] **Step 1: Add the hook entry** — append a third object to the existing `hooks.PostToolUse` array in `.claude/settings.json` (after the Vitest hook object, before the array closes). The new entry:

```json
      {
        "matcher": "Edit|Write|MultiEdit",
        "hooks": [
          {
            "type": "command",
            "description": "WARN-ONLY: print new design-token guard violations on frontend src TS/TSX edit",
            "command": "f=$(jq -r '.tool_input.file_path // empty'); case \"$f\" in *packages/HimalayaUI/frontend/src/*.ts|*packages/HimalayaUI/frontend/src/*.tsx) cd \"$CLAUDE_PROJECT_DIR/packages/HimalayaUI/frontend\" && out=$(node scripts/check-design.mjs 2>&1); code=$?; if [ \"$code\" -ne 0 ]; then printf '\\n[design-guard WARN] new appearance utility in a consumer (CI build will block; move it into src/components/ui/** or fix):\\n%s\\n' \"$out\" >&2; fi ;; esac; exit 0"
          }
        ]
      }
```

  Note three things the engineer MUST preserve: (a) the trailing `; exit 0` so the hook is warn-only even when the guard returns 2; (b) `case` matches `frontend/src/**` only (NOT `test/**` — the guard does not scan tests); (c) the existing JSON array gets a comma after the prior (Vitest) object — verify valid JSON after editing.
- [ ] **Step 2: Validate the JSON** — `jq . .claude/settings.json > /dev/null && echo OK` (run from repo root) → Expected: `OK` (no parse error).
- [ ] **Step 3: Dry-run the hook command body in warn mode** — simulate a clean edit and a dirty edit:
  - clean: `cd packages/HimalayaUI/frontend && out=$(node scripts/check-design.mjs 2>&1); echo "code=$? len=${#out}"` → Expected: `code=0` (no warning would print).
  - dirty: `printf 'export const X = () => <div className="rounded-[42px]" />;\n' > src/components/__guardprobe.tsx && (node scripts/check-design.mjs 2>&1; echo "code=$?"); rm src/components/__guardprobe.tsx` → Expected: prints the `no-arbitrary-radius` line and `code=2` — confirming the hook would WARN (the `; exit 0` in the hook then swallows the code).
- [ ] **Step 4: Commit** — `git add .claude/settings.json && git commit -m "chore(hooks): warn-only design-guard PostToolUse hook on frontend src edits"`

### Task 7: Scaffold the static catalog page (`docs/design-system.html`)

The catalog is a single hand-authored static HTML page — a visual reference and adoption
surface for the Print primitives, sibling to the existing `docs/redesign-mockups/*.html`
mockups. It is **not load-bearing**: the enforcement guard (`scripts/check-design.mjs`) and
`DESIGN.md` remain the enforced source of truth (spec §5). The HTML may drift; that drift is
knowingly accepted. There are **no tests** for this unit — it is a static doc.

This task lays the scaffold (shell, inline `<style>` mirroring the `@theme` Print tokens, page
chrome, and the section frame) plus the **Phase-0 primitive sections only**: Button, ScoreBar,
Dot, and Toast (the side-stripe-fixed look). Later primitives (SegmentedControl, PhaseChip,
PhaseStrip, Kicker, ModalShell, IconButton, Card/Plate) are **appended in their own waves** — a
visible "Appended in later waves" placeholder list at the bottom of the page marks this.

Sections are hand-authored static HTML markup (not live React) — the catalog hand-mirrors the
primitives' rendered look using the same token vars, so it carries no import dependency on the
component source and cannot break the build.

**Files:**
- Create: `packages/HimalayaUI/frontend/../../../docs/design-system.html` (repo-root `docs/design-system.html`)

> Note on path: `docs/` is at the **repo root**, not under `frontend/`. The absolute target is
> `/Users/me/projects/Himalaya.jl/.claude/worktrees/impeccable-live/docs/design-system.html`.
> Use that absolute path with the Write tool.

- [ ] **Step 1: Confirm the sibling-mockup convention and token values before authoring** —
  Run `ls docs/redesign-mockups/` (expect: `focus-workspace.html`, `sample-table.html`,
  `series-builder.html`, `series-folio.html`, `series-scoping.html`) to confirm the static-HTML
  convention this page joins, and `grep -n -- "--color-paper:\|--color-ink:\|--color-hair:\|--color-accent:\|--color-success:\|--color-warning:\|--color-error:\|--color-plate:\|--color-ink-soft:\|--color-ink-faint:\|--color-hair-strong:" packages/HimalayaUI/frontend/src/styles.css`
  to confirm the live token values reproduced inline below. Verified values at authoring time
  (reproduce these exactly): `--color-paper: oklch(0.978 0.006 85)`, `--color-paper-sunk:
  oklch(0.951 0.008 84)`, `--color-plate: oklch(0.992 0.004 90)`, `--color-ink: oklch(0.265
  0.013 68)`, `--color-ink-soft: oklch(0.467 0.012 70)`, `--color-ink-faint: oklch(0.640 0.010
  74)`, `--color-hair: oklch(0.882 0.008 80)`, `--color-hair-strong: oklch(0.806 0.010 78)`,
  `--color-accent: oklch(0.555 0.150 38)`, `--color-print-accent: oklch(0.555 0.150 38)`,
  `--color-success: oklch(0.520 0.120 162)`, `--color-warning: oklch(0.620 0.130 70)`,
  `--color-error: oklch(0.520 0.170 28)`. Phase colors (from `src/phases.ts`): `Pn3m
  oklch(0.570 0.150 58)`, `Im3m oklch(0.520 0.120 162)`, `Ia3d oklch(0.520 0.130 300)`, `Square
  oklch(0.550 0.130 132)`, `Lamellar oklch(0.505 0.150 264)`. Radius tokens per the locked
  decision: `--radius-sm` and `--radius-md` BOTH `5px`, `--radius-full 9999px` (the singular
  `--radius` is deleted — do NOT reference it).

- [ ] **Step 2: Write the full scaffold + Phase-0 sections** — Write
  `/Users/me/projects/Himalaya.jl/.claude/worktrees/impeccable-live/docs/design-system.html`
  with this exact content:

```html
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>Himalaya — Print Component Library</title>
<!--
  Hand-authored visual catalog of the Print primitives. Sibling to
  docs/redesign-mockups/*.html. NOT load-bearing: scripts/check-design.mjs + DESIGN.md
  are the enforced source of truth. This page may drift; that drift is knowingly accepted
  (spec §5). The inline <style> below MIRRORS the @theme tokens in
  packages/HimalayaUI/frontend/src/styles.css — when a token changes there, update here.
  Sections are static markup hand-mirroring each primitive's rendered look (no React import),
  so this file cannot break the build.

  Phase 0 sections present: Button, ScoreBar, Dot, Toast.
  Later primitives are appended in their own migration waves (see the placeholder at the foot).
-->
<style>
  /* ── Print tokens (mirror of src/styles.css @theme) ─────────────────── */
  :root {
    --color-paper:        oklch(0.978 0.006 85);
    --color-paper-sunk:   oklch(0.951 0.008 84);
    --color-plate:        oklch(0.992 0.004 90);
    --color-ink:          oklch(0.265 0.013 68);
    --color-ink-soft:     oklch(0.467 0.012 70);
    --color-ink-faint:    oklch(0.640 0.010 74);
    --color-hair:         oklch(0.882 0.008 80);
    --color-hair-strong:  oklch(0.806 0.010 78);
    --color-accent:       oklch(0.555 0.150 38);   /* terracotta */
    --color-print-accent: oklch(0.555 0.150 38);
    --color-success:      oklch(0.520 0.120 162);  /* sage */
    --color-warning:      oklch(0.620 0.130 70);   /* amber */
    --color-error:        oklch(0.520 0.170 28);   /* red */

    /* Phase tints (mirror of src/phases.ts phaseColor) */
    --phase-Pn3m:     oklch(0.570 0.150 58);
    --phase-Im3m:     oklch(0.520 0.120 162);
    --phase-Ia3d:     oklch(0.520 0.130 300);
    --phase-Square:   oklch(0.550 0.130 132);
    --phase-Lamellar: oklch(0.505 0.150 264);

    /* Radius (locked: sm == md == 5px; full 9999px; singular --radius deleted) */
    --radius-sm: 5px;
    --radius-md: 5px;
    --radius-full: 9999px;

    --font-sans: "Plus Jakarta Sans", ui-sans-serif, system-ui, sans-serif;
    --font-serif: "Newsreader", Georgia, ui-serif, serif;
    --font-mono: ui-monospace, SFMono-Regular, Menlo, monospace;
  }

  * { box-sizing: border-box; }
  body {
    margin: 0;
    background: var(--color-paper);
    color: var(--color-ink);
    font-family: var(--font-sans);
    font-size: 13px;
    line-height: 1.5;
  }
  .page { max-width: 920px; margin: 0 auto; padding: 48px 32px 96px; }

  /* Page header */
  header.masthead { border-bottom: 1px solid var(--color-hair-strong); padding-bottom: 20px; margin-bottom: 8px; }
  header.masthead .kicker {
    font-weight: 700; letter-spacing: 0.13em; text-transform: uppercase;
    font-size: 11px; color: var(--color-accent); margin: 0 0 8px;
  }
  header.masthead h1 { font-family: var(--font-serif); font-size: 27px; font-weight: 500; margin: 0; }
  header.masthead p { color: var(--color-ink-soft); max-width: 60ch; margin: 10px 0 0; }
  header.masthead .note { color: var(--color-ink-faint); font-size: 11.5px; font-style: italic; margin-top: 6px; }

  /* Section frame */
  section.primitive { margin-top: 48px; }
  section.primitive > .sec-kicker {
    font-weight: 700; letter-spacing: 0.09em; text-transform: uppercase;
    font-size: 11px; color: var(--color-ink-faint); margin: 0 0 4px;
  }
  section.primitive > h2 { font-family: var(--font-serif); font-size: 19px; font-weight: 500; margin: 0 0 6px; }
  section.primitive > .blurb { color: var(--color-ink-soft); max-width: 64ch; margin: 0 0 4px; }
  section.primitive > .shows { color: var(--color-ink-faint); font-size: 11.5px; margin: 0 0 18px; }
  section.primitive > .shows b { color: var(--color-ink-soft); font-weight: 600; }

  .specimen {
    background: var(--color-plate);
    border: 0.5px solid var(--color-hair);
    border-radius: var(--radius-md);
    padding: 20px;
  }
  .row { display: flex; flex-wrap: wrap; align-items: center; gap: 16px; }
  .row + .row { margin-top: 18px; }
  .cell { display: flex; flex-direction: column; gap: 6px; align-items: flex-start; }
  .cell .cap { font-size: 10.5px; color: var(--color-ink-faint); letter-spacing: 0.02em; }

  /* ── Button mirror ──────────────────────────────────────────────────── */
  .btn {
    font-family: var(--font-sans); font-size: 11.5px; line-height: 1.4;
    padding: 4px 10px; border-radius: var(--radius-md); cursor: pointer;
    transition: filter 120ms, color 120ms, background 120ms;
    border: 1px solid transparent;
  }
  .btn-solid  { background: var(--color-ink);    border-color: var(--color-ink);    color: var(--color-paper); }
  .btn-accent { background: var(--color-accent); border-color: var(--color-accent); color: var(--color-paper); }
  .btn-ghost  { background: transparent; color: var(--color-ink-soft); }
  .btn-danger { background: transparent; color: var(--color-ink-soft); }
  .btn-solid:hover, .btn-accent:hover { filter: brightness(1.1); }
  .btn-ghost:hover  { color: var(--color-ink); background: var(--color-paper-sunk); }
  .btn-danger:hover { color: var(--color-error); }
  /* simulated states (static mirror — no real interaction) */
  .is-hover.btn-solid, .is-hover.btn-accent { filter: brightness(1.1); }
  .is-hover.btn-ghost  { color: var(--color-ink); background: var(--color-paper-sunk); }
  .is-hover.btn-danger { color: var(--color-error); }
  .is-focus { outline: 2px solid var(--color-accent); outline-offset: 2px; }
  .btn:disabled, .is-disabled { opacity: 0.45; cursor: not-allowed; filter: none; }

  /* ── ScoreBar mirror ────────────────────────────────────────────────── */
  .scorebar { overflow: hidden; border-radius: var(--radius-full); background: var(--color-hair); }
  .scorebar.size-bar     { height: 4px;  width: 200px; }
  .scorebar.size-compact { height: 3.5px; width: 46px; }
  .scorebar > i { display: block; height: 100%; border-radius: var(--radius-full); }

  /* ── Dot mirror ─────────────────────────────────────────────────────── */
  .dot { display: inline-block; border-radius: var(--radius-full); flex-shrink: 0; }
  .dot.size-xs { width: 6px; height: 6px; }
  .dot.size-sm { width: 8px; height: 8px; }
  .dot.tone-accent  { background: var(--color-print-accent); }
  .dot.tone-success { background: var(--color-success); }
  .dot.tone-muted   { background: transparent; box-shadow: inset 0 0 0 1px var(--color-hair-strong); }
  .dot.tone-neutral { background: var(--color-hair-strong); }

  /* ── Toast mirror (side-stripe fix: leading status icon + word + full hairline) ─ */
  .toast {
    display: flex; align-items: center; gap: 10px;
    background: var(--color-plate);
    border: 0.5px solid var(--color-hair-strong);
    border-radius: var(--radius-md);
    padding: 8px 12px; color: var(--color-ink);
    box-shadow: 0 8px 24px oklch(0.05 0 0 / 0.12);
    min-width: 280px;
  }
  .toast .status-icon {
    font-family: var(--font-mono); font-weight: 700; font-size: 12px;
    width: 18px; text-align: center; flex-shrink: 0;
  }
  .toast .status-word { font-weight: 700; letter-spacing: 0.04em; font-size: 11.5px; text-transform: uppercase; }
  .toast .msg { flex: 1; color: var(--color-ink-soft); }
  .toast .dismiss {
    border: none; background: transparent; cursor: pointer; color: var(--color-ink-faint);
    font-size: 14px; line-height: 1; padding: 2px 4px; border-radius: var(--radius-sm);
  }
  .toast .dismiss:hover { color: var(--color-ink); background: var(--color-paper-sunk); }
  .toast.kind-info    .status-icon { color: var(--color-ink-soft); }
  .toast.kind-success .status-icon { color: var(--color-success); }
  .toast.kind-warning .status-icon { color: var(--color-warning); }
  .toast.kind-error   .status-icon { color: var(--color-error); }

  /* Foot / deferred list */
  .foot { margin-top: 64px; border-top: 1px solid var(--color-hair-strong); padding-top: 24px; }
  .foot .sec-kicker { font-weight: 700; letter-spacing: 0.09em; text-transform: uppercase; font-size: 11px; color: var(--color-ink-faint); margin: 0 0 10px; }
  .foot ul { color: var(--color-ink-soft); margin: 0; padding-left: 18px; }
  .foot li { margin: 4px 0; }
  .foot li code { font-family: var(--font-mono); font-size: 11.5px; color: var(--color-ink); }
</style>
</head>
<body>
<div class="page">

  <header class="masthead">
    <p class="kicker">Himalaya Print System</p>
    <h1>Component Library</h1>
    <p>Hand-authored visual reference for the Print primitives, built on the Print design tokens.
       Each primitive renders its variants, sizes, and states on warm paper.</p>
    <p class="note">Not load-bearing — the enforcement guard (scripts/check-design.mjs) and DESIGN.md
       are the enforced source of truth. This catalog may drift; that is accepted.</p>
  </header>

  <!-- ════════════════ Phase 0 — Button ════════════════ -->
  <section class="primitive" id="button">
    <p class="sec-kicker">Phase 0 · Action</p>
    <h2>Button</h2>
    <p class="blurb">Variant axis <code>solid | accent | ghost | danger</code>. <b>solid</b> is the
       ink-filled confirm action with paper text (not terracotta); <b>accent</b> is the reserved
       terracotta grease-pencil action; <b>ghost</b> is default chrome.</p>
    <p class="shows">Catalog must show: <b>all 4 variants × {rest, hover, focus-visible ring, disabled}</b>;
       solid reads as ink (not terracotta), accent reads as terracotta, button text reads as paper.</p>
    <div class="specimen">
      <div class="row">
        <div class="cell"><button class="btn btn-solid">Save</button><span class="cap">solid · rest</span></div>
        <div class="cell"><button class="btn btn-solid is-hover">Save</button><span class="cap">solid · hover</span></div>
        <div class="cell"><button class="btn btn-solid is-focus">Save</button><span class="cap">solid · focus-visible</span></div>
        <div class="cell"><button class="btn btn-solid" disabled>Save</button><span class="cap">solid · disabled</span></div>
      </div>
      <div class="row">
        <div class="cell"><button class="btn btn-accent">+ New series</button><span class="cap">accent · rest</span></div>
        <div class="cell"><button class="btn btn-accent is-hover">+ New series</button><span class="cap">accent · hover</span></div>
        <div class="cell"><button class="btn btn-accent is-focus">+ New series</button><span class="cap">accent · focus-visible</span></div>
        <div class="cell"><button class="btn btn-accent" disabled>+ New series</button><span class="cap">accent · disabled</span></div>
      </div>
      <div class="row">
        <div class="cell"><button class="btn btn-ghost">Cancel</button><span class="cap">ghost · rest</span></div>
        <div class="cell"><button class="btn btn-ghost is-hover">Cancel</button><span class="cap">ghost · hover</span></div>
        <div class="cell"><button class="btn btn-ghost is-focus">Cancel</button><span class="cap">ghost · focus-visible</span></div>
        <div class="cell"><button class="btn btn-ghost" disabled>Cancel</button><span class="cap">ghost · disabled</span></div>
      </div>
      <div class="row">
        <div class="cell"><button class="btn btn-danger">Remove</button><span class="cap">danger · rest</span></div>
        <div class="cell"><button class="btn btn-danger is-hover">Remove</button><span class="cap">danger · hover</span></div>
        <div class="cell"><button class="btn btn-danger is-focus">Remove</button><span class="cap">danger · focus-visible</span></div>
        <div class="cell"><button class="btn btn-danger" disabled>Remove</button><span class="cap">danger · disabled</span></div>
      </div>
    </div>
  </section>

  <!-- ════════════════ Phase 0 — ScoreBar ════════════════ -->
  <section class="primitive" id="scorebar">
    <p class="sec-kicker">Phase 0 · Data</p>
    <h2>ScoreBar</h2>
    <p class="blurb">Props <code>{ value, phase, size }</code>. Color derives from
       <code>phaseColor(phase)</code> internally; track is always <code>bg-hair</code>;
       <code>value</code> is clamped 0..1 and rendered as the fill width.</p>
    <p class="shows">Catalog must show: <b>bar and compact sizes</b>, <b>several phase values</b>
       (proving color derives from phase), and <b>clamping at value 0 / 0.5 / 1 / out-of-range</b>.</p>
    <div class="specimen">
      <div class="row">
        <div class="cell"><div class="scorebar size-bar"><i style="width:0%;  background:var(--phase-Pn3m)"></i></div><span class="cap">bar · Pn3m · value 0 → 0%</span></div>
        <div class="cell"><div class="scorebar size-bar"><i style="width:50%; background:var(--phase-Im3m)"></i></div><span class="cap">bar · Im3m · value 0.5 → 50%</span></div>
        <div class="cell"><div class="scorebar size-bar"><i style="width:100%;background:var(--phase-Square)"></i></div><span class="cap">bar · Square · value 1 → 100%</span></div>
      </div>
      <div class="row">
        <div class="cell"><div class="scorebar size-bar"><i style="width:100%;background:var(--phase-Lamellar)"></i></div><span class="cap">bar · Lamellar · value 1.4 → clamped 100%</span></div>
        <div class="cell"><div class="scorebar size-bar"><i style="width:0%;  background:var(--phase-Ia3d)"></i></div><span class="cap">bar · Ia3d · value -0.2 → clamped 0%</span></div>
      </div>
      <div class="row">
        <div class="cell"><div class="scorebar size-compact"><i style="width:80%; background:var(--phase-Pn3m)"></i></div><span class="cap">compact · Pn3m · value 0.8</span></div>
        <div class="cell"><div class="scorebar size-compact"><i style="width:35%; background:var(--phase-Im3m)"></i></div><span class="cap">compact · Im3m · value 0.35</span></div>
        <div class="cell"><div class="scorebar size-compact"><i style="width:100%;background:var(--phase-Square)"></i></div><span class="cap">compact · Square · value 1</span></div>
      </div>
    </div>
  </section>

  <!-- ════════════════ Phase 0 — Dot ════════════════ -->
  <section class="primitive" id="dot">
    <p class="sec-kicker">Phase 0 · Status</p>
    <h2>Dot</h2>
    <p class="blurb">Props <code>{ label, tone, size }</code>. Tone supplies a semantic color
       (<code>accent · success · muted · neutral</code>); keeps <code>role="img"</code> + <code>aria-label</code>,
       and a decorative dot may be <code>aria-hidden</code> to skip the label.</p>
    <p class="shows">Catalog must show: <b>every tone × size</b>, plus the <b>hollow/neutral ring</b>
       and an <b>aria-hidden decorative example</b>.</p>
    <div class="specimen">
      <div class="row">
        <div class="cell"><span class="dot size-sm tone-accent" role="img" aria-label="accent"></span><span class="cap">accent · sm</span></div>
        <div class="cell"><span class="dot size-sm tone-success" role="img" aria-label="success"></span><span class="cap">success · sm</span></div>
        <div class="cell"><span class="dot size-sm tone-muted" role="img" aria-label="muted"></span><span class="cap">muted ring · sm</span></div>
        <div class="cell"><span class="dot size-sm tone-neutral" role="img" aria-label="neutral"></span><span class="cap">neutral · sm</span></div>
      </div>
      <div class="row">
        <div class="cell"><span class="dot size-xs tone-accent" role="img" aria-label="accent"></span><span class="cap">accent · xs</span></div>
        <div class="cell"><span class="dot size-xs tone-success" role="img" aria-label="success"></span><span class="cap">success · xs</span></div>
        <div class="cell"><span class="dot size-xs tone-muted" role="img" aria-label="muted"></span><span class="cap">muted ring · xs</span></div>
        <div class="cell"><span class="dot size-xs tone-neutral" role="img" aria-label="neutral"></span><span class="cap">neutral · xs</span></div>
      </div>
      <div class="row">
        <div class="cell"><span class="dot size-sm tone-neutral" aria-hidden="true"></span><span class="cap">decorative · aria-hidden (no label)</span></div>
      </div>
    </div>
  </section>

  <!-- ════════════════ Phase 0 — Toast ════════════════ -->
  <section class="primitive" id="toast">
    <p class="sec-kicker">Phase 0 · Feedback</p>
    <h2>Toast</h2>
    <p class="blurb">Four kinds (<code>info · success · warning · error</code>) with a dismiss button.
       Kind is conveyed by a <b>leading status icon + word + full hairline border</b> — the
       side-stripe fix (no edge-hue-only severity).</p>
    <p class="shows">Catalog must show: <b>all four kinds with the leading status icon + word and the dismiss button</b>.
       (Convention note: the prior <code>border-l-4</code> hue-only stripe is replaced by the
       icon+word+hairline treatment, shipped in Phase 0 per the locked side-stripe decision.)</p>
    <div class="specimen">
      <div class="row" style="flex-direction:column; align-items:stretch; gap:12px;">
        <div class="toast kind-info" role="status">
          <span class="status-icon" aria-hidden="true">i</span><span class="status-word">Info</span>
          <span class="msg">Re-analysis queued.</span>
          <button class="dismiss" aria-label="Dismiss">×</button>
        </div>
        <div class="toast kind-success" role="status">
          <span class="status-icon" aria-hidden="true">✓</span><span class="status-word">Success</span>
          <span class="msg">Phase call saved.</span>
          <button class="dismiss" aria-label="Dismiss">×</button>
        </div>
        <div class="toast kind-warning" role="status">
          <span class="status-icon" aria-hidden="true">!</span><span class="status-word">Warning</span>
          <span class="msg">Some indices are stale.</span>
          <button class="dismiss" aria-label="Dismiss">×</button>
        </div>
        <div class="toast kind-error" role="status">
          <span class="status-icon" aria-hidden="true">×</span><span class="status-word">Error</span>
          <span class="msg">Couldn't save. Try refreshing.</span>
          <button class="dismiss" aria-label="Dismiss">×</button>
        </div>
      </div>
    </div>
  </section>

  <!-- ════════════════ Appended in later waves ════════════════ -->
  <div class="foot">
    <p class="sec-kicker">Appended in later migration waves</p>
    <ul>
      <li>Phase 1 — <code>SegmentedControl&lt;T&gt;</code>, <code>PhaseChip</code>, <code>PhaseStrip</code></li>
      <li>Phase 2 — <code>Kicker</code>, <code>ModalShell</code>, <code>IconButton</code>, <code>Card</code>/<code>Plate</code></li>
    </ul>
    <p style="color:var(--color-ink-faint); font-size:11.5px; font-style:italic; margin-top:10px;">
      Each primitive's section is added in the same wave that extracts it; this list is the running
      backlog of sections not yet authored.</p>
  </div>

</div>
</body>
</html>
```

- [ ] **Step 3: Sanity-check the page renders** — Run
  `node -e "const fs=require('fs');const h=fs.readFileSync('docs/design-system.html','utf8');if(!/<\/html>\s*$/.test(h.trim()))throw new Error('truncated');for(const id of ['button','scorebar','dot','toast'])if(!h.includes('id=\"'+id+'\"'))throw new Error('missing section '+id);console.log('catalog OK:',h.length,'bytes, 4 Phase-0 sections present');"`
  → Expected: prints `catalog OK: <N> bytes, 4 Phase-0 sections present` (well-formed, no
  truncation, all four Phase-0 sections present). This is a doc-integrity check, not a unit test.

- [ ] **Step 4: Commit** —
  `git add docs/design-system.html && git commit -m "docs(design-system): scaffold static catalog + Phase-0 primitive sections (Button/ScoreBar/Dot/Toast)"`

### Task 8: Button — rename variant axis (primary→solid|accent|ghost|danger), fix text-white→text-paper, add data-variant

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ui/Button.tsx`
- Test: `packages/HimalayaUI/frontend/test/ui/Button.test.tsx` (new)

Closed-look pattern (locked): `Record<Variant,string>` map; `className` is PLACEMENT-ONLY. `solid` = `bg-ink border border-ink text-paper` (DESIGN.md §Buttons confirm action; **terracotta→ink** sign-off §8); `accent` = `bg-accent border border-accent text-paper` (the single grease-pencil action); `ghost`/`danger` unchanged. The lone `text-white` leak in the whole src tree is `Button.tsx:11` → `text-paper`. Add a `data-variant` attribute so tests assert semantics, not classes.

- [ ] **Step 1: Write the failing test** — create `test/ui/Button.test.tsx`:
  ```tsx
  import { render, screen, fireEvent } from "@testing-library/react";
  import { describe, it, expect, vi } from "vitest";
  import { Button } from "../../src/components/ui/Button";

  describe("<Button>", () => {
    it("renders a button with its children and defaults to the ghost variant", () => {
      render(<Button>Save</Button>);
      const btn = screen.getByRole("button", { name: "Save" });
      expect(btn.getAttribute("data-variant")).toBe("ghost");
    });

    it("reflects each variant on data-variant", () => {
      const { rerender } = render(<Button variant="solid">x</Button>);
      expect(screen.getByRole("button").getAttribute("data-variant")).toBe("solid");
      rerender(<Button variant="accent">x</Button>);
      expect(screen.getByRole("button").getAttribute("data-variant")).toBe("accent");
      rerender(<Button variant="danger">x</Button>);
      expect(screen.getByRole("button").getAttribute("data-variant")).toBe("danger");
    });

    it("forwards disabled and fires onClick", () => {
      const onClick = vi.fn();
      const { rerender } = render(<Button onClick={onClick}>go</Button>);
      fireEvent.click(screen.getByRole("button"));
      expect(onClick).toHaveBeenCalledTimes(1);
      rerender(<Button disabled onClick={onClick}>go</Button>);
      expect(screen.getByRole("button")).toBeDisabled();
    });
  });
  ```
- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/Button.test.tsx` → Expected: FAIL (`data-variant` is undefined; `variant="solid"` is a tsc/type error against the current `"primary"|"ghost"|"danger"` union).
- [ ] **Step 3: Implement** — replace the whole of `src/components/ui/Button.tsx` with:
  ```tsx
  import type { ButtonHTMLAttributes } from "react";

  export type ButtonVariant = "solid" | "accent" | "ghost" | "danger";

  interface ButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
    variant?: ButtonVariant;
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

  export function Button({
    variant = "ghost",
    className = "",
    children,
    ...props
  }: ButtonProps): JSX.Element {
    return (
      <button
        data-variant={variant}
        className={`rounded-md px-2.5 py-1 transition-colors ${variantClass[variant]} ${className}`}
        {...props}
      >
        {children}
      </button>
    );
  }
  ```
- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/Button.test.tsx` → Expected: PASS.
- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/Button.tsx packages/HimalayaUI/frontend/test/ui/Button.test.tsx && git commit -m "Button: solid|accent|ghost|danger variants, text-paper, data-variant"`

---

### Task 9: Migrate the 5 `variant="primary"` Button sites to `variant="solid"`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx:72`
- Modify: `packages/HimalayaUI/frontend/src/components/SpeculativeBuilder.tsx:270`
- Modify: `packages/HimalayaUI/frontend/src/components/OnboardingFlow.tsx:292,352,356`

Transformation pattern (every site): `variant="primary"` → `variant="solid"`. All five are confirm/continue actions (Re-analyze, Save, Continue, Next, Done) → `solid` (ink), per DESIGN.md §Buttons. No import changes; `Button` is already imported from `./ui` in each file. Verify the line before editing — confirmed at authoring time via `grep -n 'variant="primary"'`.

- [ ] **Step 1: Run the build to capture the pre-migration tsc state** — `cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json` → Expected: FAIL (5 errors: `"primary"` is not assignable to `ButtonVariant` after the previous task's rename).
- [ ] **Step 2: StaleIndicesBanner.tsx:72** — change `variant="primary"` → `variant="solid"`.
- [ ] **Step 3: SpeculativeBuilder.tsx:270** — change `variant="primary"` → `variant="solid"`.
- [ ] **Step 4: OnboardingFlow.tsx:292** — change `<Button variant="primary" onClick={onSubmit} data-testid="onboarding-continue">` → `variant="solid"`.
- [ ] **Step 5: OnboardingFlow.tsx:352** — change `<Button variant="primary" onClick={onNext} data-testid="tutorial-next">` → `variant="solid"`.
- [ ] **Step 6: OnboardingFlow.tsx:356** — change `<Button variant="primary" onClick={onDone} data-testid="tutorial-done">` → `variant="solid"`.
- [ ] **Step 7: Verify tsc + grep clean** — `cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json && ! grep -rn 'variant="primary"' src` → Expected: PASS (no tsc errors, no remaining `variant="primary"` in `src`).
- [ ] **Step 8: Run the affected component suites** — `cd packages/HimalayaUI/frontend && npx vitest run test/StaleIndicesBanner.test.tsx test/OnboardingFlow.test.tsx test/SpeculativeBuilder.test.tsx` → Expected: PASS (these assert testids/behavior, not the old terracotta classes; if any pins `bg-accent` on these buttons, rewrite that assert to `data-variant="solid"` here).
- [ ] **Step 9: Commit** — `git add packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx packages/HimalayaUI/frontend/src/components/SpeculativeBuilder.tsx packages/HimalayaUI/frontend/src/components/OnboardingFlow.tsx && git commit -m "Migrate the 5 variant=primary Button sites to variant=solid (terracotta→ink)"`

---

### Task 10: Fold the two inline scoping confirm buttons into `<Button variant="solid">` and rewrite their class-string asserts to data-variant

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ScopingConfirmModal.tsx:51-55`
- Modify: `packages/HimalayaUI/frontend/src/components/ScopingFoot.tsx:37-46`
- Test: `packages/HimalayaUI/frontend/test/scoping.test.tsx:113-120`
- Test: `packages/HimalayaUI/frontend/test/ScopingFoot.test.tsx:28-34`

Both sites are hand-rolled `<button className="… border border-ink bg-ink … text-paper …">` confirm actions — exactly `Button variant="solid"`. Fold them onto the primitive (preserve `data-testid`, `type`, `onClick`, `disabled`, `title`). Their tests currently pin `bg-ink`/`text-paper`/`.not bg-accent` class strings on the inline buttons; per the locked decision these must be rewritten to `data-variant` semantics in THIS task (the primitive owns the look; consumers no longer carry those utilities directly). Verify lines before editing.

- [ ] **Step 1: Rewrite the failing asserts first (scoping.test.tsx:113-120)** — replace the test body of `it("the confirm-build button uses Print ink tokens, not the accent (S-B)", …)` so it reads:
  ```tsx
  it("the confirm-build button is the solid (ink) Button variant, not the accent (S-B)", () => {
    render(<ScopingConfirmModal open orderingKey="ratio" count={2}
      onConfirm={() => {}} onClose={() => {}} />);
    const btn = screen.getByTestId("scoping-confirm-build");
    expect(btn.getAttribute("data-variant")).toBe("solid");
  });
  ```
- [ ] **Step 2: Rewrite the failing asserts (ScopingFoot.test.tsx:28-34)** — replace the test body of `it("the build button carries Print ink tokens, not ice-blue accent (S-B/S-C)", …)` so it reads:
  ```tsx
  it("the build button is the solid (ink) Button variant, not the ice-blue accent (S-B/S-C)", () => {
    render(<ScopingFoot flagCount={0} memberCount={1} keyLabel="ratio" canBuild onBuild={() => {}} />);
    const btn = screen.getByTestId("scoping-open-confirm");
    expect(btn.getAttribute("data-variant")).toBe("solid");
  });
  ```
- [ ] **Step 3: Run both, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/scoping.test.tsx test/ScopingFoot.test.tsx` → Expected: FAIL (the inline `<button>`s carry no `data-variant`).
- [ ] **Step 4: Fold ScopingConfirmModal:51-55** — add `import { Button } from "./ui";` to the top, then replace the inline confirm `<button …>Confirm &amp; build →</button>` (lines 51-55) with:
  ```tsx
  <Button type="button" data-testid="scoping-confirm-build" variant="solid"
    onClick={onConfirm} className="px-3 py-1.5 text-sm font-semibold">
    Confirm &amp; build →
  </Button>
  ```
  (Leave the `Cancel` button alone — it is a bespoke `hover:underline` link-button, not in scope.)
- [ ] **Step 5: Fold ScopingFoot:37-46** — add `import { Button } from "./ui";` to the top, then replace the inline `<button …>Confirm &amp; build →</button>` (lines 37-46) with:
  ```tsx
  <Button
    type="button"
    data-testid="scoping-open-confirm"
    variant="solid"
    disabled={!canBuild}
    onClick={onBuild}
    title={canBuild ? undefined : "Check the flagged values above before building"}
    className="shrink-0 px-[18px] py-[11px] text-[13px] font-semibold disabled:cursor-not-allowed disabled:border-hair-strong disabled:bg-paper-sunk disabled:text-ink-faint"
  >
    Confirm &amp; build →
  </Button>
  ```
  (Placement-only `className` retained: sizing/margins + the disabled-state utilities, which the primitive does not own.)
- [ ] **Step 6: Run both, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/scoping.test.tsx test/ScopingFoot.test.tsx` → Expected: PASS (the other scoping asserts — Escape, disabled gating, success-token dot, onBuild — still pass: `Button` forwards `onClick`/`disabled`/`title`/`data-testid`).
- [ ] **Step 7: Commit** — `git add packages/HimalayaUI/frontend/src/components/ScopingConfirmModal.tsx packages/HimalayaUI/frontend/src/components/ScopingFoot.tsx packages/HimalayaUI/frontend/test/scoping.test.tsx packages/HimalayaUI/frontend/test/ScopingFoot.test.tsx && git commit -m "Fold scoping confirm buttons into Button variant=solid; assert data-variant"`

---

### Task 11: ScoreBar — rewrite to `{ value, phase, size }`, derive color via phaseColor, add data-phase

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ui/ScoreBar.tsx` (full rewrite)
- Test: `packages/HimalayaUI/frontend/test/ui/ScoreBar.test.tsx` (full rewrite — breaking API)

Breaking API per the locked decision: `score`→`value`, `color`→`phase`; color computed internally from `phaseColor(phase)`; `size: "bar" | "compact"` (`bar` = `h-1 w-full`, `compact` = `h-[3.5px] w-[46px]`, the reviewed Fixed-Scale named tokens); track = `bg-hair`; keep the `data-score-bar` hook on the fill; add `data-phase` so tests assert the phase semantic instead of a raw color. `test/ui/ScoreBar.test.tsx` is the sole importer and becomes a tsc error the instant props change — it MUST be rewritten in the same commit (assert `style.width` clamping + `data-phase`, never `style.background`/`bg-hair`/`h-1`).

- [ ] **Step 1: Rewrite the failing test** — replace the whole of `test/ui/ScoreBar.test.tsx` with:
  ```tsx
  import { render } from "@testing-library/react";
  import { describe, it, expect } from "vitest";
  import { ScoreBar } from "../../src/components/ui/ScoreBar";

  describe("<ScoreBar>", () => {
    const fill = () => document.querySelector("[data-score-bar]") as HTMLElement;

    it("renders width proportional to value", () => {
      render(<ScoreBar value={0.75} phase="Pn3m" />);
      expect(fill().style.width).toBe("75%");
    });

    it("clamps value above 1 to 100%", () => {
      render(<ScoreBar value={1.4} phase="Im3m" />);
      expect(fill().style.width).toBe("100%");
    });

    it("clamps a negative value to 0%", () => {
      render(<ScoreBar value={-0.2} phase="Ia3d" />);
      expect(fill().style.width).toBe("0%");
    });

    it("renders an empty bar for value 0", () => {
      render(<ScoreBar value={0} phase="Lamellar" />);
      expect(fill().style.width).toBe("0%");
    });

    it("exposes the phase via data-phase so color derives internally", () => {
      render(<ScoreBar value={0.5} phase="Hexagonal" />);
      expect(fill().getAttribute("data-phase")).toBe("Hexagonal");
    });
  });
  ```
- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/ScoreBar.test.tsx` → Expected: FAIL (`value`/`phase` are not valid props; `data-phase` absent — tsc + runtime).
- [ ] **Step 3: Implement** — replace the whole of `src/components/ui/ScoreBar.tsx` with:
  ```tsx
  import { phaseColor } from "../../phases";

  export type ScoreBarSize = "bar" | "compact";

  interface ScoreBarProps {
    value: number;
    phase: string;
    size?: ScoreBarSize;
    className?: string;
  }

  const sizeClass: Record<ScoreBarSize, string> = {
    bar: "h-1 w-full",
    compact: "h-[3.5px] w-[46px]",
  };

  export function ScoreBar({
    value,
    phase,
    size = "bar",
    className = "",
  }: ScoreBarProps): JSX.Element {
    const pct = `${Math.round(Math.min(1, Math.max(0, value)) * 100)}%`;
    return (
      <div className={`overflow-hidden rounded-full bg-hair ${sizeClass[size]} ${className}`}>
        <i
          data-score-bar
          data-phase={phase}
          className="block h-full"
          style={{ width: pct, background: phaseColor(phase) }}
        />
      </div>
    );
  }
  ```
  (Verify the relative import depth: `src/components/ui/ScoreBar.tsx` → `../../phases`.)
- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/ScoreBar.test.tsx` → Expected: PASS.
- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/ScoreBar.tsx packages/HimalayaUI/frontend/test/ui/ScoreBar.test.tsx && git commit -m "ScoreBar: rewrite to {value,phase,size}; color via phaseColor; data-phase"`

---

### Task 12: Adopt ScoreBar at the two inline PhasePanel phase bars

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx:58-63` (PhaseCallBlock, → `size="bar"`)
- Modify: `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx:153-155` (CandidateRow, → `size="compact"`)
- Test: `packages/HimalayaUI/frontend/test/PhasePanel.test.tsx` (run existing suite as the guard)

`PhasePanel` already imports `phaseColor` and computes `const color = phaseColor(index.phase)` + `const score = index.score ?? 0` in both `PhaseCallBlock` (line 36-37) and `CandidateRow` (line 89-91). Replace each inline track+fill `<div>` with `<ScoreBar>`; the `color` local at each site becomes dead only if unused elsewhere — it is still used by the serif title/score text, so leave it. `SamplesPage:131-139` is a progress meter (fixed `bg-print-accent`, no phase) — **EXCLUDE** it (future Meter primitive). Verify lines before editing.

- [ ] **Step 1: Add the import** — in `src/components/PhasePanel.tsx`, add `ScoreBar` to the existing `./ui` import line (it imports `HintText` from `./ui` today — confirm and extend to `import { HintText, ScoreBar } from "./ui";`).
- [ ] **Step 2: PhaseCallBlock (lines 58-63)** — replace the inline wide bar:
  ```tsx
  <div className="mt-2 h-1 overflow-hidden rounded-full bg-hair">
    <i
      className="block h-full"
      style={{ width: `${Math.round(score * 100)}%`, background: color }}
    />
  </div>
  ```
  with:
  ```tsx
  <ScoreBar value={score} phase={index.phase} size="bar" className="mt-2" />
  ```
- [ ] **Step 3: CandidateRow (lines 153-155)** — replace the inline compact bar:
  ```tsx
  <div className="mt-1 h-[3.5px] w-[46px] overflow-hidden rounded-full bg-hair">
    <i className="block h-full" style={{ width: `${Math.round(score * 100)}%`, background: color }} />
  </div>
  ```
  with:
  ```tsx
  <ScoreBar value={score} phase={index.phase} size="compact" className="mt-1" />
  ```
- [ ] **Step 4: tsc + suite guard** — `cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json && npx vitest run test/PhasePanel.test.tsx` → Expected: PASS (no unused-`color` tsc error since `color` still feeds the title/score text; PhasePanel tests assert testids/text, not the inline bar markup). If a test pins the old `<div … bg-hair>` track markup, rewrite it to query `[data-score-bar]` + `data-phase` here.
- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/components/PhasePanel.tsx && git commit -m "Adopt ScoreBar at PhasePanel wide (bar) + compact phase bars"`

---

### Task 13: Dot — add tone + size props, allow aria-hidden decorative dots, add data-tone

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ui/Dot.tsx` (rewrite)
- Test: `packages/HimalayaUI/frontend/test/ui/Dot.test.tsx` (new)

Add `tone: "accent" | "success" | "muted" | "neutral"` (accent→`bg-print-accent`, success→`bg-success`, neutral→`bg-hair-strong`, muted→hollow ring `border-[1.5px] border-hair-strong`) and `size: "xs" | "sm"` (xs=`h-1.5 w-1.5`, sm=`h-2 w-2`, default sm). Keep `label` + `role="img"` + `aria-label`, but make them OPTIONAL when the dot is decorative (`aria-hidden`): if `aria-hidden` is set OR no `label`, omit `role="img"`/`aria-label`. Add `data-tone`. `...props` still rides through (preserves `data-testid`, `title`). Note: LoupeSidebar's `h-2.5 w-2.5` (md) is normalized to `sm` per Fixed-Scale (it adopts in a later sweep wave, not Phase 0 — out of this unit's scope).

- [ ] **Step 1: Write the failing test** — create `test/ui/Dot.test.tsx`:
  ```tsx
  import { render } from "@testing-library/react";
  import { describe, it, expect } from "vitest";
  import { Dot } from "../../src/components/ui/Dot";

  describe("<Dot>", () => {
    it("is an img with the label when labelled", () => {
      const { container } = render(<Dot label="kept" tone="success" />);
      const dot = container.querySelector("span")!;
      expect(dot.getAttribute("role")).toBe("img");
      expect(dot.getAttribute("aria-label")).toBe("kept");
      expect(dot.getAttribute("data-tone")).toBe("success");
    });

    it("reflects the tone on data-tone", () => {
      const { container } = render(<Dot label="x" tone="accent" />);
      expect(container.querySelector("span")!.getAttribute("data-tone")).toBe("accent");
    });

    it("drops role=img and aria-label when decorative (aria-hidden)", () => {
      const { container } = render(<Dot tone="accent" size="xs" aria-hidden />);
      const dot = container.querySelector("span")!;
      expect(dot.getAttribute("role")).toBeNull();
      expect(dot.getAttribute("aria-label")).toBeNull();
      expect(dot.getAttribute("aria-hidden")).toBe("true");
    });

    it("passes through data-testid and title", () => {
      const { container } = render(<Dot tone="accent" data-testid="stale-dot" title="Has stale members" aria-hidden />);
      const dot = container.querySelector("span")!;
      expect(dot.getAttribute("data-testid")).toBe("stale-dot");
      expect(dot.getAttribute("title")).toBe("Has stale members");
    });
  });
  ```
- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/Dot.test.tsx` → Expected: FAIL (`tone`/`size` are not valid props; `data-tone` absent; `label` is currently required so the decorative cases are a tsc error).
- [ ] **Step 3: Implement** — replace the whole of `src/components/ui/Dot.tsx` with:
  ```tsx
  import type { HTMLAttributes } from "react";

  export type DotTone = "accent" | "success" | "muted" | "neutral";
  export type DotSize = "xs" | "sm";

  interface DotProps extends Omit<HTMLAttributes<HTMLSpanElement>, "color"> {
    label?: string;
    tone?: DotTone;
    size?: DotSize;
  }

  const toneClass: Record<DotTone, string> = {
    accent: "bg-print-accent",
    success: "bg-success",
    neutral: "bg-hair-strong",
    muted: "border-[1.5px] border-hair-strong",
  };

  const sizeClass: Record<DotSize, string> = {
    xs: "h-1.5 w-1.5",
    sm: "h-2 w-2",
  };

  export function Dot({
    label,
    tone = "neutral",
    size = "sm",
    className = "",
    ...props
  }: DotProps): JSX.Element {
    const decorative = props["aria-hidden"] === true || props["aria-hidden"] === "true" || label == null;
    return (
      <span
        className={`inline-block shrink-0 rounded-full ${sizeClass[size]} ${toneClass[tone]} ${className}`}
        data-tone={tone}
        {...(decorative ? {} : { role: "img", "aria-label": label })}
        {...props}
      />
    );
  }
  ```
- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/Dot.test.tsx` → Expected: PASS.
- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/Dot.tsx packages/HimalayaUI/frontend/test/ui/Dot.test.tsx && git commit -m "Dot: add tone+size, allow aria-hidden decorative dots, data-tone"`

---

### Task 14: Toast — add barrel export and normalize the App.tsx import

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ui/index.ts`
- Modify: `packages/HimalayaUI/frontend/src/App.tsx:7`
- Test: `packages/HimalayaUI/frontend/test/Toast.test.tsx` (existing — extend with a barrel-import smoke assertion)

Phase-0 Toast scope is ONLY the barrel export (the `border-l-4` side-stripe fix is the separate pulled-forward side-stripe unit, not this one — do not touch the stripe here). Add `export { ToastContainer } from "./Toast";` to the barrel; change `App.tsx:7` to import from the barrel (cosmetic, functional no-op). Preserve `role=status`, `data-testid="toast"`, `data-toast-kind`, and the `aria-label="Dismiss"` close button.

- [ ] **Step 1: Add the failing barrel-import smoke test** — append to `test/Toast.test.tsx` (inside the existing top-level `describe`, or a new `describe`):
  ```tsx
  it("ToastContainer is exported from the ui barrel", async () => {
    const mod = await import("../src/components/ui");
    expect(typeof mod.ToastContainer).toBe("function");
  });
  ```
- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/Toast.test.tsx` → Expected: FAIL (`mod.ToastContainer` is `undefined`; the barrel does not re-export it yet).
- [ ] **Step 3: Add the barrel export** — in `src/components/ui/index.ts`, add the line:
  ```ts
  export { ToastContainer } from "./Toast";
  ```
- [ ] **Step 4: Normalize the App.tsx import** — change `App.tsx:7` from `import { ToastContainer } from "./components/ui/Toast";` to `import { ToastContainer } from "./components/ui";` (verify the line; `<ToastContainer />` usage at :101 is unchanged).
- [ ] **Step 5: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/Toast.test.tsx` → Expected: PASS.
- [ ] **Step 6: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/index.ts packages/HimalayaUI/frontend/src/App.tsx packages/HimalayaUI/frontend/test/Toast.test.tsx && git commit -m "Toast: add ToastContainer barrel export; normalize App import"`

---

### Task 15: Delete dead SectionLabel + Input primitives and remove their barrel exports

**Files:**
- Delete: `packages/HimalayaUI/frontend/src/components/ui/SectionLabel.tsx`
- Delete: `packages/HimalayaUI/frontend/src/components/ui/Input.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/ui/index.ts`

Both have ZERO importers (confirmed at authoring time: no `<SectionLabel>`/`<Input>` JSX, no `from ".../SectionLabel"`/`Input` imports). `SectionLabel` is subsumed by the future `Kicker`; `Input`'s real text-field cases (adornment, borderless-inline, steppers) need a designed `TextField` (deferred later wave). The `text-label` utility stays in `styles.css`. Guard the deletion with a grep so a stray importer fails loudly.

- [ ] **Step 1: Confirm zero importers (must be empty)** — `cd packages/HimalayaUI/frontend && grep -rn "SectionLabel\|from .*ui/Input\|<Input\b" src test | grep -v "ui/index.ts" | grep -v "ui/Input.tsx" | grep -v "ui/SectionLabel.tsx"` → Expected: no output (zero importers).
- [ ] **Step 2: Delete the files** — `cd packages/HimalayaUI/frontend && git rm src/components/ui/SectionLabel.tsx src/components/ui/Input.tsx`
- [ ] **Step 3: Remove the barrel exports** — in `src/components/ui/index.ts`, delete the two lines `export { Input } from "./Input";` and `export { SectionLabel } from "./SectionLabel";`.
- [ ] **Step 4: tsc + suite guard** — `cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json && npx vitest run` → Expected: PASS (nothing referenced the deleted files).
- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/index.ts && git commit -m "Delete dead SectionLabel + Input primitives; drop barrel exports"`

---

### Task 16: Lock the ui barrel public surface with a surface test

**Files:**
- Test: `packages/HimalayaUI/frontend/test/ui/barrel.test.tsx` (new)

Lock the public exports so the Toast omission can't regress and Input/SectionLabel can't be re-added silently. Assert the present exports are functions/defined and the deleted ones are absent.

- [ ] **Step 1: Write the failing test** — create `test/ui/barrel.test.tsx`:
  ```tsx
  import { describe, it, expect } from "vitest";
  import * as ui from "../../src/components/ui";

  describe("ui barrel public surface", () => {
    it("exports the adopted primitives", () => {
      expect(typeof ui.Button).toBe("function");
      expect(typeof ui.Dot).toBe("function");
      expect(typeof ui.HintText).toBe("function");
      expect(typeof ui.ScoreBar).toBe("function");
      expect(typeof ui.ToastContainer).toBe("function");
    });

    it("does NOT export the deleted primitives", () => {
      expect("Input" in ui).toBe(false);
      expect("SectionLabel" in ui).toBe(false);
    });
  });
  ```
- [ ] **Step 2: Run it, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/barrel.test.tsx` → Expected: PASS (this task runs AFTER the Toast-export and delete tasks; it is the guardrail that locks their result. If sequenced before them it FAILS, which is the intended ordering signal).
- [ ] **Step 3: Commit** — `git add packages/HimalayaUI/frontend/test/ui/barrel.test.tsx && git commit -m "Lock ui barrel public surface (Toast in; Input/SectionLabel out)"`

---

### Task 17: Full unit gate — type-check + full Vitest suite green

**Files:**
- (no source changes — verification only)

Final gate for the adopt-or-delete unit before handing back to the orchestrator.

- [ ] **Step 1: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (`lint:design` + `tsc --noEmit -p tsconfig.build.json` + `vite build` all green). NOTE: `lint:design` exists only after the token-foundation unit (concern 7) lands; if this unit is sequenced before that, run `npx tsc --noEmit -p tsconfig.build.json && npx vite build` instead.
- [ ] **Step 2: Full unit suite** — `cd packages/HimalayaUI/frontend && npm test` → Expected: PASS (whole Vitest suite; capture once per `test/AGENTS.md` slow-suite guidance and grep the output if it is slow).
- [ ] **Step 3: Commit (if any incidental test fixups were needed)** — `git add -A && git commit -m "Adopt-or-delete unit: full build + Vitest suite green"` (skip if the tree is already clean).

**— Phase 0, side-stripe fix (Toast + InfrastructureBanner) —**

Pull-forward of the spec §6 Phase 0 side-stripe ban: convert `Toast.tsx:85` and
`InfrastructureBanner.tsx:55` from `border-l-4` + hue-only severity to a **leading
status icon (inline glyph) + word + full hairline border** (`border border-hair`).
Kind is conveyed by the icon and an accessible severity word, NOT by an edge hue.

> NOTE — inventory drift: the inventory JSON's `concern` text says "keep the
> `border-l-4` side-stripe as-is for Phase 0 (side-stripe ban is a Phase-3 item)".
> That note is **stale**. The spec §6 Phase 0 and the brief's locked decisions both
> pull the side-stripe fix forward into Phase 0. Follow the spec/brief; ignore the
> inventory note for this unit.

Design contract (both components):
- Severity word per kind: `error → "Error"`, `warning → "Warning"`,
  `success → "Success"`, `info → "Info"`. Banner: `stuck → "Error"`,
  `showing → "Saving"`.
- Leading glyph per kind (decorative color via `text-{tone}`, but a11y severity is
  carried by an `aria-label` on the icon span AND the visible word):
  `error → "✕"`, `warning → "!"`, `success → "✓"`, `info → "i"`.
- Container border becomes `border border-hair` (full hairline). The `border-l-4`
  and the hue-only `border-{accent|success|warning|error}` edge classes are removed.
- `bg-plate`/`shadow-lg`/`rounded-md`/text/sizing and all `data-testid`,
  `data-toast-kind`/`data-state`, `role="status"`, and the Dismiss/Refresh buttons
  are preserved exactly.
- Vitest asserts the accessible severity **label per kind** (the word text and the
  icon's `aria-label`), never a stripe/edge class string.

Independence: this unit touches only `Toast.tsx` + `InfrastructureBanner.tsx` and
their two test files. It has NO dependency on the token rename, the guard, or any
primitive wave, and nothing in later waves depends on it. Ships standalone, first.

---

### Task 18: Toast — leading status icon + word + full hairline border

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ui/Toast.tsx` (KIND_CLASS at `:17-23`; container/body at `:84-90`)
- Test: `packages/HimalayaUI/frontend/test/Toast.test.tsx`

- [ ] **Step 1: Write the failing test** — append these cases inside the existing
  `describe("Toast", …)` block in `test/Toast.test.tsx` (after the existing
  `"renders a toast…"` test). They assert the severity word + the icon's accessible
  label per kind, and that the edge-hue stripe class is gone:
  ```tsx
  it("shows a severity word label per kind", () => {
    render(<ToastContainer />);
    const cases: Array<[import("../src/lib/toast").ToastKind, string]> = [
      ["error", "Error"],
      ["warning", "Warning"],
      ["success", "Success"],
      ["info", "Info"],
    ];
    for (const [kind, word] of cases) {
      act(() => {
        showToast(`msg-${kind}`, kind);
      });
      const toast = screen.getByTestId("toast");
      expect(toast).toHaveAttribute("data-toast-kind", kind);
      // visible severity word
      expect(toast).toHaveTextContent(word);
      // accessible status icon naming the severity (second channel, not hue)
      expect(toast.querySelector(`[aria-label="${word}"]`)).not.toBeNull();
      act(() => {
        fireEvent.click(screen.getByLabelText("Dismiss"));
      });
    }
  });

  it("uses a full hairline border, not a left-edge severity stripe", () => {
    render(<ToastContainer />);
    act(() => {
      showToast("hello", "error");
    });
    const toast = screen.getByTestId("toast");
    expect(toast.className).toContain("border-hair");
    expect(toast.className).not.toContain("border-l-4");
    expect(toast.className).not.toContain("border-error");
  });
  ```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/Toast.test.tsx` → Expected: FAIL (no severity word rendered; container still has `border-l-4 border-error`, no `aria-label`'d icon).

- [ ] **Step 3: Implement** — replace `KIND_CLASS` (`Toast.tsx:17-23`) with severity word/glyph maps, and rewrite the toast body (`:84-90`) to render the leading icon + word with a full hairline border. Full replacement blocks:

  Replace lines 17-23 (`KIND_CLASS`):
  ```tsx
  // Severity is conveyed by a leading status icon + word, NOT an edge hue.
  // The icon's color is decorative; the accessible severity is the aria-label
  // (on the icon) plus the visible word.
  const KIND_WORD: Record<ToastKind, string> = {
    info: "Info",
    success: "Success",
    warning: "Warning",
    error: "Error",
  };

  const KIND_GLYPH: Record<ToastKind, string> = {
    info: "i",
    success: "✓", // ✓
    warning: "!",
    error: "✕", // ✕
  };

  const KIND_ICON_COLOR: Record<ToastKind, string> = {
    info: "text-accent",
    success: "text-success",
    warning: "text-warning",
    error: "text-error",
  };
  ```

  Replace the toast `<div>` + its children (`:79-99`, the `items.map` body) with:
  ```tsx
        <div
          key={t.id}
          data-testid="toast"
          data-toast-kind={t.kind}
          role="status"
          className={
            "pointer-events-auto flex items-start gap-2 rounded-md border border-hair " +
            "bg-plate text-ink px-3 py-2 shadow-lg text-body min-w-[220px] max-w-[360px]"
          }
        >
          <span
            aria-label={KIND_WORD[t.kind]}
            className={`flex-shrink-0 font-bold leading-tight ${KIND_ICON_COLOR[t.kind]}`}
          >
            {KIND_GLYPH[t.kind]}
          </span>
          <span className="flex-1 break-words">
            <span className="font-semibold">{KIND_WORD[t.kind]}</span>{" "}
            {t.msg}
          </span>
          <button
            type="button"
            aria-label="Dismiss"
            onClick={() => dismiss(t.id)}
            className="text-ink-soft hover:text-ink leading-none px-1"
          >
            ×
          </button>
        </div>
  ```

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/Toast.test.tsx` → Expected: PASS (all prior cases plus the two new ones; `data-toast-kind`, `data-testid`, `role`, Dismiss preserved).

- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/Toast.tsx packages/HimalayaUI/frontend/test/Toast.test.tsx && git commit -m "Toast: status-icon+word+hairline severity (kill border-l-4 hue stripe)"`

---

### Task 19: InfrastructureBanner — status icon + word + hairline border, fix em-dash

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/InfrastructureBanner.tsx` (tint/container at `:44-57`; body at `:59-78`)
- Test: `packages/HimalayaUI/frontend/test/InfrastructureBanner.test.tsx`

- [ ] **Step 1: Write the failing test** — append these cases inside the existing
  `describe("InfrastructureBanner", …)` block in
  `test/InfrastructureBanner.test.tsx` (after the existing `"upgrades to stuck…"`
  test). They assert the severity word/icon per state, the corrected prose, and that
  the edge stripe is gone:
  ```tsx
  it("uses a full hairline border, not a left-edge severity stripe (showing)", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 1000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner.className).toContain("border-hair");
    expect(banner.className).not.toContain("border-l-4");
    expect(banner.className).not.toContain("border-warning");
  });

  it("labels the showing state with a Saving severity word + icon", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 1000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner).toHaveAttribute("data-state", "showing");
    expect(banner).toHaveTextContent("Saving");
    expect(banner.querySelector('[aria-label="Saving"]')).not.toBeNull();
  });

  it("labels the stuck state with an Error severity word + icon and corrected prose", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 31000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner).toHaveAttribute("data-state", "stuck");
    expect(banner.querySelector('[aria-label="Error"]')).not.toBeNull();
    // Em-dash retired: two sentences, no " — ".
    expect(banner).toHaveTextContent("Couldn’t save. Try refreshing.");
    expect(banner.textContent ?? "").not.toContain("—");
    expect(screen.getByRole("button", { name: /refresh/i })).toBeInTheDocument();
  });
  ```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/InfrastructureBanner.test.tsx` → Expected: FAIL (no `aria-label` icon; banner still has `border-l-4 border-warning`; prose still `"Couldn’t save — try refreshing"` with an em-dash).

- [ ] **Step 3: Implement** — replace the tint/container className (`:44-57`) and the
  state body (`:59-78`). Full replacement blocks:

  Replace lines 43-57 (from `const stateAttr` through the closing `>` of the
  container `<div>` opening tag, i.e. the `tintClass` def and `className`):
  ```tsx
    const stateAttr = stuck ? "stuck" : "showing";
    // Severity is conveyed by the leading icon + word below, not an edge hue.

    return (
      <div
        data-testid="infrastructure-banner"
        data-state={stateAttr}
        role="status"
        className={
          "fixed left-1/2 -translate-x-1/2 bottom-4 z-40 flex items-center gap-3 " +
          "rounded-md border border-hair bg-plate text-ink px-4 py-2 shadow-lg text-body"
        }
      >
  ```

  Replace the `{stuck ? ( … ) : ( … )}` body (`:59-78`) with:
  ```tsx
        {stuck ? (
          <>
            <span aria-label="Error" className="flex-shrink-0 font-bold text-error">
              {"✕"}
            </span>
            <span>
              <span className="font-semibold">Error.</span> Couldn&rsquo;t save. Try
              refreshing.
            </span>
            <button
              type="button"
              onClick={() => window.location.reload()}
              className="rounded-md border border-error px-2 py-0.5 text-ink hover:bg-paper-sunk"
            >
              Refresh
            </button>
          </>
        ) : (
          <>
            <span
              aria-label="Saving"
              className="inline-block h-3 w-3 flex-shrink-0 rounded-full border-2 border-accent border-t-transparent animate-spin"
            />
            <span>
              <span className="font-semibold">Saving</span> your changes&hellip;
            </span>
          </>
        )}
  ```

  > The spinner span IS the leading status icon for the `showing` state; it now
  > carries `aria-label="Saving"` (was `aria-hidden`) so the severity has an
  > accessible channel. The `text-warning` hue on the spinner border is replaced by
  > the neutral `border-accent` so kind is no longer hue-only.

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/InfrastructureBanner.test.tsx` → Expected: PASS (new cases plus the existing `data-state`/threshold/settle tests; the existing `toHaveTextContent(/Couldn.t save/)` and `/Saving your changes/` regexes still match).

- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/components/InfrastructureBanner.tsx packages/HimalayaUI/frontend/test/InfrastructureBanner.test.tsx && git commit -m "InfrastructureBanner: status-icon+word+hairline severity; retire em-dash prose"`

---

### Task 20: Build gate for the side-stripe fix

**Files:**
- (No source change — gate only.)

- [ ] **Step 1: Run the full unit suite** — `cd packages/HimalayaUI/frontend && npx vitest run test/Toast.test.tsx test/InfrastructureBanner.test.tsx` → Expected: PASS (both files green).

- [ ] **Step 2: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (`tsc --noEmit -p tsconfig.build.json` + `vite build` clean; no type error from the edited components). If `lint:design` is already wired into `build` by an earlier Phase 0 unit, it must also pass — the converted components no longer emit a side-stripe violation.

- [ ] **Step 3: Commit (if anything regenerated)** — only if the build produced tracked artifacts: `git add -A && git commit -m "Side-stripe Phase 0: build gate green"`. Otherwise skip (no-op).


## Phase 1 — Core data primitives (SegmentedControl, PhaseChip, PhaseStrip)

Independent of each other; each migrates its call sites and lowers the guard baseline by exactly the hashes it removes.

### Task 21: Build the SegmentedControl<T> primitive

Builds the generic single-select segmented control: `bordered`/`plain` variants, `group`/`radiogroup` roles (driving `aria-pressed` vs `role=radio`+`aria-checked`), `sm`/`md` sizes, required `aria-label`, per-segment `data-active`/`data-value`/`role`, container `data-state` (+ optional `data-mode` honoring), 44px touch target gated behind `@media (pointer:coarse)`, roving arrow-key nav for `radiogroup`, and a `focus-visible` ring. Library owns appearance; consumer `className` is placement-only.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ui/SegmentedControl.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/ui/index.ts` (add export)
- Modify: `packages/HimalayaUI/frontend/src/styles.css` (add `@media (pointer:coarse)` touch rule + `.seg-touch` hook — only if a utility is needed; see Step 7)
- Test: `packages/HimalayaUI/frontend/test/ui/SegmentedControl.test.tsx`

- [ ] **Step 1: Write the failing test** — create `packages/HimalayaUI/frontend/test/ui/SegmentedControl.test.tsx`:

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { SegmentedControl } from "../../src/components/ui/SegmentedControl";

type Mode = "log" | "linear";
const TWO = [
  { value: "log" as Mode, label: "log q", testId: "scale-log" },
  { value: "linear" as Mode, label: "linear q", testId: "scale-linear" },
];

type GMode = "bySample" | "byPhase" | "distinct";
const THREE = [
  { value: "bySample" as GMode, label: "By sample" },
  { value: "byPhase" as GMode, label: "By phase" },
  { value: "distinct" as GMode, label: "Distinct" },
];

describe("SegmentedControl — semantics (group)", () => {
  it("renders one button per option with the option label as its name", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={() => {}}
        testId="scale-toggle"
      />,
    );
    expect(screen.getByRole("button", { name: "log q" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "linear q" })).toBeInTheDocument();
  });

  it("container carries the required aria-label and role=group (default)", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={() => {}}
        testId="scale-toggle"
      />,
    );
    const root = screen.getByTestId("scale-toggle");
    expect(root).toHaveAttribute("role", "group");
    expect(root).toHaveAttribute("aria-label", "q-axis scale");
  });

  it("group role drives aria-pressed; active=true, others=false", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={() => {}}
      />,
    );
    expect(screen.getByRole("button", { name: "log q" })).toHaveAttribute("aria-pressed", "true");
    expect(screen.getByRole("button", { name: "linear q" })).toHaveAttribute("aria-pressed", "false");
  });

  it("clicking an unselected segment fires onChange once with its value", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={onChange}
      />,
    );
    fireEvent.click(screen.getByRole("button", { name: "linear q" }));
    expect(onChange).toHaveBeenCalledTimes(1);
    expect(onChange).toHaveBeenCalledWith("linear");
  });

  it("clicking the already-active segment re-fires onChange with its value", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={onChange}
      />,
    );
    fireEvent.click(screen.getByRole("button", { name: "log q" }));
    expect(onChange).toHaveBeenCalledWith("log");
  });
});

describe("SegmentedControl — data contract (E2E selectors)", () => {
  it("each segment exposes data-value and data-active reflecting value", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="linear"
        onChange={() => {}}
      />,
    );
    const log = screen.getByRole("button", { name: "log q" });
    const lin = screen.getByRole("button", { name: "linear q" });
    expect(log).toHaveAttribute("data-value", "log");
    expect(log).toHaveAttribute("data-active", "false");
    expect(lin).toHaveAttribute("data-value", "linear");
    expect(lin).toHaveAttribute("data-active", "true");
  });

  it("container reflects active value via data-state", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="linear"
        onChange={() => {}}
        testId="scale-toggle"
      />,
    );
    expect(screen.getByTestId("scale-toggle")).toHaveAttribute("data-state", "linear");
  });

  it("applies per-option testId to the segment button", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={() => {}}
      />,
    );
    expect(screen.getByTestId("scale-log")).toBeInTheDocument();
    expect(screen.getByTestId("scale-linear")).toBeInTheDocument();
  });

  it("honors a stateAttr override to ALSO emit data-mode (GroupingModeToggle contract)", () => {
    render(
      <SegmentedControl<GMode>
        aria-label="Trace grouping mode"
        role="radiogroup"
        variant="plain"
        options={THREE}
        value="byPhase"
        onChange={() => {}}
        testId="grouping-mode"
        stateAttr="data-mode"
      />,
    );
    const root = screen.getByTestId("grouping-mode");
    expect(root).toHaveAttribute("data-mode", "byPhase");
    expect(root).toHaveAttribute("data-state", "byPhase");
  });
});

describe("SegmentedControl — radiogroup semantics", () => {
  it("radiogroup role makes segments role=radio with aria-checked", () => {
    render(
      <SegmentedControl<GMode>
        aria-label="Trace grouping mode"
        role="radiogroup"
        variant="plain"
        options={THREE}
        value="byPhase"
        onChange={() => {}}
      />,
    );
    expect(screen.getByRole("radiogroup", { name: "Trace grouping mode" })).toBeInTheDocument();
    expect(screen.getByRole("radio", { name: "By phase" })).toHaveAttribute("aria-checked", "true");
    expect(screen.getByRole("radio", { name: "By sample" })).toHaveAttribute("aria-checked", "false");
  });

  it("roving tabindex: only the active radio is tabbable", () => {
    render(
      <SegmentedControl<GMode>
        aria-label="Trace grouping mode"
        role="radiogroup"
        variant="plain"
        options={THREE}
        value="byPhase"
        onChange={() => {}}
      />,
    );
    expect(screen.getByRole("radio", { name: "By phase" })).toHaveAttribute("tabindex", "0");
    expect(screen.getByRole("radio", { name: "By sample" })).toHaveAttribute("tabindex", "-1");
    expect(screen.getByRole("radio", { name: "Distinct" })).toHaveAttribute("tabindex", "-1");
  });

  it("ArrowRight on the active radio selects the next option (onChange)", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<GMode>
        aria-label="Trace grouping mode"
        role="radiogroup"
        variant="plain"
        options={THREE}
        value="byPhase"
        onChange={onChange}
      />,
    );
    fireEvent.keyDown(screen.getByRole("radio", { name: "By phase" }), { key: "ArrowRight" });
    expect(onChange).toHaveBeenCalledWith("distinct");
  });

  it("ArrowLeft wraps from the first option to the last", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<GMode>
        aria-label="Trace grouping mode"
        role="radiogroup"
        variant="plain"
        options={THREE}
        value="bySample"
        onChange={onChange}
      />,
    );
    fireEvent.keyDown(screen.getByRole("radio", { name: "By sample" }), { key: "ArrowLeft" });
    expect(onChange).toHaveBeenCalledWith("distinct");
  });

  it("arrow keys do NOT move selection for role=group", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={onChange}
      />,
    );
    fireEvent.keyDown(screen.getByRole("button", { name: "log q" }), { key: "ArrowRight" });
    expect(onChange).not.toHaveBeenCalled();
  });
});

describe("SegmentedControl — disabled segment + variants/size as data attrs", () => {
  it("a disabled segment is not clickable and is marked disabled", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={[TWO[0], { ...TWO[1], disabled: true, title: "Open a sample" }]}
        value="log"
        onChange={onChange}
      />,
    );
    const lin = screen.getByRole("button", { name: "linear q" });
    expect(lin).toBeDisabled();
    expect(lin).toHaveAttribute("title", "Open a sample");
    fireEvent.click(lin);
    expect(onChange).not.toHaveBeenCalled();
  });

  it("reflects variant + size as data-variant / data-size on the container", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        variant="plain"
        size="md"
        options={TWO}
        value="log"
        onChange={() => {}}
        testId="scale-toggle"
      />,
    );
    const root = screen.getByTestId("scale-toggle");
    expect(root).toHaveAttribute("data-variant", "plain");
    expect(root).toHaveAttribute("data-size", "md");
  });

  it("defaults to data-variant=bordered data-size=sm", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={() => {}}
        testId="scale-toggle"
      />,
    );
    const root = screen.getByTestId("scale-toggle");
    expect(root).toHaveAttribute("data-variant", "bordered");
    expect(root).toHaveAttribute("data-size", "sm");
  });
});
```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/SegmentedControl.test.tsx` → Expected: FAIL (module `src/components/ui/SegmentedControl` does not exist).

- [ ] **Step 3: Implement** — create `packages/HimalayaUI/frontend/src/components/ui/SegmentedControl.tsx`:

```tsx
import type { ReactNode, KeyboardEvent } from "react";

/**
 * SegmentedControl<T> — the canonical single-select button group.
 *
 * Closed look / open placement (the Button.tsx pattern): the library owns
 * appearance via `variant`/`size` and the canonical ink-on-paper active fill
 * (DESIGN.md §211/§240 — the active segment is `bg-ink text-paper`, never the
 * recessed `bg-paper-sunk` "L-5/B-A defect"). The consumer's `className` is
 * PLACEMENT-ONLY (margin / width / flex position); the design guard bans
 * appearance utilities there.
 *
 * `role` drives child semantics: "group" -> button[aria-pressed];
 * "radiogroup" -> button[role=radio][aria-checked] with WAI-ARIA roving
 * tabindex + arrow-key navigation. The container always reflects the active
 * value via `data-state` (and, if `stateAttr` is given, an extra aliased
 * attribute e.g. `data-mode` for GroupingModeToggle's E2E contract); each
 * segment reflects `data-active`/`data-value`.
 *
 * 44px touch target is folded in via the `.seg-segment` class behind a
 * `@media (pointer:coarse)` rule in styles.css, so dense desktop toolbars
 * keep their compact density on a fine pointer.
 */

export interface SegmentOption<T extends string> {
  value: T;
  label: ReactNode;
  /** Disabled individual segment (e.g. Loupe with no sample selected). */
  disabled?: boolean;
  /** Per-segment title/tooltip. */
  title?: string;
  /** Stable test id for E2E, applied to the segment button. */
  testId?: string;
}

export type SegmentedVariant = "bordered" | "plain";
export type SegmentedSize = "sm" | "md";

export interface SegmentedControlProps<T extends string> {
  options: ReadonlyArray<SegmentOption<T>>;
  value: T;
  onChange: (next: T) => void;
  /** group | radiogroup — drives child role + keyboard model. Default "group". */
  role?: "group" | "radiogroup";
  variant?: SegmentedVariant;
  size?: SegmentedSize;
  /** Required: names the control for assistive tech. */
  "aria-label": string;
  /** Container test id. */
  testId?: string;
  /**
   * Extra container attribute name aliasing the active value, in addition to
   * the always-present `data-state` (e.g. "data-mode" for GroupingModeToggle).
   */
  stateAttr?: string;
  /** PLACEMENT-ONLY: margin / width / flex-grid position. No appearance utils. */
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

const containerClass: Record<SegmentedVariant, string> = {
  bordered: "inline-flex overflow-hidden rounded border border-hair-strong divide-x divide-hair-strong",
  plain: "inline-flex items-center gap-1",
};

const sizeClass: Record<SegmentedSize, string> = {
  sm: "text-xs px-3 py-1.5",
  md: "text-sm px-3.5 py-2",
};

const segmentBaseClass: Record<SegmentedVariant, string> = {
  bordered: "",
  plain: "rounded border border-transparent",
};

export function SegmentedControl<T extends string>({
  options,
  value,
  onChange,
  role = "group",
  variant = "bordered",
  size = "sm",
  "aria-label": ariaLabel,
  testId,
  stateAttr,
  className = "",
}: SegmentedControlProps<T>): JSX.Element {
  const isRadio = role === "radiogroup";

  const move = (delta: number, from: number): void => {
    const n = options.length;
    if (n === 0) return;
    let i = from;
    for (let step = 0; step < n; step++) {
      i = (i + delta + n) % n;
      if (!options[i].disabled) {
        onChange(options[i].value);
        return;
      }
    }
  };

  const onKeyDown = (e: KeyboardEvent<HTMLButtonElement>, idx: number): void => {
    if (!isRadio) return;
    if (e.key === "ArrowRight" || e.key === "ArrowDown") {
      e.preventDefault();
      move(1, idx);
    } else if (e.key === "ArrowLeft" || e.key === "ArrowUp") {
      e.preventDefault();
      move(-1, idx);
    }
  };

  const containerProps: Record<string, string> = {
    "data-state": value,
    "data-variant": variant,
    "data-size": size,
  };
  if (stateAttr) containerProps[stateAttr] = value;
  if (testId) containerProps["data-testid"] = testId;

  return (
    <div
      role={role}
      aria-label={ariaLabel}
      className={cx(containerClass[variant], className)}
      {...containerProps}
    >
      {options.map((opt, idx) => {
        const active = opt.value === value;
        const selectedProps = isRadio
          ? { role: "radio" as const, "aria-checked": active, tabIndex: active ? 0 : -1 }
          : { "aria-pressed": active };
        return (
          <button
            key={opt.value}
            type="button"
            disabled={opt.disabled}
            title={opt.title}
            data-testid={opt.testId}
            data-value={opt.value}
            data-active={active ? "true" : "false"}
            onClick={() => onChange(opt.value)}
            onKeyDown={(e) => onKeyDown(e, idx)}
            className={cx(
              "seg-segment transition-colors focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
              sizeClass[size],
              segmentBaseClass[variant],
              active
                ? "bg-ink text-paper"
                : "text-ink-faint hover:text-ink hover:bg-paper-sunk",
              opt.disabled && "opacity-50 cursor-not-allowed",
            )}
            {...selectedProps}
          >
            {opt.label}
          </button>
        );
      })}
    </div>
  );
}
```

- [ ] **Step 4: Add the 44px touch rule** — append to `packages/HimalayaUI/frontend/src/styles.css` (after the type-scale region; the `.seg-segment` hook keeps dense desktop density on fine pointers and only inflates on coarse pointers):

```css
/* SegmentedControl — fold in the 44px touch target only on coarse pointers
   (DESIGN.md sign-off §8) so dense desktop toolbars keep compact density. */
@media (pointer: coarse) {
  .seg-segment {
    min-height: 44px;
    display: inline-flex;
    align-items: center;
    justify-content: center;
  }
}
```

- [ ] **Step 5: Export from the barrel** — in `packages/HimalayaUI/frontend/src/components/ui/index.ts` add:

```ts
export { SegmentedControl } from "./SegmentedControl";
export type { SegmentOption, SegmentedVariant, SegmentedSize, SegmentedControlProps } from "./SegmentedControl";
```

- [ ] **Step 6: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/SegmentedControl.test.tsx` → Expected: PASS (all cases).

- [ ] **Step 7: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (`lint:design` + `tsc --noEmit -p tsconfig.build.json` + `vite build`). The consumer-`className` here is empty so the guard has nothing to flag.

- [ ] **Step 8: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/SegmentedControl.tsx packages/HimalayaUI/frontend/src/components/ui/index.ts packages/HimalayaUI/frontend/src/styles.css packages/HimalayaUI/frontend/test/ui/SegmentedControl.test.tsx && git commit -m "feat(ui): SegmentedControl<T> primitive (bordered/plain, group/radiogroup, roving nav, coarse-pointer 44px)"`

---

### Task 22: Migrate ScaleToggle + RepresentationToggle onto SegmentedControl

Thin both components down to a one-line wrapper around the primitive, preserving each component's public `value`/`onChange` props and exported union type (`ScaleMode`, `Representation`) so importers — `MultiTracePlot`, `SeriesBuilderRail`, `SeriesBuilderPage`, the figure-export adapters/marks — stay untouched. Active fill is already canonical ink-on-paper for both; this preserves every `data-testid` and `aria-pressed`.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ScaleToggle.tsx` (whole component, inventory `loc` ScaleToggle.tsx:14)
- Modify: `packages/HimalayaUI/frontend/src/components/RepresentationToggle.tsx` (whole component, inventory `loc` RepresentationToggle.tsx:22)
- Test: `packages/HimalayaUI/frontend/test/ScaleToggle.test.tsx`, `packages/HimalayaUI/frontend/test/RepresentationToggle.test.tsx` (existing — must still pass unchanged; they assert `aria-pressed` + `data-testid` + `onChange`, which the primitive preserves)

Transformation pattern (applied per component): keep the exported type + props interface; replace the hand-rolled `<div role="group">…</div>` body with `<SegmentedControl<Type> variant="bordered" role="group" aria-label=… testId=… options={…} value={value} onChange={onChange} />`, mapping each old button's `data-testid` onto `option.testId`.

- [ ] **Step 1: Migrate ScaleToggle** — replace the body of `packages/HimalayaUI/frontend/src/components/ScaleToggle.tsx`:

```tsx
import { SegmentedControl, type SegmentOption } from "./ui/SegmentedControl";

export type ScaleMode = "log" | "linear";

interface ScaleToggleProps {
  value: ScaleMode;
  onChange: (next: ScaleMode) => void;
}

const OPTIONS: ReadonlyArray<SegmentOption<ScaleMode>> = [
  { value: "log", label: "log q", testId: "scale-log" },
  { value: "linear", label: "linear q", testId: "scale-linear" },
];

/**
 * ScaleToggle — log/linear q-axis segmented control (R8 / B-F). Drives
 * MultiTracePlot's `xType`. Thin wrapper over the shared SegmentedControl
 * primitive; keeps its `ScaleMode` export + `value`/`onChange` contract so
 * importers are untouched. Active segment is the canonical ink-on-paper fill.
 */
export function ScaleToggle({ value, onChange }: ScaleToggleProps): JSX.Element {
  return (
    <SegmentedControl<ScaleMode>
      aria-label="q-axis scale"
      role="group"
      variant="bordered"
      testId="scale-toggle"
      options={OPTIONS}
      value={value}
      onChange={onChange}
    />
  );
}
```

- [ ] **Step 2: Migrate RepresentationToggle** — replace the body of `packages/HimalayaUI/frontend/src/components/RepresentationToggle.tsx` (preserve the long doc comment about waterfall/heatmap; only the render body and imports change):

```tsx
import { SegmentedControl, type SegmentOption } from "./ui/SegmentedControl";

export type Representation = "waterfall" | "heatmap";

interface RepresentationToggleProps {
  value: Representation;
  onChange: (next: Representation) => void;
}

const OPTIONS: ReadonlyArray<SegmentOption<Representation>> = [
  { value: "waterfall", label: "Waterfall", testId: "repr-waterfall" },
  { value: "heatmap", label: "Heatmap", testId: "repr-heatmap" },
];

/* …keep the existing waterfall/heatmap doc comment block here… */
export function RepresentationToggle({ value, onChange }: RepresentationToggleProps): JSX.Element {
  return (
    <SegmentedControl<Representation>
      aria-label="Plot representation"
      role="group"
      variant="bordered"
      testId="representation-toggle"
      options={OPTIONS}
      value={value}
      onChange={onChange}
    />
  );
}
```

- [ ] **Step 3: Run both existing tests, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ScaleToggle.test.tsx test/RepresentationToggle.test.tsx` → Expected: PASS (primitive preserves `data-testid` + `aria-pressed` + `onChange`).

- [ ] **Step 4: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS.

- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/components/ScaleToggle.tsx packages/HimalayaUI/frontend/src/components/RepresentationToggle.tsx && git commit -m "refactor(ui): migrate ScaleToggle + RepresentationToggle onto SegmentedControl"`

---

### Task 23: Migrate GroupingModeToggle onto SegmentedControl (radiogroup + data-mode)

Thin `GroupingModeToggle` onto the primitive's `plain`+`radiogroup` configuration. This is the highest-risk migration: it must preserve the container `data-mode` (read by `SeriesBuilderPage.test.tsx:271-275` and E2E) via the primitive's `stateAttr="data-mode"`, AND every per-segment `role=radio`/`aria-checked`/`data-active`/`data-value` (pinned by `test/GroupingModeToggle.test.tsx:51-64`). Keep the `GroupingMode` import + `value`/`onChange` props so `ComparePage`/`SeriesBuilderPage` are untouched.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/GroupingModeToggle.tsx` (whole component, inventory `loc` GroupingModeToggle.tsx:36)
- Test: `packages/HimalayaUI/frontend/test/GroupingModeToggle.test.tsx` (existing — must pass unchanged), `packages/HimalayaUI/frontend/test/SeriesBuilderPage.test.tsx:271-275` (existing — must pass unchanged)

- [ ] **Step 1: Migrate the component** — replace the body of `packages/HimalayaUI/frontend/src/components/GroupingModeToggle.tsx` (preserve the doc comment; the `OPTIONS` shape gains `SegmentOption` typing):

```tsx
/* …keep the existing doc comment block here… */
import { SegmentedControl, type SegmentOption } from "./ui/SegmentedControl";
import type { GroupingMode } from "../lib/comparison/coloring";

const OPTIONS: ReadonlyArray<SegmentOption<GroupingMode>> = [
  { value: "bySample", label: "By sample" },
  { value: "byPhase", label: "By phase" },
  { value: "distinct", label: "Distinct" },
];

export interface GroupingModeToggleProps {
  value: GroupingMode;
  onChange: (next: GroupingMode) => void;
}

export function GroupingModeToggle({
  value,
  onChange,
}: GroupingModeToggleProps): JSX.Element {
  return (
    <SegmentedControl<GroupingMode>
      aria-label="Trace grouping mode"
      role="radiogroup"
      variant="plain"
      testId="grouping-mode"
      stateAttr="data-mode"
      options={OPTIONS}
      value={value}
      onChange={onChange}
    />
  );
}
```

- [ ] **Step 2: Run existing tests, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/GroupingModeToggle.test.tsx test/SeriesBuilderPage.test.tsx` → Expected: PASS. Confirms `data-mode` reflection, `role=radio`+`aria-checked`, `data-active`/`data-value`, and the `getByRole("radio",{name:/by phase/i})` click path all hold. (Note: `SeriesBuilderPage.test.tsx:272` fires a click on the radio — the primitive's `onClick` still fires `onChange`; the new arrow-key handler does not change click behavior.)

- [ ] **Step 3: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS.

- [ ] **Step 4: Commit** — `git add packages/HimalayaUI/frontend/src/components/GroupingModeToggle.tsx && git commit -m "refactor(ui): migrate GroupingModeToggle onto SegmentedControl (radiogroup + data-mode alias + roving nav)"`

---

### Task 24: Migrate PlotCard inline XScaleToggle onto SegmentedControl

Replace the inline `XScaleToggle` (PlotCard.tsx:551, used at :500) with the primitive. This is a deliberate VISUAL change: the active fill converges from the banned recessed `bg-paper-sunk text-ink` to the canonical `bg-ink text-paper` (DESIGN.md §211/§240), and the control GAINS the `role="group"` + `aria-pressed` it lacks today. The hand-rolled `<span class="w-px bg-hair-strong"/>` divider is dropped (the bordered variant's `divide-x` supplies it). Preserve the `x-scale-log`/`x-scale-linear` testIds and the short `log`/`lin` labels verbatim (do NOT normalize "lin" to "linear").

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PlotCard.tsx:500` (call site) and `:546-578` (delete the inline `XScaleToggle` + `XScaleToggleProps`) — verify the line before editing; the inventory flagged drift, so match on the source text, not the numbers.
- Test: add a case to an existing PlotCard test file if one asserts the toggle; otherwise the build gate + the new `aria-pressed` is the regression guard. (No dedicated XScaleToggle test exists today.)

- [ ] **Step 1: Replace the call site** — at `packages/HimalayaUI/frontend/src/components/PlotCard.tsx:500`, change `<XScaleToggle xType={xType} onSetXType={onSetXType} />` to:

```tsx
<SegmentedControl<"log" | "linear">
  aria-label="x-axis scale"
  role="group"
  variant="bordered"
  options={[
    { value: "log", label: "log", testId: "x-scale-log" },
    { value: "linear", label: "lin", testId: "x-scale-linear" },
  ]}
  value={xType}
  onChange={onSetXType}
/>
```

- [ ] **Step 2: Delete the inline component** — remove the `interface XScaleToggleProps { … }` block and the `function XScaleToggle({ xType, onSetXType }) { … }` block (the `btn` helper, the `<span … w-px bg-hair-strong …>` divider, and the wrapper). Search-and-confirm there are no other `XScaleToggle` references: `cd packages/HimalayaUI/frontend && grep -n XScaleToggle src/components/PlotCard.tsx` → Expected: no matches after deletion.

- [ ] **Step 3: Add the import** — ensure `packages/HimalayaUI/frontend/src/components/PlotCard.tsx` imports `import { SegmentedControl } from "./ui/SegmentedControl";` (add to the existing import block; verify it is not already imported).

- [ ] **Step 4: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS.

- [ ] **Step 5: Run the PlotCard + any plot-strip tests, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/PlotCard.test.tsx` → Expected: PASS (if the file exists; otherwise run `npm test` once and grep for PlotCard). The testIds `x-scale-log`/`x-scale-linear` are preserved so E2E specs keep working.

- [ ] **Step 6: Commit** — `git add packages/HimalayaUI/frontend/src/components/PlotCard.tsx && git commit -m "refactor(ui): migrate PlotCard XScaleToggle onto SegmentedControl (ink-on-paper active fill + adds role/aria-pressed)"`

---

### Task 25: Migrate SeriesFolioPage sort onto SegmentedControl

Replace the inline 3-way SORT control (SeriesFolioPage.tsx:133-148) with the primitive, mapping each `SORT_OPTIONS` entry onto an option with `testId: series-folio-sort-${value}`. ADDS the missing wrapper `aria-label`. Active fill is already canonical ink-on-paper, `aria-pressed` is preserved. Leave the FILTER chips (`:111-128`, `rounded-full` pills) untouched per the spec carve-out (future `Chip` primitive). Reconcile the divergent `text-[11.5px] font-semibold` to the primitive's `sm` size (`text-xs`) — this is the intended size-token convergence (sign-off §8).

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesFolioPage.tsx:133` (sort block; verify the line before editing — inventory `loc` may have drifted)
- Test: `packages/HimalayaUI/frontend/test/SeriesFolioPage.test.tsx` (if present — must still find `series-folio-sort-*` testIds + `aria-pressed`); otherwise the build gate guards it.

- [ ] **Step 1: Replace the sort block** — at `packages/HimalayaUI/frontend/src/pages/SeriesFolioPage.tsx`, replace the `<div className="flex overflow-hidden rounded-md border border-hair-strong"> … {SORT_OPTIONS.map(…)} … </div>` block with:

```tsx
<SegmentedControl<FolioSort>
  aria-label="Sort series"
  role="group"
  variant="bordered"
  options={SORT_OPTIONS.map((s) => ({ value: s.value, label: s.label, testId: `series-folio-sort-${s.value}` }))}
  value={sort}
  onChange={setSort}
/>
```

- [ ] **Step 2: Add the import** — add `import { SegmentedControl } from "../components/ui/SegmentedControl";` to the import block of `packages/HimalayaUI/frontend/src/pages/SeriesFolioPage.tsx` (verify the relative path: pages → `../components/ui/`).

- [ ] **Step 3: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS.

- [ ] **Step 4: Run the SeriesFolioPage test, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/SeriesFolioPage.test.tsx` → Expected: PASS (sort testIds + `aria-pressed` preserved; filter chips unchanged). If the file does not exist, run `npm test` once and confirm no regressions referencing `series-folio-sort`.

- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/pages/SeriesFolioPage.tsx && git commit -m "refactor(ui): migrate SeriesFolioPage sort onto SegmentedControl (adds aria-label; filter chips deferred to Chip)"`

---

### Task 26: Fix AnnotationToggles banned active fill in place (interim canonical, baseline-noted)

`AnnotationToggles` is multi-select (two independent boolean toggles) and does NOT fit `SegmentedControl<T>` (single `value`/`onChange`) — it is a carve-out for a future `ToggleButton` primitive. But its active fill is the DESIGN.md L-5/B-A *defect* (`bg-paper-sunk text-ink`). Fix that in place to the canonical `bg-ink text-paper`, preserving every `data-testid`, `aria-pressed`, `data-active`, and the Zustand wiring. Document the in-place ink fill as the *interim canonical* multi-select treatment and add a guard-baseline note so the corrected toggle is not mistaken for un-migrated drift with nowhere to migrate to.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/AnnotationToggles.tsx:51-53` (the `ToggleButton` active branch) and the doc comment at `:18-22`
- Modify: the design-guard violation baseline file (created in the token-foundation unit; reference it, do NOT redefine it here) — add a note entry for AnnotationToggles' interim ink fill.
- Test: `packages/HimalayaUI/frontend/test/AnnotationToggles.test.tsx` (existing — must pass unchanged; it asserts `data-active`/`aria-pressed`, never class strings, so the fill change is invisible to it).

- [ ] **Step 1: Confirm the existing test is fill-agnostic** — `cd packages/HimalayaUI/frontend && npx vitest run test/AnnotationToggles.test.tsx` → Expected: PASS now (baseline). The suite asserts only `data-active`/`aria-pressed`/store wiring, so no test change is needed for the visual fix.

- [ ] **Step 2: Fix the active fill** — in `packages/HimalayaUI/frontend/src/components/AnnotationToggles.tsx`, change the `ToggleButton` active branch from:

```tsx
        active
          ? "bg-paper-sunk text-ink"
          : "text-ink-faint hover:text-ink hover:bg-paper-sunk",
```

to:

```tsx
        active
          ? "bg-ink text-paper"
          : "text-ink-faint hover:text-ink hover:bg-paper-sunk",
```

- [ ] **Step 3: Update the doc comment** — in the same file, update the "Styling" paragraph (`:18-22`) so it no longer claims the recessed fill is intended. Replace the sentence describing `bg-paper-sunk text-ink` active with:

```
 * **Styling — interim canonical multi-select toggle.** Two independent
 * on/off toggles (not a single-select SegmentedControl). Active uses the
 * canonical ink-on-paper fill (`bg-ink text-paper`, DESIGN.md §211/§240),
 * `text-ink-faint` at rest, ghost hover (`hover:text-ink hover:bg-paper-sunk`).
 * This is the interim canonical multi-select treatment until a dedicated
 * `ToggleButton` primitive lands (Phase 2); the design-guard baseline carries
 * a note so this corrected fill is not read as un-migrated drift. No native
 * checkbox — `aria-pressed` carries the toggle semantics, `data-active` backs
 * E2E selectors.
```

- [ ] **Step 4: Add the guard-baseline note** — append a note entry to the design-guard violation baseline (owned by the token-foundation unit; the baseline is keyed on (rule, normalized-violation-text) content-hash). Add an entry marking `AnnotationToggles.tsx` as an intentional interim multi-select toggle (no SegmentedControl migration target yet). If the baseline file does not yet exist at this point in the build sequence, leave a `- [ ]` follow-up note here and add it during the Phase-3 ratchet — but the in-place fill change ships now regardless.

- [ ] **Step 5: Run the test, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/AnnotationToggles.test.tsx` → Expected: PASS (unchanged — fill change is class-only, tests assert semantics).

- [ ] **Step 6: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS.

- [ ] **Step 7: Commit** — `git add packages/HimalayaUI/frontend/src/components/AnnotationToggles.tsx && git commit -m "fix(ui): AnnotationToggles active fill -> ink-on-paper (kill L-5/B-A defect; interim canonical multi-select, baseline-noted)"`

### Task 27: Build the PhaseChip primitive

The canonical monospace, phase-tinted data badge that **always renders the phase name** (DESIGN.md line 226). Closed look (variant/size own appearance), open placement (`className` is placement-only, appended last). Owns `phaseColor()` + the `color-mix` recipe internally so no consumer ever passes a raw color. Exposes `data-variant`/`data-size` for rendered-semantics assertions and defaults `data-testid="phase-chip"` (overridable via spread).

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ui/PhaseChip.tsx`
- Test: `packages/HimalayaUI/frontend/test/ui/PhaseChip.test.tsx`

- [ ] **Step 1: Write the failing test** — create `packages/HimalayaUI/frontend/test/ui/PhaseChip.test.tsx`:

```tsx
import { render, screen, cleanup } from "@testing-library/react";
import { afterEach, describe, it, expect } from "vitest";
import { PhaseChip } from "../../src/components/ui/PhaseChip";

afterEach(cleanup);

describe("<PhaseChip>", () => {
  it("always renders the phase name as text (the colour-blind second channel)", () => {
    render(<PhaseChip phase="Pn3m" />);
    expect(screen.getByText("Pn3m")).toBeInTheDocument();
  });

  it("renders the name for every variant × size combination", () => {
    const combos = [
      ["tint", "sm"],
      ["tint", "md"],
      ["solid", "sm"],
      ["solid", "md"],
    ] as const;
    for (const [variant, size] of combos) {
      const { unmount } = render(
        <PhaseChip phase="Im3m" variant={variant} size={size} />,
      );
      expect(screen.getByText("Im3m")).toBeInTheDocument();
      unmount();
    }
  });

  it("defaults to data-testid=phase-chip and exposes data-variant/data-size", () => {
    render(<PhaseChip phase="Ia3d" />);
    const chip = screen.getByTestId("phase-chip");
    expect(chip).toHaveTextContent("Ia3d");
    expect(chip).toHaveAttribute("data-variant", "tint");
    expect(chip).toHaveAttribute("data-size", "sm");
  });

  it("reflects an explicit variant and size on the data-attributes", () => {
    render(<PhaseChip phase="Square" variant="solid" size="md" />);
    const chip = screen.getByTestId("phase-chip");
    expect(chip).toHaveAttribute("data-variant", "solid");
    expect(chip).toHaveAttribute("data-size", "md");
  });

  it("lets a consumer override the testid via spread", () => {
    render(<PhaseChip phase="Hexagonal" data-testid="member-meta-phase-chip" />);
    expect(screen.getByTestId("member-meta-phase-chip")).toHaveTextContent(
      "Hexagonal",
    );
  });

  it("renders an unknown phase (phaseColor FALLBACK) without throwing", () => {
    render(<PhaseChip phase="Bogus" />);
    expect(screen.getByText("Bogus")).toBeInTheDocument();
  });

  it("appends a placement className to the element", () => {
    render(<PhaseChip phase="Lamellar" className="ml-2" />);
    expect(screen.getByTestId("phase-chip").className).toContain("ml-2");
  });
});
```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/PhaseChip.test.tsx` → Expected: FAIL (`Cannot find module '../../src/components/ui/PhaseChip'` — the primitive does not exist yet).

- [ ] **Step 3: Implement** — create `packages/HimalayaUI/frontend/src/components/ui/PhaseChip.tsx`:

```tsx
import type { CSSProperties, HTMLAttributes } from "react";
import { phaseColor } from "../../phases";

export type PhaseChipVariant = "tint" | "solid";
export type PhaseChipSize = "sm" | "md";

interface PhaseChipProps
  extends Omit<HTMLAttributes<HTMLSpanElement>, "color" | "children"> {
  /** Phase name, e.g. "Pn3m". Rendered as the chip's text (the always-on
   *  second channel) AND drives the hue via phaseColor() internally. */
  phase: string;
  /** "tint": phase-tinted fill + phase-colored text (default, the M-6 look).
   *  "solid": phase-colored fill + paper text (emphasis / dense rows). */
  variant?: PhaseChipVariant;
  /** "sm" = 11px mono (contact sheet, inline rows). "md" = 13px (candidate
   *  rows, builder). Default "sm". */
  size?: PhaseChipSize;
  /** PLACEMENT ONLY — margin / width / grid-or-flex position. Appearance
   *  utilities are banned by the lint:design guard. */
  className?: string;
}

const base = "inline-flex items-center rounded-sm border";

const sizeClass: Record<PhaseChipSize, string> = {
  sm: "font-mono text-[11px] font-bold px-2 py-0.5",
  md: "font-mono text-[13px] font-bold px-1.5 py-0.5",
};

function variantStyle(variant: PhaseChipVariant, color: string): CSSProperties {
  if (variant === "solid") {
    return {
      color: "var(--color-paper)",
      background: color,
      borderColor: "transparent",
    };
  }
  return {
    color,
    background: `color-mix(in oklab, ${color} 13%, transparent)`,
    borderColor: `color-mix(in oklab, ${color} 35%, transparent)`,
  };
}

export function PhaseChip({
  phase,
  variant = "tint",
  size = "sm",
  className = "",
  ...props
}: PhaseChipProps): JSX.Element {
  const color = phaseColor(phase);
  return (
    <span
      data-testid="phase-chip"
      data-variant={variant}
      data-size={size}
      className={`${base} ${sizeClass[size]} ${className}`}
      style={variantStyle(variant, color)}
      {...props}
    >
      {phase}
    </span>
  );
}
```

Note: `data-testid="phase-chip"` is declared BEFORE `{...props}` so a consumer-supplied `data-testid` (MemberMetaRow's `member-meta-phase-chip`) overrides the default via spread. `style` is set explicitly and is NOT in `props` (excluded from override is acceptable — no consumer passes `style`).

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/PhaseChip.test.tsx` → Expected: PASS (7 tests).

- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/PhaseChip.tsx packages/HimalayaUI/frontend/test/ui/PhaseChip.test.tsx && git commit -m "feat(ui): add PhaseChip primitive (phase-tinted mono data badge)"`

---

### Task 28: Export PhaseChip from the ui barrel

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ui/index.ts`

- [ ] **Step 1: Add the export** — in `packages/HimalayaUI/frontend/src/components/ui/index.ts`, add after the `ScoreBar` line:

```ts
export { PhaseChip } from "./PhaseChip";
export type { PhaseChipVariant, PhaseChipSize } from "./PhaseChip";
```

(Do not touch the SectionLabel/Input export lines here — those are removed by the Phase 0 adopt-or-delete unit, not this one. If a merge conflict on this file arises, keep both the PhaseChip additions and Phase 0's removals.)

- [ ] **Step 2: Verify build** — `cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json` → Expected: PASS (no type errors from the new export).

- [ ] **Step 3: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/index.ts && git commit -m "feat(ui): export PhaseChip from ui barrel"`

---

### Task 29: Migrate SampleStatusChip to delegate its chip to PhaseChip

SampleStatusChip **stays** as the M-6 domain wrapper. It keeps owning the `phase | null | undefined` branch and the hollow-dot "Not indexed" empty state (`status-dot` testid + the `ink-faint` text). It delegates ONLY the truthy-phase chip rendering to `<PhaseChip phase={phase} variant="tint" size="sm" />`. The default `phase-chip` testid flows through PhaseChip, so the ContactSheetRow test selector (`test/contact-sheet.test.tsx:278-285`) keeps working.

Migration transformation (apply at `SampleStatusChip.tsx:18-30`, verified current source):

```tsx
// BEFORE
if (phase) {
  const color = phaseColor(phase);
  return (
    <span
      data-testid="phase-chip"
      className="rounded px-2 py-0.5 font-mono text-[11px] font-bold"
      style={{
        color,
        background: `color-mix(in oklab, ${color} 13%, transparent)`,
      }}
    >
      {phase}
    </span>
  );
}

// AFTER
if (phase) {
  return <PhaseChip phase={phase} variant="tint" size="sm" />;
}
```

Deliberate visual deltas (per inventory canonical values — flag for design sign-off, will shift any ContactSheetRow pixel/snapshot baseline): radius `rounded` 4px → `rounded-sm` 5px; a 35% hairline border is added (SampleStatusChip previously had none). Tint stays 13%. No rendered-semantics test asserts on those, so the suite stays green.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/SampleStatusChip.tsx`

- [ ] **Step 1: Confirm coverage exists (no new test needed)** — the wrapper's two branches are already covered by `test/contact-sheet.test.tsx:269-286` ("Not indexed" + `status-dot`; phase chip text + `phase-chip` testid). Run the baseline now to capture green-before-change: `cd packages/HimalayaUI/frontend && npx vitest run test/contact-sheet.test.tsx` → Expected: PASS.

- [ ] **Step 2: Implement the delegation** — edit `packages/HimalayaUI/frontend/src/components/SampleStatusChip.tsx`. Replace the import line `import { phaseColor } from "../phases";` with `import { PhaseChip } from "./ui/PhaseChip";`, then replace the truthy-`phase` `<span>…</span>` block (lines 18-30) with the AFTER form above. The empty-state branch (the `status-dot` span + "Not indexed") is UNCHANGED.

- [ ] **Step 3: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/contact-sheet.test.tsx` → Expected: PASS (both status-cell tests still green; chip text "Pn3m" + `phase-chip` testid resolved through PhaseChip).

- [ ] **Step 4: Commit** — `git add packages/HimalayaUI/frontend/src/components/SampleStatusChip.tsx && git commit -m "refactor(ui): SampleStatusChip delegates chip to PhaseChip, keeps empty-state branch"`

---

### Task 30: Migrate MemberMetaRow's phase chip to PhaseChip

Replace the hand-rolled inline phase-chip span at `MemberMetaRow.tsx:336-346` (verify the line range before editing — the inventory `loc` is 336-346 and was confirmed against current source) with `<PhaseChip … data-testid="member-meta-phase-chip" className="shrink-0" />`. The `member-meta-phase-chip` testid is preserved via spread override; `shrink-0` is legal placement. This kills the third hand-rolled `color-mix(in oklab, …)` copy.

Migration transformation:

```tsx
// BEFORE (MemberMetaRow.tsx:336-346)
<span
  data-testid="member-meta-phase-chip"
  className="text-data-strong px-1.5 py-0.5 rounded-sm border shrink-0"
  style={{
    color: phaseColor(ci.phase),
    background: `color-mix(in oklab, ${phaseColor(ci.phase)} 10%, transparent)`,
    borderColor: `color-mix(in oklab, ${phaseColor(ci.phase)} 35%, transparent)`,
  }}
>
  {ci.phase}
</span>

// AFTER
<PhaseChip
  phase={ci.phase}
  variant="tint"
  size="md"
  data-testid="member-meta-phase-chip"
  className="shrink-0"
/>
```

Deliberate deltas (per the task spec, which locks `size="md"`): tint 10% → 13% (1-channel, imperceptible). **Size note for design sign-off:** the old chip used `text-data-strong`, which resolves to `--text-sm` = 11.5px / weight 600 (confirmed in `styles.css:272`). `size="md"` renders 13px / `font-bold`. This is a deliberate ~1.5px size increase in the dense member row, locked by the unit spec — flag it. (If design rejects the size bump, the one-line fallback is `size="sm"`, which renders 11px and is the closer match; do NOT change it silently.) The existing assertion `screen.getByText("Pn3m")` (`MemberMetaRow.test.tsx:101`) is text-only and stays green either way.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx` (the chip block at ~336-346; verify the line before editing) and its `phaseColor` import (drop only if no other `phaseColor` use remains in the file — see Step 2).

- [ ] **Step 1: Run the baseline** — `cd packages/HimalayaUI/frontend && npx vitest run test/MemberMetaRow.test.tsx` → Expected: PASS (captures green-before-change; `MemberMetaRow.test.tsx:99-107` asserts the "Pn3m" chip text + d / R² / NGC).

- [ ] **Step 2: Check remaining phaseColor usage** — `cd packages/HimalayaUI/frontend && grep -n "phaseColor" src/components/MemberMetaRow.tsx` → if the only matches are the three inside the chip block being replaced, the `phaseColor` import becomes dead; if there are other uses, keep the import. Add `import { PhaseChip } from "./ui/PhaseChip";` regardless.

- [ ] **Step 3: Implement** — edit `packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx`: add the `PhaseChip` import; replace the BEFORE span block with the AFTER `<PhaseChip … />`; remove the now-dead `phaseColor` import only if Step 2 confirmed it has no other use in this file.

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/MemberMetaRow.test.tsx test/MemberMetaGutter.reorder.test.tsx test/MemberMetaGutter.resize.test.tsx test/MemberMetaRow.collapse.test.tsx test/MemberMetaRow.drag.test.tsx` → Expected: PASS (chip text "Pn3m" still rendered; `member-meta-phase-chip` testid preserved via spread; no other MemberMetaRow behavior touched).

- [ ] **Step 5: Commit** — `git add packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx && git commit -m "refactor(ui): MemberMetaRow phase chip adopts PhaseChip (md, shrink-0)"`

---

### Task 31: Verify no regressions and confirm the explicitly-excluded sites are untouched

PhaseChip's contract is "always render the phase name," which is exactly why the following `phaseColor()` sites are **explicitly OUT of scope** and must remain unchanged in this unit (the inventory verified each; converting any would add a tint surface DESIGN.md does not call for, or structurally cannot represent a non-phase label):

- `PhasePanel.tsx:42-47` — the 23px **serif phase title** in PhaseCallBlock. Color-only display title, not a data badge.
- `PhasePanel.tsx:141-143` — the `font-mono text-[13px] font-bold` color-only **candidate-row phase label**. No tint surface; hue pairs with the score bar / checkbox state.
- `FocusReflectionsTable.tsx:168-179` — an 8px **dot + colored label** ("unindexed" → `ink-faint`). Semantic-Dot + colored-label pattern (DESIGN.md line 225) → belongs to a future SemanticDot/PhaseDot primitive.
- `SpeculativeBuilder.tsx:241-245` — the **`anchor` badge** renders the literal word "anchor", not a phase name. Structurally cannot be a PhaseChip. Separate one-off tag.
- `SeriesFolioCard.tsx:110-118` — color **swatches** (`h-2.5 w-6 rounded-sm`, name only in a `title` tooltip). Belongs to the PhaseStrip sibling unit.

- [ ] **Step 1: Confirm the excluded sites are unchanged** — `cd packages/HimalayaUI/frontend && git diff --name-only` → Expected: only `src/components/ui/PhaseChip.tsx`, `src/components/ui/index.ts`, `src/components/SampleStatusChip.tsx`, `src/components/MemberMetaRow.tsx`, `test/ui/PhaseChip.test.tsx` appear across this unit's commits. PhasePanel.tsx, FocusReflectionsTable.tsx, SpeculativeBuilder.tsx, SeriesFolioCard.tsx must NOT appear.

- [ ] **Step 2: Run the full unit suite** — `cd packages/HimalayaUI/frontend && npm test` → Expected: PASS (all Vitest files green; in particular `test/ui/PhaseChip.test.tsx`, `test/contact-sheet.test.tsx`, the five `MemberMetaRow*`/`MemberMetaGutter*` files).

- [ ] **Step 3: Type + build + design-lint gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (`lint:design` passes — PhaseChip's appearance utilities live inside the primitive, and the only consumer-side `className` values introduced are placement: `shrink-0` on MemberMetaRow, none on SampleStatusChip; `tsc --noEmit` clean; `vite build` succeeds). Note: `npm run build` includes `lint:design` only once the Phase 0 token/guard unit has wired `"lint:design"` into the `build` script — if that unit has not merged yet, run `npx tsc --noEmit -p tsconfig.build.json && npx vite build` for this gate instead.

- [ ] **Step 4: Commit (if any incidental fixups were needed)** — only if Step 2/3 surfaced a fixup: `git add -A && git commit -m "test(ui): confirm PhaseChip migration leaves excluded phaseColor sites intact"`. Otherwise no commit.

### Task 32: Create the PhaseStrip primitive (failing test first)

Builds the shared `PhaseStrip` primitive that both Scoping and Series folio adopt. API is `{ segments, size?, emptyLabel?, className? }` — NO `heading` prop (the Scoping kicker becomes a sibling element at the call site). Owns `segBackground` (null → ink-faint, coexist → 100deg 2-stop gradient, else `phaseColor`), the distinct-count caption, per-segment `aria-label`/`title`, decorative arrow `aria-hidden`, and `ps-seg`/`ps-cap` testids. Unindexed tone is canonical `var(--color-ink-faint)`; default size is `md` (7px), `sm` = legacy 8px Scoping bar; default `emptyLabel` is `"No clear phase"`.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ui/PhaseStrip.tsx`
- Test: `packages/HimalayaUI/frontend/test/PhaseStrip.test.tsx`

- [ ] **Step 1: Write the failing test** — create `packages/HimalayaUI/frontend/test/PhaseStrip.test.tsx`:

```tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { PhaseStrip } from "../src/components/ui/PhaseStrip";
import type { PhaseSegment } from "../src/components/ui/PhaseStrip";

const seg = (phase: string | null, coexistWith: string | null = null): PhaseSegment => ({
  phase,
  coexistWith,
});

describe("PhaseStrip", () => {
  it("renders one ps-seg per segment", () => {
    render(<PhaseStrip segments={[seg("Pn3m"), seg("Lamellar"), seg("Im3m")]} />);
    expect(screen.getAllByTestId("ps-seg")).toHaveLength(3);
  });

  it("captions a transition with both phase names and a decorative arrow", () => {
    render(<PhaseStrip segments={[seg("Pn3m"), seg("Lamellar")]} />);
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("Pn3m");
    expect(cap).toHaveTextContent("Lamellar");
    // arrow glyph present...
    expect(cap).toHaveTextContent("→");
    // ...but decorative (aria-hidden), so not part of the cap's accessible text channel.
    const arrow = cap.querySelector('[aria-hidden="true"]');
    expect(arrow).not.toBeNull();
    expect(arrow).toHaveTextContent("→");
  });

  it("captions a single distinct phase as '<phase> throughout'", () => {
    render(<PhaseStrip segments={[seg("Pn3m"), seg("Pn3m")]} />);
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("Pn3m");
    expect(cap).toHaveTextContent(/throughout/i);
  });

  it("treats a non-monotone strip as a transition (distinct-count, not first===last)", () => {
    // [Pn3m, Lamellar, Pn3m] has two distinct indexed phases → transition, NOT "throughout".
    render(<PhaseStrip segments={[seg("Pn3m"), seg("Lamellar"), seg("Pn3m")]} />);
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent(/→/);
    expect(cap).not.toHaveTextContent(/throughout/i);
  });

  it("renders the default empty label when no segment is indexed", () => {
    render(<PhaseStrip segments={[seg(null), seg(null)]} />);
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/no clear phase/i);
  });

  it("honors an emptyLabel override", () => {
    render(
      <PhaseStrip
        segments={[seg(null)]}
        emptyLabel="Members not yet indexed; phase preview unavailable."
      />,
    );
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/not yet indexed/i);
  });

  it("labels each segment with its phase via aria-label and title", () => {
    render(<PhaseStrip segments={[seg("Pn3m"), seg(null)]} />);
    const segs = screen.getAllByTestId("ps-seg");
    expect(segs[0]).toHaveAttribute("aria-label", "Pn3m");
    expect(segs[0]).toHaveAttribute("title", "Pn3m");
    expect(segs[1]).toHaveAttribute("aria-label", "Unindexed");
    expect(segs[1]).toHaveAttribute("title", "Unindexed");
  });

  it("labels a coexistence segment as 'A + B (coexistence)'", () => {
    render(<PhaseStrip segments={[seg("Pn3m", "Lamellar")]} />);
    const segment = screen.getAllByTestId("ps-seg")[0];
    expect(segment).toHaveAttribute("aria-label", "Pn3m + Lamellar (coexistence)");
    expect(segment).toHaveAttribute("title", "Pn3m + Lamellar (coexistence)");
  });

  it("exposes a data-size attribute reflecting the size prop (default md)", () => {
    const { rerender, container } = render(<PhaseStrip segments={[seg("Pn3m")]} />);
    expect(container.firstChild).toHaveAttribute("data-size", "md");
    rerender(<PhaseStrip segments={[seg("Pn3m")]} size="sm" />);
    expect(container.firstChild).toHaveAttribute("data-size", "sm");
  });

  it("applies the placement className to the root", () => {
    const { container } = render(<PhaseStrip segments={[seg("Pn3m")]} className="mt-5" />);
    expect(container.firstChild).toHaveClass("mt-5");
  });
});
```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/PhaseStrip.test.tsx` → Expected: FAIL (`src/components/ui/PhaseStrip.tsx` does not exist — module resolution error).

- [ ] **Step 3: Implement** — create `packages/HimalayaUI/frontend/src/components/ui/PhaseStrip.tsx`:

```tsx
import { phaseColor } from "../../phases";

/**
 * PhaseStrip — one cell per series member coloured by that member's confirmed
 * phase, captioned with the phase story (single "throughout", "A → B"
 * transition, or an empty label). The shared primitive behind the Series folio
 * card strip and the series-scoping preview strip.
 *
 * Closed look / open placement: appearance (size, colours, gradient) is owned
 * here; the consumer `className` is PLACEMENT ONLY (margin / width / grid).
 * Unindexed cells use the canonical `var(--color-ink-faint)` tone (matches the
 * mini-waterfall's UNINDEXED_COLOR). The "throughout vs transition" caption is
 * derived from the COUNT OF DISTINCT indexed phases (the truthful rule): a
 * non-monotone strip like [Pn3m, Lamellar, Pn3m] reads as a transition.
 *
 * Per-segment `aria-label`/`title` carry the phase name as the accessible
 * second channel (colour is never the sole signal — see phases.ts). The visual
 * glyph/pattern second channel is deferred to the plotting redesign.
 */

/** One per series member, in display (low→high) order. `coexistWith` drives the
 *  2-stop gradient; a null phase = unindexed (the neutral ink-faint cell). */
export interface PhaseSegment {
  phase: string | null;
  coexistWith?: string | null;
}

export type PhaseStripSize = "sm" | "md"; // sm = legacy 8px Scoping bar; md = 7px folio bar (default)

export interface PhaseStripProps {
  segments: PhaseSegment[];
  /** Discrete size. "md" (7px, default) | "sm" (8px legacy Scoping bar). */
  size?: PhaseStripSize;
  /** Caption when no segment is indexed. Default "No clear phase". */
  emptyLabel?: string;
  /** PLACEMENT ONLY: margin / width / grid position. No appearance utilities. */
  className?: string;
}

const UNINDEXED = "var(--color-ink-faint)";

function cx(...parts: Array<string | false | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

const sizeClass: Record<PhaseStripSize, string> = {
  sm: "h-2 gap-0.5", // 8px bar
  md: "h-[7px] gap-[2px]", // 7px bar (folio canonical)
};

function segBackground(seg: PhaseSegment): string {
  if (seg.phase === null) return UNINDEXED;
  if (seg.coexistWith) {
    return `linear-gradient(100deg, ${phaseColor(seg.phase)} 42%, ${phaseColor(seg.coexistWith)} 58%)`;
  }
  return phaseColor(seg.phase);
}

function segLabel(seg: PhaseSegment): string {
  if (seg.phase === null) return "Unindexed";
  if (seg.coexistWith) return `${seg.phase} + ${seg.coexistWith} (coexistence)`;
  return seg.phase;
}

export function PhaseStrip({
  segments,
  size = "md",
  emptyLabel = "No clear phase",
  className = "",
}: PhaseStripProps): JSX.Element {
  const indexed = segments
    .map((s) => s.phase)
    .filter((p): p is string => p !== null);
  const first = indexed.length > 0 ? indexed[0]! : null;
  const last = indexed.length > 0 ? indexed[indexed.length - 1]! : null;
  const distinct = new Set(indexed);

  return (
    <div className={className} data-size={size}>
      <div className={cx("flex", sizeClass[size])}>
        {segments.map((seg, i) => (
          <div
            key={i}
            data-testid="ps-seg"
            aria-label={segLabel(seg)}
            title={segLabel(seg)}
            className="flex-1 rounded-[1.5px]"
            style={{ background: segBackground(seg) }}
          />
        ))}
      </div>
      <div
        data-testid="ps-cap"
        className="mt-1.5 flex items-center gap-1.5 text-base text-ink-soft"
      >
        {first === null || last === null ? (
          <span className="font-semibold text-ink-faint">{emptyLabel}</span>
        ) : distinct.size > 1 ? (
          <>
            <span className="font-semibold" style={{ color: phaseColor(first) }}>
              {first}
            </span>
            <span className="text-ink-faint" aria-hidden="true">
              →
            </span>
            <span className="font-semibold" style={{ color: phaseColor(last) }}>
              {last}
            </span>
          </>
        ) : (
          <>
            <span className="font-semibold" style={{ color: phaseColor(first) }}>
              {first}
            </span>
            <span className="text-ink-faint">throughout</span>
          </>
        )}
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/PhaseStrip.test.tsx` → Expected: PASS (all 11 cases green).

- [ ] **Step 5: Commit** — `cd packages/HimalayaUI/frontend && git add src/components/ui/PhaseStrip.tsx test/PhaseStrip.test.tsx && git commit -m "feat(ui): PhaseStrip primitive (segments/size/emptyLabel, distinct-count caption, aria-labelled cells)"`

---

### Task 33: Export PhaseStrip from the ui barrel

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ui/index.ts`

- [ ] **Step 1: Add the export** — append to `packages/HimalayaUI/frontend/src/components/ui/index.ts` (exact text; keep alongside the other re-exports):

```ts
export { PhaseStrip } from "./PhaseStrip";
export type { PhaseSegment, PhaseStripSize, PhaseStripProps } from "./PhaseStrip";
```

- [ ] **Step 2: Type-check** — `cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json` → Expected: PASS (no errors).

- [ ] **Step 3: Commit** — `cd packages/HimalayaUI/frontend && git add src/components/ui/index.ts && git commit -m "feat(ui): export PhaseStrip from ui barrel"`

> Note: barrel imports are optional — call sites in this fragment import directly from `./ui/PhaseStrip` / `../components/ui/PhaseStrip` to avoid depending on barrel ordering decided by other Phase-0/Phase-1 fragments.

---

### Task 34: Migrate SeriesFolioCard to PhaseStrip and delete SeriesPhaseStrip

Migration pattern: the card now calls `buildPhaseStrip(members).segments` itself (the primitive takes `segments`, not `members`) and renders `<PhaseStrip>`. `buildPhaseStrip` STAYS in `lib/series/folioFigure.ts` (its `kind`/`first`/`last` model fields remain — only the caption RENDER moved into the primitive). The legacy `SeriesPhaseStrip` had an outer `data-testid="series-phase-strip"` container that the primitive does not reproduce; the folio test that asserted it is updated to assert the canonical `ps-seg` cells instead.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/SeriesFolioCard.tsx:5` (import) and `:184` (render — verify the line before editing; inventory `loc` cited)
- Delete: `packages/HimalayaUI/frontend/src/components/SeriesPhaseStrip.tsx`
- Delete: `packages/HimalayaUI/frontend/test/SeriesPhaseStrip.test.tsx`
- Modify: `packages/HimalayaUI/frontend/test/SeriesFolioCard.test.tsx:97`

- [ ] **Step 1: Update the folio test assertion to canonical testid** — in `packages/HimalayaUI/frontend/test/SeriesFolioCard.test.tsx`, the "renders the live miniature + phase strip once detail loads" case asserts `getByTestId("series-phase-strip")`. Replace that single line:

  - FROM: `    expect(screen.getByTestId("series-phase-strip")).toBeInTheDocument();`
  - TO: `    expect(screen.getAllByTestId("ps-seg").length).toBeGreaterThan(0);`

- [ ] **Step 2: Swap the import in SeriesFolioCard.tsx** — replace line 5:

  - FROM: `import { SeriesPhaseStrip } from "./SeriesPhaseStrip";`
  - TO:
    ```tsx
    import { PhaseStrip } from "./ui/PhaseStrip";
    import { buildPhaseStrip } from "../lib/series/folioFigure";
    ```

- [ ] **Step 3: Swap the render in SeriesFolioCard.tsx:184** — verify the line before editing (inventory `loc` cited):

  - FROM: `        {hasMiniature && <SeriesPhaseStrip members={members} />}`
  - TO: `        {hasMiniature && <PhaseStrip segments={buildPhaseStrip(members).segments} className="mt-3" />}`

- [ ] **Step 4: Delete the legacy component and its test** — `cd packages/HimalayaUI/frontend && git rm src/components/SeriesPhaseStrip.tsx test/SeriesPhaseStrip.test.tsx`

- [ ] **Step 5: Run the folio test + type-check** — `cd packages/HimalayaUI/frontend && npx vitest run test/SeriesFolioCard.test.tsx && npx tsc --noEmit -p tsconfig.build.json` → Expected: PASS (folio test green against `ps-seg`; no dangling `SeriesPhaseStrip` import).

- [ ] **Step 6: Commit** — `cd packages/HimalayaUI/frontend && git add -A && git commit -m "refactor(series): SeriesFolioCard adopts PhaseStrip; delete SeriesPhaseStrip"`

---

### Task 35: Migrate SeriesScopingPage to PhaseStrip and delete ScopingPhaseStrip

Migration pattern: Scoping maps its `PhaseRead{dominant,coexist}[]` into `PhaseSegment[]` inline (`{ phase: r.dominant, coexistWith: r.coexist }`). The legacy `ScopingPhaseStrip` owned an uppercase kicker heading and a longer empty-state sentence — the heading becomes an INTERIM inline sibling element above the strip (a real `Kicker` primitive lands in a later fragment), and the empty sentence passes through the `emptyLabel` prop. The `mt-5` wrapper moves onto the primitive via placement `className`. `dominantPhase`/`PhaseRead` are UNCHANGED (still Scoping's data source). NOTE the canonical behaviour changes that ship here: unindexed cells go from `var(--color-hair)` → `var(--color-ink-faint)` (darker, intended), and the caption switches to the distinct-count rule (a non-monotone series no longer reads "throughout").

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesScopingPage.tsx:20` (import) and `:333` (render — verify the line before editing; inventory `loc` cited)
- Delete: `packages/HimalayaUI/frontend/src/components/ScopingPhaseStrip.tsx`
- Delete: `packages/HimalayaUI/frontend/test/ScopingPhaseStrip.test.tsx`

- [ ] **Step 1: Swap the import in SeriesScopingPage.tsx:20**:

  - FROM: `import { ScopingPhaseStrip } from "../components/ScopingPhaseStrip";`
  - TO: `import { PhaseStrip } from "../components/ui/PhaseStrip";`

- [ ] **Step 2: Swap the render in SeriesScopingPage.tsx:333** — verify the line before editing (inventory `loc` cited). The kicker heading is rendered as an interim inline sibling (a `<div>` carrying the legacy kicker classes; a `Kicker` primitive replaces it in a later fragment):

  - FROM: `                <ScopingPhaseStrip reads={phaseReads} />`
  - TO:
    ```tsx
                <div className="mt-5">
                  {/* interim inline kicker — replace with <Kicker> when it lands */}
                  <div className="mb-1.5 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">
                    Preview: phase across the series
                  </div>
                  <PhaseStrip
                    size="sm"
                    emptyLabel="Members not yet indexed; phase preview unavailable."
                    segments={phaseReads.map((r) => ({ phase: r.dominant, coexistWith: r.coexist }))}
                  />
                </div>
    ```

- [ ] **Step 3: Delete the legacy component** — `cd packages/HimalayaUI/frontend && git rm src/components/ScopingPhaseStrip.tsx`

- [ ] **Step 4: Rewrite the Scoping test as PhaseStrip-backed coverage** — overwrite `packages/HimalayaUI/frontend/test/ScopingPhaseStrip.test.tsx` to exercise the Scoping mapping shape through the shared primitive (semantic assertions only, canonical `ps-seg`/`ps-cap` testids, canonical empty copy passed as `emptyLabel`):

```tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { PhaseStrip } from "../src/components/ui/PhaseStrip";
import type { PhaseRead } from "../src/lib/scoping/dominantPhase";

// Mirrors the SeriesScopingPage call site: PhaseRead[] → PhaseSegment[].
function strip(reads: PhaseRead[]) {
  return render(
    <PhaseStrip
      size="sm"
      emptyLabel="Members not yet indexed; phase preview unavailable."
      segments={reads.map((r) => ({ phase: r.dominant, coexistWith: r.coexist }))}
    />,
  );
}

describe("Scoping preview strip (PhaseStrip, size=sm)", () => {
  it("renders one segment per member", () => {
    strip([
      { dominant: "Pn3m", coexist: null },
      { dominant: "Lamellar", coexist: null },
    ]);
    expect(screen.getAllByTestId("ps-seg")).toHaveLength(2);
  });

  it("captions a transition first → last", () => {
    strip([
      { dominant: "Pn3m", coexist: null },
      { dominant: "Lamellar", coexist: null },
    ]);
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("Pn3m");
    expect(cap).toHaveTextContent("Lamellar");
  });

  it("captions 'throughout' when every indexed segment is one phase", () => {
    strip([
      { dominant: "Pn3m", coexist: null },
      { dominant: "Pn3m", coexist: null },
    ]);
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/throughout/i);
  });

  it("renders the not-yet-indexed empty label when no members are indexed", () => {
    strip([{ dominant: null, coexist: null }]);
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/not yet indexed/i);
  });
});
```

- [ ] **Step 5: Run both suites + type-check** — `cd packages/HimalayaUI/frontend && npx vitest run test/ScopingPhaseStrip.test.tsx test/PhaseStrip.test.tsx && npx tsc --noEmit -p tsconfig.build.json` → Expected: PASS (Scoping mapping green; no dangling `ScopingPhaseStrip` import).

- [ ] **Step 6: Commit** — `cd packages/HimalayaUI/frontend && git add -A && git commit -m "refactor(scoping): SeriesScopingPage adopts PhaseStrip; delete ScopingPhaseStrip; canonical ink-faint + distinct-count caption"`

---

### Task 36: Full-suite gate after the PhaseStrip migration

Confirms no other module referenced the deleted components/testids and the build passes end to end. `buildPhaseStrip` and `test/folioFigure.test.ts` are intentionally UNCHANGED — confirm they still pass.

**Files:** (none — verification only)

- [ ] **Step 1: Confirm no dangling references** — `cd packages/HimalayaUI/frontend && grep -rn "ScopingPhaseStrip\|SeriesPhaseStrip\|scoping-ps-seg\|scoping-ps-cap\|series-phase-strip" src test` → Expected: NO matches (empty output).

- [ ] **Step 2: Confirm the model builder test is intact** — `cd packages/HimalayaUI/frontend && npx vitest run test/folioFigure.test.ts` → Expected: PASS (buildPhaseStrip segments/kind/first/last/coexistResolver/empty cases unchanged).

- [ ] **Step 3: Run the full unit suite** — `cd packages/HimalayaUI/frontend && npm test` → Expected: PASS (whole Vitest suite green).

- [ ] **Step 4: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (`lint:design` + `tsc --noEmit -p tsconfig.build.json` + `vite build` all succeed).


## Phase 2 — Chrome primitives (ModalShell + focus-trap fix, Kicker, IconButton, Card/Plate)

ModalShell's first task is the global useFocusTrap fix. Kicker depends on SectionLabel having been deleted in Phase 0. Each wave lowers the guard baseline.

### Task 37: useFocusTrap FOCUSABLE fix (textarea + anchor)

Load-bearing prerequisite for the Notes-drawer migration: the trap's `FOCUSABLE`
selector omits `textarea` and `a[href]`, so a drawer whose only focusable is a
textarea computes an empty `focusable` set and the trap no-ops (early return at
`useFocusTrap.ts:24`). Add both before any drawer migration. No existing consumer
asserts the focusable-set size, so the risk is low — but ship a textarea-only trap
test to lock it.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/hooks/useFocusTrap.ts:3-4`
- Test: `packages/HimalayaUI/frontend/test/useFocusTrap.test.tsx`

- [ ] **Step 1: Write the failing test** — append this case to the existing `describe("useFocusTrap", …)` block in `test/useFocusTrap.test.tsx` (add a second container component above the `describe`):

```tsx
function TextareaTrap({ active = true }: { active?: boolean }): JSX.Element {
  const ref = useRef<HTMLDivElement>(null);
  useFocusTrap(ref, active);
  return (
    <div ref={ref}>
      <textarea data-testid="ta-only" defaultValue="" />
    </div>
  );
}
```

```tsx
  it("traps a textarea-only container (Tab on the sole textarea stays put)", () => {
    const { getByTestId } = render(<TextareaTrap />);
    const ta = getByTestId("ta-only");
    ta.focus();
    expect(document.activeElement).toBe(ta);
    // Sole focusable is both first and last → Tab wraps to itself, focus held.
    const evt = fireEvent.keyDown(ta, { key: "Tab", bubbles: true });
    expect(evt).toBe(false); // defaultPrevented: trap acted (textarea is now in FOCUSABLE)
    expect(document.activeElement).toBe(ta);
  });

  it("includes anchors with href in the focusable set", () => {
    function AnchorTrap(): JSX.Element {
      const ref = useRef<HTMLDivElement>(null);
      useFocusTrap(ref, true);
      return (
        <div ref={ref}>
          <a href="#x" data-testid="lnk">link</a>
          <button data-testid="btn-z">Z</button>
        </div>
      );
    }
    const { getByTestId } = render(<AnchorTrap />);
    const lnk = getByTestId("lnk");
    const z = getByTestId("btn-z");
    z.focus();
    fireEvent.keyDown(z, { key: "Tab", bubbles: true }); // last → wraps to first (anchor)
    expect(document.activeElement).toBe(lnk);
  });
```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/useFocusTrap.test.tsx` → Expected: FAIL (current `FOCUSABLE` excludes `textarea` and `a[href]`; the textarea case sees `focusable.length === 0` so the handler early-returns and `fireEvent` returns `true`, and the anchor case never wraps to the link).

- [ ] **Step 3: Implement** — replace the `FOCUSABLE` constant at `src/hooks/useFocusTrap.ts:3-4`:

```ts
const FOCUSABLE =
  'a[href],button:not([disabled]),input:not([disabled]),select:not([disabled]),textarea:not([disabled]),[tabindex]:not([tabindex="-1"])';
```

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/useFocusTrap.test.tsx` → Expected: PASS (all original cases plus the two new ones).

- [ ] **Step 5: Commit** — `cd packages/HimalayaUI/frontend && git add src/hooks/useFocusTrap.ts test/useFocusTrap.test.tsx && git commit -m "fix(useFocusTrap): include textarea + a[href] in FOCUSABLE (unblocks textarea-only Notes drawer trap)"`


### Task 38: ModalShell primitive

Generalize the proven `ConflictModalShell`/`NavModal` chrome into one primitive
that owns scrim + frame + Esc + focus-trap + outside-click. Closed-look: appearance
comes from `variant`/`size`/`align` internally (the `Record<…,string>` map pattern
from `Button.tsx`); consumer `className` is placement-only and appended LAST to the
frame. Returns `null` when `!open` (matches every site's `if (!open) return null`).
Crucially it does **not** impose an initial-focus target — children own focus
(NavModal synchronously focuses its input on open; stealing focus to the frame would
regress that and the jsdom-synchronous focus test).

Canonical look (from the inventory): scrim `bg-scrim backdrop-blur-sm anim-pal-in`;
frame `bg-plate border border-hair-strong rounded-md shadow-2xl anim-pal-scale`
(radius 5px = `rounded-md` per the locked card/plate decision, NOT the shipped
`rounded-xl` 12px). `bg-scrim` is the Phase-0 `--color-scrim` token
(`oklch(0.05 0 0 / 0.65)`).

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ui/ModalShell.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/ui/index.ts`
- Test: `packages/HimalayaUI/frontend/test/ui/ModalShell.test.tsx`

- [ ] **Step 1: Write the failing test** — create `test/ui/ModalShell.test.tsx`:

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ModalShell } from "../../src/components/ui/ModalShell";

describe("ModalShell", () => {
  it("renders nothing when open=false", () => {
    render(
      <ModalShell open={false} onClose={() => {}} aria-label="X" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(screen.queryByTestId("m")).toBeNull();
  });

  it("renders a role=dialog frame with aria-modal and the passed aria-label when open", () => {
    render(
      <ModalShell open onClose={() => {}} aria-label="Build thing" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    const dialog = screen.getByRole("dialog");
    expect(dialog).toHaveAttribute("aria-modal", "true");
    expect(dialog).toHaveAttribute("aria-label", "Build thing");
    expect(dialog).toHaveAttribute("data-testid", "m");
  });

  it("forwards aria-labelledby / aria-describedby to the frame", () => {
    render(
      <ModalShell open onClose={() => {}} aria-labelledby="t" aria-describedby="s" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    const dialog = screen.getByRole("dialog");
    expect(dialog).toHaveAttribute("aria-labelledby", "t");
    expect(dialog).toHaveAttribute("aria-describedby", "s");
  });

  it("role=dialog is on the frame, NOT on the scrim (regression guard for the scrim-role bug)", () => {
    render(
      <ModalShell open onClose={() => {}} aria-label="X" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    const scrim = screen.getByTestId("m-scrim");
    expect(scrim).not.toHaveAttribute("role", "dialog");
    expect(scrim).toHaveAttribute("role", "presentation");
  });

  it("exposes data-variant / data-size / data-align on the frame", () => {
    render(
      <ModalShell open onClose={() => {}} size="lg" align="top" aria-label="X" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    const dialog = screen.getByRole("dialog");
    expect(dialog).toHaveAttribute("data-variant", "dialog");
    expect(dialog).toHaveAttribute("data-size", "lg");
    expect(dialog).toHaveAttribute("data-align", "top");
  });

  it("Esc fires onClose by default; closeOnEsc=false suppresses it", () => {
    const onClose = vi.fn();
    const { rerender } = render(
      <ModalShell open onClose={onClose} aria-label="X" testId="m"><button>hi</button></ModalShell>,
    );
    fireEvent.keyDown(document, { key: "Escape" });
    expect(onClose).toHaveBeenCalledTimes(1);

    onClose.mockClear();
    rerender(
      <ModalShell open onClose={onClose} closeOnEsc={false} aria-label="X" testId="m"><button>hi</button></ModalShell>,
    );
    fireEvent.keyDown(document, { key: "Escape" });
    expect(onClose).not.toHaveBeenCalled();
  });

  it("clicking the scrim fires onClose; clicking inside the frame does not", () => {
    const onClose = vi.fn();
    render(
      <ModalShell open onClose={onClose} aria-label="X" testId="m">
        <button data-testid="inner">hi</button>
      </ModalShell>,
    );
    fireEvent.click(screen.getByTestId("inner"));
    expect(onClose).not.toHaveBeenCalled();
    fireEvent.click(screen.getByTestId("m-scrim"));
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("closeOnOutsideClick=false suppresses scrim-click dismiss", () => {
    const onClose = vi.fn();
    render(
      <ModalShell open onClose={onClose} closeOnOutsideClick={false} aria-label="X" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    fireEvent.click(screen.getByTestId("m-scrim"));
    expect(onClose).not.toHaveBeenCalled();
  });

  it("does NOT steal initial focus to the frame (children own focus)", () => {
    const trigger = document.createElement("button");
    document.body.appendChild(trigger);
    trigger.focus();
    render(
      <ModalShell open onClose={() => {}} aria-label="X" testId="m"><button>hi</button></ModalShell>,
    );
    // ModalShell imposes no autofocus → focus stays where the caller left it.
    expect(document.activeElement).toBe(trigger);
    document.body.removeChild(trigger);
  });

  it("appends placement className to the frame", () => {
    render(
      <ModalShell open onClose={() => {}} aria-label="X" testId="m" className="max-h-[80vh]">
        <button>hi</button>
      </ModalShell>,
    );
    expect(screen.getByRole("dialog").className).toContain("max-h-[80vh]");
  });

  it("drawer variant: data-variant=drawer on the frame, scrim still present + dismissable", () => {
    const onClose = vi.fn();
    render(
      <ModalShell open onClose={onClose} variant="drawer" aria-label="Notes" testId="d">
        <textarea />
      </ModalShell>,
    );
    expect(screen.getByRole("dialog")).toHaveAttribute("data-variant", "drawer");
    fireEvent.click(screen.getByTestId("d-scrim"));
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("traps focus inside a textarea-only frame (the Notes-drawer case)", () => {
    render(
      <ModalShell open onClose={() => {}} variant="drawer" aria-label="Notes" testId="d">
        <textarea data-testid="ta" />
      </ModalShell>,
    );
    const ta = screen.getByTestId("ta");
    ta.focus();
    const evt = fireEvent.keyDown(ta, { key: "Tab", bubbles: true });
    expect(evt).toBe(false); // trap acted (textarea is the only focusable, wraps to itself)
    expect(document.activeElement).toBe(ta);
  });
});
```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/ModalShell.test.tsx` → Expected: FAIL (module `src/components/ui/ModalShell` does not exist).

- [ ] **Step 3: Implement** — create `src/components/ui/ModalShell.tsx`:

```tsx
import { useEffect, useRef } from "react";
import type { ReactNode } from "react";
import { useFocusTrap } from "../../hooks/useFocusTrap";

export type ModalSize = "sm" | "md" | "lg";
export type ModalAlign = "center" | "top";
export type ModalVariant = "dialog" | "drawer";

export interface ModalShellProps {
  /** Render gate. Returns null when false (matches every call site's `if (!open) return null`). */
  open: boolean;
  onClose: () => void;
  size?: ModalSize;            // default "md"; ignored for variant="drawer"
  align?: ModalAlign;          // default "center"; ignored for variant="drawer"
  /** Esc dismiss. Default true. Set false where the parent owns Escape (e.g. NavModal). */
  closeOnEsc?: boolean;
  /** Backdrop-click dismiss. Default true. */
  closeOnOutsideClick?: boolean;
  /** "drawer" → right-edge sheet with a lower-z scrim. Default "dialog". */
  variant?: ModalVariant;
  /** a11y: the caller must supply at least one of these. */
  "aria-label"?: string;
  "aria-labelledby"?: string;
  "aria-describedby"?: string;
  /** Forwarded to the FRAME (role=dialog). Scrim gets `${testId}-scrim`. */
  testId?: string;
  /** PLACEMENT-ONLY on the frame (max-height, flex/grid, one-off width). No appearance utilities. */
  className?: string;
  children: ReactNode;
}

const sizeClass: Record<ModalSize, string> = {
  sm: "w-[min(440px,calc(100vw-48px))]",
  md: "w-[min(640px,calc(100vw-48px))]",
  lg: "w-[min(820px,calc(100vw-48px))]",
};

const alignClass: Record<ModalAlign, string> = {
  center: "items-center justify-center",
  top: "items-start justify-center pt-[12vh]",
};

/** Tiny placement-only class joiner (brief-sanctioned; no cva/clsx/tailwind-merge). */
function cx(...parts: (string | false | undefined)[]): string {
  return parts.filter(Boolean).join(" ");
}

export function ModalShell({
  open,
  onClose,
  size = "md",
  align = "center",
  closeOnEsc = true,
  closeOnOutsideClick = true,
  variant = "dialog",
  testId,
  className = "",
  children,
  ...aria
}: ModalShellProps): JSX.Element | null {
  const frameRef = useRef<HTMLDivElement>(null);
  useFocusTrap(frameRef, open);

  useEffect(() => {
    if (!open || !closeOnEsc) return;
    const onKey = (e: KeyboardEvent): void => {
      if (e.key === "Escape") {
        e.preventDefault();
        onClose();
      }
    };
    document.addEventListener("keydown", onKey);
    return () => document.removeEventListener("keydown", onKey);
  }, [open, closeOnEsc, onClose]);

  if (!open) return null;

  const onScrimClick = (e: React.MouseEvent): void => {
    if (closeOnOutsideClick && e.target === e.currentTarget) onClose();
  };

  const ariaProps = {
    "aria-label": aria["aria-label"],
    "aria-labelledby": aria["aria-labelledby"],
    "aria-describedby": aria["aria-describedby"],
  };

  if (variant === "drawer") {
    return (
      <div data-testid={testId} className="contents">
        <div
          data-testid={testId ? `${testId}-scrim` : undefined}
          role="presentation"
          onClick={onScrimClick}
          className="fixed inset-0 z-40 bg-scrim anim-pal-in"
        />
        <div
          ref={frameRef}
          role="dialog"
          aria-modal="true"
          data-variant="drawer"
          {...ariaProps}
          className={cx(
            "fixed right-0 top-14 bottom-0 z-50 w-[300px] max-w-[85vw] overflow-y-auto",
            "bg-plate border-l border-hair-strong shadow-2xl",
            className,
          )}
        >
          {children}
        </div>
      </div>
    );
  }

  return (
    <div
      data-testid={testId ? `${testId}-scrim` : undefined}
      role="presentation"
      onClick={onScrimClick}
      className={cx(
        "fixed inset-0 z-50 flex bg-scrim backdrop-blur-sm anim-pal-in",
        alignClass[align],
      )}
    >
      <div
        ref={frameRef}
        data-testid={testId}
        role="dialog"
        aria-modal="true"
        data-variant="dialog"
        data-size={size}
        data-align={align}
        {...ariaProps}
        className={cx(
          sizeClass[size],
          "bg-plate border border-hair-strong rounded-md shadow-2xl anim-pal-scale",
          "flex flex-col overflow-hidden",
          className,
        )}
      >
        {children}
      </div>
    </div>
  );
}
```

NOTE — `bg-scrim` depends on the Phase-0 `--color-scrim` token
(`oklch(0.05 0 0 / 0.65)`). If that token has not yet landed in `styles.css` when
this task runs, the `lint:design`/build will not fail (Tailwind emits the class),
but the scrim renders transparent. Confirm the token concern is merged first, or
temporarily inline `bg-[oklch(0.05_0_0/0.65)]` and leave a `// TODO bg-scrim once
--color-scrim lands` comment. Tests assert nodes/roles, not the class, so they pass
either way.

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/ModalShell.test.tsx` → Expected: PASS (all 12 cases).

- [ ] **Step 5: Export from the ui barrel** — add to `src/components/ui/index.ts`:

```ts
export { ModalShell } from "./ModalShell";
export type { ModalSize, ModalAlign, ModalVariant, ModalShellProps } from "./ModalShell";
```

- [ ] **Step 6: Type-gate + commit** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (tsc clean). Then `git add src/components/ui/ModalShell.tsx src/components/ui/index.ts test/ui/ModalShell.test.tsx && git commit -m "feat(ui): ModalShell primitive (scrim+frame+Esc+trap+outside-click, dialog|drawer)"`


### Task 39: Migrate ConflictModalShell to compose ModalShell

Keep `ConflictModalShell` as a thin conflict-specific composition over `ModalShell`
(preserves its test surface and its sole consumer `SeriesCommitConflictModal`). Drop
the hand-rolled scrim + frame + Esc effect + `useFocusTrap`; keep the header/body/
footer JSX as children. Preserves `data-testid="conflict-modal"`,
`role=dialog`, `aria-labelledby="conflict-title"`, `aria-describedby="conflict-subtitle"`.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ConflictModalShell.tsx:18-128`
- Test (must stay green): `packages/HimalayaUI/frontend/test/ConflictModalShell.test.tsx`, `packages/HimalayaUI/frontend/test/SeriesCommitConflictModal.test.tsx`

- [ ] **Step 1: Run the existing tests first (baseline green)** — `cd packages/HimalayaUI/frontend && npx vitest run test/ConflictModalShell.test.tsx test/SeriesCommitConflictModal.test.tsx` → Expected: PASS (capture as the regression baseline; these tests are NOT modified — they assert role/aria/testid/text on the frame).

- [ ] **Step 2: Replace imports + the outer scrim/frame/trap/Esc machinery** — in `ConflictModalShell.tsx`:

Change the imports at the top (`:18-20`):
```tsx
import type { ReactNode } from "react";
import { ModalShell } from "./ui";
```
(Remove `useEffect`, `useRef`, and the `useFocusTrap` import — the shell owns them.)

Replace the component body from the `dialogRef`/`useFocusTrap`/`useEffect`/`if (!open) return null` block (`:52-87`) and the matching closing `</div></div>` of the scrim+frame (`:125-126`) so the function reads:

```tsx
export function ConflictModalShell({
  open, heading, subtitle, serverPanel, localPanel,
  onClose, onDiscard, discardLabel, onOverwrite, overwriteBusy, extraAction,
}: ConflictModalShellProps): JSX.Element | null {
  return (
    <ModalShell
      open={open}
      onClose={onClose}
      size="lg"
      testId="conflict-modal"
      aria-labelledby="conflict-title"
      aria-describedby="conflict-subtitle"
      className="max-h-[80vh]"
    >
      <header className="px-5 py-4 border-b border-hair-strong">
        <h2 id="conflict-title" className="text-ink text-lg font-medium">
          {heading}
        </h2>
        <p id="conflict-subtitle" className="text-ink-soft text-sm mt-1">
          {subtitle}
        </p>
      </header>

      <div className="flex-1 min-h-0 overflow-y-auto grid grid-cols-2 gap-3 p-5">
        <Panel {...serverPanel} />
        <Panel {...localPanel} />
      </div>

      <footer className="flex items-center gap-2 px-5 py-3 border-t border-hair-strong">
        <button
          type="button"
          data-testid="conflict-discard"
          onClick={onDiscard}
          className="px-3 py-1.5 rounded border border-hair-strong text-ink text-sm
                     hover:bg-paper-sunk"
        >
          {discardLabel}
        </button>
        {extraAction}
        <span className="flex-1" />
        <button
          type="button"
          data-testid="conflict-overwrite"
          onClick={onOverwrite}
          disabled={overwriteBusy}
          className="px-3 py-1.5 rounded border border-accent bg-accent
                     text-paper text-sm disabled:opacity-60"
        >
          {overwriteBusy ? "Saving…" : "Overwrite with mine"}
        </button>
      </footer>
    </ModalShell>
  );
}
```

Leave the `Panel` helper (`:130-162`) and the `ConflictPanelData`/`ConflictModalShellProps` interfaces unchanged.

- [ ] **Step 3: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ConflictModalShell.test.tsx test/SeriesCommitConflictModal.test.tsx` → Expected: PASS (frame still `role=dialog`, `data-testid="conflict-modal"`, `aria-labelledby`/`aria-describedby` preserved; Esc now via document listener inside ModalShell; outside-click preserved; the parent `overwriteInFlightRef` double-click guard is untouched — ModalShell does not touch the footer handlers).

- [ ] **Step 4: Commit** — `cd packages/HimalayaUI/frontend && git add src/components/ConflictModalShell.tsx && git commit -m "refactor(ConflictModalShell): compose ModalShell (drop hand-rolled scrim/frame/Esc/trap)"`


### Task 40: Migrate NavModal to ModalShell

NavModal owns Escape inside its input keydown (`:166-199`, with `preventDefault` +
Backspace-pops-chips + arrow nav + Enter/Tab commit), so it passes
`closeOnEsc={false}`. It synchronously focuses its input on open (`:88`) — ModalShell
must not steal that, which is why ModalShell imposes no initial focus. Uses
`align="top"` (command-palette placement) and `size="md"`. Outside-click + trap still
come from the shell.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/NavModal.tsx:1-7,77-78,117,213-230,314-316`
- Test (must stay green): `packages/HimalayaUI/frontend/test/NavModal.test.tsx`

- [ ] **Step 1: Run the existing test first (baseline green)** — `cd packages/HimalayaUI/frontend && npx vitest run test/NavModal.test.tsx` → Expected: PASS (baseline; tests assert `nav-modal`/`nav-modal-input`/`nav-modal-results` testids, Esc-via-input, Backspace chip-pop, arrow nav — none modified).

- [ ] **Step 2: Replace imports + remove the hand-rolled trap** — at `NavModal.tsx`:
  - Remove `useRef` of `dialogRef` (`:77`) and the `useFocusTrap(dialogRef, open)` line (`:78`). Keep `inputRef` (`:76`).
  - Remove the `import { useFocusTrap } from "../hooks/useFocusTrap";` line (`:7`).
  - Add `import { ModalShell } from "./ui";` near the other imports.
  - Remove the standalone `if (!open) return null;` at `:117` (ModalShell now gates on `open`). Verify the line before editing.

- [ ] **Step 3: Replace the scrim+frame wrapper** — replace the outer scrim `<div data-testid="nav-modal" …>` and the inner frame `<div ref={dialogRef} role="dialog" …>` (`:213-230`) with a single `<ModalShell>` opener, and replace the two matching closing `</div>` tags at `:314-315` with `</ModalShell>`:

Opening (replaces `:213-230`):
```tsx
  return (
    <ModalShell
      open={open}
      onClose={closeModal}
      size="md"
      align="top"
      closeOnEsc={false}
      testId="nav-modal"
      aria-label="Navigate to experiment or sample"
      className="max-h-[72vh]"
    >
```

Closing (replaces the frame `</div>` and scrim `</div>` at `:314-315`):
```tsx
    </ModalShell>
  );
}
```

The chip/input header `<div>`, the `<Skeleton>` results region, and the footer
keyboard-hint row stay as `ModalShell` children unchanged.

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/NavModal.test.tsx` → Expected: PASS (input still synchronously focused on open since ModalShell imposes no autofocus; Esc still handled by the input keydown because `closeOnEsc={false}`; `nav-modal` testid now on the frame; `nav-modal` was previously on the SCRIM div — verify the NavModal test queries that survive: `data-testid="nav-modal"` null when closed is still satisfied because ModalShell returns null, and the input/results testids are unchanged).

- [ ] **Step 5: Commit** — `cd packages/HimalayaUI/frontend && git add src/components/NavModal.tsx && git commit -m "refactor(NavModal): adopt ModalShell (align=top, closeOnEsc=false, input keeps Esc/Backspace/focus)"`


### Task 41: Migrate OnboardingFlow to ModalShell (collapse double-dialog)

Today: a presentation scrim wraps EITHER `NameStep` OR `TutorialStep`, and each step
renders its OWN `role=dialog` frame + its OWN `useFocusTrap` (a double-dialog). After
migration: ONE `ModalShell` frame whose children swap between the two step bodies.
Non-dismissable → `closeOnOutsideClick={false}`; Esc/arrows live in the step keydown
handlers → `closeOnEsc={false}`. Verify NameStep's Enter-submit still reaches its
handler via bubbling once the per-step `role=dialog` is gone (it attaches `onKeyDown`
to a content wrapper that contains the input, so bubbling carries Enter up).

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/OnboardingFlow.tsx:1-6,144-176,200-201,220-229,313-314,317-327`
- Test (must stay green): any `test/*Onboarding*` spec if present (search first)

- [ ] **Step 1: Find + run the existing Onboarding tests** — `cd packages/HimalayaUI/frontend && ls test | grep -i onboard; npx vitest run $(ls test | grep -i onboard | sed 's#^#test/#' | tr '\n' ' ')` → Expected: PASS if a spec exists (baseline). If none exists, note "no Onboarding unit test — rely on `npm run build` + e2e" and proceed. (The e2e `onboarding-overlay`/`onboarding-name`/`onboarding-tutorial`/`onboarding-continue`/`tutorial-next` testids must all survive — preserve them exactly.)

- [ ] **Step 2: Wrap the shared scrim in ModalShell** — at `OnboardingFlow.tsx`:
  - Add `import { ModalShell } from "./ui";` (keep the existing `import { Button } from "./ui";` — or merge into one line `import { Button, ModalShell } from "./ui";`).
  - Replace the outer presentation scrim `<div data-testid="onboarding-overlay" … role="presentation">` (`:144-150`) and its closing `</div>` (`:175`) with:

Opening (replaces `:144-150`):
```tsx
  return (
    <ModalShell
      open
      onClose={closeTutorial}
      size="sm"
      closeOnEsc={false}
      closeOnOutsideClick={false}
      testId="onboarding-overlay"
      aria-label="Onboarding"
    >
```
Closing (replaces `:175`):
```tsx
    </ModalShell>
  );
```
The `{phase === "name" && <NameStep …/>}` and `{phase === "tutorial" && <TutorialStep …/>}` blocks stay as children.

- [ ] **Step 3: Strip NameStep's frame chrome + per-step trap** — in `NameStep` (`:200-229`):
  - Remove `const dialogRef = useRef<HTMLDivElement>(null);` (`:200`) and `useFocusTrap(dialogRef, true);` (`:201`).
  - Change the outer `<div ref={dialogRef} data-testid="onboarding-name" role="dialog" aria-modal="true" className="bg-plate border border-hair-strong rounded-lg p-6 min-w-[360px] max-w-[480px] flex flex-col gap-4" onKeyDown={…}>` (`:220-228`) to a plain content wrapper that keeps the testid + the Enter-submit handler but drops `role`/`aria-modal`/`ref` and the surface chrome (ModalShell's `size="sm"` + plate surface supply width + background):

```tsx
    <div
      data-testid="onboarding-name"
      className="p-6 flex flex-col gap-4"
      onKeyDown={(e) => { if (e.key === "Enter") { e.preventDefault(); onSubmit(); } }}
    >
```
(`useRef`/`useFocusTrap` imports stay in the file — TutorialStep still uses neither after this task, so remove the now-unused `useFocusTrap` import only after Step 4; see Step 5.)

- [ ] **Step 4: Strip TutorialStep's frame chrome + per-step trap** — in `TutorialStep` (`:313-327`):
  - Remove `const dialogRef = useRef<HTMLDivElement>(null);` (`:313`) and `useFocusTrap(dialogRef, true);` (`:314`).
  - Change the outer `<div ref={dialogRef} data-testid="onboarding-tutorial" role="dialog" aria-modal="true" tabIndex={-1} onKeyDown={onKeyDown} className="bg-plate border border-hair-strong rounded-lg p-7 min-w-[420px] max-w-[520px] flex flex-col gap-4 outline-0">` (`:317-327`) to a focus-receiving content wrapper (keep `tabIndex={-1}` so the arrow/Esc `onKeyDown` still receives keys; drop `role`/`aria-modal`/`ref`/surface chrome):

```tsx
    <div
      data-testid="onboarding-tutorial"
      tabIndex={-1}
      onKeyDown={onKeyDown}
      className="p-7 flex flex-col gap-4 outline-0"
    >
```

- [ ] **Step 5: Remove the now-unused useFocusTrap import** — both step traps are gone; delete `import { useFocusTrap } from "../hooks/useFocusTrap";` (`:6`) and remove `useRef` from the React import (`:1`) if no other `useRef` remains in the file (grep `useRef` first; `OnboardingFlow` uses none after this change). Verify with the build in Step 6.

- [ ] **Step 6: Type-gate + run** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (no unused-import / tsc errors). Then run any Onboarding unit spec from Step 1 → Expected: PASS. RISK to verify in the run/build: only ONE `role=dialog` now exists (was two nested); Enter-submit still fires (the `onKeyDown` is on the wrapper containing the input → bubbling carries it); arrow/Esc still fire (TutorialStep wrapper keeps `tabIndex=-1`).

- [ ] **Step 7: Commit** — `cd packages/HimalayaUI/frontend && git add src/components/OnboardingFlow.tsx && git commit -m "refactor(OnboardingFlow): single ModalShell frame (collapse double-dialog, closeOnEsc/OutsideClick=false)"`


### Task 42: Migrate SpeculativeBuilder to ModalShell (fixes role=dialog-on-scrim a11y bug)

The current scrim div carries `role=dialog aria-modal aria-label` (`:114-119`) — an
a11y bug (AT treats the backdrop as the dialog). ModalShell puts `role=dialog` on the
frame. Drops the `window` Esc keydown effect (`:100-104`), the manual
`useFocusTrap` (`:38-39`), `bg-black/40`→scrim token, and `bg-paper`→plate. The
builder is mount-gated by the parent (`PhasePanel:402`), so pass `open` (always true
while mounted). The success-close path (`closedOnSuccessRef`) and the
`variant="primary"` Save button are untouched by THIS task (the Save-button variant
rename is the Button unit's job — see cross-unit notes).

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/SpeculativeBuilder.tsx:1-7,38-39,99-104,113-126,279-281`
- Test (must stay green): `packages/HimalayaUI/frontend/test/SpeculativeBuilder.test.tsx`

- [ ] **Step 1: Run the existing test first (baseline green)** — `cd packages/HimalayaUI/frontend && npx vitest run test/SpeculativeBuilder.test.tsx` → Expected: PASS (baseline). Verify the spec queries the FRAME (`speculative-builder` testid / `getByRole("dialog")`) and not the scrim — it queries `data-testid="speculative-builder"` (the frame) and `getByRole("button", {name:/cancel/i})`, so the role-move is safe.

- [ ] **Step 2: Replace imports + remove the manual trap/Esc** — at `SpeculativeBuilder.tsx`:
  - Add `ModalShell` to the ui import: change `import { Button } from "./ui";` (`:7`) to `import { Button, ModalShell } from "./ui";`.
  - Remove `import { useFocusTrap } from "../hooks/useFocusTrap";` (`:4`).
  - Remove `const dialogRef = useRef<HTMLDivElement>(null);` and `useFocusTrap(dialogRef, true);` (`:38-39`). Keep the other `useRef` calls (`lastGoodSnapRef`, `closedOnSuccessRef`) — so leave `useRef` in the React import.
  - Remove the Esc `useEffect` block (`:99-104`).

- [ ] **Step 3: Replace the scrim+frame wrapper** — replace the outer scrim `<div role="dialog" aria-modal … bg-black/40 onClick={onClose}>` and the inner frame `<div ref={dialogRef} data-testid="speculative-builder" … bg-paper … onClick={stopPropagation}>` (`:113-126`) with a single `ModalShell` opener; replace the two closing `</div>` tags (`:279-280`) with `</ModalShell>`:

Opening (replaces `:113-126`):
```tsx
  return (
    <ModalShell
      open
      onClose={onClose}
      size="sm"
      aria-label="Build speculative index"
      testId="speculative-builder"
      className="max-h-[90vh] overflow-y-auto gap-4 p-5"
    >
      <div className="flex items-center justify-between">
        <h2 className="text-lg font-semibold text-ink">Speculative index</h2>
        <button
          onClick={onClose}
          aria-label="Close"
          className="text-ink-faint hover:text-ink text-xl px-2 leading-none"
        >×</button>
      </div>
```
Closing (replaces `:279-280`):
```tsx
    </ModalShell>
  );
}
```
The `<p>` description, phase/anchor/ratio pickers, snap preview, error, and footer
button row stay as children unchanged.

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/SpeculativeBuilder.test.tsx` → Expected: PASS (`speculative-builder` testid now on the frame which IS `role=dialog`; Cancel button still calls `onClose`; save→onClose path unchanged; outside-click + Esc now come from ModalShell — the old `onClick={onClose}` on the scrim is replaced by ModalShell's `closeOnOutsideClick` default true).

- [ ] **Step 5: Commit** — `cd packages/HimalayaUI/frontend && git add src/components/SpeculativeBuilder.tsx && git commit -m "refactor(SpeculativeBuilder): adopt ModalShell (fixes role=dialog-on-scrim a11y bug, drops window Esc + manual trap)"`


### Task 43: Migrate ScopingConfirmModal to ModalShell

Replaces the divergent `bg-black/40` scrim + `bg-paper rounded-lg` frame with the
canonical scrim token + plate surface. The current frame-scoped `onKeyDown` Esc
(`:36`) is a latent bug (only fires when the frame has focus); ModalShell's
document-level Esc fixes it (Esc now works even when a child input has focus) — a
behavior change to note, and the existing `fireEvent.keyDown(dialog, {key:"Escape"})`
test still passes because keydown bubbles to document. The confirm-build button's
class-string assertion + its fold into `<Button variant="solid">` is the BUTTON
unit's task — preserve the button markup and its `scoping-confirm-build` testid here.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ScopingConfirmModal.tsx:1-2,20-37,57-59`
- Test (must stay green): `packages/HimalayaUI/frontend/test/scoping.test.tsx` (the modal cases)

- [ ] **Step 1: Run the existing test first (baseline green)** — `cd packages/HimalayaUI/frontend && npx vitest run test/scoping.test.tsx` → Expected: PASS (baseline; cases at `:85-119` assert `scoping-confirm-modal` role/text, closed→null, Esc→onClose, and the confirm-build class strings).

- [ ] **Step 2: Replace imports + the scrim/frame** — rewrite `ScopingConfirmModal.tsx` so the wrapper is `ModalShell` (keep the `Props`/heading/body/buttons identical):

Replace `:1-2` (imports):
```tsx
import { ModalShell } from "./ui";
```
(Remove the `useRef` + `useFocusTrap` imports — ModalShell owns the trap and the open-gate.)

Replace the body `:20-37` (the `dialogRef`/`useFocusTrap`/`if (!open) return null`/
scrim `<div>`/frame `<div ref=… role=dialog …>` opener) with:
```tsx
}: Props): JSX.Element | null {
  return (
    <ModalShell
      open={open}
      onClose={onClose}
      size="sm"
      testId="scoping-confirm-modal"
      aria-label="Confirm series scoping"
      className="p-6"
    >
```
Replace the two closing `</div>` tags at `:57-58` with:
```tsx
    </ModalShell>
  );
}
```
The `<h2>`, `<p>`, and the Cancel/Confirm `<div>` button row stay as children
unchanged (preserving `scoping-confirm-cancel` and `scoping-confirm-build` testids).

- [ ] **Step 3: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/scoping.test.tsx` → Expected: PASS (frame still `role=dialog`, `scoping-confirm-modal` testid preserved, closed→null via ModalShell, `fireEvent.keyDown(dialog, {Escape})` still calls onClose via document-level listener and bubbling; confirm-build class-string test untouched here).

- [ ] **Step 4: Commit** — `cd packages/HimalayaUI/frontend && git add src/components/ScopingConfirmModal.tsx && git commit -m "refactor(ScopingConfirmModal): adopt ModalShell (scrim+plate tokens, document-level Esc)"`


### Task 44: Migrate the Notes drawer (FocusWorkspaceLayout) to ModalShell variant=drawer

This is the site the `useFocusTrap` textarea fix unblocks: the drawer body
(`FocusNotesMargin`, textarea at `:108`) is its sole/last focusable, so the
pre-fix trap leaked. Use `variant="drawer"` (right-edge sheet + lower-z-40 scrim).
Drops the manual scrim + `useFocusTrap` (`:43-44`) + `window` Esc effect (`:45-50`).
The `${testId}-scrim` convention must match the existing
`focus-notes-drawer-scrim` testid (`:152`), so pass `testId="focus-notes-drawer"`.
DO THIS AFTER the useFocusTrap fix task (hard ordering dependency).

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/FocusWorkspaceLayout.tsx:1-7,43-50,149-182`
- Test (must stay green): `packages/HimalayaUI/frontend/test/FocusWorkspaceLayout.header.test.tsx` (+ search for any drawer-specific spec)

- [ ] **Step 1: Find + run the existing FocusWorkspace tests** — `cd packages/HimalayaUI/frontend && npx vitest run test/FocusWorkspaceLayout.header.test.tsx test/FocusWorkspacePage.layout.test.tsx` → Expected: PASS (baseline). Grep for drawer testids: `grep -rn "focus-notes-drawer" test` → preserve every match (`focus-notes-drawer`, `focus-notes-drawer-scrim`, `focus-notes-drawer-close`).

- [ ] **Step 2: Remove the manual trap + Esc effect** — at `FocusWorkspaceLayout.tsx`:
  - Remove the `const notesDrawerRef = useRef<HTMLDivElement>(null);` + `useFocusTrap(notesDrawerRef, notesDrawerOpen);` (`:43-44`) and the Esc `useEffect` (`:45-50`).
  - Remove `import { useFocusTrap } from "../hooks/useFocusTrap";` (`:7`).
  - Add `import { ModalShell } from "./ui";` near the component imports.
  - Remove `useRef` from the React import (`:1`) only if no other `useRef` remains (grep first).

- [ ] **Step 3: Replace the drawer scrim+frame** — replace the `{notesSample !== undefined && notesDrawerOpen && ( … )}` block's inner scrim `<div data-testid="focus-notes-drawer-scrim" …>` and frame `<div ref={notesDrawerRef} role="dialog" … bg-paper …>` (`:149-182`) with a ModalShell. Keep the outer `xl:hidden` gate so the drawer never shows when the margin column is visible:

```tsx
      {notesSample !== undefined && (
        <div className="xl:hidden">
          <ModalShell
            open={notesDrawerOpen}
            onClose={closeNotesDrawer}
            variant="drawer"
            testId="focus-notes-drawer"
            aria-label="Notes"
          >
            <div className="flex items-center justify-between border-b border-hair px-4 py-2">
              <span className="text-meta uppercase tracking-wider text-ink-faint">Notes</span>
              <button
                type="button"
                data-testid="focus-notes-drawer-close"
                onClick={closeNotesDrawer}
                aria-label="Close notes"
                className="rounded px-1.5 text-base leading-none text-ink-faint hover:text-ink"
              >
                &#215;
              </button>
            </div>
            <FocusNotesMargin
              sample={notesSample}
              onSaveNotes={(notes) => updateSample.mutate({ notes })}
            />
          </ModalShell>
        </div>
      )}
```
NOTE: ModalShell now owns the `open` gate, so the `&& notesDrawerOpen` is dropped
from the outer condition and moved to `open={notesDrawerOpen}` (the outer `xl:hidden`
wrapper renders an empty `<div className="contents">`-equivalent when closed because
ModalShell returns null). The `focus-notes-drawer` testid moves from the wrapper to
the ModalShell frame, and `focus-notes-drawer-scrim` is now emitted by ModalShell as
`${testId}-scrim` — identical id, so existing tests/E2E match.

- [ ] **Step 4: Type-gate + run** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS. Then `npx vitest run test/FocusWorkspaceLayout.header.test.tsx test/FocusWorkspacePage.layout.test.tsx` → Expected: PASS. Verify the drawer trap now holds (the textarea fix makes `FocusNotesMargin`'s textarea a valid trap boundary) and the `focus-notes-drawer-scrim` testid still dismisses.

- [ ] **Step 5: Commit** — `cd packages/HimalayaUI/frontend && git add src/components/FocusWorkspaceLayout.tsx && git commit -m "refactor(Notes drawer): ModalShell variant=drawer (relies on useFocusTrap textarea fix; lower-z-40 scrim preserved)"`


### Task 45: Full-suite + build gate for the ModalShell unit

Final convergence check across every migrated surface + their indirect consumers
(`SeriesCommitConflictModal` via `ConflictModalShell`; `App.tsx` top-level mounts of
`OnboardingFlow`/`NavModal`/`SeriesCommitConflictModal`; `PhasePanel:402` builder
mount-gate; `SeriesScopingPage:347` `ScopingConfirmModal` consumer — none change
prop shape).

**Files:**
- Test: the full Vitest suite + the design lint + tsc + vite build.

- [ ] **Step 1: Run the full unit suite** — `cd packages/HimalayaUI/frontend && npm test` → Expected: PASS (all specs, especially `useFocusTrap`, `ModalShell`, `ConflictModalShell`, `SeriesCommitConflictModal`, `NavModal`, `SpeculativeBuilder`, `scoping`, `FocusWorkspaceLayout.*`). If any fail, the failure is a missed testid/role/aria forward — fix the migration, not the test (the only intentionally-rewritten tests in this unit are the two added `useFocusTrap` cases; the Button-unit owns the scoping/ScopingFoot class-string rewrites).

- [ ] **Step 2: Build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (`lint:design` + `tsc --noEmit -p tsconfig.build.json` + `vite build` all clean; no unused-import errors from the removed `useFocusTrap`/`useRef` imports).

- [ ] **Step 3: Commit (if any fixups were needed)** — `cd packages/HimalayaUI/frontend && git add -A && git commit -m "test(ModalShell): full-suite + build gate green across all 6 migrated overlays"` (skip if the tree is clean).

### Task 46: Add the `.text-kicker` typography role to styles.css

The Kicker is its OWN typography role (NOT `.text-label`): weight **700**, uppercase, `~11px` (we use the existing `--text-sm` = 11.5px scale step to avoid reintroducing `text-[Npx]` drift), with **per-tone tracking** — accent runs wide at `0.13em`, faint tighter at `0.09em`. Geometry lives in `styles.css` so no call site carries appearance utilities. We define a base `.text-kicker` (size/weight/uppercase/line-height) plus two tone classes `.text-kicker-accent` (color `print-accent` + `letter-spacing: 0.13em`) and `.text-kicker-faint` (color `ink-faint` + `letter-spacing: 0.09em`).

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/styles.css:259-263` (add new role inside the `@layer`/`@theme` utilities block, after `.text-label`)

- [ ] **Step 1: Add the kicker role classes** — In `src/styles.css`, immediately after the `.text-label` rule (currently lines 261-263), insert the three rules below. The base sets everything except color + tracking; the tone classes set color + tracking so the two divergent tracking families collapse into named roles, not inline `tracking-[…]`:

```css
  .text-kicker        { font-size: var(--text-sm);  line-height: 1.4;
                        font-weight: 700; text-transform: uppercase; }
  .text-kicker-accent { color: var(--color-print-accent); letter-spacing: 0.13em; }
  .text-kicker-faint  { color: var(--color-ink-faint);    letter-spacing: 0.09em; }
```

- [ ] **Step 2: Document the role in the scale comment** — In the `text-*` role legend comment block above (near line 239, the `text-label` legend line), add a sibling line so the catalog/legend stays truthful:

```
 *   text-kicker       eyebrow / section label — uppercase, 700, tracked, tone-colored (accent|faint)
```

- [ ] **Step 3: Build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (CSS-only addition; `lint:design` + `tsc` + `vite build` all clean — no call sites reference the new class yet so the guard baseline is unaffected).
- [ ] **Step 4: Commit** — `git add packages/HimalayaUI/frontend/src/styles.css && git commit -m "feat(ui): add .text-kicker typography role (700, per-tone tracking accent 0.13em / faint 0.09em)"`

---

### Task 47: Create the Kicker primitive (TDD)

Closed-look / open-placement primitive following the `Button.tsx` `Record<Variant,string>` pattern. Tone selects color + tracking via the `.text-kicker-*` classes; `as` selects the element (`div` default for inline labels, `span`, or `h2`/`h3` for true section eyebrows so headings stay in the a11y tree); `className` is **placement-only**. Emits `data-tone` so Vitest asserts semantics, never Tailwind strings. Spreads `{...props}` so `data-testid` / `aria-hidden` / `aria-*` round-trip.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ui/Kicker.tsx`
- Test: `packages/HimalayaUI/frontend/test/ui/Kicker.test.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/ui/index.ts`

- [ ] **Step 1: Write the failing test** — Create `test/ui/Kicker.test.tsx`:

```tsx
import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Kicker } from "../../src/components/ui/Kicker";

describe("Kicker", () => {
  it("renders its children text", () => {
    render(<Kicker>Integration</Kicker>);
    expect(screen.getByText("Integration")).toBeInTheDocument();
  });

  it("exposes the tone via data-tone (accent)", () => {
    render(<Kicker tone="accent">Folio</Kicker>);
    expect(screen.getByText("Folio")).toHaveAttribute("data-tone", "accent");
  });

  it("defaults tone to faint", () => {
    render(<Kicker>Notes</Kicker>);
    expect(screen.getByText("Notes")).toHaveAttribute("data-tone", "faint");
  });

  it("renders a heading when as='h3'", () => {
    render(<Kicker as="h3" tone="accent">Folio</Kicker>);
    expect(screen.getByRole("heading", { name: "Folio" })).toBeInTheDocument();
  });

  it("renders no heading by default (inline label)", () => {
    render(<Kicker>Sort</Kicker>);
    expect(screen.queryByRole("heading")).not.toBeInTheDocument();
  });

  it("renders a span when as='span'", () => {
    render(<Kicker as="span">Offset</Kicker>);
    expect(screen.getByText("Offset").tagName).toBe("SPAN");
  });

  it("forwards data-testid", () => {
    render(<Kicker data-testid="focus-plot-kicker">Integration</Kicker>);
    expect(screen.getByTestId("focus-plot-kicker")).toBeInTheDocument();
  });

  it("forwards aria-hidden", () => {
    render(<Kicker aria-hidden="true" data-testid="heatmap-axis-title">x →</Kicker>);
    expect(screen.getByTestId("heatmap-axis-title")).toHaveAttribute("aria-hidden", "true");
  });

  it("accepts placement-only className (margin) and applies it", () => {
    render(<Kicker className="mb-2">samples screened</Kicker>);
    expect(screen.getByText("samples screened")).toHaveClass("mb-2");
  });
});
```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/Kicker.test.tsx` → Expected: FAIL (`src/components/ui/Kicker.tsx` does not exist — import error).

- [ ] **Step 3: Implement** — Create `src/components/ui/Kicker.tsx`:

```tsx
import type { HTMLAttributes } from "react";

export type KickerTone = "accent" | "faint";

interface KickerProps extends HTMLAttributes<HTMLElement> {
  /** terracotta page/section eyebrow vs. dim metric/field label. Default: "faint". */
  tone?: KickerTone;
  /**
   * Element to render. Section headings that introduce a region should be a
   * heading ("h2"/"h3") for the a11y tree; inline/aside labels stay "div"/"span".
   * Default: "div" (matches the prevailing inline usage).
   */
  as?: "div" | "span" | "h2" | "h3";
}

// tone selects color + per-tone tracking; base geometry (700/uppercase/size) is .text-kicker.
const toneClass: Record<KickerTone, string> = {
  accent: "text-kicker-accent",
  faint: "text-kicker-faint",
};

export function Kicker({
  tone = "faint",
  as: Tag = "div",
  className = "",
  children,
  ...props
}: KickerProps): JSX.Element {
  return (
    <Tag
      data-tone={tone}
      className={`text-kicker ${toneClass[tone]} ${className}`}
      {...props}
    >
      {children}
    </Tag>
  );
}
```

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/Kicker.test.tsx` → Expected: PASS (9 tests).

- [ ] **Step 5: Export from the barrel** — In `src/components/ui/index.ts`, add after the `HintText` export line:

```ts
export { Kicker } from "./Kicker";
export type { KickerTone } from "./Kicker";
```

(Note: the `SectionLabel` export removal is owned by the Phase 0 adopt-or-delete task — see Dependencies. Do not touch the `SectionLabel` line here.)

- [ ] **Step 6: Build + commit** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS. Then `git add packages/HimalayaUI/frontend/src/components/ui/Kicker.tsx packages/HimalayaUI/frontend/test/ui/Kicker.test.tsx packages/HimalayaUI/frontend/src/components/ui/index.ts && git commit -m "feat(ui): add Kicker primitive (tone accent|faint, as div|span|h2|h3, data-tone)"`

---

### Task 48: Migrate accent-tone page/section eyebrows to Kicker

These are the terracotta `text-print-accent` eyebrows. Transformation pattern (shown once):

```tsx
// BEFORE — inline appearance utilities (banned by the guard)
<div className="text-xs font-bold uppercase tracking-[0.13em] text-print-accent">Integration</div>
// AFTER — appearance flows through tone; only placement classes (if any) stay in className
<Kicker tone="accent">Integration</Kicker>
```

Drop at every site: `text-xs`/`text-sm`/`text-[Npx]`, `font-bold`/`font-semibold`, `uppercase`, `tracking-[…]`/`tracking-wide(r)(st)`, `text-print-accent`. Keep ONLY placement utilities (`mb-*`, `mt-*`, `ml-auto`, `block`). Add `import { Kicker } from "@/components/ui"` (or the file's existing relative `ui` barrel path — match sibling imports in each file) where not already present. Preserve element semantics: convert `<span>` sources with `as="span"`.

**Files (verify each line before editing — inventory `loc` may have drifted):**
- Modify: `src/components/FocusPlotHeader.tsx:44`
- Modify: `src/pages/SeriesFolioPage.tsx:71`
- Modify: `src/pages/SamplesPage.tsx:94`
- Modify: `src/components/SeriesFolioCard.tsx:129`
- Modify: `src/pages/SeriesScopingPage.tsx:238`
- Modify: `src/pages/SeriesBuilderPage.tsx:293`

- [ ] **`FocusPlotHeader.tsx:44`** — `<div data-testid="focus-plot-kicker" className="text-xs font-bold uppercase tracking-[0.13em] text-print-accent">Integration</div>` → `<Kicker tone="accent" data-testid="focus-plot-kicker">Integration</Kicker>`. **PRESERVE `data-testid="focus-plot-kicker"`** (asserted in tests/e2e).
- [ ] **`SeriesFolioPage.tsx:71`** — `<div className="text-[11px] font-bold uppercase tracking-[0.14em] text-print-accent">Folio</div>` → `<Kicker tone="accent">Folio</Kicker>` (page eyebrow; 11px → text-sm 11.5px; tracking 0.14em → role 0.13em).
- [ ] **`SamplesPage.tsx:94`** — `<div className="text-xs font-semibold uppercase tracking-[0.14em] text-print-accent">Contact sheet</div>` → `<Kicker tone="accent">Contact sheet</Kicker>` (font-semibold drift dropped).
- [ ] **`SeriesFolioCard.tsx:129`** — element with `className="text-[10.5px] font-bold uppercase tracking-[0.11em] text-print-accent"` → `<Kicker tone="accent">…</Kicker>` (card eyebrow; preserve the original children expression). If the source element is a `<span>`, use `as="span"`.
- [ ] **`SeriesScopingPage.tsx:238`** — `<div className="mb-1.5 text-[11px] font-bold uppercase tracking-[0.14em] text-print-accent">…</div>` → `<Kicker tone="accent" className="mb-1.5">…</Kicker>` (keep `mb-1.5` placement).
- [ ] **`SeriesBuilderPage.tsx:293`** — `<span className="text-[11px] font-bold uppercase tracking-[0.14em] text-print-accent">…</span>` → `<Kicker as="span" tone="accent">…</Kicker>`.

- [ ] **Run affected unit tests** — `cd packages/HimalayaUI/frontend && npx vitest run test/` (or scope to the touched components' specs if a full run is slow) → Expected: PASS, including any spec that does `getByTestId("focus-plot-kicker")`.
- [ ] **Build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS. The `lint:design` baseline shrinks (six fewer banned-utility occurrences) — this is allowed (baseline is a ceiling; not-in-baseline is what blocks).
- [ ] **Commit** — `git add -A && git commit -m "refactor(ui): migrate accent page/section eyebrows to Kicker (6 sites)"`

---

### Task 49: Migrate faint-tone metric/field/aside labels to Kicker

These are the dim `text-ink-faint` labels (metric sublabels, scoping field labels, aside labels, section labels). Transformation pattern (shown once):

```tsx
// BEFORE
<div className="mt-0.5 text-[10.5px] font-bold uppercase tracking-[0.08em] text-ink-faint">samples screened</div>
// AFTER — placement (mt-0.5) stays; all appearance drops; tone=faint is the default but pass it explicitly for grep-ability
<Kicker tone="faint" className="mt-0.5">samples screened</Kicker>
```

Several sources set NO font-weight (effectively 400) and now become 700 — this is the deliberate role unification, accepted per the brief. Sources spanning 10px–11.5px collapse to text-sm (11.5px). Convert `<span>` sources with `as="span"`. Verify each line before editing.

**Files (verify lines before editing):**
- Modify: `src/pages/SeriesFolioPage.tsx:87,130`
- Modify: `src/pages/SamplesPage.tsx:128`
- Modify: `src/components/PhasePanel.tsx:185,318`
- Modify: `src/components/LoupeSidebar.tsx:61`
- Modify: `src/components/ScopingLooseMatches.tsx:31`
- Modify: `src/components/ScopingOrderField.tsx:16`
- Modify: `src/components/ScopingPhaseStrip.tsx:26`
- Modify: `src/components/OffsetDock.tsx:22`
- Modify: `src/components/OnboardingFlow.tsx:328`
- Modify: `src/components/SeriesBuilderRail.tsx:87`
- Modify: `src/pages/SeriesScopingPage.tsx:293`
- Modify: `src/pages/SeriesBuilderPage.tsx:318`
- Modify: `src/components/CorpusTopbar.tsx:270`
- Modify: `src/components/FocusNotesMargin.tsx:80`
- Modify: `src/components/FocusWorkspaceLayout.tsx:165`

- [ ] **`SeriesFolioPage.tsx:87`** — `<div className="mt-1 text-[10.5px] font-bold uppercase tracking-[0.08em] text-ink-faint">{isFiltered ? \`of ${series.length} shown\` : "series in the folio"}</div>` → `<Kicker tone="faint" className="mt-1">{isFiltered ? \`of ${series.length} shown\` : "series in the folio"}</Kicker>` (keep `mt-1`).
- [ ] **`SeriesFolioPage.tsx:130`** — `<span className="ml-auto text-[10.5px] font-bold uppercase tracking-[0.07em] text-ink-faint">Sort</span>` → `<Kicker as="span" tone="faint" className="ml-auto">Sort</Kicker>` (keep `ml-auto`).
- [ ] **`SamplesPage.tsx:128`** — `<div className="mt-0.5 text-[10.5px] font-bold uppercase tracking-[0.08em] text-ink-faint">samples screened</div>` → `<Kicker tone="faint" className="mt-0.5">samples screened</Kicker>`.
- [ ] **`PhasePanel.tsx:185`** — `<div className="text-[10.5px] font-bold uppercase tracking-[0.09em] text-ink-faint">…</div>` → `<Kicker tone="faint">…</Kicker>` (preserve children).
- [ ] **`PhasePanel.tsx:318`** — element with `className="px-4 pt-2.5 text-[10px] font-bold uppercase tracking-[0.08em] text-ink-faint"` → `<Kicker tone="faint" className="px-4 pt-2.5">…</Kicker>` (keep padding placement).
- [ ] **`LoupeSidebar.tsx:61`** — `<div className="mb-2 text-[10px] font-bold uppercase tracking-wide text-ink-faint">…</div>` → `<Kicker tone="faint" className="mb-2">…</Kicker>`.
- [ ] **`ScopingLooseMatches.tsx:31`** — `<div className="mb-2 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">…</div>` → `<Kicker tone="faint" className="mb-2">…</Kicker>`.
- [ ] **`ScopingOrderField.tsx:16`** — `<div className="mb-1.5 mt-5 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">…</div>` → `<Kicker tone="faint" className="mb-1.5 mt-5">…</Kicker>`.
- [ ] **`ScopingPhaseStrip.tsx:26`** — `<div className="mb-1.5 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">…</div>` → `<Kicker tone="faint" className="mb-1.5">…</Kicker>`.
- [ ] **`OffsetDock.tsx:22`** — `<span className="text-[10px] font-bold uppercase tracking-wide text-ink-faint">Offset</span>` → `<Kicker as="span" tone="faint">Offset</Kicker>`.
- [ ] **`OnboardingFlow.tsx:328`** — `<div className="text-xs uppercase tracking-widest text-ink-faint">…</div>` → `<Kicker tone="faint">…</Kicker>` (was unweighted 400 → 700).
- [ ] **`SeriesBuilderRail.tsx:87`** — element with `className="text-xs uppercase tracking-[0.14em] text-ink-faint"` → `<Kicker tone="faint">…</Kicker>` (was 400 → 700).
- [ ] **`SeriesScopingPage.tsx:293`** — `<span className="text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">…</span>` → `<Kicker as="span" tone="faint">…</Kicker>`.
- [ ] **`SeriesBuilderPage.tsx:318`** — `<span data-testid="series-builder-editing-badge" className="text-[10.5px] uppercase tracking-wide text-ink-faint">…</span>` → `<Kicker as="span" tone="faint" data-testid="series-builder-editing-badge">…</Kicker>`. **PRESERVE `data-testid="series-builder-editing-badge"`**.
- [ ] **`CorpusTopbar.tsx:270`** — `<span className="text-[10px] uppercase tracking-wide text-ink-faint">sample {n} of {m}</span>` → `<Kicker as="span" tone="faint">sample {n} of {m}</Kicker>` (preserve the exact children expression; was 400 → 700).
- [ ] **`FocusNotesMargin.tsx:80`** — `<span className="text-meta uppercase tracking-wider text-ink-faint">Notes</span>` → `<Kicker as="span" tone="faint">Notes</Kicker>`.
- [ ] **`FocusWorkspaceLayout.tsx:165`** — `<span className="text-meta uppercase tracking-wider text-ink-faint">Notes</span>` → `<Kicker as="span" tone="faint">Notes</Kicker>`.

- [ ] **Run affected unit tests** — `cd packages/HimalayaUI/frontend && npx vitest run test/` → Expected: PASS, including any spec asserting `getByTestId("series-builder-editing-badge")`.
- [ ] **Build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (baseline shrinks).
- [ ] **Commit** — `git add -A && git commit -m "refactor(ui): migrate faint metric/field/aside labels to Kicker (17 sites)"`

---

### Task 50: Migrate the rotated heatmap axis title (placement-only transform stack)

The rotated decorative axis title is faint-tone; its `className` is a heavy but **placement/transform-only** stack (position + rotate + pointer-events) which is allowed under the placement-only contract. `aria-hidden="true"` and `data-testid="heatmap-axis-title"` must round-trip via `{...props}`.

```tsx
// BEFORE
<div data-testid="heatmap-axis-title" aria-hidden="true"
  className="pointer-events-none absolute left-0 top-1/2 z-10 origin-left -translate-y-1/2 -rotate-90 whitespace-nowrap text-[10.5px] uppercase tracking-[0.14em] text-ink-faint">{ordering_variable} →</div>
// AFTER — drop appearance (text-[10.5px]/uppercase/tracking/color); keep the entire transform+position stack as placement
<Kicker tone="faint" aria-hidden="true" data-testid="heatmap-axis-title"
  className="pointer-events-none absolute left-0 top-1/2 z-10 origin-left -translate-y-1/2 -rotate-90 whitespace-nowrap">{ordering_variable} →</Kicker>
```

**Files:**
- Modify: `src/pages/SeriesBuilderPage.tsx:376` (verify line before editing)

- [ ] **`SeriesBuilderPage.tsx:376`** — apply the transform above. **PRESERVE `aria-hidden="true"` and `data-testid="heatmap-axis-title"`**; keep the exact `{ordering_variable} →` children. Confirm the placement className contains NO appearance utilities (`text-[…]`, `uppercase`, `tracking-`, `text-ink-faint` all removed).
- [ ] **Run affected unit/e2e-touching spec** — `cd packages/HimalayaUI/frontend && npx vitest run test/` → Expected: PASS (any `getByTestId("heatmap-axis-title")` assertion still resolves and remains `aria-hidden`).
- [ ] **Build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS. Confirm `lint:design` accepts the transform/position utilities as placement (per the migrationRisk note — if the guard flags any of the transform classes, that is a guard-whitelist bug owned by the Token/guard task, NOT a reason to keep appearance classes here).
- [ ] **Commit** — `git add packages/HimalayaUI/frontend/src/pages/SeriesBuilderPage.tsx && git commit -m "refactor(ui): migrate rotated heatmap axis title to Kicker (placement-only transform stack)"`

---

### Task 51: Resolve borderline kicker-look sites (leave-as-is decisions, documented)

These sites LOOK like kickers but are out of Kicker's `accent|faint` two-tone scope per the locked decision. This task records the deliberate decisions so the guard baseline (Phase 0) keeps their inline utilities as accepted entries rather than churned. **No code change is made here** beyond optionally normalizing the `text-meta` faint-exact ones already handled in the faint task. The deliverable is the decision record below, applied by NOT migrating these and ensuring the guard baseline includes them.

**Decisions (do NOT migrate to Kicker):**

- [ ] **Table column-header strips STAY inline.** `SamplesPage.tsx:144` (grid-row `CONTACT_SHEET_COLS` header with border/padding) and `FocusReflectionsTable.tsx:110,115,120,125` (table column headers; `:120,:125` are `font-mono`). These carry table/grid-cell layout semantics and (for the mono pair) `font-mono`, which is appearance and cannot live on a plain Kicker. Leave as-is. Their banned-utility occurrences are recorded in the Phase 0 guard baseline.
- [ ] **Full-ink card-header labels STAY inline.** `SeriesRecipeEditor.tsx:131` ("Recipe", `text-ink`), `FocusDetectorPanel.tsx:98` ("Detector image", `text-meta`/ink), `FocusReflectionsTable.tsx:230` ("Reflections", `text-meta`/ink). The locked tone set is `accent|faint` only — there is no `ink` tone. Per the brief, do NOT silently recolor these to faint (it loses emphasis); leave them inline. Baseline records them.
- [ ] **`ConflictModalShell.tsx:138` `<header>` STAYS inline.** Element is `<header>`, which is not in the `as` union (`div|span|h2|h3`) and the locked `as` set is not to be widened. Leave inline (faint look preserved). Baseline records it.
- [ ] **`ScopingValueCell.tsx:74` 9px accent flag STAYS inline.** `text-[9px]` inline annotation in a dense table cell; the locked Fixed-Scale (no `size` prop) would jump it to 11.5px (+2.5px) and risk cell overflow. Leave as a one-off. Baseline records it.
- [ ] **Confirmed exclusions (NOT kickers at all):** `LoupeFrame.tsx:56` "Dropped" badge (badge/chip, `text-paper` on `bg-print-accent`), `FocusDetectorPanel.tsx:151` "Set rep" button label (Button concern), `CorpusTopbar.tsx:99` HIMALAYA wordmark (brand link), `CorpusTopbar.tsx:133,153` nav-tab buttons (SegmentedControl/Button concern). No action.

- [ ] **No-op verification** — `cd packages/HimalayaUI/frontend && git status --short` → Expected: no changes from this task (decisions only). If the Phase 0 guard task has already flipped Kicker-rule enforcement, run `npm run build` → Expected: PASS with these sites present in the baseline (not net-new violations).

### Task 52: Build the IconButton primitive (TDD)

The canonical icon-only / dismiss button. Standardizes 11+ hand-rolled buttons onto one
44px touch target (WCAG 2.5.5) + the canonical accent focus-visible ring (seeded from
`PhasePanel.tsx:163`, the only site that has it today) + canonical dismiss glyph `×`
(U+00D7).

**Locked decisions for this primitive (resolving the inventory's open questions):**
- **Hit area:** the 44px target is delivered by a centered `::before` pseudo-element, NOT
  a layout box, so it does not blow out dense chip rows (resolves migrationRisks #1). The
  visible button stays compact (`p-1`). The pseudo-element lives in a tiny `.icon-button`
  CSS rule in `styles.css` (component layer); tone/focus/hover stay Tailwind in the
  component so tests assert `data-tone`, never classes.
- **Dismiss-hover convention** (resolves openQuestion / migrationRisks #2): `tone="ghost"`
  hovers to `ink`; `tone="danger"` hovers to `error` (destructive delete). Chip-remove ×
  uses `tone="accent"` → hovers to the terracotta accent (DESIGN.md:121 makes terracotta
  the reject/grease-pencil mark). So the tone axis is **`ghost | accent | danger`** — three
  tones, not two. (The brief names `ghost|danger`; we add `accent` to faithfully preserve
  the three existing `hover:text-print-accent` tag-removes without relitigating their look.)
- **`label` is required** and is the accessible name; `aria-label` is `Omit`ted from the
  spread so a caller cannot bypass it. `title` defaults to `label`.
- **`dismiss`** renders the canonical `×` (U+00D7). The grease-pencil reject **✕** (U+2715)
  on detector frames (`CrossTraceTrackingLayer:124`, `ContactSheetRow` frame mark,
  `LoupeFrame`) is a DOMAIN mark and stays untouched — NOT an IconButton.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ui/IconButton.tsx`
- Modify: `packages/HimalayaUI/frontend/src/styles.css` (append `.icon-button` rule to the `@layer components` block, after `.text-data-strong`)
- Modify: `packages/HimalayaUI/frontend/src/components/ui/index.ts` (add the export)
- Test: `packages/HimalayaUI/frontend/test/ui/IconButton.test.tsx`

- [ ] **Step 1: Write the failing test** — create `test/ui/IconButton.test.tsx`:

```tsx
import { render, screen } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import userEvent from "@testing-library/user-event";
import { IconButton } from "../../src/components/ui/IconButton";

describe("IconButton", () => {
  it("uses label as the accessible name", () => {
    render(<IconButton label="Close notes" dismiss />);
    expect(screen.getByRole("button", { name: "Close notes" })).toBeInTheDocument();
  });

  it("dismiss renders the canonical × glyph (U+00D7)", () => {
    render(<IconButton label="Dismiss" dismiss />);
    expect(screen.getByRole("button", { name: "Dismiss" })).toHaveTextContent("×");
  });

  it("renders children instead of the glyph when not dismiss", () => {
    render(<IconButton label="Move up">{"▲"}</IconButton>);
    expect(screen.getByRole("button", { name: "Move up" })).toHaveTextContent("▲");
  });

  it("defaults tone to ghost and exposes data-tone", () => {
    render(<IconButton label="x" dismiss />);
    expect(screen.getByRole("button", { name: "x" })).toHaveAttribute("data-tone", "ghost");
  });

  it("exposes data-tone for accent and danger", () => {
    const { rerender } = render(<IconButton label="remove" dismiss tone="accent" />);
    expect(screen.getByRole("button", { name: "remove" })).toHaveAttribute("data-tone", "accent");
    rerender(<IconButton label="delete" dismiss tone="danger" />);
    expect(screen.getByRole("button", { name: "delete" })).toHaveAttribute("data-tone", "danger");
  });

  it("title defaults to label but a passed title wins", () => {
    const { rerender } = render(<IconButton label="Remove from recipe" dismiss />);
    expect(screen.getByRole("button", { name: "Remove from recipe" })).toHaveAttribute("title", "Remove from recipe");
    rerender(<IconButton label="Remove from recipe" dismiss title="Custom" />);
    expect(screen.getByRole("button", { name: "Remove from recipe" })).toHaveAttribute("title", "Custom");
  });

  it("defaults type to button (does not submit forms)", () => {
    render(<IconButton label="x" dismiss />);
    expect(screen.getByRole("button", { name: "x" })).toHaveAttribute("type", "button");
  });

  it("is disabled and does not fire onClick when disabled", async () => {
    const onClick = vi.fn();
    render(<IconButton label="Move up" disabled onClick={onClick}>{"▲"}</IconButton>);
    const btn = screen.getByRole("button", { name: "Move up" });
    expect(btn).toBeDisabled();
    await userEvent.click(btn);
    expect(onClick).not.toHaveBeenCalled();
  });

  it("fires onClick when enabled", async () => {
    const onClick = vi.fn();
    render(<IconButton label="Dismiss" dismiss onClick={onClick} />);
    await userEvent.click(screen.getByRole("button", { name: "Dismiss" }));
    expect(onClick).toHaveBeenCalledTimes(1);
  });

  it("forwards data-testid and placement className to the button", () => {
    render(<IconButton label="x" dismiss data-testid="my-id" className="shrink-0" />);
    const btn = screen.getByTestId("my-id");
    expect(btn).toHaveClass("shrink-0");
  });
});
```

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/IconButton.test.tsx` → Expected: FAIL (module `../../src/components/ui/IconButton` does not exist).

- [ ] **Step 3: Implement** — create `src/components/ui/IconButton.tsx`:

```tsx
import type { ButtonHTMLAttributes, ReactNode } from "react";

export type IconButtonTone = "ghost" | "accent" | "danger";

interface IconButtonProps
  extends Omit<ButtonHTMLAttributes<HTMLButtonElement>, "aria-label"> {
  /** REQUIRED — icon-only buttons have no text; this is the accessible name. */
  label: string;
  /** hover intent. ghost → hover:text-ink (default); accent → hover:text-print-accent
   *  (chip-remove, the terracotta reject mark); danger → hover:text-error (destructive). */
  tone?: IconButtonTone;
  /** render the canonical dismiss glyph (×, U+00D7) as the content. */
  dismiss?: boolean;
  /** glyph or SVG when not a dismiss button (chevron, trash, reorder arrow). */
  children?: ReactNode;
}

const toneClass: Record<IconButtonTone, string> = {
  ghost: "text-ink-faint hover:text-ink",
  accent: "text-ink-faint hover:text-print-accent",
  danger: "text-ink-faint hover:text-error",
};

function cx(...parts: Array<string | false | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export function IconButton({
  label,
  tone = "ghost",
  dismiss = false,
  className = "",
  children,
  type = "button",
  title,
  ...props
}: IconButtonProps): JSX.Element {
  return (
    <button
      type={type}
      aria-label={label}
      title={title ?? label}
      data-tone={tone}
      className={cx(
        // `.icon-button` (styles.css) supplies the >=44px hit-area pseudo-element
        // so the visible box stays compact and never balloons dense chip rows.
        "icon-button relative inline-flex items-center justify-center",
        "rounded p-1 leading-none transition-colors",
        "focus-visible:outline focus-visible:outline-2",
        "focus-visible:outline-offset-2 focus-visible:outline-accent",
        "disabled:cursor-not-allowed disabled:opacity-30",
        toneClass[tone],
        className,
      )}
      {...props}
    >
      {dismiss ? "×" : children}
    </button>
  );
}
```

  Then append the hit-area rule to `src/styles.css` inside the existing `@layer components { … }`
  block (after the `.text-data-strong` rule, before the block's closing `}`):

```css
  /* Icon-only buttons keep a compact visible box but expand the pointer/touch
     target to >=44px (WCAG 2.5.5) via a centered overlay, so dense chip rows
     do not balloon. The overlay is transparent and pointer-events:auto so it
     forwards clicks to the button. */
  .icon-button::before {
    content: "";
    position: absolute;
    top: 50%;
    left: 50%;
    width: 44px;
    height: 44px;
    transform: translate(-50%, -50%);
  }
```

  Then add the export to `src/components/ui/index.ts`:

```ts
export { IconButton } from "./IconButton";
export type { IconButtonTone } from "./IconButton";
```

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/IconButton.test.tsx` → Expected: PASS (all cases green).

- [ ] **Step 5: Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (`lint:design` + `tsc --noEmit -p tsconfig.build.json` + `vite build` all clean; new primitive is unused so far but compiles).

- [ ] **Step 6: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/IconButton.tsx packages/HimalayaUI/frontend/src/components/ui/index.ts packages/HimalayaUI/frontend/src/styles.css packages/HimalayaUI/frontend/test/ui/IconButton.test.tsx && git commit -m "feat(ui): add IconButton primitive — 44px hit area, accent focus ring, canonical × dismiss"`

---

### Task 53: Migrate dismiss / icon-only buttons to IconButton

Convert every icon-only / dismiss `<button>` in the inventory's IconButton callSites to the
new primitive, preserving each `data-testid` exactly (they back E2E + unit specs) and adding
`aria-label` to the three title-only buttons (`recipe-row-up/down/remove`). Each gains the
accent focus-visible ring + the 44px hit area it lacked.

**Transformation pattern (apply at every site):**

```tsx
// BEFORE (a dismiss × with no focus ring):
<button type="button" aria-label="Dismiss" onClick={onClose}
  className="text-ink-soft hover:text-ink leading-none px-1">×</button>

// AFTER:
<IconButton label="Dismiss" dismiss onClick={onClose} />

// BEFORE (a glyph/SVG icon button, e.g. a stepper):
<button type="button" data-testid="sample-stepper-prev" disabled={…}
  onClick={…} aria-label="Previous sample"
  className="rounded px-1.5 py-0.5 text-base leading-none text-ink-faint
             hover:text-ink disabled:cursor-not-allowed disabled:opacity-30">
  &#8249;
</button>

// AFTER (tone/focus/disabled/hit-area owned by the primitive; only placement stays):
<IconButton label="Previous sample" tone="ghost" disabled={…}
  onClick={…} data-testid="sample-stepper-prev">{"‹"}</IconButton>
```

Rules applied uniformly: drop `type="button"` (primitive defaults it), drop the color/hover/
focus/disabled/padding/leading utilities (primitive owns them), keep ONLY placement utilities
(`shrink-0`, etc.) in `className`, normalize every glyph (`&#215;`, `✕` U+2715) to the
canonical `×` via `dismiss`. Each call site needs `import { IconButton } from "./ui/IconButton";`
(or `"../components/ui/IconButton"` per file depth — verify the relative path per file).

Convert these sites (verify the line before editing — inventory `loc` may have drifted):

- [ ] **`src/components/ui/Toast.tsx:93-98`** — dismiss button. `aria-label="Dismiss"` glyph `×`. → `<IconButton label="Dismiss" dismiss tone="ghost" onClick={() => dismiss(t.id)} />`. (Toast lives in `ui/`, so import is `./IconButton`.) Preserve the `role="status"`, `data-testid="toast"`, `data-toast-kind` on the *container* — they are not on this button, leave them. Gains focus ring + 44px (was missing).
- [ ] **`src/components/NavModal.tsx:333-341`** — Chip remove `×`, `aria-label={`Remove ${label}`}`, `data-testid={testId ? `${testId}-remove` : undefined}`, hover:text-error → chip-remove ⇒ `tone="accent"` per the locked convention. → `<IconButton label={`Remove ${label}`} dismiss tone="accent" onClick={onRemove} data-testid={testId ? `${testId}-remove` : undefined} />`. **Preserve the `${testId}-remove` testid exactly.** Lives inside NavModal's local `Chip` component — convert in place (do NOT extract a Chip primitive; that is a deferred future wave per the inventory).
- [ ] **`src/components/FocusWorkspaceLayout.tsx:168-174`** — `data-testid="focus-notes-drawer-close"`, `aria-label="Close notes"`, glyph `&#215;`. → `<IconButton label="Close notes" dismiss tone="ghost" onClick={closeNotesDrawer} data-testid="focus-notes-drawer-close" />`. Normalizes `&#215;` → `×`. **Preserve `focus-notes-drawer-close` testid.**
- [ ] **`src/components/SpeculativeBuilder.tsx:130-133`** — close `×`, `aria-label="Close"`, `text-xl` glyph sizing. → `<IconButton label="Close" dismiss tone="ghost" onClick={onClose} />`. (Glyph size normalizes to the primitive's compact size — accept the minor shrink from `text-xl`; flagged as visual-parity risk in the inventory, no testid here.)
- [ ] **`src/components/CorpusTopbar.tsx:255-264`** — stepper-prev. `data-testid="sample-stepper-prev"`, `disabled`, `aria-label="Previous sample"`, glyph `&#8249;` (‹). This is a glyph icon button, **NOT dismiss**. → `<IconButton label="Previous sample" tone="ghost" disabled={prevSample === undefined} onClick={() => prevSample && navigate(`/sample/${prevSample.id}`)} data-testid="sample-stepper-prev">{"‹"}</IconButton>`. **Preserve testid + disabled.**
- [ ] **`src/components/CorpusTopbar.tsx:275-284`** — stepper-next. Same shape, glyph `&#8250;` (›). → `<IconButton label="Next sample" tone="ghost" disabled={nextSample === undefined} onClick={() => nextSample && navigate(`/sample/${nextSample.id}`)} data-testid="sample-stepper-next">{"›"}</IconButton>`. **Preserve `sample-stepper-next` testid + disabled.**
- [ ] **`src/components/PhasePanel.tsx:160-176`** — the canonical seed (`data-testid={`spec-delete-${index.id}`}`, trash SVG, full focus ring already present, `hover:text-error`, `title="Delete this speculative index"`, `onClick` with `stopPropagation`). → `<IconButton label={`Delete speculative index ${index.id}`} tone="danger" title="Delete this speculative index" data-testid={`spec-delete-${index.id}`} onClick={(e) => { e.stopPropagation(); onDelete(); }} className="shrink-0">` wrapping the existing `<svg viewBox="0 0 16 16" className="h-3.5 w-3.5" aria-hidden="true">…</svg>` child unchanged `</IconButton>`. **Preserve testid + title + the SVG glyph + stopPropagation; `shrink-0` stays in placement className.** This is a behavior-preserving swap (already had the ring) — verify the rendered focus utilities are unchanged.
- [ ] **`src/components/LoupeSidebar.tsx:121-128`** — tag remove `×`, `aria-label={`Remove ${tag.key || tag.value} tag`}`, `hover:text-print-accent`. Chip-remove ⇒ `tone="accent"` (preserves the print-accent hover). → `<IconButton label={`Remove ${tag.key || tag.value} tag`} dismiss tone="accent" onClick={() => remove.mutate(tag.id)} />`. Gains focus ring + 44px.
- [ ] **`src/components/LoupeSidebar.tsx:173-180`** — cancel-add `×`, `aria-label="cancel adding tag"`, `hover:text-print-accent`. → `<IconButton label="cancel adding tag" dismiss tone="accent" onClick={reset} />`.
- [ ] **`src/components/ContactSheetRow.tsx:425-432`** — tag remove `×`, `aria-label={`Remove ${t.key || t.value} tag`}`, `hover:text-print-accent`. → `<IconButton label={`Remove ${t.key || t.value} tag`} dismiss tone="accent" onClick={() => removeTag.mutate(t.id)} />`.
- [ ] **`src/components/ContactSheetRow.tsx:477-483`** (the second `×`, the tag-form cancel) — `aria-label="cancel adding tag"`, `hover:text-print-accent`, `onClick={resetTagForm}`. → `<IconButton label="cancel adding tag" dismiss tone="accent" onClick={resetTagForm} />`. (Note: the "Add" button just above at ~470 is a TEXT button — leave it alone; only the `×` cancel migrates.)
- [ ] **`src/components/SeriesRecipeEditor.tsx:144-153`** — `recipe-row-up`. `disabled`, `title="Move up"`, **NO aria-label today** (a11y fix), glyph `▲`, NOT dismiss. → `<IconButton label="Move up" tone="ghost" disabled={i === 0} onClick={() => reorderSample(i, i - 1)} title="Move up" data-testid="recipe-row-up">{"▲"}</IconButton>`. **Adds aria-label; preserves testid + disabled.**
- [ ] **`src/components/SeriesRecipeEditor.tsx:154-163`** — `recipe-row-down`. `disabled`, `title="Move down"`, no aria-label, glyph `▼`. → `<IconButton label="Move down" tone="ghost" disabled={i === draft.recipe.length - 1} onClick={() => reorderSample(i, i + 1)} title="Move down" data-testid="recipe-row-down">{"▼"}</IconButton>`. **Adds aria-label; preserves testid + disabled.**
- [ ] **`src/components/SeriesRecipeEditor.tsx:164-172`** — `recipe-row-remove`. `title="Remove from recipe"`, no aria-label, glyph `✕` (U+2715). → `<IconButton label="Remove from recipe" dismiss tone="ghost" onClick={() => removeSample(row.id)} title="Remove from recipe" data-testid="recipe-row-remove" />`. **Adds aria-label; normalizes ✕ → ×; preserves testid.** (Recipe-row remove is a list-row action, ghost not accent — it is not a chip-remove pill.)

- [ ] **Leave untouched (domain marks, NOT IconButton):** the grease-pencil reject **✕** on detector frames at `CrossTraceTrackingLayer.tsx:124`, the `ContactSheetRow` frame reject mark (the inventory's `ContactSheetRow:50`), and `LoupeFrame`. Do not convert these — verify by grepping `✕`/`✕` in those files and confirming each is a frame overlay glyph, not an action button.

- [ ] **Run the affected specs** — `cd packages/HimalayaUI/frontend && npx vitest run test/Toast.test.tsx test/NavModal.test.tsx test/FocusWorkspaceLayout.test.tsx test/SpeculativeBuilder.test.tsx test/CorpusTopbar.test.tsx test/PhasePanel.test.tsx test/LoupeSidebar.test.tsx test/ContactSheetRow.test.tsx test/SeriesRecipeEditor.test.tsx` → Expected: PASS (every existing `data-testid` still resolves; queries by accessible name still match because `label` reproduces each prior `aria-label`). If a spec file name differs, discover it with `ls test | grep -i <component>` first; if any test queried by `aria-label` for the formerly-title-only buttons, it now passes where it previously could not have.

- [ ] **Full unit suite** — `cd packages/HimalayaUI/frontend && npm test` → Expected: PASS (no regressions across the suite).

- [ ] **Type + build gate** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS. The migrated files dropped color/border/focus/padding utilities from these buttons' `className`, which lowers the `lint:design` violation backlog — if `lint:design` reports baseline entries are now *unused*/over-counted for these sites, drain them from `scripts`/`design-baseline.json` per the guard's "lower the baseline by deleted call sites" rule (Phase 0 owns the baseline mechanism; here just remove the now-fixed entries so CI stays green).

- [ ] **Commit** — `git add -A && git commit -m "refactor(ui): migrate 14 dismiss/icon buttons to IconButton — focus rings, 44px targets, +3 aria-labels"`

**Scope.** This unit ships the `Card` primitive (`{ elevated?, border?, as? }`, default flat) and migrates
the genuine card surfaces from the inventory: the **4 lifted plates** → `elevated`, and the **flat hairline
cards** → `Card` default. `elevated` reuses the existing `.card` CSS recipe so the Plate-Lift shadow stays a
single source of truth. Every over-captured `bg-plate` grep hit that is *not* a card (inputs, buttons,
chips, banners, modals, segmented container, floaters/tooltips/ink-toolbar) is explicitly routed elsewhere
and **NOT touched here** (enumerated in the final "Explicitly-not-Card routing" task).

**Hard dependency (cross-unit).** The radius decision (`--radius-sm` = `--radius-md` = 5px; `.card`'s
hardcoded `border-radius: 12px` → 5px; `--radius` singular deleted) is owned by the **Token-foundation /
Phase 0** unit, NOT this one. `Card`'s flat branch emits `rounded-md` and the `elevated` branch applies
`.card`; both render 5px **only after** the Phase 0 token+`.card` change has landed. Land Phase 0 first.
This unit must NOT re-edit `.card`'s box-shadow or radius.

**Locked design facts used below (from the brief + spec §8):**
- `.card` rule (`src/styles.css:172`): `background: var(--color-plate)` + `border: 0.5px solid var(--color-hair)` + radius (5px post-Phase-0) + the Plate-Lift double shadow. `elevated` = apply this class verbatim; do not re-inline the shadow.
- Flat default = `rounded-md bg-plate` + a hairline (`border-hair` default, `border-hair-strong` when `border="strong"`) + **NO shadow** (Flat-Except-the-Plate).
- `data-elevated="true"` is added when elevated; absent when flat. Tests assert this `data-*`, never class strings.
- `className` is **placement-only** (margin, width, max-w, overflow, grid/flex position, sticky/fixed offsets, motion like `hover:-translate-y-0.5` / `transition-[max-width]`). Appearance (`bg-*`/`border-*`/`shadow-*`/`rounded-*`) is owned by the primitive.
- Polymorphic `as` is required: `SeriesFolioCard` is a `<button>`; the masonry tile is a `<div>`. Support `as?: "div" | "button" | "section" | "li"` (default `"div"`), spreading the right HTML attrs.

---

### Task 54: Scaffold the `Card` primitive (failing test first)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ui/Card.tsx`
- Test: `packages/HimalayaUI/frontend/test/ui/Card.test.tsx`

- [ ] **Step 1: Write the failing test** — create `test/ui/Card.test.tsx`:
```tsx
import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Card } from "../../src/components/ui/Card";

describe("Card", () => {
  it("renders children inside a div by default", () => {
    render(<Card data-testid="c">hello</Card>);
    const el = screen.getByTestId("c");
    expect(el.tagName).toBe("DIV");
    expect(el).toHaveTextContent("hello");
  });

  it("is flat by default — no data-elevated attribute", () => {
    render(<Card data-testid="c">x</Card>);
    expect(screen.getByTestId("c")).not.toHaveAttribute("data-elevated");
  });

  it("elevated sets data-elevated=\"true\"", () => {
    render(<Card data-testid="c" elevated>x</Card>);
    expect(screen.getByTestId("c")).toHaveAttribute("data-elevated", "true");
  });

  it("passes placement className through to the root", () => {
    render(<Card data-testid="c" className="mb-5 max-w-[760px]">x</Card>);
    expect(screen.getByTestId("c").className).toContain("mb-5");
    expect(screen.getByTestId("c").className).toContain("max-w-[760px]");
  });

  it("renders as a button when as=\"button\" and forwards button props", () => {
    const onClick = vi.fn();
    render(<Card as="button" data-testid="c" onClick={onClick}>go</Card>);
    const el = screen.getByTestId("c");
    expect(el.tagName).toBe("BUTTON");
    el.click();
    expect(onClick).toHaveBeenCalledOnce();
  });

  it("renders as li / section when requested", () => {
    render(<Card as="li" data-testid="li">x</Card>);
    render(<Card as="section" data-testid="sec">y</Card>);
    expect(screen.getByTestId("li").tagName).toBe("LI");
    expect(screen.getByTestId("sec").tagName).toBe("SECTION");
  });
});
```
(Add `import { vi } from "vitest";` to the import line — fold `vi` into the existing `{ describe, it, expect }` import: `import { describe, it, expect, vi } from "vitest";`.)

- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/Card.test.tsx` → Expected: FAIL (module `../../src/components/ui/Card` does not exist).

- [ ] **Step 3: Implement** — create `src/components/ui/Card.tsx`:
```tsx
import type {
  ButtonHTMLAttributes,
  HTMLAttributes,
  LiHTMLAttributes,
} from "react";

type CardElement = "div" | "button" | "section" | "li";

type ElementProps = {
  div: HTMLAttributes<HTMLDivElement>;
  section: HTMLAttributes<HTMLElement>;
  li: LiHTMLAttributes<HTMLLIElement>;
  button: ButtonHTMLAttributes<HTMLButtonElement>;
};

type CardOwnProps<T extends CardElement> = {
  /** Render element. Default "div". The folio card is a <button>; masonry
   *  tiles are <li>/<div>. Placement only — appearance is fixed. */
  as?: T;
  /** true → the single lifted "Plate Lift" object (applies the `.card` CSS
   *  rule: bg-plate + hairline + radius + inset highlight + soft warm drop
   *  shadow). false (default) → flat hairline card: bg-plate + hairline +
   *  radius, NO shadow (Flat-Except-the-Plate). */
  elevated?: boolean;
  /** hairline weight for the FLAT variant. "hair" (default) for inner cards;
   *  "strong" for outer/standalone surfaces. Ignored when elevated (`.card`
   *  owns the border). */
  border?: "hair" | "strong";
};

export type CardProps<T extends CardElement = "div"> = CardOwnProps<T> &
  Omit<ElementProps[T], keyof CardOwnProps<T>>;

const borderClass: Record<NonNullable<CardOwnProps<CardElement>["border"]>, string> = {
  hair: "border border-hair",
  strong: "border border-hair-strong",
};

export function Card<T extends CardElement = "div">({
  as,
  elevated = false,
  border = "hair",
  className = "",
  children,
  ...rest
}: CardProps<T>): JSX.Element {
  const Tag = (as ?? "div") as CardElement;
  // elevated → `.card` (single source of the shadow). flat → tonal hairline,
  // never a shadow. Radius token (rounded-md = 5px) comes from Phase 0.
  const appearance = elevated
    ? "card"
    : `rounded-md bg-plate ${borderClass[border]}`;
  const props = {
    className: `${appearance} ${className}`.trim(),
    ...(elevated ? { "data-elevated": "true" } : {}),
    ...rest,
  } as Record<string, unknown>;
  // @ts-expect-error — runtime element name; props are validated by CardProps<T>.
  return <Tag {...props}>{children}</Tag>;
}
```

- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/Card.test.tsx` → Expected: PASS (6 tests).

- [ ] **Step 5: Type-gate** — `cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json` → Expected: no errors from `Card.tsx`.

- [ ] **Step 6: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/Card.tsx packages/HimalayaUI/frontend/test/ui/Card.test.tsx && git commit -m "Add Card/Plate primitive (elevated? + border + polymorphic as)"`

---

### Task 55: Export `Card` from the `ui` barrel

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ui/index.ts`

- [ ] **Step 1: Add the export** — append to `src/components/ui/index.ts` (place after the `Button` export):
```ts
export { Card } from "./Card";
export type { CardProps } from "./Card";
```
(Note: removal of the `Input` and `SectionLabel` exports + adding `ToastContainer` is owned by the **Phase 0 adopt/delete** unit; do not touch those lines here.)

- [ ] **Step 2: Type-gate** — `cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json` → Expected: no errors.

- [ ] **Step 3: Commit** — `git add packages/HimalayaUI/frontend/src/components/ui/index.ts && git commit -m "Export Card from ui barrel"`

---

### Task 56: Migrate the 4 elevated plates → `<Card elevated>`

The four genuine lifted surfaces named in DESIGN.md §Elevation: the figure plate, the focus-workspace
panels (detector panel + reflections table both apply `.card` today), the folio card, and the scoping +
builder worksheets. **Transformation pattern** (shown once): the element that currently carries `.card`
(or a hand-rolled `bg-plate + border + inline Plate-Lift shadow`) becomes `<Card elevated …>`; the
appearance utilities (`bg-plate`, `border-*`, `rounded-*`, `shadow-[…]`, the literal `card` class) are
**removed** from `className`, leaving only placement/layout/motion utilities. Preserve every `data-testid`
and `data-*`.

> **Visual change shipping (signed off, spec §8):** the two worksheet inline shadows
> (`SeriesScopingPage`, `SeriesBuilderPage`) use a *warmer/shallower* recipe than `.card`; folding them into
> `elevated` adopts `.card`'s cooler recipe. This is intentional (single source of truth). `SeriesFolioCard`
> currently has **no** shadow; `elevated` **adds** the Plate-Lift — also intentional (DESIGN.md names folio
> cards as lifted plates).

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PlotCard.tsx:362-366` (verify lines before editing — inventory `loc` is `:359`, drift flagged)
- Modify: `packages/HimalayaUI/frontend/src/components/FocusDetectorPanel.tsx:96`
- Modify: `packages/HimalayaUI/frontend/src/components/FocusReflectionsTable.tsx:227`
- Modify: `packages/HimalayaUI/frontend/src/components/SeriesFolioCard.tsx:80-92`
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesScopingPage.tsx:236`
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesBuilderPage.tsx:277`

- [ ] **Step 1: PlotCard figure plate** — the outer wrapper conditionally applies `.card` only in the
  focus variant (`headerSlot` present). To preserve the *non-focus* contract (no plate, no bg) while making
  the figure plate a `Card elevated`, keep the conditional and route only the elevated branch through
  `Card`. Replace the `outerClass`/return at `PlotCard.tsx:362-366`:
```tsx
  // R3-N3: the trace plate is the hero. In the focus variant it is the single
  // elevated object (Card elevated). Non-focus consumers stay an unstyled flow
  // container (no plate, no bg) — gated on headerSlot exactly as before.
  const layout = "flex flex-col h-full min-h-0 overflow-hidden";
  return (
    headerSlot ? (
      <Card elevated data-testid="plot-card" className={layout}>
        <TitleStrip
          headerSlot={headerSlot}
          experimentName={experimentName}
```
  …and the matching closing tag for this branch becomes `</Card>`, with an `else` branch rendering the
  original `<div data-testid="plot-card" className={layout}>…</div>`. **Verify the JSX block before
  editing** — preserve the full `TitleStrip`/children subtree and the `experimentName` etc. props exactly;
  only the wrapper element and the `.card` class move. Add `import { Card } from "./ui";` at the top.
  > RISK: this is the only conditional `.card` site. Do NOT collapse the two branches into one `<Card>` — the
  > non-focus branch must stay shadow/plate-free. If splitting the return is awkward, hoist the children
  > into a `const body = (…)` and render `headerSlot ? <Card elevated …>{body}</Card> : <div …>{body}</div>`.

- [ ] **Step 2: FocusDetectorPanel** — `FocusDetectorPanel.tsx:96`: replace the wrapper
  `<… className="card flex min-h-0 flex-col p-4">` with `<Card elevated className="flex min-h-0 flex-col p-4">`
  (drop `card`; keep `flex min-h-0 flex-col p-4` placement; preserve any existing `data-testid`/`ref`).
  Add `import { Card } from "./ui";`. Verify the matching close tag is updated to `</Card>`.

- [ ] **Step 3: FocusReflectionsTable** — `FocusReflectionsTable.tsx:227`: same transform —
  `className="card flex min-h-0 flex-col p-4"` → `<Card elevated className="flex min-h-0 flex-col p-4">`,
  close tag `</Card>`. Add `import { Card } from "./ui";`.

- [ ] **Step 4: SeriesFolioCard** — `SeriesFolioCard.tsx:80-92` is a `<button>`; migrate to
  `<Card as="button" elevated border="strong" …>`. The draft state uses a **dashed** `border-hair-strong`;
  the primitive does not model dashed borders, so the dashed-when-draft border stays a placement-className
  override (`border-dashed`) layered on top, and `hover:-translate-y-0.5` (motion) stays in `className`.
  Replace lines 80-92 with:
```tsx
    <Card
      as="button"
      elevated
      type="button"
      data-testid={`series-card-${series.id}`}
      data-member-count={series.member_count}
      data-draft={isDraft ? "true" : "false"}
      data-stale={series.has_stale_members ? "true" : "false"}
      onClick={() => onOpen(series.id)}
      className={[
        "mb-5 block w-full break-inside-avoid overflow-hidden text-left",
        "transition-transform hover:-translate-y-0.5",
        isDraft ? "border-dashed" : "",
      ].join(" ")}
    >
```
  Update the matching closing `</button>` → `</Card>`. Add `import { Card } from "./ui";`.
  > NOTE: `elevated` applies `.card` which sets `border-hair`; the draft dashed override only changes the
  > style, not the color. If a draft card must read as `hair-strong` dashed, append `border-hair-strong` to
  > the draft branch string. Confirm against the bone/E2E snapshot for `series-card-*`.

- [ ] **Step 5: SeriesScopingPage worksheet** — `SeriesScopingPage.tsx:236`: replace
  `className="w-full max-w-[760px] rounded-md border border-hair bg-plate px-8 py-7 shadow-[0_1px_0_rgba(255,255,255,.6)_inset,0_1px_1px_rgba(60,52,40,.04),0_18px_40px_-22px_rgba(60,52,40,.22)]"`
  with `<Card elevated className="w-full max-w-[760px] px-8 py-7">` (drop the inline shadow + bg + border +
  rounded; keep width + padding). Update the close tag to `</Card>`. Add `import { Card } from "../components/ui";`.

- [ ] **Step 6: SeriesBuilderPage worksheet** — `SeriesBuilderPage.tsx:277`: replace
  `className={\`w-full ${collapsed ? "max-w-[1336px]" : "max-w-[1180px]"} rounded border border-hair bg-plate p-8 shadow-[…] transition-[max-width] duration-200\`}`
  with
  `<Card elevated className={\`w-full ${collapsed ? "max-w-[1336px]" : "max-w-[1180px]"} p-8 transition-[max-width] duration-200\`}>`
  (drop shadow + bg + border + rounded; keep the conditional max-w, padding, and the `transition-[max-width]`
  motion). Update the close tag to `</Card>`. Add `import { Card } from "../components/ui";`.

- [ ] **Step 7: Type-gate** — `cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json` → Expected: no errors.

- [ ] **Step 8: Run affected component suites** — `cd packages/HimalayaUI/frontend && npx vitest run test/PlotCard.test.tsx test/SeriesFolioCard.test.tsx test/scoping.test.tsx` (run only the files that exist; use `npx vitest run -t "PlotCard"` / file globs as available) → Expected: PASS, every `series-card-*` / `plot-card` `data-testid` still resolves. If a test pins the literal `card` class string, that is a class-string assertion that should be relaxed to a `data-elevated`/`data-testid` assertion in the same commit (cite the failing line).

- [ ] **Step 9: Commit** — `git add packages/HimalayaUI/frontend/src/components/PlotCard.tsx packages/HimalayaUI/frontend/src/components/FocusDetectorPanel.tsx packages/HimalayaUI/frontend/src/components/FocusReflectionsTable.tsx packages/HimalayaUI/frontend/src/components/SeriesFolioCard.tsx packages/HimalayaUI/frontend/src/pages/SeriesScopingPage.tsx packages/HimalayaUI/frontend/src/pages/SeriesBuilderPage.tsx && git commit -m "Migrate the 4 lifted plates to <Card elevated> (single shadow source)"`

---

### Task 57: Migrate the flat hairline cards → `<Card>` (default)

The genuine *flat* inner/standalone cards — `bg-plate + hairline + radius`, no shadow. **Transformation
pattern** (shown once): drop `bg-plate`, `border …`, and the radius utility from `className`; render
`<Card …>` (no `elevated`), choosing `border="strong"` only if the original used `border-hair-strong`
(none in this set do — all are `border-hair`). Normalize the two off-token `rounded-lg` PhasePanel
wrappers onto the primitive's `rounded-md`. Preserve `data-testid`s and layout utilities (`overflow-hidden`,
`p-3.5`, etc.).

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesFolioPage.tsx:18` (bone fixture)
- Modify: `packages/HimalayaUI/frontend/src/components/AutogroupCard.tsx:20`
- Modify: `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx:210` (bone fixture) and `:308` (live)
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesBuilderPage.tsx:30` (skeleton block)

- [ ] **Step 1: AutogroupCard** — `AutogroupCard.tsx:20`: replace
  `className="rounded-md border border-hair bg-plate p-3.5"` with `<Card data-testid="autogroup-card" className="p-3.5">`
  (move `data-testid` onto `Card`; keep `p-3.5`; drop appearance). Update the close tag to `</Card>`. Add
  `import { Card } from "./ui";`. Verify the `data-testid="autogroup-card"` lands on the root (it currently does).

- [ ] **Step 2: PhasePanel live phase-call wrapper** — `PhasePanel.tsx:308`: replace
  `<div className="overflow-hidden rounded-lg border border-hair bg-plate">` with
  `<Card className="overflow-hidden">` (drop `rounded-lg` → primitive `rounded-md`; keep `overflow-hidden`).
  Update the matching `</div>` → `</Card>`. Add `import { Card } from "./ui";` if not already importing from `./ui`.

- [ ] **Step 3: PhasePanel bone fixture wrapper** — `PhasePanel.tsx:210` (inside `PHASE_PANEL_FIXTURE`):
  same transform as Step 2 — `<div className="overflow-hidden rounded-lg border border-hair bg-plate">` →
  `<Card className="overflow-hidden">`, close `</Card>`. (Bone fixtures render through the same primitive so
  the captured skeleton matches the live layout.)

- [ ] **Step 4: SeriesFolioPage bone fixture tile** — `SeriesFolioPage.tsx:18`: replace
  `<div key={i} className="mb-5 break-inside-avoid rounded border border-hair bg-plate">` with
  `<Card key={i} className="mb-5 break-inside-avoid">` (drop `rounded`/`border`/`bg-plate`; keep masonry
  placement). Update close tag `</div>` → `</Card>`. Add `import { Card } from "../components/ui";`.

- [ ] **Step 5: SeriesBuilderPage skeleton block** — `SeriesBuilderPage.tsx:30`: replace
  `<div className="m-4 rounded border border-hair bg-plate" style={{ aspectRatio: "10 / 3" }} />` with
  `<Card className="m-4" style={{ aspectRatio: "10 / 3" }} />` (keep `m-4` placement + the inline
  `aspectRatio` style; drop appearance). The `import { Card } from "../components/ui";` is already added in
  the elevated-plates task for this file — reuse it.

- [ ] **Step 6: Type-gate** — `cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json` → Expected: no errors.

- [ ] **Step 7: Run affected suites** — `cd packages/HimalayaUI/frontend && npx vitest run test/AutogroupCard.test.tsx test/PhasePanel.test.tsx` (run only files that exist) → Expected: PASS, `autogroup-card` / `phase-call-empty` `data-testid`s still resolve.

- [ ] **Step 8: Commit** — `git add packages/HimalayaUI/frontend/src/components/AutogroupCard.tsx packages/HimalayaUI/frontend/src/components/PhasePanel.tsx packages/HimalayaUI/frontend/src/pages/SeriesFolioPage.tsx packages/HimalayaUI/frontend/src/pages/SeriesBuilderPage.tsx && git commit -m "Migrate flat hairline cards to <Card> default (normalize rounded-lg → rounded-md)"`

---

### Task 58: Document the explicitly-not-Card routing (no code change here, but assert the boundary)

The `bg-plate` grep over-captures. The following sites match `bg-plate (+ border + rounded)` but are **other
primitives** and are deliberately **NOT migrated by this unit**. This task records the routing so a reviewer
can confirm nothing was force-fit into `Card`, and adds one guard-side test that `Card`'s consumer
`className` never reintroduces appearance utilities (so the primitive can't silently regress
closed-look/open-placement).

**Routing table (do NOT migrate to Card — owned by the named unit):**
- **Inputs / value wells → Input primitive:** `ScopingOrderField.tsx:21`, `ScopingValueCell.tsx:50`, `CorpusTopbar.tsx:177` (select pill).
- **Buttons / toggles → Button / SegmentedControl:** `ScopingLooseMatches.tsx:60`, `PlotCard.tsx:511`, `CorpusTopbar.tsx:299` (notes-toggle).
- **Tag pills / chips → Chip/Tag primitive:** `LoupeSidebar.tsx:118,135`, `ContactSheetRow.tsx:422,439`, `SamplesPage.tsx:215`.
- **Banners → Phase 3 status-icon+label banner:** `StaleIndicesBanner.tsx:66`, `InfrastructureBanner.tsx:45-46`.
- **Modal/dialog shells → ModalShell (Phase 2 sibling):** `NavModal.tsx:227`, `ConflictModalShell.tsx:85`, `OnboardingFlow.tsx:226,325`, `SpeculativeBuilder.tsx:124` (note: these carry `shadow-2xl`/`shadow-xl` Tailwind defaults — Flat-Except-the-Plate violations fixed *inside* ModalShell, not here).
- **Segmented container → SegmentedControl (Phase 1):** `CorpusTopbar.tsx:213`.
- **Floaters / tooltips / ink toolbar → their own primitive or sanctioned floating shadow (out of Card scope):** `OffsetDock.tsx:20`, `FigureExportControls.tsx:172-176` (currently doubles `.card` + redundant `shadow-lg` + `bg-plate` — its de-dup is owned by whichever unit takes floaters, not Card), `MemberMetaRow.tsx:449` (dropdown), `MultiTracePlot.tsx:679` (tooltip), `CullBar.tsx:26` (`bg-ink` floating action bar — explicitly not a plate).

**Files:**
- Test: `packages/HimalayaUI/frontend/test/ui/Card.test.tsx` (extend)

- [ ] **Step 1: Write the boundary/closed-look test** — append to `test/ui/Card.test.tsx`:
```tsx
describe("Card closed-look contract", () => {
  it("never emits its own appearance via consumer className — flat carries no shadow utility", () => {
    render(<Card data-testid="flat">x</Card>);
    const cls = screen.getByTestId("flat").className;
    expect(cls).toContain("bg-plate");
    expect(cls).toContain("rounded-md");
    expect(cls).not.toMatch(/shadow/);
    expect(screen.getByTestId("flat")).not.toHaveAttribute("data-elevated");
  });

  it("elevated delegates appearance to the `.card` recipe (no inline shadow, no bg-plate util, no rounded util)", () => {
    render(<Card data-testid="lift" elevated>x</Card>);
    const cls = screen.getByTestId("lift").className;
    expect(cls).toContain("card");
    expect(cls).not.toMatch(/shadow-\[/);
    expect(cls).not.toContain("bg-plate");
    expect(cls).not.toContain("rounded-md");
    expect(screen.getByTestId("lift")).toHaveAttribute("data-elevated", "true");
  });
});
```
> This is the ONE permitted class-string assertion in this unit — it pins the *primitive's own internal*
> appearance contract (single source of the shadow), not a consumer's look. It guards that `elevated` keeps
> appearance inside `.card` and flat keeps it off the shadow path.

- [ ] **Step 2: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run test/ui/Card.test.tsx` → Expected: PASS (all Card cases).

- [ ] **Step 3: Commit** — `git add packages/HimalayaUI/frontend/test/ui/Card.test.tsx && git commit -m "Lock Card closed-look contract (flat=no shadow, elevated=.card single source)"`

---

### Task 59: Full-suite + build gate for the Card unit

**Files:** (none — verification only)

- [ ] **Step 1: Run the full unit suite** — `cd packages/HimalayaUI/frontend && npm test` → Expected: PASS (no regressions; all migrated `data-testid`s resolve). Capture output to a file and grep for failures rather than scrolling (suite is large).
- [ ] **Step 2: Build gate** — `cd packages/HimalayaUI/frontend && npm run build` (= `lint:design` + `tsc --noEmit -p tsconfig.build.json` + `vite build`) → Expected: PASS. The `lint:design` baseline must show **no new** `bg-plate`/`rounded-*` violations introduced (the migrations should *reduce* the count); if `check-design.mjs` reports newly-cleared violations still pinned in `design-baseline.json`, that baseline drain is owned by the Phase 3 ratchet unit — do not edit the baseline here unless `lint:design` *fails* on a now-absent violation, in which case remove only those entries and note it in the commit.
- [ ] **Step 3: Final commit (if baseline touched)** — `git add packages/HimalayaUI/frontend/scripts/design-baseline.json && git commit -m "Drop design-baseline entries cleared by Card migration"` (skip if untouched).


## Phase 3 — Sweep & per-rule ratchet to enforcing

Depends on all prior waves. Drains the remaining `text-[`/`rounded-[` backlog onto the scale and hardens each guard rule to absolute as its baseline reaches zero.

> **Scope.** Drain the two largest guard backlogs — `text-[` (93 sites across 24 files) and `rounded-[`
> (8 sites across 5 files) — onto the named scale roles and the real `--radius-*` scale, then flip each
> guard rule to absolute as its baseline empties, and finally delete `design-baseline.json` + the
> allowance code. This unit assumes Phase 0 has already landed: the `@theme` tokens (`--radius-sm`/`-md`
> both **5px**, `--radius-full`, `--color-scrim`, `--color-print-accent: var(--color-accent)`), the
> `scripts/check-design.mjs` guard, the seed `scripts/design-baseline.json`, the `lint:design` build step,
> and the side-stripe fix (Toast + InfrastructureBanner) with rule #4 already at zero. It also assumes the
> rule #3 (raw-color-utility) and rule #5 (raw-color-literal/scrim) backlogs were drained by the
> token-adoption call sites in Phase 0/2 (NavModal/ConflictModalShell/OnboardingFlow scrim, FocusNotesMargin
> color-via-`text-[`). **Verify the rule #3/#5 baseline is already empty at the start of this unit (Step 0
> of the first ratchet-flip task); if not, those residual sites must be drained first — see the
> cross-unit dependency note.**

> **Locked values used here (from the brief, authoritative over the inventory):**
> scale = `text-xs` 10.5 / `text-sm` 11.5 / `text-base` 13 / `text-lg` 15 / `text-xl` 19 / `text-display`
> 27 (px). Radius scale `rounded-sm` = `rounded-md` = **5px**, `rounded-full` = 9999px. The inventory's
> `--radius-md: 7px` is superseded; every `rounded-[7px]` therefore snaps to `rounded-md` (5px) — a
> deliberate, reviewed −2px change.

> **Off-scale policy (default = SNAP).** Each off-scale `text-[Npx]` is mapped to the nearest named role
> (round to closest step; ties round up). Snapping is the default because the brief locks the 5-step scale
> and authorizes no new size tokens. The ONE genuinely distinct typographic size — `PhasePanel.tsx:43`
> `text-[23px]` (the serif phase title, between `text-xl` 19 and `text-display` 27) — is the single
> reviewed scale candidate: snap it to `text-xl` unless the redesign reviewer elects to add a
> `--text-title: 23px` step (a one-line `@theme` add + `.text-title` role). Default in this plan: **snap to
> `text-xl`**; the extend option is called out inline at that site for the reviewer.

> **The migration mechanic — closed-look, placement-only className.** Every edit replaces an arbitrary
> bracket utility in a *consumer* className with a named scale class. No primitive is touched (those live
> under `src/components/ui/**`, which the guard excludes). Vitest asserts rendered semantics (text, role,
> `data-*`), never class strings, so these pure-className swaps need **no new tests** — the regression net
> is (a) the existing unit + E2E suites still passing, (b) `npm run build`'s `lint:design` step proving the
> baseline shrank, and (c) the guard's own subset-diff CI assertion (`baseline_new ⊆ baseline_old`). Run
> `npm test` once at the end of each file group, and `npm run build` after every baseline edit.

> **Snap table (memorize once, apply everywhere).** `8.5/9/9.5/10/10.5 → text-xs` (10.5) · `11/11.5 →
> text-sm` (11.5) · `12.5/13/13.5 → text-base` (13) · `15 → text-lg` · `23 → text-xl` (see policy above).
> Always preserve every other token on the line and the element's `data-testid`.

---

### Task 60: Sweep `text-[` — Scoping cluster (10 files, 27 sites)

This task drains the Scoping* component family. The transform pattern, applied once and repeated per site:

```
- find the call site at `<file>:<line>` (verify the line — the inventory `loc` predates this unit's
  edits and adjacent edits shift lines; grep the file for the exact `text-[Npx]` token instead of trusting
  the number).
- replace ONLY the `text-[Npx]` token with its snap-table role; leave the rest of the className verbatim.
  e.g.  `className="text-[10.5px] font-medium text-ink-soft"`  →  `className="text-xs font-medium text-ink-soft"`
- after each FILE is fully converted, grep the file for `text-[` → expect zero matches.
```

**Files:**
- Modify: `src/components/ScopingFoot.tsx` (:25 `12.5`→`text-base`, :33 `10.5`→`text-xs`)
- Modify: `src/components/ScopingLooseMatches.tsx` (:31 `10.5`→xs, :35 `11.5`→sm, :50 `12.5`→base, :51 `11`→sm, :60 `11.5`→sm, :71 `11.5`→sm)
- Modify: `src/components/ScopingOrderField.tsx` (:16 `10.5`→xs, :23 `15`→`text-lg`, :25 `11`→sm)
- Modify: `src/components/ScopingPhaseStrip.tsx` (:26 `10.5`→xs, :39 `11.5`→sm)
- Modify: `src/components/ScopingRow.tsx` (:32 `13`→base, :33 `10.5`→xs)
- Modify: `src/components/ScopingValueCell.tsx` (:50 `13`→base, :65 `13`→base, :74 `9`→`text-xs`)
- Modify: `src/components/SeriesScopingPage.tsx` (:238 `11`→sm, :293 `10.5`→xs, :303 `11`→sm, :308 `10.5`→xs)
- Test: none new — assert via existing suites + `lint:design`.

- [ ] **Step 1: Convert `ScopingFoot.tsx`** — grep `rg -n 'text-\['` it (expect 2), snap both per the table, re-grep → expect 0.
- [ ] **Step 2: Convert `ScopingLooseMatches.tsx`** — 6 sites; snap each per the table; re-grep → 0.
- [ ] **Step 3: Convert `ScopingOrderField.tsx`** — 3 sites (note :23 → `text-lg`); re-grep → 0.
- [ ] **Step 4: Convert `ScopingPhaseStrip.tsx`** — 2 sites; re-grep → 0. (NOTE: if Phase-1 PhaseStrip primitive already absorbed this file's phase-cell markup, the `text-[` may be gone already — grep first; if 0, skip and remove its baseline hashes in Step 8 anyway.)
- [ ] **Step 5: Convert `ScopingRow.tsx`** — 2 sites; re-grep → 0.
- [ ] **Step 6: Convert `ScopingValueCell.tsx`** — 3 sites; re-grep → 0.
- [ ] **Step 7: Convert `SeriesScopingPage.tsx`** — 4 sites; re-grep → 0.
- [ ] **Step 8: Lower the baseline** — remove exactly the 22 corresponding `(rule="text-[", text=…)` hashes from `scripts/design-baseline.json`. Regenerate the hashes the same way the guard does: for each converted line's original violation text, recompute its normalized hash (run `node scripts/check-design.mjs --print-hashes` if the guard exposes it, else delete the entries whose normalized text matches the removed `text-[Npx]` tokens). Then `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (`lint:design` confirms `baseline_new ⊆ baseline_old` and no new violations).
- [ ] **Step 9: Run the affected unit tests** — `cd packages/HimalayaUI/frontend && npx vitest run test/scoping.test.tsx test/ScopingFoot.test.tsx` → Expected: PASS (these assert text/`data-*`, unaffected by class swaps).
- [ ] **Step 10: Commit** — `git add packages/HimalayaUI/frontend/src/components/Scoping*.tsx packages/HimalayaUI/frontend/src/components/SeriesScopingPage.tsx packages/HimalayaUI/frontend/scripts/design-baseline.json && git commit -m "sweep(text-scale): Scoping cluster off text-[ onto named roles; shrink guard baseline"`

---

### Task 61: Sweep `text-[` — Loupe + Contact-sheet cluster (5 files, 31 sites)

Same transform pattern as the Scoping task. This is the heaviest single group (LoupeSidebar alone has 17).

**Files:**
- Modify: `src/components/LoupeSidebar.tsx` (17 sites — see snap notes below)
- Modify: `src/components/LoupeFrame.tsx` (:56 `10`→`text-xs`, :62 `11`→`text-sm`)
- Modify: `src/components/LoupePage.tsx` (:184 `11.5`→sm)
- Modify: `src/components/ContactSheetRow.tsx` (:126 `8.5`→xs, :348 `13.5`→base, :351 `10.5`→xs, :400 `10`→xs, :422 `10.5`→xs, :450 `10.5`→xs, :465 `10.5`→xs, :472 `10.5`→xs, :494 `10.5`→xs)
- Modify: `src/components/LoupeSidebar.tsx` rounded swatch is handled in the rounded task, not here.
- Test: none new.

LoupeSidebar snap map (all from the table): :30 `11.5`→sm, :61 `10`→xs, :107 `10.5`→xs, :118 `10.5`→xs,
:146 `10.5`→xs, :161 `10.5`→xs, :168 `10.5`→xs, :189 `10.5`→xs, :249 `11.5`→sm, :270 `13`→base,
:273 `10.5`→xs, :283 `11.5`→sm, :291 `10`→xs, :315 `11.5`→sm, :323 `11.5`→sm, :337 `11.5`→sm, :356 `11`→sm.

- [ ] **Step 1: Convert `LoupeSidebar.tsx`** — 17 sites; grep `text-\[` first (expect 17), snap each per the map, re-grep → 0. (Many lines repeat `text-[10.5px]`; do them top-down so line numbers don't drift under you, or grep-and-replace each unique surrounding context.)
- [ ] **Step 2: Convert `LoupeFrame.tsx`** — 2 sites; re-grep → 0.
- [ ] **Step 3: Convert `LoupePage.tsx`** — 1 site; re-grep → 0.
- [ ] **Step 4: Convert `ContactSheetRow.tsx` (text only)** — 9 `text-[` sites; re-grep `text-\[` → 0 (a `rounded-[` will remain at :108 — that is the rounded task; do NOT touch it here).
- [ ] **Step 5: Lower the baseline** — remove the 29 corresponding `text-[` hashes from `scripts/design-baseline.json`; `npm run build` → Expected: PASS.
- [ ] **Step 6: Run affected unit + E2E-backing tests** — `cd packages/HimalayaUI/frontend && npx vitest run test/` (filter to loupe/contact specs if present) → Expected: PASS.
- [ ] **Step 7: Commit** — `git add packages/HimalayaUI/frontend/src/components/Loupe*.tsx packages/HimalayaUI/frontend/src/components/ContactSheetRow.tsx packages/HimalayaUI/frontend/scripts/design-baseline.json && git commit -m "sweep(text-scale): Loupe + ContactSheetRow off text-[ onto named roles; shrink baseline"`

---

### Task 62: Sweep `text-[` — Series + Phase + misc cluster (9 files, 35 sites)

Same transform pattern.

**Files:**
- Modify: `src/components/PhasePanel.tsx` (9 sites — incl. the reviewed serif title)
- Modify: `src/components/SeriesBuilderPage.tsx` (8 sites)
- Modify: `src/components/SeriesFolioCard.tsx` (5 sites)
- Modify: `src/components/SeriesFolioPage.tsx` (5 sites)
- Modify: `src/components/SeriesPhaseStrip.tsx` (:46 `11`→sm)
- Modify: `src/components/SamplesPage.tsx` (:128 `10.5`→xs, :207 `11.5`→sm, :215 `10.5`→xs)
- Modify: `src/components/CorpusTopbar.tsx` (:207 `11.5`→sm, :270 `10`→xs)
- Modify: `src/components/CullBar.tsx` (:28 `12.5`→base, :42 `9.5`→xs, :52 `9.5`→xs)
- Modify: `src/components/OffsetDock.tsx` (:22 `10`→xs)
- Modify: `src/components/DetectorImage.tsx` (:217 `11`→sm) — NOTE: DetectorImage is rule-#3 allowlisted, but the allowlist is for *color* utilities only; `text-[11px]` is a **size** (rule #1), which has NO allowlist, so this MUST be converted.
- Modify: `src/components/FocusReflectionsTable.tsx` (:175 `9.5`→xs, :208 `11`→sm; leave the adjacent `mt-[11px] pt-[10px]` spacing utilities — out of scope, not guard-banned)
- Modify: `src/components/SampleStatusChip.tsx` (:22 `11`→sm)
- Test: none new.

PhasePanel snap map: :43 `23`→**`text-xl`** (reviewed scale candidate — see off-scale policy; if reviewer
elects `--text-title: 23px`, swap to `text-title` here and add the `@theme` step + `.text-title` role in
`styles.css`), :52 `11`→sm, :65 `10.5`→xs, :141 `13`→base, :144 `10.5`→xs, :152 `13`→base, :185 `10.5`→xs,
:318 `10`→xs, :352 `11`→sm. (The `:65`/`:152` inline phase score bars themselves are migrated to `ScoreBar
size` in the ScoreBar-adoption task of another unit; here only their *text* labels' `text-[` are snapped.)

SeriesBuilderPage snap: :293 `11`→sm, :297/300/303/310/318 `10.5`→xs, :326 `11`→sm, :376 `10.5`→xs.
SeriesFolioCard snap: :119 `10`→xs, :129 `10.5`→xs, :137 `10`→xs, :170 `11.5`→sm, :189 `10.5`→xs.
SeriesFolioPage snap: :71 `11`→sm, :87 `10.5`→xs, :120 `11.5`→sm, :130 `10.5`→xs, :141 `11.5`→sm.

- [ ] **Step 1: Convert `PhasePanel.tsx`** — 9 sites per the map (apply the `:43` reviewed decision); grep `text-\[` → 0.
- [ ] **Step 2: Convert `SeriesBuilderPage.tsx`** — 8 sites; re-grep → 0.
- [ ] **Step 3: Convert `SeriesFolioCard.tsx`** — 5 sites; re-grep → 0.
- [ ] **Step 4: Convert `SeriesFolioPage.tsx`** — 5 sites; re-grep → 0.
- [ ] **Step 5: Convert `SeriesPhaseStrip.tsx`** — 1 site; re-grep `text-\[` → 0 (a `rounded-[1.5px]` remains at :39 — handled in the rounded task; if Phase-1 PhaseStrip primitive already absorbed this file, grep first and skip).
- [ ] **Step 6: Convert the misc singles** — `SamplesPage.tsx` (3), `CorpusTopbar.tsx` (2), `CullBar.tsx` (3 text; leave the `rounded-[` at :25/:38/:48 for the rounded task), `OffsetDock.tsx` (1), `DetectorImage.tsx` (1), `FocusReflectionsTable.tsx` (2), `SampleStatusChip.tsx` (1). Re-grep each `text-\[` → 0.
- [ ] **Step 7: Lower the baseline** — remove the remaining ~42 `text-[` hashes from `scripts/design-baseline.json`. After this, `node scripts/check-design.mjs` should report **zero** `rule="text-["` baseline entries. `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS.
- [ ] **Step 8: Confirm rule #1 backlog is drained** — `cd packages/HimalayaUI/frontend && rg -n 'text-\[' src --glob '!src/components/ui/**'` → Expected: zero matches. This is the precondition for the next task's rule flip.
- [ ] **Step 9: Run the unit suite** — `cd packages/HimalayaUI/frontend && npm test` → Expected: PASS (all class-only swaps; semantics unchanged).
- [ ] **Step 10: Commit** — `git add packages/HimalayaUI/frontend/src/components/PhasePanel.tsx packages/HimalayaUI/frontend/src/components/Series*.tsx packages/HimalayaUI/frontend/src/components/SamplesPage.tsx packages/HimalayaUI/frontend/src/components/CorpusTopbar.tsx packages/HimalayaUI/frontend/src/components/CullBar.tsx packages/HimalayaUI/frontend/src/components/OffsetDock.tsx packages/HimalayaUI/frontend/src/components/DetectorImage.tsx packages/HimalayaUI/frontend/src/components/FocusReflectionsTable.tsx packages/HimalayaUI/frontend/src/components/SampleStatusChip.tsx packages/HimalayaUI/frontend/scripts/design-baseline.json && git commit -m "sweep(text-scale): Series/Phase/misc off text-[; rule #1 backlog drained to zero"`

---

### Task 63: Flip rule #1 (`text-[`) to absolute

With the `text-[` baseline empty (verified in the previous task's Step 8), harden the rule so any *future*
`text-[` is blocked unconditionally — even though other rules may still carry backlog.

**Files:**
- Modify: `scripts/check-design.mjs` (the per-rule allowance branch for rule #1)
- Modify: `scripts/design-baseline.json` (assert no `text-[` entries remain)
- Test: `scripts/check-design.test.mjs` (the guard's own fixture suite from Phase 0)

- [ ] **Step 1: Write the failing test** — add a fixture case asserting a `text-[10px]` violation is reported as a **hard error (exit 2)** regardless of baseline contents:
  ```js
  // scripts/check-design.test.mjs (append)
  test("rule #1 (text-[) is absolute: a text-[ violation errors even if baseline-listed", () => {
    const baseline = { entries: [{ rule: "text-[", text: "text-[10px]" }] }; // stale baseline entry
    const result = runGuardOnSource('<div className="text-[10px]" />', { baseline });
    // after hardening, rule #1 ignores the baseline allowance entirely:
    expect(result.exitCode).toBe(2);
    expect(result.violations.some(v => v.rule === "text-[")).toBe(true);
  });
  ```
  (Match the existing test harness's actual helper names — read `scripts/check-design.test.mjs` first and reuse its `runGuardOnSource`/fixture shape; the above is the intent.)
- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run scripts/check-design.test.mjs` → Expected: FAIL (rule #1 currently consults the baseline allowance, so a baseline-listed hash is suppressed → exit 0).
- [ ] **Step 3: Implement** — in `check-design.mjs`, change rule #1's evaluation so it never consults the baseline allowance: emit every `text-[` hit straight to the hard-error channel (the same channel net-new violations use), bypassing the `isInBaseline(hash)` suppression that the still-soft rules use. Concretely, gate the suppression on a per-rule `absolute` flag and set `absolute: true` on rule #1's descriptor; the dispatch becomes `if (!rule.absolute && isInBaseline(hash)) continue;`.
- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run scripts/check-design.test.mjs` → Expected: PASS.
- [ ] **Step 5: Sanity-build** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (no `text-[` exists in `src`, so absolute rule #1 finds nothing).
- [ ] **Step 6: Commit** — `git add packages/HimalayaUI/frontend/scripts/check-design.mjs packages/HimalayaUI/frontend/scripts/check-design.test.mjs && git commit -m "guard(harden): rule #1 text-[ is now absolute (baseline drained)"`

---

### Task 64: Sweep `rounded-[` onto the radius scale (5 files, 8 sites)

Transform pattern: replace each `rounded-[Npx]` with a scale class. With `--radius-sm` = `--radius-md` =
**5px**, both `rounded-sm` and `rounded-md` render 5px; choose `rounded-sm` for chips/swatches/buttons and
`rounded-md` for plates/cards by *intent* (they are pixel-identical today but keep semantic meaning for any
future scale split). The two decorative sub-px phase-strip cells (`rounded-[1.5px]`) and the concentric
inset thumbnails (`rounded-[3px]`) are NOT allowlisted by the spec — they snap to `rounded-sm` (the brief
adds no rounded allowlist; rule #2 takes no file allowlist). The nested-radius reminder in MEMORY applies
to `ThumbnailGallery` chip recompute, but here the inventory's `ThumbnailGallery.tsx:88 rounded-[3px]` simply
snaps to `rounded-sm` (5px) per the no-allowlist decision — note this is a +2px change to that inset frame.

**Files:**
- Modify: `src/components/ThumbnailGallery.tsx` (:88 `rounded-[3px]`→`rounded-sm`)
- Modify: `src/components/LoupeSidebar.tsx` (:50 `rounded-[1px]`→`rounded-sm`)
- Modify: `src/components/CullBar.tsx` (:25 `rounded-[10px]`→`rounded-md`, :38 `rounded-[7px]`→`rounded-md`, :48 `rounded-[7px]`→`rounded-md`)
- Modify: `src/components/ScopingPhaseStrip.tsx` (:34 `rounded-[1.5px]`→`rounded-sm`) — skip if Phase-1 PhaseStrip primitive already absorbed it (grep first).
- Modify: `src/components/ContactSheetRow.tsx` (:108 `rounded-[3px]`→`rounded-sm`)
- Modify: `src/components/SeriesPhaseStrip.tsx` (:39 `rounded-[1.5px]`→`rounded-sm`) — skip if absorbed by the PhaseStrip primitive.
- Test: none new — assert via existing suites + `lint:design`.

- [ ] **Step 1: Convert `ThumbnailGallery.tsx:88`** — grep `rounded-\[` (expect 1), `rounded-[3px]`→`rounded-sm`; re-grep → 0.
- [ ] **Step 2: Convert `LoupeSidebar.tsx:50`** — `rounded-[1px]`→`rounded-sm`; re-grep `rounded-\[` → 0.
- [ ] **Step 3: Convert `CullBar.tsx`** — 3 sites (`:25`→`rounded-md`, `:38`/`:48`→`rounded-md`); re-grep → 0.
- [ ] **Step 4: Convert `ScopingPhaseStrip.tsx:34`** — grep first; if present, `rounded-[1.5px]`→`rounded-sm`; re-grep → 0.
- [ ] **Step 5: Convert `ContactSheetRow.tsx:108`** — `rounded-[3px]`→`rounded-sm`; re-grep `rounded-\[` → 0.
- [ ] **Step 6: Convert `SeriesPhaseStrip.tsx:39`** — grep first; if present, `rounded-[1.5px]`→`rounded-sm`; re-grep → 0.
- [ ] **Step 7: Confirm rule #2 backlog drained** — `cd packages/HimalayaUI/frontend && rg -n 'rounded-\[' src --glob '!src/components/ui/**'` → Expected: zero matches.
- [ ] **Step 8: Lower the baseline** — remove all 8 `(rule="rounded-[", …)` hashes from `scripts/design-baseline.json`; `npm run build` → Expected: PASS.
- [ ] **Step 9: Run the unit suite** — `cd packages/HimalayaUI/frontend && npm test` → Expected: PASS.
- [ ] **Step 10: Commit** — `git add packages/HimalayaUI/frontend/src/components/ThumbnailGallery.tsx packages/HimalayaUI/frontend/src/components/LoupeSidebar.tsx packages/HimalayaUI/frontend/src/components/CullBar.tsx packages/HimalayaUI/frontend/src/components/ScopingPhaseStrip.tsx packages/HimalayaUI/frontend/src/components/ContactSheetRow.tsx packages/HimalayaUI/frontend/src/components/SeriesPhaseStrip.tsx packages/HimalayaUI/frontend/scripts/design-baseline.json && git commit -m "sweep(radius): drain rounded-[ onto rounded-sm/md scale; rule #2 backlog zero"`

---

### Task 65: Flip rule #2 (`rounded-[`) to absolute

**Files:**
- Modify: `scripts/check-design.mjs` (set `absolute: true` on rule #2's descriptor)
- Test: `scripts/check-design.test.mjs`

- [ ] **Step 1: Write the failing test** — append a fixture asserting `rounded-[3px]` errors hard (exit 2) even if a stale baseline entry lists it:
  ```js
  test("rule #2 (rounded-[) is absolute after backlog drain", () => {
    const baseline = { entries: [{ rule: "rounded-[", text: "rounded-[3px]" }] };
    const result = runGuardOnSource('<div className="rounded-[3px]" />', { baseline });
    expect(result.exitCode).toBe(2);
    expect(result.violations.some(v => v.rule === "rounded-[")).toBe(true);
  });
  ```
- [ ] **Step 2: Run it, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run scripts/check-design.test.mjs` → Expected: FAIL (rule #2 still soft; baseline entry suppresses it).
- [ ] **Step 3: Implement** — set `absolute: true` on rule #2's descriptor in `check-design.mjs` (the `!rule.absolute && isInBaseline(hash)` gate added in the rule-#1 flip already handles the rest).
- [ ] **Step 4: Run, expect pass** — same command → Expected: PASS.
- [ ] **Step 5: Sanity-build** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS.
- [ ] **Step 6: Commit** — `git add packages/HimalayaUI/frontend/scripts/check-design.mjs packages/HimalayaUI/frontend/scripts/check-design.test.mjs && git commit -m "guard(harden): rule #2 rounded-[ is now absolute"`

---

### Task 66: Flip rules #3, #4, #5 to absolute (close out remaining soft rules)

Rule #4 (side-stripe) was drained in Phase 0; rules #3 (raw-color-utility) and #5 (raw-color-literal/scrim)
were drained by the token-adoption call sites in Phases 0/2. This task confirms each is empty and flips it
absolute. **Run Step 0 first; if any rule still carries baseline entries, that rule's residual sites must
be drained before flipping it (see cross-unit dependency note) — flip only the empty ones and leave the
others soft until their owning unit lands.**

**Files:**
- Modify: `scripts/check-design.mjs` (set `absolute: true` on rules #3, #4, #5)
- Test: `scripts/check-design.test.mjs`

- [ ] **Step 0: Confirm each rule's baseline is empty** — `cd packages/HimalayaUI/frontend && node scripts/check-design.mjs --print-baseline-rules` (or read `design-baseline.json`) → Expected: zero entries for rules #3, #4, #5. Also `rg -n 'border-(l|r)-(4|[2-9]|\[)' src --glob '!src/components/ui/**'` → zero (rule #4 surface clean).
- [ ] **Step 1: Write the failing tests** — append three fixtures: `border-l-4` (rule #4), a raw `bg-[oklch(0.05_0_0/0.65)]` color utility (rule #3), and a `style={{ background: "oklch(0.05 0 0 / 0.65)" }}` literal (rule #5) each error hard (exit 2) with a stale baseline entry present:
  ```js
  test("rule #4 (side-stripe) absolute", () => {
    const baseline = { entries: [{ rule: "border-l", text: "border-l-4" }] };
    expect(runGuardOnSource('<div className="border-l-4" />', { baseline }).exitCode).toBe(2);
  });
  test("rule #3 (raw color utility) absolute", () => {
    const baseline = { entries: [{ rule: "color-util", text: "bg-[oklch(0.05_0_0/0.65)]" }] };
    expect(runGuardOnSource('<div className="bg-[oklch(0.05_0_0/0.65)]" />', { baseline }).exitCode).toBe(2);
  });
  test("rule #5 (raw color literal) absolute", () => {
    const baseline = { entries: [{ rule: "color-literal", text: 'background: "oklch(0.05 0 0 / 0.65)"' }] };
    expect(runGuardOnSource('<div style={{ background: "oklch(0.05 0 0 / 0.65)" }} />', { baseline }).exitCode).toBe(2);
  });
  ```
  (Match the guard's actual rule identifiers — read the descriptor table in `check-design.mjs` and use its real `rule` keys instead of the placeholders above.)
- [ ] **Step 2: Run, expect fail** — `cd packages/HimalayaUI/frontend && npx vitest run scripts/check-design.test.mjs` → Expected: FAIL (rules still soft).
- [ ] **Step 3: Implement** — set `absolute: true` on the descriptors for rules #3, #4, #5 in `check-design.mjs`.
- [ ] **Step 4: Run, expect pass** — same command → Expected: PASS.
- [ ] **Step 5: Sanity-build** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS.
- [ ] **Step 6: Commit** — `git add packages/HimalayaUI/frontend/scripts/check-design.mjs packages/HimalayaUI/frontend/scripts/check-design.test.mjs && git commit -m "guard(harden): rules #3/#4/#5 now absolute — all rules enforcing"`

---

### Task 67: Retire `design-baseline.json` and the baseline-allowance machinery

With all five rules absolute, the baseline is dead weight. Delete the file and excise the allowance code so
the guard is a pure absolute linter.

**Files:**
- Delete: `scripts/design-baseline.json`
- Modify: `scripts/check-design.mjs` (remove baseline load, `isInBaseline`, the `--print-*baseline*` flags, the `absolute` flag plumbing, and the subset-diff CI assertion)
- Modify: `scripts/check-design.test.mjs` (drop the baseline-subset and stale-baseline fixtures; keep the pure pattern fixtures — `text-[10px]` flagged, `text-base` clean, `rounded-[3px]` flagged, `rounded-md` clean, `style={{ color }}` clean, var-only clean, `border-l-4` flagged, `shadow-[…rgba…]` clean, `ui/` excluded, allowlisted `bg-[oklch]` clean)
- Modify: none in `package.json` — `lint:design` and the `build` wiring stay (the guard still runs; it just no longer reads a baseline)

- [ ] **Step 1: Confirm a clean tree** — `cd packages/HimalayaUI/frontend && npm run build` → Expected: PASS (proves nothing depends on the baseline allowance for a pass).
- [ ] **Step 2: Update the guard tests first (TDD)** — edit `scripts/check-design.test.mjs` to remove baseline-coupled fixtures and to call the guard with no `baseline` arg; the remaining fixtures must still pin pure-pattern behavior. Run `cd packages/HimalayaUI/frontend && npx vitest run scripts/check-design.test.mjs` → Expected: FAIL (the guard still requires/loads a baseline that the tests no longer pass).
- [ ] **Step 3: Implement** — in `check-design.mjs`: delete the baseline file read, `isInBaseline`, the `absolute` flag (now every rule is unconditionally hard), the subset-diff assertion, and the `--print-*` baseline flags. Every matched pattern → exit 2. Keep the `src/components/ui/**` exclusion and rule #3's file allowlist (`src/phases.ts`, `src/lib/comparison/coloring.ts`, `src/lib/figure-export/**`, `MemberHeatmapLayer.tsx`, `DetectorImage.tsx`, `FocusDetectorPanel.tsx`, `src/main.tsx`).
- [ ] **Step 4: Run, expect pass** — `cd packages/HimalayaUI/frontend && npx vitest run scripts/check-design.test.mjs` → Expected: PASS.
- [ ] **Step 5: Delete the baseline file** — `git rm packages/HimalayaUI/frontend/scripts/design-baseline.json`.
- [ ] **Step 6: Full build + suite gate** — `cd packages/HimalayaUI/frontend && npm run build && npm test` → Expected: PASS (guard runs absolute, finds nothing; all units green).
- [ ] **Step 7: Final commit** — `git add packages/HimalayaUI/frontend/scripts/check-design.mjs packages/HimalayaUI/frontend/scripts/check-design.test.mjs && git commit -m "guard(final): retire design-baseline.json + allowance machinery; check-design is now pure absolute"`

---
