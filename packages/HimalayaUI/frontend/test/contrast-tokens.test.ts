import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";
import { contrastRatio } from "./helpers/contrast";

/*
 * Token-contrast regression guard — pins the F-CONTRAST settled colour roles.
 * ---------------------------------------------------------------------------
 * F-CONTRAST fixed how the de-emphasis text tones map onto WCAG context:
 *
 *   ink        primary text             AAA          (>= 7.0:1)
 *   ink-soft   small / normal TEXT      AA-normal    (>= 4.5:1)  ← the key guard
 *   ink-faint  large / decorative ONLY  AA-large/non (>= 3.0:1)  intentionally < 4.5
 *   accent     accent text + non-text   AA-normal where it carries small text
 *   success    status                   AA-normal where it carries small text
 *   warning    status, large/non-text   AA-large/non ONLY        never normal text
 *
 * This test reads the REAL token oklch() values out of src/styles.css at test
 * time, so darkening/lightening any token there re-runs the check. Thresholds
 * are ROLE FLOORS, not the measured values: AA-normal = 4.5 for small text,
 * AA-large/non-text = 3.0. The measured contrast today is noted inline beside
 * each row so a future tweak can see the headroom it is spending.
 *
 * Assertions compute contrast from the parsed oklch() values — they do NOT
 * inspect Tailwind class strings. This is TEST-ONLY: it guards styles.css, it
 * does not change it.
 */

const __dirname = dirname(fileURLToPath(import.meta.url));
const STYLES_CSS = resolve(__dirname, "../src/styles.css");

// Build name -> oklch-string from the real stylesheet. Regex each
//   --color-<name>: oklch(...)
// declaration. Only matches oklch() values (the design tokens we guard are all
// authored as oklch); colour tokens authored any other way simply won't appear.
function readColorTokens(): Record<string, string> {
  const css = readFileSync(STYLES_CSS, "utf8");
  const re = /--color-([a-z0-9-]+):\s*(oklch\([^)]*\))/gi;
  const out: Record<string, string> = {};
  let m: RegExpExecArray | null;
  while ((m = re.exec(css)) !== null) {
    out[m[1]!] = m[2]!;
  }
  return out;
}

const TOKENS = readColorTokens();

// The names every assertion below depends on. Asserting they were all found
// means a rename in styles.css can't silently skip a contrast check.
const EXPECTED_TOKEN_NAMES = [
  "paper",
  "paper-sunk",
  "plate",
  "ink",
  "ink-soft",
  "ink-faint",
  "accent",
  "success",
  "warning",
] as const;

// foreground token on background token, with the role floor it must clear.
// `measured` is informational only (the value today) — the assertion uses `min`.
interface Pair {
  fg: string;
  bg: string;
  min: number;
  role: string;
  measured: number;
}

const PAIRS: Pair[] = [
  // ink — primary text, AAA (>= 7.0) on every surface.
  { fg: "ink", bg: "paper", min: 7.0, role: "AAA primary text", measured: 14.4 },
  { fg: "ink", bg: "paper-sunk", min: 7.0, role: "AAA primary text", measured: 13.3 },
  { fg: "ink", bg: "plate", min: 7.0, role: "AAA primary text", measured: 15.0 },

  // ink-soft — small/normal informational TEXT, AA-normal (>= 4.5). THE key guard.
  { fg: "ink-soft", bg: "paper", min: 4.5, role: "AA-normal small text", measured: 6.5 },
  { fg: "ink-soft", bg: "paper-sunk", min: 4.5, role: "AA-normal small text", measured: 6.0 },
  { fg: "ink-soft", bg: "plate", min: 4.5, role: "AA-normal small text", measured: 6.8 },

  // ink-faint — LARGE / decorative / non-text ONLY, AA-large/non-text floor (>= 3.0).
  // ink-faint is INTENTIONALLY below AA-normal (4.5): small text uses ink-soft.
  // NOTE: ink-faint on paper-sunk is ~2.92 — BELOW even the 3.0 floor, so faint
  // must be purely decorative on the sunk rail bg; that pair is deliberately NOT
  // asserted here.
  { fg: "ink-faint", bg: "paper", min: 3.0, role: "AA-large/non-text only", measured: 3.16 },
  { fg: "ink-faint", bg: "plate", min: 3.0, role: "AA-large/non-text only", measured: 3.29 },

  // accent (terracotta) — carries small accent text on paper/plate, AA-normal (>= 4.5).
  { fg: "accent", bg: "paper", min: 4.5, role: "AA-normal accent text", measured: 4.77 },
  { fg: "accent", bg: "plate", min: 4.5, role: "AA-normal accent text", measured: 4.97 },
  // accent on the sunk rail bg is ~4.40 — just under AA-normal — so on paper-sunk
  // accent must be non-text/large only; floor here is the AA-large/non-text 3.0.
  { fg: "accent", bg: "paper-sunk", min: 3.0, role: "non-text/large on sunk rail", measured: 4.40 },

  // success (sage) — carries small status text on paper/plate, AA-normal (>= 4.5).
  { fg: "success", bg: "paper", min: 4.5, role: "AA-normal status text", measured: 4.83 },
  { fg: "success", bg: "plate", min: 4.5, role: "AA-normal status text", measured: 5.03 },
  // success on the sunk rail bg is ~4.46 — just under AA-normal — non-text/large only.
  { fg: "success", bg: "paper-sunk", min: 3.0, role: "non-text/large on sunk rail", measured: 4.46 },

  // warning (amber) — AA-large/non-text ONLY everywhere (max ~3.66 < 4.5), so
  // warning must NEVER be small/normal-weight text. Floor is 3.0 on every surface.
  { fg: "warning", bg: "paper", min: 3.0, role: "AA-large/non-text ONLY", measured: 3.51 },
  { fg: "warning", bg: "paper-sunk", min: 3.0, role: "AA-large/non-text ONLY", measured: 3.24 },
  { fg: "warning", bg: "plate", min: 3.0, role: "AA-large/non-text ONLY", measured: 3.66 },
];

describe("contrast tokens — F-CONTRAST settled roles (styles.css)", () => {
  it("found every expected --color token in styles.css", () => {
    for (const name of EXPECTED_TOKEN_NAMES) {
      expect(TOKENS[name], `--color-${name} should be present (as an oklch() value) in styles.css`).toBeTruthy();
    }
  });

  it("every role pair clears its WCAG floor on the real token values", () => {
    for (const { fg, bg, min, role } of PAIRS) {
      const fgVal = TOKENS[fg];
      const bgVal = TOKENS[bg];
      expect(fgVal, `--color-${fg} present`).toBeTruthy();
      expect(bgVal, `--color-${bg} present`).toBeTruthy();
      const ratio = contrastRatio(fgVal!, bgVal!);
      expect(
        ratio,
        `--color-${fg} on --color-${bg} (${role}) must clear ${min}:1 — got ${ratio.toFixed(2)}:1`,
      ).toBeGreaterThanOrEqual(min);
    }
  });
});
