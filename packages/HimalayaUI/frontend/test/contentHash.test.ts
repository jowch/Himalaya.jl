/**
 * Cross-language hash parity test (Plan §Phase 3, Task 3.4).
 *
 * Pins that the JS `contentHash` produces byte-identical SHA-256 output
 * to the Julia `compute_content_hash` for a fixed input. The fixture
 * lives in `test/fixtures/contentHash.fixture.json` and is loaded by
 * BOTH this Vitest test (asserting the JS hash) AND the Julia test
 * `test_comparisons.jl::cross-language hash parity` (asserting the
 * Julia hash) — so a drift in either canonical-serialization step
 * fails the suite on whichever side moved.
 *
 * Why this matters: chat citations like `@comparison:42@<hash8>`
 * eagerly write the current content_hash at insertion time. If client
 * and server canonicalize differently, the @-resolver mismatches and
 * shows "(changed)" indicators on freshly-cited frozen snapshots.
 */
import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join } from "node:path";
import { contentHash, canonicalJson } from "../src/lib/comparison/contentHash";
import type { ContentHashInput } from "../src/lib/comparison/contentHash";
import { webcrypto } from "node:crypto";

// JSDOM doesn't ship `crypto.subtle`. Polyfill from Node's webcrypto so the
// browser-side digest path is exercised in unit tests. (Vite/jsdom env in
// production loads `crypto.subtle` natively in the browser.)
if (typeof globalThis.crypto === "undefined" || !globalThis.crypto.subtle) {
  Object.defineProperty(globalThis, "crypto", {
    value: webcrypto, writable: false, configurable: true,
  });
}

interface Fixture {
  description: string;
  input: {
    title: string;
    description: string | null;
    members: Array<Record<string, unknown>>;
  };
  expected_hash: string;
  expected_canonical: string;
}

function loadFixture(): Fixture {
  const path = join(__dirname, "fixtures", "contentHash.fixture.json");
  return JSON.parse(readFileSync(path, "utf-8"));
}

describe("canonicalJson — alphabetical key ordering at every nesting level", () => {
  it("sorts object keys alphabetically", () => {
    expect(canonicalJson({ b: 2, a: 1 })).toBe('{"a":1,"b":2}');
  });

  it("recurses into nested objects", () => {
    expect(canonicalJson({ z: { b: 2, a: 1 }, a: 1 })).toBe('{"a":1,"z":{"a":1,"b":2}}');
  });

  it("preserves array order", () => {
    expect(canonicalJson([3, 1, 2])).toBe("[3,1,2]");
  });

  it("recurses into objects inside arrays", () => {
    expect(canonicalJson([{ b: 2, a: 1 }])).toBe('[{"a":1,"b":2}]');
  });

  it("emits null for null and undefined", () => {
    expect(canonicalJson(null)).toBe("null");
    expect(canonicalJson(undefined)).toBe("null");
  });

  it("emits scalars via JSON.stringify (matches JSON3.write on the Julia side)", () => {
    expect(canonicalJson(42)).toBe("42");
    expect(canonicalJson(1.5)).toBe("1.5");
    expect(canonicalJson("hi")).toBe('"hi"');
    expect(canonicalJson(true)).toBe("true");
    expect(canonicalJson(false)).toBe("false");
  });

  it("escapes strings the same way JSON.stringify does", () => {
    expect(canonicalJson('a"b')).toBe('"a\\"b"');
  });
});

describe("contentHash — SHA-256 over canonical-JSON, sha256: prefixed lowercase hex", () => {
  it("returns a sha256:<hex> string", async () => {
    const h = await contentHash({ title: "x", description: null, members: [] });
    expect(h).toMatch(/^sha256:[0-9a-f]{64}$/);
  });

  it("ordering of input keys does not matter (canonicalization)", async () => {
    const h1 = await contentHash({ title: "x", description: null, members: [] });
    const h2 = await contentHash({ description: null, members: [], title: "x" } as any);
    expect(h1).toBe(h2);
  });

  it("ordering of nested member fields does not matter", async () => {
    const member1 = {
      exposure_id: 1, display_order: 0, band_height: 1, y_offset: 0,
      normalization: "none", color_override: null, label_override: null,
      q_window_min: null, q_window_max: null, peak_display: null,
      snapshot: { effective_peaks: [], confirmed_index: null,
                  analysis_inputs_hash: "h" },
    };
    const member2 = {
      // all the same fields but in a wildly different insertion order
      snapshot: { analysis_inputs_hash: "h", effective_peaks: [],
                  confirmed_index: null },
      peak_display: null, q_window_max: null, q_window_min: null,
      label_override: null, color_override: null,
      normalization: "none", y_offset: 0, band_height: 1,
      display_order: 0, exposure_id: 1,
    };
    const h1 = await contentHash({ title: "x", description: null, members: [member1] });
    const h2 = await contentHash({ title: "x", description: null, members: [member2] });
    expect(h1).toBe(h2);
  });

  it("hash changes when title changes", async () => {
    const h1 = await contentHash({ title: "A", description: null, members: [] });
    const h2 = await contentHash({ title: "B", description: null, members: [] });
    expect(h1).not.toBe(h2);
  });

  it("hash changes when a member's display_order changes (display_order is part of the canonical form)", async () => {
    const baseMember = {
      exposure_id: 1, band_height: 1, y_offset: 0,
      normalization: "none", color_override: null, label_override: null,
      q_window_min: null, q_window_max: null, peak_display: null,
      snapshot: { effective_peaks: [], confirmed_index: null,
                  analysis_inputs_hash: "h" },
    };
    const h1 = await contentHash({
      title: "x", description: null,
      members: [{ ...baseMember, display_order: 0 }],
    });
    const h2 = await contentHash({
      title: "x", description: null,
      members: [{ ...baseMember, display_order: 5 }],
    });
    expect(h1).not.toBe(h2);
  });
});

describe("contentHash cross-language fixture parity", () => {
  it("matches the expected canonical-JSON byte stream from the fixture", () => {
    const f = loadFixture();
    expect(canonicalJson(f.input)).toBe(f.expected_canonical);
  });

  it("matches the expected SHA-256 hash from the fixture", async () => {
    const f = loadFixture();
    // The fixture's `members` is typed loosely as `Record<string, unknown>[]`
    // (the JSON loader can't infer the precise per-member shape). The fixture
    // data IS structurally a valid ContentHashInput at runtime — assert that
    // contract via the cast so `tsc --noEmit` accepts the call without
    // weakening `ContentHashInput.members` itself.
    const h = await contentHash(f.input as ContentHashInput);
    expect(h).toBe(f.expected_hash);
  });
});
