/**
 * Client-side `contentHash` (Plan §Phase 3, Task 3.4).
 *
 * SHA-256 over a canonical JSON serialization of the comparison's identity:
 * `(title, description, members)` with members in `display_order ASC`
 * order. Each member's tuple includes `display_order` itself, so the hash
 * is unambiguous when display_order ties exist (spec §Hash memoization).
 *
 * **Canonical form invariant:** every nesting level emits keys in
 * alphabetical order. Array order is preserved (display_order is the
 * sortable field for member ordering; the producer is responsible for
 * passing members in the desired order). Scalars use JSON.stringify
 * which matches JSON3.write on the Julia side for IEEE-754 doubles,
 * integers, strings, and booleans. `null` and `undefined` both serialize
 * to `null`.
 *
 * Mirrors `canonical_json` + `compute_content_hash` in
 * `packages/HimalayaUI/src/comparisons.jl`. The fixture in
 * `test/fixtures/contentHash.fixture.json` pins both implementations to
 * byte-identical output for the same input.
 */

/**
 * Serialize a JS value to a canonical JSON string. Object keys are emitted
 * in alphabetical order at every level; arrays preserve insertion order;
 * `null` and `undefined` both map to `null`.
 *
 * Exported for testability — the cross-language parity test asserts
 * `canonicalJson(input)` equals the bytes the Julia side hashes.
 */
export function canonicalJson(x: unknown): string {
  if (x === null || x === undefined) return "null";
  if (typeof x === "string") return JSON.stringify(x);
  if (typeof x === "number") {
    // JSON.stringify handles finite doubles correctly and produces
    // canonical short forms (e.g. 1.5, 1e21). Match JSON3.write on the
    // Julia side. NaN and Infinity → "null" (JSON.stringify default).
    return JSON.stringify(x);
  }
  if (typeof x === "boolean") return x ? "true" : "false";
  if (Array.isArray(x)) {
    return "[" + x.map(canonicalJson).join(",") + "]";
  }
  if (typeof x === "object") {
    const keys = Object.keys(x as object).sort();
    const parts: string[] = [];
    for (const k of keys) {
      const v = (x as Record<string, unknown>)[k];
      parts.push(JSON.stringify(k) + ":" + canonicalJson(v));
    }
    return "{" + parts.join(",") + "}";
  }
  // Fallback for anything exotic (BigInt, function): JSON.stringify
  // produces "undefined"/"null"; the input shape for content_hash should
  // never include these. If you land here, you passed something the spec
  // doesn't allow.
  return JSON.stringify(x);
}

/**
 * The hash input shape. Mirrors the server's `compute_content_hash`
 * payload Dict — title + description + members. Members must already be
 * sorted by `display_order` ASC (the producer is the
 * `useSaveComparison` flow which iterates Zustand draft slots in
 * order; the dispatcher then reads from DB sorted by
 * `display_order ASC, id ASC`).
 */
export interface ContentHashInput {
  title: string;
  description: string | null;
  members: Array<{
    exposure_id: number | null;
    display_order: number;
    band_height: number;
    y_offset: number;
    normalization: string;
    color_override: string | null;
    label_override: string | null;
    q_window_min: number | null;
    q_window_max: number | null;
    peak_display: unknown;
    snapshot: unknown;
  }>;
}

/**
 * Compute the canonical content_hash. Returns a `sha256:<hex>` string
 * matching the Julia side. Uses `crypto.subtle.digest` (browser-native);
 * Node tests inject `node:crypto.webcrypto` for parity.
 */
export async function contentHash(input: ContentHashInput): Promise<string> {
  const canonical = canonicalJson(input);
  const bytes = new TextEncoder().encode(canonical);
  const digest = await crypto.subtle.digest("SHA-256", bytes);
  const hex = bytesToHex(new Uint8Array(digest));
  return "sha256:" + hex;
}

function bytesToHex(buf: Uint8Array): string {
  let s = "";
  for (let i = 0; i < buf.length; i++) {
    s += buf[i]!.toString(16).padStart(2, "0");
  }
  return s;
}

/**
 * Short form of a comparison's `content_hash`: first 8 hex chars after
 * stripping the `sha256:` prefix, lowercased. The eager-hash fragment of
 * `[[comparison:N@hash8]]` mention tokens.
 *
 * RETAINED through the 2026-05-29 chat retirement.
 * Its only callers (the now-deleted MentionPicker emit path and MentionChip
 * drift detector) are gone, but this is the canonical fix for issue #62 — where
 * two diverging slice conventions made every well-formed mention render a
 * permanent "(changed)" annotation. It is kept as the single anti-drift
 * primitive for the parked comparison-mention vocabulary, so a chat revival
 * cannot silently reintroduce #62.
 */
export function comparisonHash8(content_hash: string): string {
  return content_hash.replace(/^sha256:/, "").slice(0, 8).toLowerCase();
}
