import { describe, it, expect } from "vitest";
import { parseMentions } from "../src/lib/renderMentions";

describe("parseMentions", () => {
  it("returns single string segment for plain text", () => {
    const result = parseMentions("hello world");
    expect(result).toEqual([{ kind: "text", text: "hello world" }]);
  });

  it("returns a mention token for a single [[type:id]]", () => {
    const result = parseMentions("[[peak:42]]");
    expect(result).toEqual([{ kind: "mention", type: "peak", id: 42 }]);
  });

  it("splits mixed text and tokens correctly", () => {
    const result = parseMentions("shoulder at [[peak:42]] matches [[index:17]]");
    expect(result).toEqual([
      { kind: "text", text: "shoulder at " },
      { kind: "mention", type: "peak", id: 42 },
      { kind: "text", text: " matches " },
      { kind: "mention", type: "index", id: 17 },
    ]);
  });

  it("handles token at start and end", () => {
    const result = parseMentions("[[sample:8]] and [[exposure:33]]");
    expect(result).toEqual([
      { kind: "mention", type: "sample", id: 8 },
      { kind: "text", text: " and " },
      { kind: "mention", type: "exposure", id: 33 },
    ]);
  });

  it("ignores empty text segments between adjacent tokens", () => {
    const result = parseMentions("[[peak:1]][[peak:2]]");
    expect(result).toEqual([
      { kind: "mention", type: "peak", id: 1 },
      { kind: "mention", type: "peak", id: 2 },
    ]);
  });

  it("returns plain text for malformed tokens (unknown type)", () => {
    const result = parseMentions("[[notatype:42]]");
    expect(result.every((s) => s.kind === "text")).toBe(true);
  });

  it("handles empty string", () => {
    expect(parseMentions("")).toEqual([]);
  });

  // ─── Comparison mention syntax (Phase 10) ───────────────────────────────
  // The eager-hash-insertion model means comparison mentions land in chat
  // bodies as `[[comparison:42@a1b2c3d4]]`, where the 8-char hex suffix is
  // the truncated content_hash at compose time. The hash is optional —
  // legacy `[[comparison:42]]` (no hash) must still parse to support
  // bodies authored before the eager-hash policy or for non-comparison
  // types that don't carry a hash.

  it("parses comparison mention with hash suffix", () => {
    const result = parseMentions("[[comparison:42@a1b2c3d4]]");
    expect(result).toEqual([
      { kind: "mention", type: "comparison", id: 42, hash: "a1b2c3d4" },
    ]);
  });

  it("parses comparison mention without hash suffix", () => {
    const result = parseMentions("[[comparison:42]]");
    expect(result).toEqual([
      { kind: "mention", type: "comparison", id: 42 },
    ]);
  });

  it("hash field is undefined for non-comparison types", () => {
    const result = parseMentions("[[peak:1]]");
    // No `hash` key on peak/index/exposure/sample mentions.
    expect(result).toEqual([{ kind: "mention", type: "peak", id: 1 }]);
  });

  it("handles mixed comparison + other mentions in one body", () => {
    const result = parseMentions(
      "see [[comparison:42@a1b2c3d4]] vs [[peak:8]]",
    );
    expect(result).toEqual([
      { kind: "text", text: "see " },
      { kind: "mention", type: "comparison", id: 42, hash: "a1b2c3d4" },
      { kind: "text", text: " vs " },
      { kind: "mention", type: "peak", id: 8 },
    ]);
  });

  it("rejects an obviously-malformed hash (too long)", () => {
    // The lexer is strict — only [0-9a-f]{8} qualifies. Anything else falls
    // through to the legacy literal-text fallback.
    const result = parseMentions("[[comparison:42@abc]]");
    expect(result.every((s) => s.kind === "text")).toBe(true);
  });
});
