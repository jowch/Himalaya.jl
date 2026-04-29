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
});
