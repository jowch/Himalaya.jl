import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { slugifyForFilename, buildFilename } from "../../src/lib/figure-export/filename";

describe("slugifyForFilename", () => {
  it("lowercases and replaces spaces with dashes", () => {
    expect(slugifyForFilename("My Sample 04")).toBe("my-sample-04");
  });
  it("collapses runs of non-alphanumerics into a single dash", () => {
    expect(slugifyForFilename("a / b // c")).toBe("a-b-c");
  });
  it("trims leading and trailing dashes", () => {
    expect(slugifyForFilename("--foo--")).toBe("foo");
  });
  it("returns sentinel 'figure' for empty input", () => {
    expect(slugifyForFilename("")).toBe("figure");
  });
  it("returns sentinel 'figure' for input that slugifies to empty", () => {
    expect(slugifyForFilename("///")).toBe("figure");
  });
  it("preserves digits and basic latin letters", () => {
    expect(slugifyForFilename("AgBe-2026-001")).toBe("agbe-2026-001");
  });
  it("strips diacritics by replacing with dashes (acceptable for filename)", () => {
    // Non-ASCII alphabetic chars are not preserved; the contract is "safe filename",
    // not "Unicode-fidelity". Crash-safe, predictable.
    const out = slugifyForFilename("café");
    expect(out).toMatch(/^[a-z0-9-]+$/);
    expect(out.length).toBeGreaterThan(0);
  });
});

describe("buildFilename", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    // Pin to a known local date. Spec mandates en-CA YYYY-MM-DD output.
    vi.setSystemTime(new Date(2026, 4, 8, 14, 30)); // 2026-05-08 local
  });
  afterEach(() => { vi.useRealTimers(); });

  it("appends '-{YYYY-MM-DD}.{ext}' (en-CA)", () => {
    expect(buildFilename("himalaya-trace-jc23", "png")).toBe("himalaya-trace-jc23-2026-05-08.png");
  });
  it("supports svg extension", () => {
    expect(buildFilename("himalaya-trace-jc23", "svg")).toBe("himalaya-trace-jc23-2026-05-08.svg");
  });

  // Regression: the raw title was interpolated unsanitized — spaces/colons/dots
  // leaked into the filename and could defeat extension detection (and the
  // anchor `download` attribute). buildFilename must slugify the stem.
  it("slugifies a spaced title so the extension is unambiguous", () => {
    expect(buildFilename("himalaya-series-Series by dose", "svg")).toBe(
      "himalaya-series-series-by-dose-2026-05-08.svg",
    );
  });
  it("strips dots and colons from the stem (no mid-name '.' before the ext)", () => {
    const out = buildFilename("himalaya-series-1:2.5 ratio", "png");
    expect(out).toBe("himalaya-series-1-2-5-ratio-2026-05-08.png");
    // Exactly one dot, and it precedes the extension.
    expect(out.split(".").length).toBe(2);
    expect(out.endsWith(".png")).toBe(true);
  });
});
