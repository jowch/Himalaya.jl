import { describe, it, expect } from "vitest";
import { replacePlaceholder } from "../../src/lib/queue/replacePlaceholder";

interface Row { id: number; q: number; label?: string }

describe("replacePlaceholder", () => {
  it("replaces a single negative-id placeholder matching the predicate", () => {
    const list: Row[] = [
      { id: -1, q: 0.1, label: "optimistic" },
      { id:  5, q: 0.2, label: "existing" },
    ];
    const server: Row = { id: 42, q: 0.1, label: "server" };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    expect(out).toEqual([
      { id: 42, q: 0.1, label: "server" },
      { id:  5, q: 0.2, label: "existing" },
    ]);
  });

  it("appends server item when no placeholder matches", () => {
    const list: Row[] = [{ id: 5, q: 0.2 }];
    const server: Row = { id: 42, q: 0.1 };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    expect(out).toEqual([{ id: 5, q: 0.2 }, { id: 42, q: 0.1 }]);
  });

  it("dedupes against a concurrent SSE that already inserted the server id", () => {
    const list: Row[] = [
      { id: -1, q: 0.1, label: "optimistic" },
      { id: 42, q: 0.1, label: "from-sse" },
    ];
    const server: Row = { id: 42, q: 0.1, label: "server" };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    // Placeholder replaced with server; the pre-existing positive-id 42 is dropped.
    expect(out).toEqual([{ id: 42, q: 0.1, label: "server" }]);
  });

  it("replaces only the first placeholder when multiple match", () => {
    const list: Row[] = [
      { id: -1, q: 0.1, label: "a" },
      { id: -2, q: 0.1, label: "b" },
      { id:  5, q: 0.2, label: "keep" },
    ];
    const server: Row = { id: 42, q: 0.1, label: "server" };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    expect(out).toEqual([
      { id: 42, q: 0.1, label: "server" },
      { id: -2, q: 0.1, label: "b" },
      { id:  5, q: 0.2, label: "keep" },
    ]);
  });

  it("handles an empty list (appends server item)", () => {
    const out = replacePlaceholder<Row>([], { id: 42, q: 0.1 }, () => true);
    expect(out).toEqual([{ id: 42, q: 0.1 }]);
  });

  it("ignores positive-id rows even when matches() returns true", () => {
    const list: Row[] = [
      { id: 7, q: 0.1, label: "real" },
    ];
    const server: Row = { id: 42, q: 0.1, label: "server" };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    // No placeholder; server appended; pre-existing positive row preserved.
    expect(out).toEqual([
      { id: 7, q: 0.1, label: "real" },
      { id: 42, q: 0.1, label: "server" },
    ]);
  });

  it("preserves placeholder position (does not move the row to the end)", () => {
    const list: Row[] = [
      { id:  1, q: 0.2 },
      { id: -1, q: 0.1, label: "opt" },
      { id:  2, q: 0.3 },
    ];
    const server: Row = { id: 42, q: 0.1, label: "server" };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    expect(out).toEqual([
      { id:  1, q: 0.2 },
      { id: 42, q: 0.1, label: "server" },
      { id:  2, q: 0.3 },
    ]);
  });

  describe("isDuplicate override", () => {
    interface Tagged { id: number; source: "manual" | "auto" }

    it("preserves rows the override classifies as non-duplicates", () => {
      const list: Tagged[] = [
        { id: 42, source: "auto" },       // same id as server, but different source
        { id: -1, source: "manual" },     // placeholder
      ];
      const server: Tagged = { id: 42, source: "manual" };
      const out = replacePlaceholder(
        list,
        server,
        (r) => r.source === "manual",
        { isDuplicate: (r) => r.source === "manual" && r.id === server.id },
      );
      // Auto row with id=42 survives; placeholder replaced with server row.
      expect(out).toEqual([
        { id: 42, source: "auto" },
        { id: 42, source: "manual" },
      ]);
    });

    it("drops rows the override classifies as duplicates", () => {
      const list: Tagged[] = [
        { id: 42, source: "manual" },     // SSE-race duplicate
        { id: -1, source: "manual" },     // placeholder
      ];
      const server: Tagged = { id: 42, source: "manual" };
      const out = replacePlaceholder(
        list,
        server,
        (r) => r.source === "manual",
        { isDuplicate: (r) => r.source === "manual" && r.id === server.id },
      );
      expect(out).toEqual([{ id: 42, source: "manual" }]);
    });

    it("preserves peakAdd's exact three-row race: placeholder + auto-same-id + manual-different-id", () => {
      // Mirrors the original peakAdd loop's defended scenario:
      // placeholder gets replaced; auto peak sharing id namespace survives;
      // unrelated manual peak survives.
      const list: Tagged[] = [
        { id: -1, source: "manual" },     // placeholder for the q the user clicked
        { id: 42, source: "auto"   },     // auto peak that happens to have id=42
        { id: 15, source: "manual" },     // unrelated existing manual peak
      ];
      const server: Tagged = { id: 42, source: "manual" };
      const out = replacePlaceholder(
        list,
        server,
        (r) => r.source === "manual",
        { isDuplicate: (r) => r.source === "manual" && r.id === server.id },
      );
      expect(out).toEqual([
        { id: 42, source: "manual" },     // placeholder replaced
        { id: 42, source: "auto"   },     // auto survives — id collision is allowed across sources
        { id: 15, source: "manual" },     // unrelated manual survives
      ]);
    });
  });
});
