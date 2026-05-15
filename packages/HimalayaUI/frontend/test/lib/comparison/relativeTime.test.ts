import { describe, it, expect } from "vitest";
import { relativeTime } from "../../../src/lib/comparison/relativeTime";

describe("relativeTime — Compare UX C-2", () => {
  const NOW = new Date("2026-05-14T12:00:00Z").getTime();
  const at = (s: string) => relativeTime(s, NOW);

  it("renders 'just now' for <60s", () => {
    expect(at("2026-05-14T11:59:30Z")).toBe("just now");
  });
  it("renders minutes", () => {
    expect(at("2026-05-14T11:55:00Z")).toBe("5m ago");
  });
  it("renders hours", () => {
    expect(at("2026-05-14T10:00:00Z")).toBe("2h ago");
  });
  it("renders days up to 30d", () => {
    expect(at("2026-05-10T12:00:00Z")).toBe("4d ago");
  });
  it("still renders '30d ago' exactly at the 30-day boundary", () => {
    expect(at("2026-04-14T12:00:00Z")).toBe("30d ago");
  });
  it("falls back to a date once past 30 days", () => {
    expect(at("2026-04-13T12:00:00Z")).toMatch(/^\d{4}-\d{2}-\d{2}$/);
  });
  it("falls back to a date for >30d", () => {
    expect(at("2026-01-01T12:00:00Z")).toMatch(/\d{4}-\d{2}-\d{2}/);
  });
  it("returns null for null input", () => {
    expect(relativeTime(null, NOW)).toBeNull();
  });
});
