import { describe, it, expect } from "vitest";
import { formatRelativeTime } from "../src/lib/relativeTime";

describe("formatRelativeTime", () => {
  const now = new Date("2026-06-22T12:00:00Z").getTime();
  const ago = (sec: number) => new Date(now - sec * 1000).toISOString();

  it("buckets recency coarsely", () => {
    expect(formatRelativeTime(ago(10), now)).toBe("just now");
    expect(formatRelativeTime(ago(120), now)).toBe("2 min ago");
    expect(formatRelativeTime(ago(2 * 3600), now)).toBe("2 hr ago");
    expect(formatRelativeTime(ago(26 * 3600), now)).toBe("1 day ago");
    expect(formatRelativeTime(ago(3 * 86400), now)).toBe("3 days ago");
  });

  it("clamps a future timestamp to just now and rejects garbage", () => {
    expect(formatRelativeTime(ago(-60), now)).toBe("just now");
    expect(formatRelativeTime("not-a-date", now)).toBe("");
  });
});
