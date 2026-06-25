import { describe, it, expect } from "vitest";
import { effectiveIngestStatus } from "../src/lib/ingestStatus";

describe("effectiveIngestStatus", () => {
  it("terminal persisted state wins over a stale overlay — both scan kinds (8c self-heal)", () => {
    // The missed-terminal-frame case: a dropped ingest_complete strands the
    // overlay at "scanning" (initial scan) or "analyzing" (rescan); the polled
    // terminal row wins either way.
    expect(effectiveIngestStatus("scanning", "complete")).toBe("complete");
    expect(effectiveIngestStatus("scanning", "failed")).toBe("failed");
    expect(effectiveIngestStatus("analyzing", "complete")).toBe("complete");
    expect(effectiveIngestStatus("analyzing", "failed")).toBe("failed");
  });

  it("the live overlay leads while persisted is non-terminal (carries the phase)", () => {
    // During a live rescan the persisted row is "scanning"; the overlay's
    // "analyzing" phase is what distinguishes it from an initial scan.
    expect(effectiveIngestStatus("analyzing", "scanning")).toBe("analyzing");
    expect(effectiveIngestStatus("scanning", "scanning")).toBe("scanning");
    expect(effectiveIngestStatus("scanning", "idle")).toBe("scanning");
  });

  it("falls back to persisted then idle when no overlay", () => {
    expect(effectiveIngestStatus(undefined, "scanning")).toBe("scanning");
    expect(effectiveIngestStatus(undefined, undefined)).toBe("idle");
  });
});
