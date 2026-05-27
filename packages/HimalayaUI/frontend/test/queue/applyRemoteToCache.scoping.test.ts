import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import type { SseEvent } from "../../src/lib/queue/types";

describe("applyRemoteToCache — add_tag(sample) scoping fan-out", () => {
  it("invalidates the corpus tag + picker caches so /series/new refreshes", () => {
    const qc = new QueryClient();
    const inv = vi.spyOn(qc, "invalidateQueries");
    // A foreign scoping write arrives as an add_tag/sample frame.
    const remote: SseEvent = {
      id: 1, kind: "add_tag", entity_type: "sample", entity_id: 10,
      actor: "peer", client_id: "peer", client_op_id: "peer-op",
      ts: "2026-05-14T10:00:00Z",
      payload: { key: "ratio", value: "1:1", experiment_id: 1, sample_id: 10 },
    };
    applyRemoteToCache(remote, qc);
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus-sample-tags"] });
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus-picker-samples"] });
    // Existing fan-out still fires (regression guard). corpusSamples key is
    // ["corpus","samples"] (queries.ts), not ["corpus-samples"].
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus", "samples"] });
  });

  it("leaves the exposure-tag branch unchanged (tri-scope safety)", () => {
    const qc = new QueryClient();
    const inv = vi.spyOn(qc, "invalidateQueries");
    const remote: SseEvent = {
      id: 2, kind: "add_tag", entity_type: "exposure", entity_id: 5,
      actor: "x", client_id: "y", client_op_id: "z",
      ts: "2026-05-14T10:00:00Z",
      payload: { sample_id: 10 },
    };
    applyRemoteToCache(remote, qc);
    // Exposure branch: exposures() only, never the corpus-scoping keys.
    expect(inv).toHaveBeenCalledWith({ queryKey: ["sample", 10, "exposures"] });
    expect(inv).not.toHaveBeenCalledWith({ queryKey: ["corpus-sample-tags"] });
    expect(inv).not.toHaveBeenCalledWith({ queryKey: ["corpus-picker-samples"] });
  });
});
