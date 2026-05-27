import { describe, it, expect, vi, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { scopeSeriesMutator } from "../../src/lib/queue/mutators/scopeSeries";
import * as api from "../../src/api";

describe("scopeSeriesMutator", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("kind is add_tag (reuses the existing event kind)", () => {
    expect(scopeSeriesMutator.kind).toBe("add_tag");
  });

  it("onMutate is a no-op (no optimistic write; user navigates away)", () => {
    const qc = new QueryClient();
    const ctx = scopeSeriesMutator.onMutate(
      { key: "ratio", tags: [{ sampleId: 10, value: "1:1" }],
        username: "a", clientId: "c", clientOpId: "op1" } as never,
      qc,
    );
    expect(typeof ctx.restore).toBe("function");
    expect(() => ctx.restore()).not.toThrow(); // no-op
  });

  it("request POSTs the batch with source='scoping'", async () => {
    const spy = vi.spyOn(api, "batchSampleTags").mockResolvedValue([]);
    await scopeSeriesMutator.request(
      { key: "ratio", tags: [{ sampleId: 10, value: "1:1" }, { sampleId: 11, value: "2:1" }],
        username: "a", clientId: "c", clientOpId: "op1" } as never,
      new AbortController().signal,
    );
    expect(spy).toHaveBeenCalledWith(
      "ratio",
      [{ sample_id: 10, value: "1:1" }, { sample_id: 11, value: "2:1" }],
      "scoping",
      expect.anything(), // AuthOpts incl. X-Client-Op-Id
    );
  });

  it("onSuccess invalidates only the corpus tag + corpus picker caches (no broad blast)", () => {
    const qc = new QueryClient();
    const inv = vi.spyOn(qc, "invalidateQueries");
    scopeSeriesMutator.onSuccess(
      { key: "ratio", tags: [], username: "a", clientId: "c", clientOpId: "op1" } as never, [], qc);
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus-sample-tags"] });
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus-picker-samples"] });
    // Narrowed: must NOT blast the whole experiment subtree (queue anti-pattern).
    expect(inv).not.toHaveBeenCalledWith({ queryKey: ["experiment"] });
    expect(inv).toHaveBeenCalledTimes(2);
  });

  it("defines no synthesizeFromSse (never invoked; foreign add_tag routes to addSampleTagMutator)", () => {
    expect(scopeSeriesMutator.synthesizeFromSse).toBeUndefined();
  });
});
