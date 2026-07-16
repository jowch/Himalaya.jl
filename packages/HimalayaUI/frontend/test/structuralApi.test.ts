import { describe, it, expect, vi, afterEach } from "vitest";
import { renameSample, moveExposure, mergeSamples, splitSample, dismissGroupingFlag } from "../src/api";

afterEach(() => vi.restoreAllMocks());

function mockOnce(status = 200, body: unknown = {}) {
  return vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(body), { status, headers: { "Content-Type": "application/json" } }),
  );
}

describe("structural-edit fetchers (use request<T>, thread X-Client-Op-Id)", () => {
  it("renameSample PATCHes /api/samples/:id/name with { name }", async () => {
    const f = mockOnce(200, { id: 1, name: "X", notes: null });
    await renameSample(1, "X", { username: "alice", clientOpId: "op1" });
    const [url, init] = f.mock.calls[0]!;
    expect(String(url)).toContain("/api/samples/1/name");
    expect((init as RequestInit).method).toBe("PATCH");
    expect(JSON.parse(String((init as RequestInit).body))).toEqual({ name: "X" });
    expect((init as RequestInit).headers).toMatchObject({ "X-Client-Op-Id": "op1" });
  });

  it("moveExposure POSTs /api/exposures/:id/move with { sample_id }", async () => {
    const f = mockOnce();
    await moveExposure(100, 20, { clientOpId: "op2" });
    const [url, init] = f.mock.calls[0]!;
    expect(String(url)).toContain("/api/exposures/100/move");
    expect((init as RequestInit).method).toBe("POST");
    expect(JSON.parse(String((init as RequestInit).body))).toEqual({ sample_id: 20 });
  });

  it("mergeSamples POSTs /api/samples/:loser/merge with { survivor_id }", async () => {
    const f = mockOnce(200, { loser_id: 5, survivor_id: 3 });
    await mergeSamples(5, 3, { clientOpId: "op3" });
    const [url, init] = f.mock.calls[0]!;
    expect(String(url)).toContain("/api/samples/5/merge");
    expect(JSON.parse(String((init as RequestInit).body))).toEqual({ survivor_id: 3 });
  });

  it("splitSample POSTs /api/samples/:id/split with { exposure_ids, name }", async () => {
    const f = mockOnce(201, { new_sample_id: 9 });
    await splitSample(1, [101, 102], "HA85 (S01P01b)", { clientOpId: "op4" });
    const [url, init] = f.mock.calls[0]!;
    expect(String(url)).toContain("/api/samples/1/split");
    expect(JSON.parse(String((init as RequestInit).body))).toEqual({
      exposure_ids: [101, 102], name: "HA85 (S01P01b)",
    });
  });

  it("dismissGroupingFlag POSTs the dismiss route with { flag_kind, merge_with_sample_id? }", async () => {
    const f = mockOnce(200, {});
    await dismissGroupingFlag(20, { flag_kind: "merge", merge_with_sample_id: 10 }, { clientOpId: "op5" });
    const [url, init] = f.mock.calls[0]!;
    expect(String(url)).toContain("/api/samples/20/dismiss-flag"); // confirm exact path vs Phase D
    expect((init as RequestInit).method).toBe("POST");
    expect(JSON.parse(String((init as RequestInit).body))).toEqual({ flag_kind: "merge", merge_with_sample_id: 10 });
  });
});
