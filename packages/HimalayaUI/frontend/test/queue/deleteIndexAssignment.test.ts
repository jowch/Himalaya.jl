/**
 * deleteIndexMutator ↔ assignment cache coherence.
 *
 * `assignment_members.index_id` is `ON DELETE CASCADE` (db.jl:223), so
 * DELETE /api/indices/:id silently removes the index from the durable
 * assignment as well. The mutator's optimistic write must mirror that or the
 * cart renders a member id that no longer resolves to any index — a phantom
 * PhaseBlock until something else refetches the assignment.
 *
 * Layer 5 of the six-layer contract (docs/contract-testing.md); the symmetric
 * layer-3 case lives in sseEventPayload.contract.test.ts, and rollback
 * symmetry for both caches in rollbackSymmetry.test.ts.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { deleteIndexMutator } from "../../src/lib/queue/mutators/indexGroup";
import { queryKeys } from "../../src/queries";
import type { Assignment } from "../../src/api";

const EXPOSURE = 5;

function mutate(qc: QueryClient, indexId: number) {
  return deleteIndexMutator.onMutate({
    kind: "delete_index", clientOpId: "op",
    payload: { indexId },
    exposureId: EXPOSURE, username: "alice", clientId: "tab",
    indexId,
  } as never, qc);
}

function readAssignment(qc: QueryClient): Assignment | undefined {
  return qc.getQueryData<Assignment>(queryKeys.assignment(EXPOSURE));
}

describe("deleteIndexMutator optimistic assignment drop", () => {
  let qc: QueryClient;

  beforeEach(() => {
    qc = new QueryClient();
    qc.setQueryData(queryKeys.assignment(EXPOSURE),
      { exposure_id: EXPOSURE, state: "indexed", members: [10, 11] });
  });

  it("drops the deleted index from members", () => {
    mutate(qc, 10);
    expect(readAssignment(qc)!.members).toEqual([11]);
  });

  it("leaves sibling members and the state untouched", () => {
    mutate(qc, 10);
    const a = readAssignment(qc)!;
    expect(a.members).not.toContain(10);
    expect(a.state).toBe("indexed");
    expect(a.exposure_id).toBe(EXPOSURE);
  });

  it("is a no-op when the index was not in the call", () => {
    mutate(qc, 99);
    expect(readAssignment(qc)!.members).toEqual([10, 11]);
  });

  it("does not create an assignment cache entry when none was seeded", () => {
    // A Focus page that never read the assignment must not gain a synthetic
    // one from a delete — the next real read would then skip the fetch.
    const bare = new QueryClient();
    mutate(bare, 10);
    expect(bare.getQueryData(queryKeys.assignment(EXPOSURE))).toBeUndefined();
  });

  it("rollback restores the member", () => {
    const ctx = mutate(qc, 10);
    ctx.restore();
    expect(readAssignment(qc)!.members).toEqual([10, 11]);
  });
});
