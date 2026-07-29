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

/**
 * Cross-language constant pin.
 *
 * `CUSTOM_SNAP_TOL` (speculative.jl) and `landsOn`'s default `relTol`
 * (lib/customIndex.ts) MUST be equal: the modal counts "N of M land" with one
 * and the backend claims peaks with the other, so a drift between them means
 * the commit claims a different set than the user was shown — the exact bug
 * this PR fixes. Both docstrings say MUST; nothing enforced it.
 */
describe("CUSTOM_SNAP_TOL ↔ landsOn relTol", () => {
  it("the Julia constant equals the TS default", async () => {
    const fs = await import("node:fs");
    const path = await import("node:path");
    const here = path.dirname(new URL(import.meta.url).pathname);
    const root = path.resolve(here, "../../../../..");

    const jl = fs.readFileSync(
      path.join(root, "packages/HimalayaUI/src/speculative.jl"), "utf8");
    const ts = fs.readFileSync(
      path.join(root, "packages/HimalayaUI/frontend/src/lib/customIndex.ts"), "utf8");

    const jlTol = jl.match(/^const CUSTOM_SNAP_TOL = ([0-9.]+)/m)?.[1];
    const tsTol = ts.match(/export function landsOn\([^)]*relTol = ([0-9.]+)/s)?.[1];

    // A null here means the constant was renamed or reformatted — fix the
    // regex AND re-check the pairing, don't delete the test.
    expect(jlTol, "CUSTOM_SNAP_TOL not found in speculative.jl").not.toBeUndefined();
    expect(tsTol, "landsOn relTol default not found in customIndex.ts").not.toBeUndefined();
    expect(Number(jlTol)).toBe(Number(tsTol));
  });
});
