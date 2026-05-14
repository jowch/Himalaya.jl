import { afterEach, describe, expect, it } from "vitest";
import { __resetOptimisticIdForTest, nextOptimisticId } from "../../src/lib/queue/optimisticId";

// Simulates a fresh tab: clears the per-tab client_id, resets the
// optimisticId module's cached fingerprint and counter, and installs the
// given uuid as the new client_id.
function asTab(clientId: string): void {
  sessionStorage.clear();
  sessionStorage.setItem("himalaya.client_id", clientId);
  __resetOptimisticIdForTest();
}

afterEach(() => {
  sessionStorage.clear();
  __resetOptimisticIdForTest();
});

describe("nextOptimisticId", () => {
  it("returns strictly-negative monotonic ids within a single tab", () => {
    asTab("aaaaaaaa-aaaa-4aaa-8aaa-aaaaaaaaaaaa");
    const ids = Array.from({ length: 32 }, () => nextOptimisticId());
    for (const id of ids) {
      expect(id).toBeLessThan(0);
      expect(Number.isSafeInteger(id)).toBe(true);
    }
    // Strictly decreasing magnitude — counter increments by one each call.
    for (let i = 1; i < ids.length; i++) {
      expect(ids[i]).toBeLessThan(ids[i - 1]!);
    }
  });

  it("two tabs minted in the same tick produce disjoint id streams", () => {
    asTab("11111111-1111-4111-8111-111111111111");
    const tabA = Array.from({ length: 64 }, () => nextOptimisticId());

    asTab("22222222-2222-4222-8222-222222222222");
    const tabB = Array.from({ length: 64 }, () => nextOptimisticId());

    const setA = new Set(tabA);
    for (const id of tabB) {
      expect(setA.has(id)).toBe(false);
    }
    // The old encoding (`-(Date.now() * 1000 + counter)`) would have made
    // these arrays bit-for-bit identical, since both tabs increment from
    // counter=0 against the same Date.now() in a tight test loop.
    expect(tabA).not.toEqual(tabB);
  });

  it("the same tab is deterministic across module-state resets", () => {
    asTab("33333333-3333-4333-8333-333333333333");
    const first = Array.from({ length: 8 }, () => nextOptimisticId());

    asTab("33333333-3333-4333-8333-333333333333");
    const second = Array.from({ length: 8 }, () => nextOptimisticId());

    expect(first).toEqual(second);
  });
});
