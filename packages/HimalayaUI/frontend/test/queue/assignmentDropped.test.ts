/**
 * F-WIPE W2 — assignment_dropped aggregation copy + toast gating.
 *
 * Backend reanalysis (W1) re-attaches assignment members semantically; a
 * member only drops when its index genuinely no longer exists in the new
 * candidate set. When that happens the peak_* SSE frame carries
 * `post_state.assignment_dropped: string[]` (one phase name PER MEMBER,
 * may repeat). The frontend aggregates the list into one toast.
 */
import { describe, it, expect, beforeEach, afterEach } from "vitest";
import {
  buildAssignmentDroppedMessage,
  announceAssignmentDropped,
} from "../../src/lib/queue/assignmentDropped";
import { setToastImpl } from "../../src/lib/toast";

describe("buildAssignmentDroppedMessage", () => {
  it("one member: phase + singular noun + singular tail", () => {
    expect(buildAssignmentDroppedMessage(["Pn3m"])).toBe(
      "Pn3m index dropped from the call. Its index no longer fits the peaks.",
    );
  });

  it("two members of one phase: count-prefixed, plural noun + tail", () => {
    expect(buildAssignmentDroppedMessage(["Pn3m", "Pn3m"])).toBe(
      "2 Pn3m indices dropped from the call. They no longer fit the peaks.",
    );
  });

  it("mixed phases: 'and'-joined, plural", () => {
    expect(buildAssignmentDroppedMessage(["Pn3m", "Lamellar"])).toBe(
      "Pn3m and Lamellar indices dropped from the call. They no longer fit the peaks.",
    );
  });

  it("three phase groups: comma list with trailing 'and'", () => {
    expect(buildAssignmentDroppedMessage(["Pn3m", "Im3m", "Lamellar"])).toBe(
      "Pn3m, Im3m and Lamellar indices dropped from the call. They no longer fit the peaks.",
    );
  });

  it("repeats aggregate within a mixed list, first-seen phase order preserved", () => {
    expect(buildAssignmentDroppedMessage(["Lamellar", "Pn3m", "Lamellar"])).toBe(
      "2 Lamellar and Pn3m indices dropped from the call. They no longer fit the peaks.",
    );
  });
});

describe("announceAssignmentDropped", () => {
  let toastCalls: Array<{ msg: string; kind: string }> = [];
  beforeEach(() => {
    toastCalls = [];
    setToastImpl((msg, kind) => { toastCalls.push({ msg, kind }); });
  });
  afterEach(() => { setToastImpl(null); });

  it("fires one warning toast with the aggregated copy", () => {
    announceAssignmentDropped(["Pn3m"]);
    expect(toastCalls).toEqual([{
      msg: "Pn3m index dropped from the call. Its index no longer fits the peaks.",
      kind: "warning",
    }]);
  });

  it("no toast on an empty list", () => {
    announceAssignmentDropped([]);
    expect(toastCalls).toEqual([]);
  });

  it("no toast when the field is absent (undefined)", () => {
    announceAssignmentDropped(undefined);
    expect(toastCalls).toEqual([]);
  });
});
