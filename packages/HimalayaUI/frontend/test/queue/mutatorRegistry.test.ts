import { describe, it, expect } from "vitest";
import { resolveMutator, resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import {
  updateSampleMutator,
  setExposureStatusMutator,
} from "../../src/lib/queue/mutators/trivial";
import { peakAddMutator } from "../../src/lib/queue/mutators/peakAdd";
import { peakRemoveMutator } from "../../src/lib/queue/mutators/peakRemove";
import {
  peakExcludeMutator,
  peakUnexcludeMutator,
} from "../../src/lib/queue/mutators/peakSetExcluded";
import { deleteIndexMutator } from "../../src/lib/queue/mutators/indexGroup";
import { createSpeculativeMutator } from "../../src/lib/queue/mutators/createSpeculative";
import { reanalyzeExposureMutator } from "../../src/lib/queue/mutators/reanalyzeExposure";

describe("resolveMutator", () => {
  it("returns the canonical mutator for non-dual kinds", () => {
    expect(resolveMutator({ kind: "peak_added", payload: {} })).toBe(
      peakAddMutator,
    );
    expect(resolveMutator({ kind: "peak_removed", payload: {} })).toBe(
      peakRemoveMutator,
    );
    expect(resolveMutator({ kind: "peak_excluded", payload: {} })).toBe(
      peakExcludeMutator,
    );
    expect(resolveMutator({ kind: "peak_unexcluded", payload: {} })).toBe(
      peakUnexcludeMutator,
    );
    expect(resolveMutator({ kind: "delete_index", payload: {} })).toBe(
      deleteIndexMutator,
    );
    expect(
      resolveMutator({ kind: "speculative_created", payload: {} }),
    ).toBe(createSpeculativeMutator);
    expect(
      resolveMutator({ kind: "reanalyze_exposure", payload: {} }),
    ).toBe(reanalyzeExposureMutator);
    expect(resolveMutator({ kind: "update_sample", payload: {} })).toBe(
      updateSampleMutator,
    );
    expect(
      resolveMutator({ kind: "set_exposure_status", payload: {} }),
    ).toBe(setExposureStatusMutator);
  });

  it("returns undefined for speculative_deleted (no outbound op)", () => {
    expect(
      resolveMutator({ kind: "speculative_deleted", payload: {} }),
    ).toBeUndefined();
  });

});

describe("resolveMutatorForEvent", () => {
  it("routes peak_added to peakAddMutator", () => {
    expect(resolveMutatorForEvent("peak_added")).toBe(peakAddMutator);
  });

  it("routes peak_excluded and peak_unexcluded to their dedicated mutators", () => {
    expect(resolveMutatorForEvent("peak_excluded")).toBe(peakExcludeMutator);
    expect(resolveMutatorForEvent("peak_unexcluded")).toBe(peakUnexcludeMutator);
  });

  it("routes analyze_run to reanalyzeExposure", () => {
    expect(resolveMutatorForEvent("analyze_run")).toBe(reanalyzeExposureMutator);
  });

  it("returns undefined for unknown event kinds", () => {
    expect(resolveMutatorForEvent("not_a_kind")).toBeUndefined();
  });
});

describe("resolveMutator ↔ resolveMutatorForEvent consistency", () => {
  const cases = [
    { mutator: peakAddMutator,          eventKind: "peak_added" },
    { mutator: peakRemoveMutator,       eventKind: "peak_removed" },
    { mutator: peakExcludeMutator,      eventKind: "peak_excluded" },
    { mutator: peakUnexcludeMutator,    eventKind: "peak_unexcluded" },
    { mutator: deleteIndexMutator,      eventKind: "speculative_deleted" },
    { mutator: createSpeculativeMutator, eventKind: "speculative_created" },
    { mutator: reanalyzeExposureMutator, eventKind: "analyze_run" },
    { mutator: updateSampleMutator,     eventKind: "update_sample" },
    { mutator: setExposureStatusMutator, eventKind: "set_exposure_status" },
  ];

  cases.forEach(({ mutator, eventKind }) => {
    it(`${eventKind} resolves to ${mutator.kind}`, () => {
      expect(resolveMutatorForEvent(eventKind)).toBe(mutator);
    });
  });
});
