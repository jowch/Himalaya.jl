import { describe, it, expect } from "vitest";
import { resolveMutator } from "../../src/lib/queue/mutatorRegistry";
import {
  addSampleTagMutator,
  removeSampleTagMutator,
  addExposureTagMutator,
  removeExposureTagMutator,
  updateSampleMutator,
  postSampleMessageMutator,
  setExposureStatusMutator,
  selectExposureMutator,
} from "../../src/lib/queue/mutators/trivial";
import { peakAddMutator } from "../../src/lib/queue/mutators/peakAdd";
import { peakRemoveMutator } from "../../src/lib/queue/mutators/peakRemove";
import {
  peakExcludeMutator,
  peakUnexcludeMutator,
} from "../../src/lib/queue/mutators/peakSetExcluded";
import {
  addIndexToGroupMutator,
  removeIndexFromGroupMutator,
  deleteIndexMutator,
} from "../../src/lib/queue/mutators/indexGroup";
import { createSpeculativeMutator } from "../../src/lib/queue/mutators/createSpeculative";
import { reanalyzeExposureMutator } from "../../src/lib/queue/mutators/reanalyzeExposure";

describe("resolveMutator", () => {
  it("dispatches add_tag based on payload experimentId presence", () => {
    expect(
      resolveMutator({
        kind: "add_tag",
        payload: { experimentId: 1, sampleId: 10, key: "k", value: "v" },
      }),
    ).toBe(addSampleTagMutator);
    expect(
      resolveMutator({
        kind: "add_tag",
        payload: { sampleId: 10, exposureId: 99, key: "k", value: "v" },
      }),
    ).toBe(addExposureTagMutator);
  });

  it("dispatches remove_tag the same way", () => {
    expect(
      resolveMutator({
        kind: "remove_tag",
        payload: { experimentId: 1, sampleId: 10, tagId: 7 },
      }),
    ).toBe(removeSampleTagMutator);
    expect(
      resolveMutator({
        kind: "remove_tag",
        payload: { sampleId: 10, exposureId: 99, tagId: 7 },
      }),
    ).toBe(removeExposureTagMutator);
  });

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
    expect(resolveMutator({ kind: "index_confirmed", payload: {} })).toBe(
      addIndexToGroupMutator,
    );
    expect(resolveMutator({ kind: "index_unconfirmed", payload: {} })).toBe(
      removeIndexFromGroupMutator,
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
    expect(resolveMutator({ kind: "post_message", payload: {} })).toBe(
      postSampleMessageMutator,
    );
    expect(
      resolveMutator({ kind: "set_exposure_status", payload: {} }),
    ).toBe(setExposureStatusMutator);
    expect(
      resolveMutator({ kind: "select_exposure", payload: {} }),
    ).toBe(selectExposureMutator);
  });

  it("returns undefined for speculative_deleted (no outbound op)", () => {
    expect(
      resolveMutator({ kind: "speculative_deleted", payload: {} }),
    ).toBeUndefined();
  });

  it("tolerates undefined/null payload by treating add_tag as exposure-scoped", () => {
    // No experimentId means exposure-scoped fallback. This should not throw.
    expect(
      resolveMutator({ kind: "add_tag", payload: undefined }),
    ).toBe(addExposureTagMutator);
    expect(resolveMutator({ kind: "add_tag", payload: null })).toBe(
      addExposureTagMutator,
    );
  });
});
