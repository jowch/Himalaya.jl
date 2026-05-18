import { describe, it, expect } from "vitest";
import { resolveMutator, resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import {
  addSampleTagMutator,
  removeSampleTagMutator,
  addExposureTagMutator,
  removeExposureTagMutator,
  addCorpusSampleTagMutator,
  removeCorpusSampleTagMutator,
  updateSampleMutator,
  postSampleMessageMutator,
  postComparisonMessageMutator,
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
import { saveComparisonMutator } from "../../src/lib/queue/mutators/saveComparison";
import { deleteComparisonMutator } from "../../src/lib/queue/mutators/deleteComparison";

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

  it("dispatches a corpus sample-tag op (no experimentId, no exposureId) to the corpus mutator", () => {
    expect(
      resolveMutator({
        kind: "add_tag",
        payload: { sampleId: 10, key: "k", value: "v" },
      }),
    ).toBe(addCorpusSampleTagMutator);
    expect(
      resolveMutator({
        kind: "remove_tag",
        payload: { sampleId: 10, tagId: 7 },
      }),
    ).toBe(removeCorpusSampleTagMutator);
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

  it("tolerates undefined/null payload by treating add_tag as corpus-scoped", () => {
    // No exposureId and no experimentId → corpus fallthrough. Must not throw.
    expect(
      resolveMutator({ kind: "add_tag", payload: undefined }),
    ).toBe(addCorpusSampleTagMutator);
    expect(resolveMutator({ kind: "add_tag", payload: null })).toBe(
      addCorpusSampleTagMutator,
    );
  });
});

describe("resolveMutatorForEvent", () => {
  it("routes peak_added to peakAddMutator", () => {
    expect(resolveMutatorForEvent("peak_added", "exposure")).toBe(peakAddMutator);
  });

  it("routes peak_excluded and peak_unexcluded to their dedicated mutators", () => {
    expect(resolveMutatorForEvent("peak_excluded",   "exposure")).toBe(peakExcludeMutator);
    expect(resolveMutatorForEvent("peak_unexcluded", "exposure")).toBe(peakUnexcludeMutator);
  });

  it("splits add_tag by entity_type", () => {
    expect(resolveMutatorForEvent("add_tag", "sample"))  .toBe(addSampleTagMutator);
    expect(resolveMutatorForEvent("add_tag", "exposure")).toBe(addExposureTagMutator);
  });

  it("splits remove_tag by entity_type", () => {
    expect(resolveMutatorForEvent("remove_tag", "sample"))  .toBe(removeSampleTagMutator);
    expect(resolveMutatorForEvent("remove_tag", "exposure")).toBe(removeExposureTagMutator);
  });

  it("splits post_message by entity_type (wire-string entity types)", () => {
    expect(resolveMutatorForEvent("post_message", "sample_message"))    .toBe(postSampleMessageMutator);
    expect(resolveMutatorForEvent("post_message", "comparison_message")).toBe(postComparisonMessageMutator);
  });

  it("routes both comparison_created and comparison_submitted to saveComparison", () => {
    expect(resolveMutatorForEvent("comparison_created",   "comparison")).toBe(saveComparisonMutator);
    expect(resolveMutatorForEvent("comparison_submitted", "comparison")).toBe(saveComparisonMutator);
  });

  it("routes comparison_deleted to deleteComparison", () => {
    expect(resolveMutatorForEvent("comparison_deleted", "comparison")).toBe(deleteComparisonMutator);
  });

  it("routes analyze_run to reanalyzeExposure", () => {
    expect(resolveMutatorForEvent("analyze_run", "exposure")).toBe(reanalyzeExposureMutator);
  });

  it("returns undefined for unknown event kinds", () => {
    expect(resolveMutatorForEvent("not_a_kind", "exposure")).toBeUndefined();
  });
});

describe("resolveMutator ↔ resolveMutatorForEvent consistency", () => {
  // Cross-check: for each mutator, every (eventKind, entityType) pair that
  // resolveMutatorForEvent's switch maps to it must round-trip back to the
  // same mutator instance. Catches drift if someone adds a mutator + updates
  // resolveMutator but forgets resolveMutatorForEvent (or vice versa).
  const cases = [
    { mutator: peakAddMutator,          eventKind: "peak_added",          entityType: "exposure"   },
    { mutator: peakRemoveMutator,       eventKind: "peak_removed",        entityType: "exposure"   },
    { mutator: peakExcludeMutator,      eventKind: "peak_excluded",       entityType: "exposure"   },
    { mutator: peakUnexcludeMutator,    eventKind: "peak_unexcluded",     entityType: "exposure"   },
    { mutator: addIndexToGroupMutator,  eventKind: "index_confirmed",     entityType: "exposure"   },
    { mutator: removeIndexFromGroupMutator, eventKind: "index_unconfirmed", entityType: "exposure" },
    { mutator: deleteIndexMutator,      eventKind: "speculative_deleted", entityType: "exposure"   },
    { mutator: createSpeculativeMutator, eventKind: "speculative_created", entityType: "exposure"  },
    { mutator: reanalyzeExposureMutator, eventKind: "analyze_run",        entityType: "exposure"   },
    { mutator: saveComparisonMutator,   eventKind: "comparison_created",  entityType: "comparison" },
    { mutator: saveComparisonMutator,   eventKind: "comparison_submitted", entityType: "comparison" },
    { mutator: deleteComparisonMutator, eventKind: "comparison_deleted",  entityType: "comparison" },
    { mutator: addSampleTagMutator,     eventKind: "add_tag",             entityType: "sample"     },
    { mutator: addExposureTagMutator,   eventKind: "add_tag",             entityType: "exposure"   },
    { mutator: removeSampleTagMutator,  eventKind: "remove_tag",          entityType: "sample"     },
    { mutator: removeExposureTagMutator, eventKind: "remove_tag",         entityType: "exposure"   },
    { mutator: postSampleMessageMutator, eventKind: "post_message",       entityType: "sample_message"     },
    { mutator: postComparisonMessageMutator, eventKind: "post_message",   entityType: "comparison_message" },
    { mutator: updateSampleMutator,     eventKind: "update_sample",       entityType: "sample"     },
    { mutator: setExposureStatusMutator, eventKind: "set_exposure_status", entityType: "exposure"  },
    { mutator: selectExposureMutator,   eventKind: "select_exposure",     entityType: "exposure"   },
  ];

  // Use forEach + it() instead of it.each($mutator.kind) — dot-property
  // interpolation in it.each titles is Vitest-2.1+ but the codebase has no
  // prior example; a template-literal title is unambiguous.
  cases.forEach(({ mutator, eventKind, entityType }) => {
    it(`${eventKind} / ${entityType} resolves to ${mutator.kind}`, () => {
      expect(resolveMutatorForEvent(eventKind, entityType)).toBe(mutator);
    });
  });
});
