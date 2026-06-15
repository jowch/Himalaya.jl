/**
 * Contract tests for the edit_tag queue chain (Slice 2 — LO-TAGDUP).
 *
 * Covers the three legs of the six-layer contract-testing rule:
 *   (a) optimistic onMutate patches by id (id 7 edits, id 3 untouched)
 *   (b) own-op SSE confirmation — synthesizeFromSse returns the synthesized
 *       partial for an edit_tag SSE frame (the SSE-wins-partial half that
 *       mocked tests usually hide)
 *   (c) foreign-event replay — resolveMutatorForEvent routes "edit_tag" to
 *       the corpus edit mutator
 *
 * Paired with: packages/HimalayaUI/test/test_routes_samples.jl
 * (the edit_tag backend tests — two sides of the six-layer rule).
 */
import { describe, it, expect } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { editCorpusSampleTagMutator } from "../../src/lib/queue/mutators/trivial";
import { resolveMutator, resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import { queryKeys } from "../../src/queries";

describe("editCorpusSampleTagMutator", () => {
  // -------------------------------------------------------------------------
  // (a) optimistic onMutate
  // -------------------------------------------------------------------------
  it("optimistically patches the tag in place by id (NOT key+value)", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.corpusSamples, [
      {
        id: 1,
        tags: [
          { id: 3, key: "dose", value: "10", source: "manual" },
          { id: 7, key: "dose", value: "10", source: "scoping" },
        ],
      },
    ]);

    editCorpusSampleTagMutator.onMutate(
      { sampleId: 1, tagId: 7, key: "dose", value: "12" } as any,
      qc,
    );

    const tags = (qc.getQueryData(queryKeys.corpusSamples) as any)[0].tags;
    // id 7 was edited
    expect(tags.find((t: any) => t.id === 7).value).toBe("12");
    // id 3 is untouched — id-exact, not key+value match
    expect(tags.find((t: any) => t.id === 3).value).toBe("10");
  });

  it("onMutate returns a restore closure that reverts the optimistic change", () => {
    const qc = new QueryClient();
    const original = [
      {
        id: 1,
        tags: [{ id: 7, key: "dose", value: "10", source: "manual" }],
      },
    ];
    qc.setQueryData(queryKeys.corpusSamples, original);

    const ctx = editCorpusSampleTagMutator.onMutate(
      { sampleId: 1, tagId: 7, key: "dose", value: "99" } as any,
      qc,
    );
    // Optimistic patch applied
    expect((qc.getQueryData(queryKeys.corpusSamples) as any)[0].tags[0].value).toBe("99");

    // Rollback
    ctx.restore();
    expect((qc.getQueryData(queryKeys.corpusSamples) as any)[0].tags[0].value).toBe("10");
  });

  // -------------------------------------------------------------------------
  // (b) own-op SSE confirmation — synthesizeFromSse (the SSE-wins-partial leg)
  // -------------------------------------------------------------------------
  it("synthesizeFromSse returns a synthesized SampleTag for a valid edit_tag SSE frame", () => {
    const remote = {
      kind: "edit_tag",
      entity_type: "sample",
      entity_id: 1,
      payload: { tag_id: 7, key: "dose", value: "12", experiment_id: 42 },
    } as any;
    const base = { id: -999 } as any; // placeholder id (not used in edit — but ...base should spread)

    const result = editCorpusSampleTagMutator.synthesizeFromSse!(remote, base);

    expect(result).not.toBeUndefined();
    expect(result!.id).toBe(7);
    expect(result!.key).toBe("dose");
    expect(result!.value).toBe("12");
    // source is hardcoded "manual" (faint-marker inaccuracy — acceptable, see comment in mutator)
    expect(result!.source).toBe("manual");
  });

  it("synthesizeFromSse returns undefined when tag_id is missing (guard)", () => {
    const remote = {
      kind: "edit_tag",
      entity_type: "sample",
      entity_id: 1,
      payload: { key: "dose", value: "12" }, // tag_id missing
    } as any;
    const base = {} as any;

    const result = editCorpusSampleTagMutator.synthesizeFromSse!(remote, base);
    expect(result).toBeUndefined();
  });

  it("synthesizeFromSse handles missing payload gracefully", () => {
    const remote = { kind: "edit_tag", entity_type: "sample", entity_id: 1 } as any;
    const base = {} as any;

    const result = editCorpusSampleTagMutator.synthesizeFromSse!(remote, base);
    expect(result).toBeUndefined();
  });
});

// ---------------------------------------------------------------------------
// (c) foreign-event replay — registry routing
// ---------------------------------------------------------------------------
describe("resolveMutator + resolveMutatorForEvent for edit_tag", () => {
  it("resolveMutator routes a sampleId-only edit_tag op to the corpus edit mutator", () => {
    const mutator = resolveMutator({
      kind: "edit_tag",
      payload: { sampleId: 10, tagId: 7, key: "dose", value: "12" },
    });
    expect(mutator).toBe(editCorpusSampleTagMutator);
  });

  it("resolveMutatorForEvent routes edit_tag/sample to the corpus edit mutator (foreign replay)", () => {
    const mutator = resolveMutatorForEvent("edit_tag", "sample");
    expect(mutator).toBe(editCorpusSampleTagMutator);
  });

  it("resolveMutatorForEvent routes edit_tag with any entity_type to the corpus edit mutator", () => {
    // edit_tag is sample-only — no exposure variant exists
    expect(resolveMutatorForEvent("edit_tag", "exposure")).toBe(editCorpusSampleTagMutator);
  });
});
