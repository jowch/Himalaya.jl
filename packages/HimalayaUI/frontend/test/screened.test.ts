import { describe, it, expect } from "vitest";
import { isSampleScreened } from "../src/lib/sample/screened";
import type { CorpusSample, Exposure } from "../src/api";

function sample(over: Partial<CorpusSample> = {}): CorpusSample {
  return {
    id: 1,
    experiment_id: 1,
    name: "s",
    display_name: null,
    notes: null,
    tags: [],
    q_units: "A-1",
    ...over,
  };
}

function exp(status: Exposure["status"]): Exposure {
  return {
    id: 1,
    sample_id: 1,
    filename: "f",
    kind: "file",
    selected: false,
    status,
    image_path: null,
    image_version: "",
    tags: [],
    sources: [],
    trace_hash: null,
    analysis_inputs_hash: null,
  };
}

describe("isSampleScreened", () => {
  it("trusts an explicit screened flag (the #162 field) over derivation", () => {
    expect(isSampleScreened(sample({ screened: true }), [exp(null)])).toBe(true);
    expect(
      isSampleScreened(sample({ screened: false }), [exp("rejected")]),
    ).toBe(false);
  });

  it("derives screened when every exposure has a non-null status", () => {
    expect(
      isSampleScreened(sample(), [exp("accepted"), exp("rejected")]),
    ).toBe(true);
  });

  it("is not screened when any exposure is untriaged (null status)", () => {
    expect(isSampleScreened(sample(), [exp("accepted"), exp(null)])).toBe(false);
  });

  it("is not screened with no or undefined exposures", () => {
    expect(isSampleScreened(sample(), [])).toBe(false);
    expect(isSampleScreened(sample(), undefined)).toBe(false);
  });
});
