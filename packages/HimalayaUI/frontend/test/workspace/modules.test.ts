import { describe, expect, it } from "vitest";
import { modulesFor } from "../../src/workspace/modules";

describe("modulesFor", () => {
  it("gives a capillary setup every module", () => {
    expect(modulesFor("capillary")).toEqual({ peakEditing: true, calls: true, indexing: true, imageSelection: true, review: true });
  });

  it("treats a setup whose mode is not decided, or unknown, like a capillary", () => {
    expect(modulesFor(null)).toEqual(modulesFor("capillary"));
    expect(modulesFor(undefined)).toEqual(modulesFor("capillary"));
  });

  it("gives an autosampler setup none of them", () => {
    expect(modulesFor("autosampler")).toEqual({ peakEditing: false, calls: false, indexing: false, imageSelection: false, review: false });
  });
});
