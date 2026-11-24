import { describe, expect, it } from "vitest";
import { modulesFor, NO_MODULES } from "../../src/workspace/modules";

describe("modulesFor", () => {
  it("gives a capillary setup every module", () => {
    expect(modulesFor("capillary")).toEqual({ peakEditing: true, calls: true, indexing: true, imageSelection: true, review: true, background: false });
  });

  it("treats a setup whose mode is not decided, or unknown, like a capillary", () => {
    expect(modulesFor(null)).toEqual(modulesFor("capillary"));
    expect(modulesFor(undefined)).toEqual(modulesFor("capillary"));
  });

  it("gives an autosampler setup only the background module", () => {
    expect(modulesFor("autosampler")).toEqual({ peakEditing: false, calls: false, indexing: false, imageSelection: false, review: false, background: true });
  });

  it("shows no module while the setup modes are still loading", () => {
    expect(NO_MODULES.background).toBe(false);
  });
});
