import { describe, it, expect } from "vitest";
import cases from "../fixtures/representative-exposure.json";
import { keptExposures, representativeExposure } from "../../src/lib/exposure";
import type { Exposure } from "../../src/api";

type Case = { case: string; frames: Pick<Exposure, "id" | "selected" | "status">[]; representative: number | null };

function exposures(c: Case): Exposure[] {
  return c.frames.map((f) => ({
    ...f, sample_id: 1, filename: `f${f.id}`, kind: "file", image_path: null, tags: [],
  }) as unknown as Exposure);
}

describe.each(cases as Case[])("representative frame: $case", (c) => {
  it("is what the shared rule picks", () => {
    expect(representativeExposure(exposures(c))?.id ?? null).toBe(c.representative);
  });

});

describe("keptExposures", () => {
  it("drops rejected frames and keeps accepted and unmarked ones", () => {
    const out = keptExposures([{ status: "accepted" }, { status: "rejected" }, { status: null }]);
    expect(out.map((e) => e.status)).toEqual(["accepted", null]);
  });

  it("treats a not-yet-loaded list as empty", () => {
    expect(keptExposures(undefined)).toEqual([]);
  });
});
