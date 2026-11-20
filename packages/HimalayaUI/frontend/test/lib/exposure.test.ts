import { describe, it, expect } from "vitest";
import cases from "../fixtures/representative-exposure.json";
import { keptExposures, representativeExposure } from "../../src/lib/exposure";
import { resolveActiveExposure } from "../../src/hooks/useAutoPickExposure";
import { defaultExposureId } from "../../src/print/pages/loupeAdapters";
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

  it("is the frame Focus opens on a cold arrival", () => {
    expect(resolveActiveExposure(undefined, exposures(c)) ?? null).toBe(c.representative);
  });

  it("is the frame Focus moves to when the stored frame is from another sample", () => {
    expect(resolveActiveExposure(9999, exposures(c)) ?? null).toBe(c.representative);
  });

  it("is the frame Loupe opens on", () => {
    const list = exposures(c);
    const expected = c.representative ?? (list.length ? Math.max(...list.map((e) => e.id)) : null);
    expect(defaultExposureId(list) ?? null).toBe(expected);
    expect(defaultExposureId([...list].reverse()) ?? null).toBe(expected);
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
