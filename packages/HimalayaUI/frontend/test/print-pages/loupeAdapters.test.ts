import { describe, it, expect } from "vitest";
import type { Exposure, SampleTag } from "../../src/api";
import {
  buildExposureImageUrl,
  toGalleryExposures,
  toLoupeTags,
} from "../../src/print/pages/loupeAdapters";

function exp(over: Partial<Exposure>): Exposure {
  return {
    id: 1, sample_id: 1, filename: "JC000-001.dat", kind: "file",
    selected: false, status: null, image_path: "/x.tif", image_version: "v9",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null, ...over,
  };
}

describe("buildExposureImageUrl", () => {
  it("builds a full image url with the version param", () => {
    expect(buildExposureImageUrl(exp({ id: 5, image_version: "v9" }))).toBe("/api/exposures/5/image?v=v9");
  });
  it("builds a thumb url", () => {
    expect(buildExposureImageUrl(exp({ id: 5, image_version: "v9" }), { thumb: true })).toBe("/api/exposures/5/image?thumb=1&v=v9");
  });
  it("omits ?v when image_version is empty", () => {
    expect(buildExposureImageUrl(exp({ id: 5, image_version: "" }))).toBe("/api/exposures/5/image");
  });
  it("returns null when there is no image", () => {
    expect(buildExposureImageUrl(exp({ id: 5, image_path: null }))).toBeNull();
  });
});

describe("toGalleryExposures", () => {
  it("maps id/src(thumb)/frameNo/rejected/representative (binary: kept needs no marker)", () => {
    const out = toGalleryExposures([
      exp({ id: 10, selected: true, status: null }),
      exp({ id: 11, status: "rejected" }),
    ]);
    expect(out[0]).toEqual({ id: 10, src: "/api/exposures/10/image?thumb=1&v=v9", frameNo: 1, representative: true, rejected: false });
    expect(out[1]).toEqual({ id: 11, src: "/api/exposures/11/image?thumb=1&v=v9", frameNo: 2, representative: false, rejected: true });
  });

  it("marks the representative the views open on, not a dropped selected frame", () => {
    const out = toGalleryExposures([
      exp({ id: 10, status: null }),
      exp({ id: 11, selected: true, status: "rejected" }),
      exp({ id: 12, status: null }),
    ]);
    expect(out.map((g) => g.representative)).toEqual([false, false, true]);
  });

  it("a kept (null-status) frame is simply not-rejected", () => {
    const out = toGalleryExposures([exp({ id: 12, status: null })]);
    expect(out[0]).toMatchObject({ rejected: false });
  });
});

describe("toLoupeTags (LO-TAGDUP)", () => {
  it("carries the tag id and source so byte-identical duplicates are distinguishable", () => {
    const tags: SampleTag[] = [
      { id: 3, key: "dose", value: "10", source: "manual" },
      { id: 7, key: "dose", value: "10", source: "scoping" },
    ];
    const out = toLoupeTags(tags);
    expect(out).toEqual([
      { id: 3, key: "dose", value: "10", source: "manual" },
      { id: 7, key: "dose", value: "10", source: "scoping" },
    ]);
    // The two pills are now distinguishable by id even though key+value match.
    expect(out[0]!.id).not.toBe(out[1]!.id);
  });

  it("omits an empty value (keeps the exactOptionalPropertyTypes contract)", () => {
    const out = toLoupeTags([{ id: 1, key: "k", value: "", source: "manual" }]);
    // value key is OMITTED (not set to undefined) per exactOptionalPropertyTypes.
    expect(out[0]).toEqual({ id: 1, key: "k", source: "manual" });
    expect("value" in out[0]!).toBe(false);
  });
});
