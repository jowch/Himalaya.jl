import { describe, it, expect } from "vitest";
import type { Exposure, SampleTag } from "../../src/api";
import {
  defaultExposureId,
  buildExposureImageUrl,
  toGalleryExposures,
  toMetaEntries,
  toLoupeTags,
  findSampleTagId,
} from "../../src/print/pages/loupeAdapters";

function exp(over: Partial<Exposure>): Exposure {
  return {
    id: 1, sample_id: 1, filename: "JC000-001.dat", kind: "file",
    selected: false, status: "accepted", image_path: "/x.tif", image_version: "v9",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null, ...over,
  };
}

describe("defaultExposureId", () => {
  it("prefers the representative (selected)", () => {
    expect(defaultExposureId([exp({ id: 1 }), exp({ id: 2, selected: true })])).toBe(2);
  });
  it("falls back to first accepted", () => {
    expect(defaultExposureId([exp({ id: 1, status: "rejected" }), exp({ id: 2, status: "accepted" })])).toBe(2);
  });
  it("falls back to first exposure", () => {
    expect(defaultExposureId([exp({ id: 7, status: null }), exp({ id: 8, status: null })])).toBe(7);
  });
  it("returns undefined for an empty list", () => {
    expect(defaultExposureId([])).toBeUndefined();
  });
});

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
  it("maps id/src(thumb)/frameNo/rejected/representative", () => {
    const out = toGalleryExposures([
      exp({ id: 10, selected: true, status: "accepted" }),
      exp({ id: 11, status: "rejected" }),
    ]);
    expect(out[0]).toEqual({ id: 10, src: "/api/exposures/10/image?thumb=1&v=v9", frameNo: 1, representative: true, rejected: false });
    expect(out[1]).toEqual({ id: 11, src: "/api/exposures/11/image?thumb=1&v=v9", frameNo: 2, representative: false, rejected: true });
  });
});

describe("toMetaEntries", () => {
  it("builds frame position + the deferred-field placeholders, no signal row", () => {
    const exposures = [exp({ id: 1 }), exp({ id: 2 }), exp({ id: 3 })];
    expect(toMetaEntries(exposures[1]!, exposures)).toEqual([
      { key: "frame", value: "2 of 3" },
      { key: "integration", value: "—" },
      { key: "collected", value: "—" },
    ]);
  });
});

describe("toLoupeTags / findSampleTagId", () => {
  const tags: SampleTag[] = [
    { id: 100, key: "LL37", value: "", source: "user" },
    { id: 101, key: "temp", value: "37C", source: "user" },
  ];
  it("maps SampleTag -> Tag, omitting empty value (exactOptionalPropertyTypes)", () => {
    expect(toLoupeTags(tags)).toEqual([{ key: "LL37" }, { key: "temp", value: "37C" }]);
  });
  it("finds the tag id by key+value for removal", () => {
    expect(findSampleTagId(tags, { key: "temp", value: "37C" })).toBe(101);
    expect(findSampleTagId(tags, { key: "LL37" })).toBe(100);
    expect(findSampleTagId(tags, { key: "missing" })).toBeUndefined();
  });
});
