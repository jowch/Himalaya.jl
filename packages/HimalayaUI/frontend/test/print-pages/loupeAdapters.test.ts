import { describe, it, expect } from "vitest";
import type { Exposure, SampleTag } from "../../src/api";
import {
  defaultExposureId,
  buildExposureImageUrl,
  toGalleryExposures,
  toMetaEntries,
  toLoupeTags,
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
  it("maps id/src(thumb)/frameNo/rejected/kept/representative", () => {
    const out = toGalleryExposures([
      exp({ id: 10, selected: true, status: "accepted" }),
      exp({ id: 11, status: "rejected" }),
    ]);
    expect(out[0]).toEqual({ id: 10, src: "/api/exposures/10/image?thumb=1&v=v9", frameNo: 1, representative: true, rejected: false, kept: true });
    expect(out[1]).toEqual({ id: 11, src: "/api/exposures/11/image?thumb=1&v=v9", frameNo: 2, representative: false, rejected: true, kept: false });
  });

  it("keeps the SA-SCREENED tri-state: an unscreened (null) frame is neither kept nor rejected", () => {
    const out = toGalleryExposures([exp({ id: 12, status: null })]);
    expect(out[0]).toMatchObject({ kept: false, rejected: false });
  });
});

describe("toMetaEntries", () => {
  it("builds the frame position and OMITS the integration/collected rows (no API fields, controls-don't-lie)", () => {
    const exposures = [exp({ id: 1 }), exp({ id: 2 }), exp({ id: 3 })];
    expect(toMetaEntries(exposures[1]!, exposures)).toEqual([
      { key: "frame", value: "2 of 3" },
    ]);
  });

  it("shows the rejection reason row on a dropped frame that carries a rejection_reason tag", () => {
    const dropped = exp({
      id: 2,
      status: "rejected",
      tags: [{ id: 9, key: "rejection_reason", value: "beam flare", source: "manual" }],
    });
    expect(toMetaEntries(dropped, [exp({ id: 1 }), dropped])).toEqual([
      { key: "frame", value: "2 of 2" },
      { key: "reason", value: "beam flare" },
    ]);
  });

  it("joins multiple rejection_reason tags", () => {
    const dropped = exp({
      id: 1,
      status: "rejected",
      tags: [
        { id: 9, key: "rejection_reason", value: "flare", source: "manual" },
        { id: 10, key: "rejection_reason", value: "no signal", source: "manual" },
      ],
    });
    expect(toMetaEntries(dropped, [dropped])).toEqual([
      { key: "frame", value: "1 of 1" },
      { key: "reason", value: "flare; no signal" },
    ]);
  });

  it("no reason row on a dropped frame without the tag, with an empty-value tag, or on a kept frame carrying one", () => {
    const droppedBare = exp({ id: 1, status: "rejected" });
    expect(toMetaEntries(droppedBare, [droppedBare])).toEqual([
      { key: "frame", value: "1 of 1" },
    ]);
    const droppedEmpty = exp({
      id: 1,
      status: "rejected",
      tags: [{ id: 9, key: "rejection_reason", value: "  ", source: "manual" }],
    });
    expect(toMetaEntries(droppedEmpty, [droppedEmpty])).toEqual([
      { key: "frame", value: "1 of 1" },
    ]);
    const kept = exp({
      id: 1,
      status: "accepted",
      tags: [{ id: 9, key: "rejection_reason", value: "stale", source: "manual" }],
    });
    expect(toMetaEntries(kept, [kept])).toEqual([
      { key: "frame", value: "1 of 1" },
    ]);
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
