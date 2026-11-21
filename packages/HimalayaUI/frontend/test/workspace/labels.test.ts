import { describe, it, expect } from "vitest";
import type { SampleLabel } from "../../src/api";
import {
  applyAddLabel,
  bufferRole,
  buffersFirst,
  applyRemoveLabel,
  applyTagColor,
  defaultTagColor,
  cleanTag,
  displayTag,
  labelsFor,
  sampleMatchesFilter,
  tagColor,
  tagColorIsDefault,
  tagOptions,
  tagVocabulary,
} from "../../src/lib/labels";

function label(id: number, sample_id: number, text: string, state: SampleLabel["state"] = "confirmed"): SampleLabel {
  return { id, sample_id, label: text, state, source: "user" };
}

describe("tag vocabulary", () => {
  it("lists tags used in the experiment, most used first, then the standard tags not yet used", () => {
    const list = [label(1, 1, "37C"), label(2, 2, "retake"), label(3, 3, "37c"), label(4, 4, "Buffer")];
    expect(tagVocabulary(list)).toEqual(["37C", "Buffer", "retake", "Calibrant", "Empty capillary"]);
  });

  it("shows a built-in tag capitalized however it was typed", () => {
    expect(tagVocabulary([label(1, 1, "buffer"), label(2, 2, "EMPTY capillary")])).toEqual(["Buffer", "Empty capillary", "Calibrant"]);
    expect(displayTag("calibrant")).toBe("Calibrant");
    expect(displayTag("37C")).toBe("37C");
  });

  it("leaves suggestions out of the vocabulary", () => {
    expect(tagVocabulary([label(1, 1, "odd", "suggested")])).toEqual(["Buffer", "Calibrant", "Empty capillary"]);
  });
});

describe("tag options", () => {
  const vocab = ["37C", "retake", "buffer", "calibrant", "empty capillary"];

  it("offers every tag the sample lacks when nothing is typed", () => {
    expect(tagOptions(vocab, "", ["retake"])).toEqual({ matches: ["37C", "buffer", "calibrant", "empty capillary"], create: undefined });
  });

  it("puts prefix matches first and offers no Create row for a known tag", () => {
    expect(tagOptions(vocab, "Buf", [])).toEqual({ matches: ["buffer"], create: "Buf" });
    expect(tagOptions(vocab, "BUFFER", [])).toEqual({ matches: ["buffer"], create: undefined });
  });

  it("offers buffer for a misspelling, ahead of creating the misspelling", () => {
    const o = tagOptions(vocab, "bufer", []);
    expect(o.matches[0]).toBe("buffer");
    expect(o.create).toBe("bufer");
    expect(tagOptions(vocab, "calbrant", []).matches).toEqual(["calibrant"]);
  });

  it("matches inside a tag and ignores tags already on the sample", () => {
    expect(tagOptions(vocab, "capill", []).matches).toEqual(["empty capillary"]);
    expect(tagOptions(vocab, "buffer", ["Buffer"])).toEqual({ matches: [], create: undefined });
  });

  it("offers no Create row for text that cannot be a tag", () => {
    expect(tagOptions(vocab, "   ", []).create).toBeUndefined();
    expect(tagOptions(vocab, "x".repeat(65), []).create).toBeUndefined();
    expect(cleanTag("  hot  ")).toBe("hot");
    expect(cleanTag("a\nb")).toBeUndefined();
  });
});

describe("label list edits", () => {
  it("adding shows the tag at once and confirms an existing suggestion instead of duplicating it", () => {
    const list = [label(1, 1, "buffer", "suggested")];
    expect(applyAddLabel(list, 1, "Buffer")).toEqual([label(1, 1, "buffer", "confirmed")]);
    const added = applyAddLabel(list, 2, " 37C ");
    expect(labelsFor(added, 2)).toEqual([{ id: -1, sample_id: 2, label: "37C", state: "confirmed", source: "user" }]);
  });

  it("removing drops the tag regardless of case and leaves other samples alone", () => {
    const list = [label(1, 1, "retake"), label(2, 2, "retake")];
    expect(applyRemoveLabel(list, 1, "RETAKE")).toEqual([label(2, 2, "retake")]);
  });
});

describe("sample filter", () => {
  it("matches the sample name or any of its tags", () => {
    const tags = [label(1, 1, "empty capillary")];
    expect(sampleMatchesFilter("MO_H2O", tags, "capill")).toBe(true);
    expect(sampleMatchesFilter("MO_H2O", tags, "h2o")).toBe(true);
    expect(sampleMatchesFilter("MO_H2O", tags, "buffer")).toBe(false);
    expect(sampleMatchesFilter("MO_H2O", [], "  ")).toBe(true);
  });
});

describe("tag colours", () => {
  it("built-ins default to blue and other tags to violet, and a stored colour wins without case", () => {
    expect(tagColor([], "BUFFER")).toBe("blue");
    expect(tagColor([], "37C")).toBe("violet");
    expect(tagColor([{ label: "buffer", color: "stone" }], "Buffer")).toBe("stone");
    expect(tagColorIsDefault([{ label: "37c", color: "teal" }], "37C")).toBe(false);
    expect(tagColorIsDefault([], "37C")).toBe(true);
  });

  it("matches names ignoring ASCII case only, as the database does", () => {
    expect(tagColor([{ label: "abc", color: "rose" }], "ABC")).toBe("rose");
    expect(tagColor([{ label: "über", color: "rose" }], "Über")).toBe("violet");
    expect(defaultTagColor("Buffer")).toBe("blue");
  });

  it("setting replaces the entry for that name and null removes it", () => {
    const list = [{ label: "37C", color: "teal" as const }];
    expect(applyTagColor(list, " 37c ", "rose")).toEqual([{ label: "37c", color: "rose" }]);
    expect(applyTagColor(list, "37c", null)).toEqual([]);
  });
});

describe("buffers first", () => {
  const label = (sample_id: number, l: string, state: "confirmed" | "suggested"): SampleLabel =>
    ({ id: sample_id * 10, sample_id, label: l, state, source: state === "suggested" ? "ingest" : "user" });
  const samples = [1, 2, 3, 4, 5].map((id) => ({ id, name: `S${id}` }));
  const labels = [label(2, "buffer", "suggested"), label(4, "Buffer", "confirmed"), label(5, "Calibrant", "confirmed"), label(3, "37C", "confirmed")];

  it("puts confirmed buffers, then suggested buffers, ahead of every other sample, keeping list order within each", () => {
    expect(buffersFirst(samples, labels).map((s) => s.id)).toEqual([4, 2, 1, 3, 5]);
  });

  it("leaves the order alone when no sample is tagged Buffer", () => {
    expect(buffersFirst(samples, [label(5, "Calibrant", "confirmed")]).map((s) => s.id)).toEqual([1, 2, 3, 4, 5]);
  });

  it("says whether a sample's Buffer tag is confirmed or only suggested", () => {
    expect([1, 2, 4].map((id) => bufferRole(labels, id))).toEqual([undefined, "suggested", "confirmed"]);
  });
});
