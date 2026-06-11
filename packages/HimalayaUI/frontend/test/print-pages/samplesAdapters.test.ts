import { toSampleRowModel } from "../../src/print/pages/samplesAdapters";

it("derives kept/total/dropped/screened/tags/phase from a sample + exposures", () => {
  const sample = {
    id: 9,
    experiment_id: 1,
    name: "JC009",
    display_name: "JC009 — LL37",
    notes: null,
    q_units: "nm^-1",
    tags: [{ id: 1, key: "LL37", value: "", source: "manual" }],
    phase: null,
  } as any;
  const exposures = [
    { id: 101, status: "accepted", selected: true, image_path: null, image_version: "" },
    { id: 102, status: "rejected", selected: false, image_path: null, image_version: "" },
  ] as any;
  const m = toSampleRowModel(sample, exposures);
  expect(m.sampleId).toBe("#9");
  expect(m.kept).toBe(1);
  expect(m.total).toBe(2);
  expect(m.dropped).toBe(1);
  expect(m.exposures).toHaveLength(2);
  expect(m.tags).toEqual([{ key: "LL37" }]);
  expect(m.phase).toBeNull();
});

// SA-SCREENED: "Kept" on the sheet means EXPLICITLY accepted — the same truth
// the loupe caption and the hero metric tell. An unscreened (null) frame is
// neither kept nor dropped.
it("a null-status frame counts neither kept nor dropped; accepted counts kept", () => {
  const sample = {
    id: 9, name: "JC009", display_name: null, notes: null,
    q_units: "nm^-1", tags: [], phase: null,
  } as any;
  const exposures = [
    { id: 101, status: "accepted", selected: true, image_path: null, image_version: "" },
    { id: 102, status: null, selected: false, image_path: null, image_version: "" },
    { id: 103, status: "rejected", selected: false, image_path: null, image_version: "" },
  ] as any;
  const m = toSampleRowModel(sample, exposures);
  expect(m.total).toBe(3);
  expect(m.kept).toBe(1);    // only the accepted frame
  expect(m.dropped).toBe(1); // only the rejected frame — null is NOT dropped
});

it("treats undefined exposures as not-yet-loaded (empty derivation)", () => {
  const sample = {
    id: 9,
    name: null,
    display_name: null,
    notes: null,
    q_units: "nm^-1",
    tags: [],
    phase: undefined,
  } as any;
  const m = toSampleRowModel(sample, undefined);
  expect(m.total).toBe(0);
  expect(m.kept).toBe(0);
  expect(m.dropped).toBe(0);
});
