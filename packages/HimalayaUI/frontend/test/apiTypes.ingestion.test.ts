// test/apiTypes.ingestion.test.ts
import { describe, it, expect } from "vitest";
import type { Experiment, Sample, Load, GroupingFlag } from "../src/api";

describe("ingestion api types (Phase E1)", () => {
  it("Experiment carries typed geometry source tags + scan columns + E1 additive fields", () => {
    const e: Experiment = {
      id: 1, name: "SSRL", description: "AgBe SAXS run", path: "/d",
      data_dir: "/d/data", analysis_dir: "/d/an",
      manifest_path: null, created_at: "2026-04-12", q_units: "A^-1",
      beam_center_x: 421.4, beam_center_y: 836.9, pixel_size_um: 172, energy_kev: 9,
      flight_path_m: 1.8095,
      energy_kev_source: "prp", flight_path_m_source: "setup",
      beam_center_x_source: "setup", beam_center_y_source: "setup",
      pixel_size_um_source: "prp", q_units_source: "prp",
      last_scanned_at: "2026-04-12T10:00:00", scan_signature: "sig", ingest_status: "idle",
      image_pattern: "*.tif", metadata_pattern: null, integration_pattern: null,
    };
    expect(e.flight_path_m_source).toBe("setup");
    expect(e.ingest_status).toBe("idle");
    expect(e.description).toBe("AgBe SAXS run");
    expect(e.image_pattern).toBe("*.tif");
  });

  it("Sample has a single non-null `name` label and no `display_name`", () => {
    const s: Sample = {
      id: 1, experiment_id: 1, name: "HA85 (S01P15)", notes: null, tags: [],
    };
    expect(s.name).toBe("HA85 (S01P15)");
    // @ts-expect-error display_name is removed from Sample
    expect(s.display_name).toBeUndefined();
    // @ts-expect-error name is non-null after the collapse — null is not assignable
    const _nn: Sample = { id: 2, experiment_id: 1, name: null, notes: null, tags: [] };
    void _nn;
  });

  it("Load is the NESTED roll-up (Load ▸ Sample ▸ Exposures)", () => {
    const l: Load = {
      load_id: 3, load_index: 1, session_id: 1,
      start_time: "2026-04-12T10:00:00", end_time: "2026-04-12T10:10:00",
      frame_count: 64, note: null,
      samples: [
        {
          sample_id: 9, name: "HA85 (S01P15)", slot_index: 15,
          grouping_source: "computed", name_source: "computed",
          merged_into_id: null, flag: null,
          exposures: [
            { id: 51, filename: "ha85_001.tif", horizontal_position: 12.4, timestamp: "2026-04-12T10:00:01" },
          ],
        },
      ],
    };
    expect(l.load_index).toBe(1);
    expect(l.samples[0]!.exposures[0]!.id).toBe(51);
    expect(l.samples[0]!.flag).toBeNull();
  });

  it("a Load sample can carry a merge or split GroupingFlag", () => {
    const merge: GroupingFlag = { kind: "merge", merge_with_sample_id: 4, merge_with_label: "HA85 (S01P14)" };
    const split: GroupingFlag = { kind: "split", split_at_index: 32, jump_from: 12.4, jump_to: 48.1 };
    expect(merge.kind).toBe("merge");
    expect(split.kind).toBe("split");
  });
});
