// test/apiTypes.ingestion.test.ts
import { describe, it, expect } from "vitest";
import type { Experiment, Sample } from "../src/api";

describe("ingestion api types (Phase E1)", () => {
  it("Experiment keeps its folder configuration on its setups", () => {
    const e: Experiment = {
      id: 1, name: "SSRL", description: "AgBe SAXS run", path: "/d",
      created_at: "2026-04-12", last_scanned_at: "2026-04-12T10:00:00", ingest_status: "idle",
      setups: [{
        id: 3, experiment_id: 1, folder: "setup1", display_name: null, name: "setup1",
        position: 0, included: true, root_dir: "/d/setup1",
        data_dir: "/d/setup1/data", analysis_dir: "/d/setup1/analysis",
        image_pattern: "*.tif", metadata_pattern: null, integration_pattern: null,
        q_units: "A^-1", beam_center_x: 421.4, beam_center_y: 836.9, pixel_size_um: 172,
        energy_kev: 9, flight_path_m: 1.8095,
        energy_kev_source: "prp", flight_path_m_source: "setup",
        beam_center_x_source: "setup", beam_center_y_source: "setup",
        pixel_size_um_source: "prp", q_units_source: "prp",
      }],
    };
    expect(e.setups?.[0].flight_path_m_source).toBe("setup");
    expect(e.ingest_status).toBe("idle");
    expect(e.description).toBe("AgBe SAXS run");
    // @ts-expect-error folder configuration moved to setups
    expect(e.image_pattern).toBeUndefined();
  });

  it("Sample has a single non-null `name` label and no `display_name`", () => {
    const s: Sample = {
      id: 1, experiment_id: 1, name: "HA85 (S01P15)", notes: null,
    };
    expect(s.name).toBe("HA85 (S01P15)");
    // @ts-expect-error display_name is removed from Sample
    expect(s.display_name).toBeUndefined();
    // @ts-expect-error name is non-null after the collapse — null is not assignable
    const _nn: Sample = { id: 2, experiment_id: 1, name: null, notes: null };
    void _nn;
  });
});
