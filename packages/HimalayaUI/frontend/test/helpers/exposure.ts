import type { Exposure } from "../../src/api";

export function makeExposure(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 1,
    sample_id: 1,
    filename: null,
    kind: "file",
    selected: false,
    status: null,
    image_path: null,
    image_version: "",
    sources: [],
    trace_hash: null,
    analysis_inputs_hash: null,
    exposure_time: null,
    counting_time: null,
    transmitted_counts: null,
    incident_monitor: null,
    timestamp: null,
    beam_energy_ev: null,
    pipe_length_m: null,
    detector: null,
    setup_id: null,
    ...over,
  };
}
