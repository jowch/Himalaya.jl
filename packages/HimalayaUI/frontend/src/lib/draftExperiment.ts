import { create } from "zustand";
import type { ResolveLayoutResponse } from "../api";

/** File-pattern overrides for the Configuration first-run step.
 *  Shape mirrors `CreateExperimentBody.patterns` (api.ts). */
export interface DraftPatterns {
  image?: string;
  metadata?: string;
  integration?: string;
}

/** First-run geometry overrides. A field present here = the user edited it on the
 *  review screen; it is sent to create as source='user'. Shape mirrors
 *  `CreateExperimentBody.geometry`. Absent fields are auto-derived by the scan. */
export interface DraftGeometry {
  beam_center_x?: number;
  beam_center_y?: number;
  flight_path_m?: number;
  pixel_size_um?: number;
  energy_kev?: number;
}

/**
 * useDraftExperiment — client-side-only ephemeral store for the two-phase
 * ingest funnel. The picker (NewExperimentPage) commits the experiment ROOT via
 * `setRoot`; ConfigurationPage first-run resolves that root structurally
 * (`/api/fs/resolve` → `applyResolved`) into auto-discovered, user-correctable
 * defaults (name / data_dir / analysis_dir / setup_file), then the user edits
 * them with `patch` before Approve. No DB row is created until Approve.
 */
export interface DraftExperimentState {
  /** The picked experiment root directory ("" = no draft). */
  root: string;
  /** Whether resolveLayout has populated the fields for the current root. */
  resolved: boolean;
  /** User-supplied experiment name (prefilled from the resolver). */
  name: string;
  /** Resolved/corrected data directory (raw images). */
  data_dir: string;
  /** Resolved/corrected analysis directory (integration .dat). "" = none. */
  analysis_dir: string;
  /** Resolved/corrected geometry setup file. "" = none. */
  setup_file: string;
  /** True when setup discovery found none or multiple (the user should confirm). */
  setup_ambiguous: boolean;
  /** File-pattern overrides (image/metadata/integration globs). */
  patterns: DraftPatterns;
  /** Geometry overrides the user edited on the review screen (empty = all auto). */
  geometry: DraftGeometry;

  /** Picker: commit a root; resets the resolved fields (config re-resolves). */
  setRoot: (root: string) => void;
  /** Config: fill the correctable fields from a structural resolve. */
  applyResolved: (r: ResolveLayoutResponse) => void;
  /** Config: the user corrects a field in place. */
  patch: (
    fields: Partial<
      Pick<DraftExperimentState, "name" | "data_dir" | "analysis_dir" | "setup_file" | "patterns" | "geometry">
    >,
  ) => void;
  /** Clear the draft (on Cancel or after successful creation). */
  clear: () => void;
}

const EMPTY = {
  root: "",
  resolved: false,
  name: "",
  data_dir: "",
  analysis_dir: "",
  setup_file: "",
  setup_ambiguous: false,
  patterns: {} as DraftPatterns,
  geometry: {} as DraftGeometry,
};

export const useDraftExperiment = create<DraftExperimentState>((set) => ({
  ...EMPTY,
  setRoot: (root) => set({ ...EMPTY, root }),
  applyResolved: (r) =>
    set({
      resolved: true,
      name: r.name,
      data_dir: r.data_dir,
      analysis_dir: r.analysis_dir ?? "",
      setup_file: r.setup_file ?? "",
      setup_ambiguous: r.setup_ambiguous,
      geometry: {},   // geometry overrides are manifest-time; start clean per root
      // Seed detected patterns (e.g. the SSRL tot_files convention). When the
      // resolver returns null (undetected), keep patterns empty so the funnel's
      // `{name}.*` defaults apply.
      patterns: {
        ...(r.image_pattern ? { image: r.image_pattern } : {}),
        ...(r.metadata_pattern ? { metadata: r.metadata_pattern } : {}),
        ...(r.integration_pattern ? { integration: r.integration_pattern } : {}),
      },
    }),
  patch: (fields) => set((s) => ({ ...s, ...fields })),
  clear: () => set({ ...EMPTY }),
}));
