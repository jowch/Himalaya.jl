import { create } from "zustand";

/** File-pattern overrides for the Configuration first-run step.
 *  Shape mirrors `CreateExperimentBody.patterns` (api.ts). */
export interface DraftPatterns {
  image?: string;
  metadata?: string;
  integration?: string;
}

/**
 * useDraftExperiment — client-side-only ephemeral store for the two-phase
 * ingest funnel (spec §2/§6.1). The picker (NewExperimentPage) commits a
 * validated path here; ConfigurationPage in first-run mode reads it.
 * No DB row is created until the user approves at the Configuration step.
 */
export interface DraftExperimentState {
  /** Validated directory path chosen in the picker, or "" if no draft. */
  path: string;
  /** File-pattern overrides (image/metadata/integration globs). */
  patterns: DraftPatterns;
  /** Commit a validated path (and optionally patterns) into the draft slot. */
  setDraft: (draft: { path: string; patterns?: DraftPatterns }) => void;
  /** Clear the draft (on Cancel or after successful experiment creation). */
  clear: () => void;
}

export const useDraftExperiment = create<DraftExperimentState>((set) => ({
  path: "",
  patterns: {},
  setDraft: ({ path, patterns }) =>
    set({ path, patterns: patterns ?? {} }),
  clear: () => set({ path: "", patterns: {} }),
}));
