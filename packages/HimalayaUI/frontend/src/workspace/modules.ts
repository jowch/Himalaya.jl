import type { SetupMode } from "../api";

export interface WorkspaceModules {
  peakEditing: boolean;
  calls: boolean;
  indexing: boolean;
  imageSelection: boolean;
  review: boolean;
}

const ALL: WorkspaceModules = { peakEditing: true, calls: true, indexing: true, imageSelection: true, review: true };
const NONE: WorkspaceModules = { peakEditing: false, calls: false, indexing: false, imageSelection: false, review: false };

export const NO_MODULES: WorkspaceModules = NONE;

export function modulesFor(mode: SetupMode | null | undefined): WorkspaceModules {
  return mode === "autosampler" ? NONE : ALL;
}
