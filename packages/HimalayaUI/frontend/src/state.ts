import { create } from "zustand";
import { persist } from "zustand/middleware";
import type { IngestStatus } from "./api";

export const LS_KEY = "himalaya-ui:state";

const RETIRED_KEYS = [
  "activePage", "theme", "activeExperimentId", "activeSampleId", "activeExposureId",
] as const;

export type IngestProgressStatus = Exclude<IngestStatus, "idle">;
export interface IngestProgress {
  processed: number;
  total: number;
  status: IngestProgressStatus;
  stage?: string | undefined;
}

export interface AppState {
  username: string | undefined;
  firstName: string | undefined;
  lastName: string | undefined;
  tutorialSeen: boolean;
  tutorialPendingFor: string | undefined;

  ingestInFlight: Record<number, IngestProgress> | null;
  paletteOpen: boolean;

  setUsername: (name: string) => void;
  setUser: (u: { username: string; firstName?: string | undefined; lastName?: string | undefined }) => void;
  setIngestProgress: (experimentId: number, progress: IngestProgress) => void;
  clearIngestProgress: (experimentId: number) => void;
  setTutorialSeen: (seen: boolean) => void;
  setTutorialPendingFor: (handle: string | undefined) => void;
  openPalette: () => void;
  closePalette: () => void;
  clearUsername: () => void;
}

export const useAppState = create<AppState>()(
  persist(
    (set) => ({
      username: undefined,
      firstName: undefined,
      lastName: undefined,
      tutorialSeen: false,
      tutorialPendingFor: undefined,
      ingestInFlight: null,
      paletteOpen: false,

      setUsername: (username) => set({ username }),
      setUser: ({ username, firstName, lastName }) =>
        set({ username, firstName, lastName }),
      setIngestProgress: (experimentId, progress) =>
        set((s) => ({
          ingestInFlight: { ...(s.ingestInFlight ?? {}), [experimentId]: progress },
        })),
      clearIngestProgress: (experimentId) =>
        set((s) => {
          if (s.ingestInFlight === null) return {};
          const next = { ...s.ingestInFlight };
          delete next[experimentId];
          return { ingestInFlight: Object.keys(next).length === 0 ? null : next };
        }),
      setTutorialSeen: (tutorialSeen) => set({ tutorialSeen }),
      setTutorialPendingFor: (tutorialPendingFor) => set({ tutorialPendingFor }),
      openPalette: () => set({ paletteOpen: true }),
      closePalette: () => set({ paletteOpen: false }),
      clearUsername: () => set({ username: undefined, firstName: undefined, lastName: undefined }),
    }),
    {
      name: LS_KEY,
      version: 6,
      partialize: (s) => ({
        username: s.username,
        firstName: s.firstName,
        lastName: s.lastName,
        tutorialSeen: s.tutorialSeen,
        tutorialPendingFor: s.tutorialPendingFor,
      }),
      migrate: (persisted) => {
        if (persisted && typeof persisted === "object") {
          const rest = { ...(persisted as Record<string, unknown>) };
          for (const key of RETIRED_KEYS) delete rest[key];
          return rest as unknown as AppState;
        }
        return persisted as AppState;
      },
      merge: mergePersistedState,
    },
  ),
);

export function mergePersistedState(
  persisted: unknown,
  current: AppState,
): AppState {
  return { ...current, ...(persisted as Partial<AppState> | undefined) };
}
