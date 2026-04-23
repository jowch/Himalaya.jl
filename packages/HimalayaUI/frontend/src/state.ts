import { create } from "zustand";
import { persist } from "zustand/middleware";

export const LS_KEY = "himalaya-ui:state";

export interface AppState {
  username: string | undefined;
  activeSampleId: number | undefined;
  activeExposureId: number | undefined;
  hoveredIndexId: number | undefined;
  setUsername: (name: string) => void;
  setActiveSample: (id: number | undefined) => void;
  setActiveExposure: (id: number | undefined) => void;
  setHoveredIndex: (id: number | undefined) => void;
}

export const useAppState = create<AppState>()(
  persist(
    (set) => ({
      username: undefined,
      activeSampleId: undefined,
      activeExposureId: undefined,
      hoveredIndexId: undefined,
      setUsername: (username) => set({ username }),
      setActiveSample: (activeSampleId) =>
        set({ activeSampleId, activeExposureId: undefined }),
      setActiveExposure: (activeExposureId) => set({ activeExposureId }),
      setHoveredIndex: (hoveredIndexId) => set({ hoveredIndexId }),
    }),
    {
      name: LS_KEY,
      version: 1,
      partialize: (s) => ({
        username: s.username,
        activeSampleId: s.activeSampleId,
        activeExposureId: s.activeExposureId,
      }),
    },
  ),
);
