import { create } from "zustand";
import { persist } from "zustand/middleware";

export const LS_KEY = "himalaya-ui:state";

export type PageId = "index" | "compare";
export type ThemeId = "dark" | "light";
export type NavModalStep = "experiment" | "sample";

export interface AppState {
  // persisted
  username: string | undefined;
  activeExperimentId: number | undefined;
  activeSampleId: number | undefined;
  activeExposureId: number | undefined;
  activePage: PageId;
  tutorialSeen: boolean;
  theme: ThemeId;

  // ephemeral (not persisted)
  hoveredIndexId: number | undefined;
  navModalOpen: boolean;
  navModalStep: NavModalStep;

  // setters
  setUsername: (name: string) => void;
  setActiveExperiment: (id: number | undefined) => void;
  setActiveSample: (id: number | undefined) => void;
  setActiveExposure: (id: number | undefined) => void;
  setHoveredIndex: (id: number | undefined) => void;
  setActivePage: (page: PageId) => void;
  setTutorialSeen: (seen: boolean) => void;
  setTheme: (theme: ThemeId) => void;
  openNavModal: (step?: NavModalStep) => void;
  closeNavModal: () => void;
  setNavModalStep: (step: NavModalStep) => void;
  clearUsername: () => void;
}

export const useAppState = create<AppState>()(
  persist(
    (set) => ({
      username: undefined,
      activeExperimentId: undefined,
      activeSampleId: undefined,
      activeExposureId: undefined,
      activePage: "index",
      tutorialSeen: false,
      theme: "dark",

      hoveredIndexId: undefined,
      navModalOpen: false,
      navModalStep: "experiment",

      setUsername: (username) => set({ username }),
      setActiveExperiment: (activeExperimentId) =>
        set({
          activeExperimentId,
          activeSampleId: undefined,
          activeExposureId: undefined,
        }),
      setActiveSample: (activeSampleId) =>
        set({ activeSampleId, activeExposureId: undefined }),
      setActiveExposure: (activeExposureId) => set({ activeExposureId }),
      setHoveredIndex: (hoveredIndexId) => set({ hoveredIndexId }),
      setActivePage: (activePage) => set({ activePage }),
      setTutorialSeen: (tutorialSeen) => set({ tutorialSeen }),
      setTheme: (theme) => set({ theme }),
      openNavModal: (step) =>
        set(step ? { navModalOpen: true, navModalStep: step } : { navModalOpen: true }),
      closeNavModal: () => set({ navModalOpen: false }),
      setNavModalStep: (navModalStep) => set({ navModalStep }),
      clearUsername: () => set({ username: undefined }),
    }),
    {
      name: LS_KEY,
      version: 2,
      partialize: (s) => ({
        username: s.username,
        activeExperimentId: s.activeExperimentId,
        activeSampleId: s.activeSampleId,
        activeExposureId: s.activeExposureId,
        activePage: s.activePage,
        tutorialSeen: s.tutorialSeen,
        theme: s.theme,
      }),
    },
  ),
);
