import { create } from "zustand";
import { persist } from "zustand/middleware";

export const LS_KEY = "himalaya-ui:state";

export type PageId = "index" | "compare" | "inspect";
export type ThemeId = "dark" | "light";
export type NavModalStep = "experiment" | "sample";

export interface AppState {
  // persisted
  username: string | undefined;
  firstName: string | undefined;
  lastName: string | undefined;
  activeExperimentId: number | undefined;
  activeSampleId: number | undefined;
  activeExposureId: number | undefined;
  activePage: PageId;
  tutorialSeen: boolean;
  theme: ThemeId;

  // ephemeral (not persisted)
  hoveredIndexId: number | undefined;
  hoveredPeakId: number | undefined;
  navModalOpen: boolean;
  navModalStep: NavModalStep;

  // setters
  setUsername: (name: string) => void;
  setUser: (u: { username: string; firstName?: string | undefined; lastName?: string | undefined }) => void;
  setActiveExperiment: (id: number | undefined) => void;
  setActiveSample: (id: number | undefined) => void;
  setActiveExposure: (id: number | undefined) => void;
  setHoveredIndex: (id: number | undefined) => void;
  setHoveredPeak: (id: number | undefined) => void;
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
      firstName: undefined,
      lastName: undefined,
      activeExperimentId: undefined,
      activeSampleId: undefined,
      activeExposureId: undefined,
      activePage: "index",
      tutorialSeen: false,
      theme: "dark",

      hoveredIndexId: undefined,
      hoveredPeakId: undefined,
      navModalOpen: false,
      navModalStep: "experiment",

      setUsername: (username) => set({ username }),
      setUser: ({ username, firstName, lastName }) =>
        set({ username, firstName, lastName }),
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
      setHoveredPeak: (hoveredPeakId) => set({ hoveredPeakId }),
      setActivePage: (activePage) => set({ activePage }),
      setTutorialSeen: (tutorialSeen) => set({ tutorialSeen }),
      setTheme: (theme) => set({ theme }),
      openNavModal: (step) =>
        set(step ? { navModalOpen: true, navModalStep: step } : { navModalOpen: true }),
      closeNavModal: () => set({ navModalOpen: false }),
      setNavModalStep: (navModalStep) => set({ navModalStep }),
      clearUsername: () => set({ username: undefined, firstName: undefined, lastName: undefined }),
    }),
    {
      name: LS_KEY,
      version: 3,
      partialize: (s) => ({
        username: s.username,
        firstName: s.firstName,
        lastName: s.lastName,
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
