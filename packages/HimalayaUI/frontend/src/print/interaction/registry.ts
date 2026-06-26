import { create } from "zustand";
import type { Action, CursorStepperProps, ListCursor } from "./types";

interface InteractionState {
  cursor: ListCursor | null;
  actions: Action[];
  extraSteppers: CursorStepperProps[];
  setPage: (cursor: ListCursor | null, actions: Action[], extraSteppers?: CursorStepperProps[]) => void;
  clearPage: () => void;
}

/** The single source the Dock and the keyboard layer both derive from. No
 *  middleware (mirrors floatingDock.ts) — purely in-memory, never persisted. */
export const useInteraction = create<InteractionState>((set) => ({
  cursor: null,
  actions: [],
  extraSteppers: [],
  setPage: (cursor, actions, extraSteppers = []) => set({ cursor, actions, extraSteppers }),
  clearPage: () => set({ cursor: null, actions: [], extraSteppers: [] }),
}));
