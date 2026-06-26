import { create } from "zustand";
import type { ReactNode } from "react";
import type { Action, CursorStepperProps, ListCursor } from "./types";

interface InteractionState {
  cursor: ListCursor | null;
  actions: Action[];
  extraSteppers: CursorStepperProps[];
  dockExtra: ReactNode | null;
  setPage: (cursor: ListCursor | null, actions: Action[], extraSteppers?: CursorStepperProps[], dockExtra?: ReactNode | null) => void;
  clearPage: () => void;
}

/** The single source the Dock and the keyboard layer both derive from. No
 *  middleware (mirrors floatingDock.ts) — purely in-memory, never persisted. */
export const useInteraction = create<InteractionState>((set) => ({
  cursor: null,
  actions: [],
  extraSteppers: [],
  dockExtra: null,
  setPage: (cursor, actions, extraSteppers = [], dockExtra = null) => set({ cursor, actions, extraSteppers, dockExtra }),
  clearPage: () => set({ cursor: null, actions: [], extraSteppers: [], dockExtra: null }),
}));
