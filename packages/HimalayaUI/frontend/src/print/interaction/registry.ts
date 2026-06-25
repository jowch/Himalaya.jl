import { create } from "zustand";
import type { Action, ListCursor } from "./types";

interface InteractionState {
  cursor: ListCursor | null;
  actions: Action[];
  setPage: (cursor: ListCursor | null, actions: Action[]) => void;
  clearPage: () => void;
}

/** The single source the Dock and the keyboard layer both derive from. No
 *  middleware (mirrors floatingDock.ts) — purely in-memory, never persisted. */
export const useInteraction = create<InteractionState>((set) => ({
  cursor: null,
  actions: [],
  setPage: (cursor, actions) => set({ cursor, actions }),
  clearPage: () => set({ cursor: null, actions: [] }),
}));
