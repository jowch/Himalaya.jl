import { create } from "zustand";
import type { ReactNode } from "react";
import type { Action, CursorStepperProps, ListCursor } from "./types";

interface InteractionState {
  cursor: ListCursor | null;
  actions: Action[];
  shellActions: Action[];
  extraSteppers: CursorStepperProps[];
  dockExtra: ReactNode | null;
  /** Page-declared arrow-key handler. Invoked by useKeyboardLayer for Arrow*
   *  keys SCOPE-EXEMPT (arrows aren't WCAG-2.1.4 character shortcuts), so the
   *  active surface's cursor responds no matter where focus sits — never the
   *  page scroll. The page preventDefaults the arrows it claims. */
  arrowHandler: ((e: KeyboardEvent) => void) | null;
  setPage: (cursor: ListCursor | null, actions: Action[], extraSteppers?: CursorStepperProps[], dockExtra?: ReactNode | null, arrowHandler?: ((e: KeyboardEvent) => void) | null) => void;
  clearPage: () => void;
  setShellActions: (a: Action[]) => void;
}

/** The single source the Dock and the keyboard layer both derive from. No
 *  middleware (mirrors floatingDock.ts) — purely in-memory, never persisted. */
export const useInteraction = create<InteractionState>((set) => ({
  cursor: null,
  actions: [],
  shellActions: [],
  extraSteppers: [],
  dockExtra: null,
  arrowHandler: null,
  setPage: (cursor, actions, extraSteppers = [], dockExtra = null, arrowHandler = null) => set({ cursor, actions, extraSteppers, dockExtra, arrowHandler }),
  clearPage: () => set({ cursor: null, actions: [], extraSteppers: [], dockExtra: null, arrowHandler: null }),
  setShellActions: (a) => set({ shellActions: a }),
}));
