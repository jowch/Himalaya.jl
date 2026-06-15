import { create } from "zustand";

/**
 * floatingDock — coordinates the bottom-centre "dock" lane so the global
 * InfrastructureBanner never collides with a page's floating action bars
 * (the contact sheet's CullBar / ComposeBar). LA-COLLIDE.
 *
 * The banner is fixed bottom-centre by default (most prominent for a "Saving…"
 * / "Couldn't save" status). But a page that mounts a bottom-centre action bar
 * overlaps that lane with an OPAQUE, higher-z bar (`bg-ink z-50`) that paints
 * straight over the banner (`z-40`), occluding it entirely. A page declares the
 * lane occupied while any such bar is shown; the banner reads this and steps
 * aside to the bottom-right corner (free — toasts own top-right) until the lane
 * clears, so both stay visible.
 *
 * Deliberately a small dedicated store rather than a field on the big
 * `AppState`: single responsibility, and it keeps the publisher (a page) and
 * the reader (the app shell) decoupled from app/draft/url state. The signal is
 * a plain boolean, NOT a pixel offset — no brittle coupling to bar heights.
 */
interface FloatingDockState {
  /** True while a bottom-centre action bar occupies the dock lane. */
  centerLaneOccupied: boolean;
  setCenterLaneOccupied: (occupied: boolean) => void;
}

export const useFloatingDock = create<FloatingDockState>((set) => ({
  centerLaneOccupied: false,
  setCenterLaneOccupied: (centerLaneOccupied) => set({ centerLaneOccupied }),
}));
