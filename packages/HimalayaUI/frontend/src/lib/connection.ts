import { create } from "zustand";

interface ConnectionState {
  streamDownSince: number | null;
  setStreamDown: () => void;
  setStreamOpen: () => void;
}

export const useConnection = create<ConnectionState>((set, get) => ({
  streamDownSince: null,
  setStreamDown: () => {
    if (get().streamDownSince === null) set({ streamDownSince: Date.now() });
  },
  setStreamOpen: () => set({ streamDownSince: null }),
}));
