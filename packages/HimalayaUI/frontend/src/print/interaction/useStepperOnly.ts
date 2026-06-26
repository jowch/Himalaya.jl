import { useMemo } from "react";
import type { CursorStepperProps } from "./types";

interface Opts {
  ids: number[];
  currentId: number | undefined;
  onGo: (id: number) => void;
  label: string;
  testIdBase: string;
  axis?: "vertical" | "horizontal";
}

/** A cursor-less stepper for an axis whose state lives elsewhere (e.g. the URL).
 *  Same CursorStepperProps shape as useListCursor.stepperProps so the dock
 *  renders it identically — no roving focus, no selection, no in-page rows. */
export function useStepperOnly(opts: Opts): CursorStepperProps {
  const { ids, currentId, onGo, label, testIdBase, axis = "vertical" } = opts;
  return useMemo<CursorStepperProps>(() => {
    const i = currentId === undefined ? -1 : ids.indexOf(currentId);
    return {
      label, axis, testIdBase,
      count: ids.length > 0 ? `${Math.max(0, i) + 1} / ${ids.length}` : "0 / 0",
      onPrev: () => { if (i > 0) onGo(ids[i - 1]!); },
      onNext: () => { if (i >= 0 && i < ids.length - 1) onGo(ids[i + 1]!); },
      prevDisabled: i <= 0,
      nextDisabled: i < 0 || i >= ids.length - 1,
    };
  }, [ids, currentId, onGo, label, testIdBase, axis]);
}
