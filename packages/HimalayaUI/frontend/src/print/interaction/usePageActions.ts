import { useEffect } from "react";
import type { ReactNode } from "react";
import type { Action, CursorStepperProps, ListCursor } from "./types";
import { useInteraction } from "./registry";
import { assertNoCoreCollision } from "./core";

/** A page declares its cursor + actions here. Written to the store in a
 *  dependency-less effect so it runs every commit — the closures captured in
 *  enabled()/run are therefore always the latest (no manual stateRef). The
 *  store clears on unmount so an un-migrated next page shows no stale dock. */
export function usePageActions(decl: {
  cursor?: ListCursor | null;
  actions: Action[];
  extraSteppers?: CursorStepperProps[];
  dockExtra?: ReactNode;
}): void {
  const setPage = useInteraction((s) => s.setPage);
  const clearPage = useInteraction((s) => s.clearPage);

  assertNoCoreCollision(decl.actions); // dev/test guard; cheap, pure

  // No dependency array: refresh the registry on every commit.
  useEffect(() => {
    setPage(decl.cursor ?? null, decl.actions, decl.extraSteppers ?? [], decl.dockExtra ?? null);
  });

  useEffect(() => clearPage, [clearPage]); // unmount only
}
