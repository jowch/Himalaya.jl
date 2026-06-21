/**
 * Dock — the single fixed light bottom bar (plate background + top hairline +
 * soft upward shadow). Accepts placement CHILDREN only; ALL appearance lives
 * here in the ui/ primitive (closed-look / open-placement contract). Consumers
 * compose their own buttons from existing primitives and pass them as children.
 *
 * Lane registration: claims the floatingDock `centerLaneOccupied` slot on
 * mount and releases it on unmount, so the global InfrastructureBanner steps
 * aside to the bottom-right corner while the dock is visible (LA-COLLIDE).
 *
 * Radius note: uses `rounded-md` (5px) — the single design-system radius step.
 */

import { useEffect } from "react";
import { useFloatingDock } from "../shell/floatingDock";

export interface DockProps {
  children?: React.ReactNode;
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

export function Dock({ children, className = "" }: DockProps): JSX.Element {
  const setCenterLaneOccupied = useFloatingDock((s) => s.setCenterLaneOccupied);

  // Claim the dock lane on mount so InfrastructureBanner steps aside.
  // Release it on unmount so the banner re-centres when the page exits.
  useEffect(() => {
    setCenterLaneOccupied(true);
    return () => setCenterLaneOccupied(false);
  }, [setCenterLaneOccupied]);

  return (
    <div
      data-testid="dock"
      className={`fixed bottom-0 left-0 right-0 z-40 flex items-center gap-2 border-t border-hair bg-plate px-4 py-2 shadow-[0_-2px_8px_0_rgba(0,0,0,0.06)] ${className}`}
    >
      {children}
    </div>
  );
}
