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
  /** Override the test id (default "dock"). Lets a consumer keep a semantic
   *  hook like "grouping-footer" while still reusing the one dock primitive so
   *  the bar height stays identical across pages. */
  testId?: string;
}

export function Dock({ children, className = "", testId = "dock" }: DockProps): JSX.Element {
  const setCenterLaneOccupied = useFloatingDock((s) => s.setCenterLaneOccupied);

  // Claim the dock lane on mount so InfrastructureBanner steps aside.
  // Release it on unmount so the banner re-centres when the page exits.
  useEffect(() => {
    setCenterLaneOccupied(true);
    return () => setCenterLaneOccupied(false);
  }, [setCenterLaneOccupied]);

  return (
    <div
      data-testid={testId}
      // The dock is a command bar, not a focus destination. Preventing the
      // mousedown default keeps keyboard focus on the page's interaction scope
      // when a control is clicked: without it, clicking a dock button moves focus
      // ONTO the button (outside [data-interaction-scope]), and the scope-gated
      // keyboard layer then ignores every subsequent bare key until the user
      // clicks back into the list. The click still fires; focus simply never moves.
      onMouseDown={(e) => e.preventDefault()}
      className={`fixed bottom-0 left-0 right-0 z-40 flex items-center gap-2 border-t border-hair-strong bg-plate px-4 py-2 shadow-[0_-7px_24px_-10px_rgba(40,30,20,0.22)] ${className}`}
    >
      {children}
    </div>
  );
}
