import { useEffect, useState } from "react";
import { useMutationState } from "@tanstack/react-query";

const SHOW_THRESHOLD_MS = 500;
const STUCK_THRESHOLD_MS = 30000;

/**
 * InfrastructureBanner — persistent status strip for long-pending mutations.
 *
 * Hidden when no mutation has been pending for >500ms (avoids flicker on
 * fast successes). Upgrades to a "stuck" state with a Refresh button after
 * any pending mutation has been in-flight >30s.
 */
export function InfrastructureBanner(): JSX.Element | null {
  const [, setTick] = useState(0);

  const pendingSubmittedAts = useMutationState({
    filters: { status: "pending" },
    select: (m) => m.state.submittedAt,
  });

  // Re-tick once a second so time-since-submit comparisons update without
  // depending on mutation cache events. Gated on whether any mutation is
  // pending — the steady-state idle case (no pending writes) doesn't need
  // a wakeup every second. Effect only re-runs when crossing the empty
  // boundary, not on every count change.
  const hasPending = pendingSubmittedAts.length > 0;
  useEffect(() => {
    if (!hasPending) return;
    const handle = window.setInterval(() => setTick((t) => t + 1), 1000);
    return () => window.clearInterval(handle);
  }, [hasPending]);

  const now = Date.now();
  const visible = pendingSubmittedAts.filter(
    (t) => typeof t === "number" && now - t > SHOW_THRESHOLD_MS,
  );

  if (visible.length === 0) return null;

  const stuck = visible.some((t) => now - t > STUCK_THRESHOLD_MS);

  const stateAttr = stuck ? "stuck" : "showing";
  // Severity is conveyed by the leading icon + word below, not an edge hue.

  return (
    <div
      data-testid="infrastructure-banner"
      data-state={stateAttr}
      role="status"
      className={
        "fixed left-1/2 -translate-x-1/2 bottom-4 z-40 flex items-center gap-3 " +
        "rounded-md border border-hair bg-plate text-ink px-4 py-2 shadow-lg text-body"
      }
    >
      {stuck ? (
        <>
          <span aria-label="Error" className="flex-shrink-0 font-bold text-error">
            {"✕"}
          </span>
          <span>
            <span className="font-semibold">Error.</span> Couldn&rsquo;t save. Try
            refreshing.
          </span>
          <button
            type="button"
            onClick={() => window.location.reload()}
            className="rounded-md border border-error px-2 py-0.5 text-ink hover:bg-paper-sunk"
          >
            Refresh
          </button>
        </>
      ) : (
        <>
          <span
            aria-label="Saving"
            className="inline-block h-3 w-3 flex-shrink-0 rounded-full border-2 border-accent border-t-transparent animate-spin"
          />
          <span>
            <span className="font-semibold">Saving</span> your changes&hellip;
          </span>
        </>
      )}
    </div>
  );
}
