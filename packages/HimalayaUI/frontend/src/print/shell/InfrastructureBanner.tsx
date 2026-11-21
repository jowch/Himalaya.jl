import { useEffect, useState } from "react";
import { useMutationState } from "@tanstack/react-query";
import { useConnection } from "../../lib/connection";

const SHOW_THRESHOLD_MS = 500;
const STUCK_THRESHOLD_MS = 30000;
const OFFLINE_THRESHOLD_MS = 5000;

function changes(n: number): string {
  return n === 1 ? "1 change" : `${n} changes`;
}

/**
 * InfrastructureBanner — persistent status strip for long-pending mutations.
 *
 * Hidden when no mutation has been pending for >500ms (avoids flicker on
 * fast successes). After 30s pending it counts the changes not saved yet.
 * When the live stream has been down for 5s it says the connection is lost,
 * pending edits or not; queued edits keep retrying until the server answers.
 */
export function InfrastructureBanner(): JSX.Element | null {
  const [, setTick] = useState(0);

  const streamDownSince = useConnection((s) => s.streamDownSince);

  const pendingSubmittedAts = useMutationState({
    filters: { status: "pending", predicate: (m) => m.options.meta?.savesWork !== false },
    select: (m) => m.state.submittedAt,
  });

  const needsTick = pendingSubmittedAts.length > 0 || streamDownSince !== null;
  useEffect(() => {
    if (!needsTick) return;
    const handle = window.setInterval(() => setTick((t) => t + 1), 1000);
    return () => window.clearInterval(handle);
  }, [needsTick]);

  const now = Date.now();
  const visible = pendingSubmittedAts.filter(
    (t) => typeof t === "number" && now - t > SHOW_THRESHOLD_MS,
  );
  const offline =
    streamDownSince !== null && now - streamDownSince > OFFLINE_THRESHOLD_MS;

  if (visible.length === 0 && !offline) return null;

  const stuck = visible.some((t) => now - t > STUCK_THRESHOLD_MS);
  const stateAttr = offline ? "offline" : stuck ? "stuck" : "showing";
  const pending = pendingSubmittedAts.length;

  return (
    <div
      data-testid="infrastructure-banner"
      data-state={stateAttr}
      role="status"
      className={
        "fixed bottom-4 left-1/2 z-40 flex -translate-x-1/2 items-center gap-3 " +
        "rounded-md border border-hair bg-plate text-ink px-4 py-2 shadow-lg text-body"
      }
    >
      {offline || stuck ? (
        <>
          <span aria-label={offline ? "Offline" : "Not saved"} className="flex-shrink-0 font-bold text-warning">
            {"✕"}
          </span>
          {offline ? (
            <span>
              <span className="font-semibold">Connection lost.</span> Retrying&hellip;
              {pending > 0
                ? ` ${changes(pending)} will save when it\u2019s back.`
                : null}
            </span>
          ) : (
            <span>
              <span className="font-semibold">{changes(pending)} not saved yet.</span>{" "}
              Still trying.
            </span>
          )}
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
