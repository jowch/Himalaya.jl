import { useEffect, useRef } from "react";
import { useAppState } from "../state";
import { ChatCard } from "../components/ChatCard";
import { PlotCard } from "../components/PlotCard";
import { IndicesCard } from "../components/IndicesCard";

/**
 * IndexPage — the three-card workspace.
 *
 * Vertical layout:
 *   [ ~28px breathing room ]
 *   Title row (centered over the plot column via the same 22/56/22 grid)
 *   [ small gap ]
 *   Three-card grid  (chat | plot | indices)
 *
 * On narrow viewports (<1100px) the chat card reflows below the others.
 *
 * On mount, if there's no experiment/sample selected (and a username is set),
 * we auto-open the nav modal at the relevant step. We do this once per mount
 * to avoid an annoying re-open loop if the user dismisses with Esc.
 */
export function IndexPage(): JSX.Element {
  const username     = useAppState((s) => s.username);
  const experimentId = useAppState((s) => s.activeExperimentId);
  const sampleId     = useAppState((s) => s.activeSampleId);
  const openModal    = useAppState((s) => s.openNavModal);

  const autoOpenedRef = useRef(false);
  useEffect(() => {
    if (autoOpenedRef.current) return;
    if (username === undefined) return; // wait for onboarding
    if (experimentId === undefined) {
      autoOpenedRef.current = true;
      openModal("experiment");
    } else if (sampleId === undefined) {
      autoOpenedRef.current = true;
      openModal("sample");
    }
  }, [username, experimentId, sampleId, openModal]);

  return (
    <div
      data-testid="index-page"
      className="flex-1 min-h-0 flex flex-col gap-4 px-4 pb-4 pt-2"
    >
      <div
        data-testid="three-card-grid"
        className="min-h-0 grid gap-3
                   grid-cols-1
                   min-[1100px]:grid-cols-[22fr_56fr_22fr]
                   h-auto min-[1100px]:flex-1
                   min-[1100px]:max-h-[min(700px,calc(100dvh-var(--chrome-h)-1.5rem))]"
      >
        <section className="card order-3 min-[1100px]:order-none min-h-[280px] min-[1100px]:min-h-0 overflow-hidden">
          <ChatCard />
        </section>
        <section className="card order-1 min-[1100px]:order-none min-h-[420px] min-[1100px]:min-h-0 overflow-hidden">
          <PlotCard />
        </section>
        <section className="card order-2 min-[1100px]:order-none min-h-[360px] min-[1100px]:min-h-0 overflow-hidden">
          <IndicesCard />
        </section>
      </div>
    </div>
  );
}
