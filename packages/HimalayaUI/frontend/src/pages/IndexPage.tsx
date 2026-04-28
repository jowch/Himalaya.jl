import { useEffect, useRef } from "react";
import { useAppState } from "../state";
import { ChatCard } from "../components/ChatCard";
import { PlotCard } from "../components/PlotCard";
import { IndicesCard } from "../components/IndicesCard";
import { WorkspaceGrid } from "../components/WorkspaceGrid";

/**
 * IndexPage — three-card workspace (chat | plot | indices).
 *
 * Layout is driven by `WorkspaceGrid` (shared with InspectPage):
 *   < 1400px  — single column stacked: plot → indices → chat
 *   ≥ 1400px  — three columns: minmax(320px,22fr) | 56fr | minmax(320px,22fr)
 *
 * On mount, if there's no experiment/sample selected (and a username is set),
 * we auto-open the nav modal at the relevant step. Once per mount, so dismissing
 * with Esc doesn't re-open.
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
      className="flex-1 min-h-0 flex flex-col gap-4 px-4 pb-6 pt-2"
    >
      <WorkspaceGrid
        left={<ChatCard />}
        center={<PlotCard />}
        right={<IndicesCard />}
        slotClassName={{
          left:   "min-h-[280px]",
          center: "min-h-[420px]",
          right:  "min-h-[360px]",
        }}
      />
    </div>
  );
}
