import { useEffect, useRef } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import { useQueryClient } from "@tanstack/react-query";
import { AppRoutes } from "./shell/AppRoutes";
import { OnboardingFlow } from "./shell/OnboardingFlow";
import { ToastContainer, LiveRegion } from "./ui";
import { InfrastructureBanner } from "./shell/InfrastructureBanner";
import { DetectorLutFilters } from "./detector";
import { handleRemoteEvent } from "./lib/queue/replayCoordinator";
import { invalidateIngestFrameCache } from "./lib/queue/applyRemoteToCache";
import { attachPersistence, rehydrate } from "./lib/queue/persistence";
import { resolveMutator } from "./lib/queue/mutatorRegistry";
import { exposeTestHelpers } from "./lib/queue/testHelpers";
import { showToast } from "./lib/toast";
import { openLiveStream } from "./lib/liveStream";
import { useConnection } from "./lib/connection";
import { createReconnectRefresh } from "./lib/queue/reconnectRefresh";
import { useAppState } from "./state";
import { forgetExperiment } from "./queries";
import { getClientId } from "./lib/clientId";
import type { SseEvent } from "./lib/queue/types";
import { useDocumentTheme } from "./lib/theme";

export function App(): JSX.Element {
  useDocumentTheme();
  const qc = useQueryClient();
  const mc = qc.getMutationCache();
  const setIngestProgress = useAppState((s) => s.setIngestProgress);
  const clearIngestProgress = useAppState((s) => s.clearIngestProgress);
  const navigate = useNavigate();
  const location = useLocation();
  const here = useRef({ navigate, pathname: location.pathname });
  here.current = { navigate, pathname: location.pathname };

  useEffect(() => {
    exposeTestHelpers(qc, mc);
  }, [qc, mc]);

  useEffect(() => {
    const onCuration = (e: MessageEvent) => {
      try {
        const parsed = JSON.parse((e as MessageEvent).data as string) as SseEvent;
        if (parsed.kind === "experiment_deleted") {
          const expId = parsed.entity_id;
          const { navigate: go, pathname } = here.current;
          const open = pathname === `/experiments/${expId}` || pathname.startsWith(`/experiments/${expId}/`);
          if (open && parsed.client_id !== getClientId()) {
            go("/experiments", { replace: true });
            showToast("This experiment was removed in another window.");
          }
          forgetExperiment(qc, expId);
          return;
        }
        if (
          parsed.kind === "ingest_started" || parsed.kind === "ingest_progress" ||
          parsed.kind === "ingest_complete" || parsed.kind === "ingest_failed"
        ) {
          const p = (parsed as { payload?: { experiment_id?: number; processed?: number; total?: number; phase?: string; stage?: string } }).payload;
          const expId = p?.experiment_id;
          if (expId !== undefined) {
            if (parsed.kind === "ingest_started" || parsed.kind === "ingest_progress") {
              setIngestProgress(expId, {
                processed: p?.processed ?? 0,
                total: p?.total ?? 0,
                status: p?.phase === "rescan" ? "analyzing" : "scanning",
                stage: p?.stage,
              });
            } else {
              clearIngestProgress(expId);
            }
            invalidateIngestFrameCache(qc, expId, parsed.kind === "ingest_complete" || parsed.kind === "ingest_failed", p?.stage);
          }
          return;
        }
        handleRemoteEvent(parsed, qc, mc);
      } catch {
      }
    };
    return openLiveStream("/api/events", {
      listeners: { curation: onCuration },
      onReconnect: createReconnectRefresh(qc, mc),
      onStatus: (status) => {
        const c = useConnection.getState();
        if (status === "open") c.setStreamOpen();
        else c.setStreamDown();
      },
    });
  }, [qc, mc, setIngestProgress, clearIngestProgress]);

  useEffect(() => {
    void rehydrate(qc, resolveMutator).then(({ dropped, lost }) => {
      if (dropped + lost > 0) {
        showToast(
          `${dropped + lost} edits from a previous session couldn't be restored`,
          "warning",
        );
      }
    });
    return attachPersistence(mc);
  }, [qc, mc]);

  return (
    <>
      <AppRoutes />
      <OnboardingFlow />
      <ToastContainer />
      <LiveRegion />
      <InfrastructureBanner />
      <DetectorLutFilters />
    </>
  );
}
