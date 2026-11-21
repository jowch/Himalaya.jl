import { useEffect } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { AppRoutes } from "./shell/AppRoutes";
import { OnboardingFlow } from "./shell/OnboardingFlow";
import { NavModal } from "./shell/NavModal";
import { KbdOverlay } from "./shell/KbdOverlay";
import { ToastContainer, LiveRegion } from "./ui";
import { InfrastructureBanner } from "./shell/InfrastructureBanner";
import { DetectorLutFilters } from "./detector";
import { handleRemoteEvent } from "../lib/queue/replayCoordinator";
import { invalidateIngestFrameCache } from "../lib/queue/applyRemoteToCache";
import { attachPersistence, rehydrate } from "../lib/queue/persistence";
import { resolveMutator } from "../lib/queue/mutatorRegistry";
import { exposeTestHelpers } from "../lib/queue/testHelpers";
import { showToast } from "../lib/toast";
import { openLiveStream } from "../lib/liveStream";
import { useConnection } from "../lib/connection";
import { createReconnectRefresh } from "../lib/queue/reconnectRefresh";
import { useAppState } from "../state";
import type { SseEvent } from "../lib/queue/types";

/**
 * PrintApp — root of "The Print", the sole production app
 * (index.html → src/print/main.tsx → here). Two concerns: the persistent
 * workspace shell and the onboarding overlay (rendered only when username is
 * unset), plus the app-wide SSE + mutation-queue wiring.
 */
export function PrintApp(): JSX.Element {
  const qc = useQueryClient();
  const mc = qc.getMutationCache();
  const setIngestProgress = useAppState((s) => s.setIngestProgress);
  const clearIngestProgress = useAppState((s) => s.clearIngestProgress);

  // Expose minimal test helpers on `window.__himalayaTest` in DEV only.
  // Production bundles tree-shake this out (Vite + DEV gate).
  useEffect(() => {
    exposeTestHelpers(qc, mc);
  }, [qc, mc]);

  useEffect(() => {
    const onCuration = (e: MessageEvent) => {
      try {
        const parsed = JSON.parse((e as MessageEvent).data as string) as SseEvent;
        // Ingest progress is broadcast-only (spec §9.3): never an own-op. Handle
        // it HERE and RETURN — do NOT feed it to handleRemoteEvent, whose queue
        // reconciliation would roll back + re-run every pending edit on each
        // progress tick. The store write lives here (Zustand is NOT imported in
        // applyRemoteToCache, which stays pure); the cache invalidation delegates
        // to the shared invalidateIngestFrameCache helper (defined + tested in
        // applyRemoteToCache.ts, Task 18) — no duplication.
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
              // Terminal (complete/failed): drop the in-flight entry; the
              // experiment's ingest_status (refetched below) is the resting truth.
              clearIngestProgress(expId);
            }
            // Delegate cache invalidation to the shared helper (defined in
            // applyRemoteToCache.ts) — the single source of truth for which
            // query keys the ingest frames affect.
            invalidateIngestFrameCache(qc, expId, parsed.kind === "ingest_complete" || parsed.kind === "ingest_failed", p?.stage);
          }
          return; // do NOT run the queue reconciler for a broadcast-only frame
        }
        handleRemoteEvent(parsed, qc, mc);
      } catch {
        // malformed frame, ignore
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
  }, [qc, mc, setIngestProgress, clearIngestProgress]); // stable selector results; EventSource lifetime unchanged

  useEffect(() => {
    void rehydrate(qc, resolveMutator).then(({ dropped }) => {
      if (dropped > 0) {
        showToast(
          `${dropped} edits from a previous session couldn't be restored`,
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
      {/* I5.1 (#182): NavModal was mounted in the retired AppShell. It must
          stay mounted app-wide so the `/` + ⌘K shortcuts and StaleUrlPage
          recovery can open it from any surface. */}
      <NavModal />
      <KbdOverlay />
      <ToastContainer />
      <LiveRegion />
      <InfrastructureBanner />
      {/* Shared GPU colormap filters for every DetectorImage (Focus/Loupe/thumbs).
          Mounted once here so `filter: url(#detector-lut-…)` resolves app-wide. */}
      <DetectorLutFilters />
    </>
  );
}
