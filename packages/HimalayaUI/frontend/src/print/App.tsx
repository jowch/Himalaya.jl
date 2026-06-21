import { useEffect } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { AppRoutes } from "./shell/AppRoutes";
import { OnboardingFlow } from "./shell/OnboardingFlow";
import { NavModal } from "./shell/NavModal";
import { KbdOverlay } from "./shell/KbdOverlay";
import { ToastContainer, LiveRegion } from "./ui";
import { InfrastructureBanner } from "./shell/InfrastructureBanner";
import { handleRemoteEvent } from "../lib/queue/replayCoordinator";
import { invalidateIngestFrameCache } from "../lib/queue/applyRemoteToCache";
import { attachPersistence, rehydrate } from "../lib/queue/persistence";
import { resolveMutator } from "../lib/queue/mutatorRegistry";
import { exposeTestHelpers } from "../lib/queue/testHelpers";
import { showToast } from "../lib/toast";
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
    const es = new EventSource("/api/events");
    es.addEventListener("curation", (e) => {
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
          const p = (parsed as { payload?: { experiment_id?: number; processed?: number; total?: number; phase?: string } }).payload;
          const expId = p?.experiment_id;
          if (expId !== undefined) {
            if (parsed.kind === "ingest_started" || parsed.kind === "ingest_progress") {
              // The backend tags a rescan's frames phase:"rescan" (the /{id}/scan
              // route) → "analyzing" → the inline ProgressBar, since the
              // experiment's table data is already present. An initial scan (the
              // create route) sends no phase → "scanning" → GroupingReviewPage.
              setIngestProgress(expId, {
                processed: p?.processed ?? 0,
                total: p?.total ?? 0,
                status: p?.phase === "rescan" ? "analyzing" : "scanning",
              });
            } else {
              // Terminal (complete/failed): drop the in-flight entry; the
              // experiment's ingest_status (refetched below) is the resting truth.
              clearIngestProgress(expId);
            }
            // Delegate cache invalidation to the shared helper (defined in
            // applyRemoteToCache.ts) — the single source of truth for which
            // query keys the ingest frames affect.
            invalidateIngestFrameCache(qc, expId, parsed.kind === "ingest_complete");
          }
          return; // do NOT run the queue reconciler for a broadcast-only frame
        }
        handleRemoteEvent(parsed, qc, mc);
      } catch {
        // malformed frame, ignore
      }
    });
    return () => es.close();
  }, [qc, mc, setIngestProgress, clearIngestProgress]); // stable selector results; EventSource lifetime unchanged

  // Mirror pending mutation queue to sessionStorage so a tab reload can
  // rehydrate it. attach is safe to mount unconditionally; subscription is
  // cheap.
  useEffect(() => {
    return attachPersistence(mc);
  }, [mc]);

  // On mount, replay any persisted ops left over from a previous tab
  // session through their matching mutators. Server-side request-level
  // idempotency (X-Client-Op-Id) makes this safe even if the original op
  // already landed. Surfacing dropped/failed counts as toasts so the user
  // knows when edits couldn't be restored.
  useEffect(() => {
    void rehydrate(qc, resolveMutator).then(({ dropped, failed }) => {
      if (dropped > 0) {
        showToast(
          `${dropped} edits from a previous session couldn't be restored`,
          "warning",
        );
      }
      if (failed > 0) {
        showToast(
          `${failed} edits failed to replay; please retry`,
          "error",
        );
      }
    });
    // Intentionally empty deps — rehydrate runs once at mount; qc is stable.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

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
    </>
  );
}
