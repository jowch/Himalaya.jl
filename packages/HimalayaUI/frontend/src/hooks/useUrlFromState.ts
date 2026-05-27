import { useEffect, useRef } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import { useQueryClient, useQuery } from "@tanstack/react-query";
import { useAppState, type PageId } from "../state";
import { useExperiments, queryKeys } from "../queries";
import * as api from "../api";
import type { Experiment, Sample, Exposure } from "../api";
import { consumeEmitMode } from "../lib/url/emitMode";

// Spec §4.3 — Zustand → URL via react-router useNavigate. Subscribes to
// experiments + samples queries via TanStack so SSE-driven cache rewrites
// trigger a re-render of this hook (§7 invalidation).

function nameForExperiment(list: Experiment[] | undefined, id: number | undefined):
  string | undefined
{
  if (id === undefined || list === undefined) return undefined;
  const found = list.find((e) => e.id === id);
  if (found === undefined) return undefined;
  return found.name === null ? undefined : found.name;
}

function nameForSample(list: Sample[] | undefined, sId: number | undefined):
  string | undefined
{
  if (sId === undefined || list === undefined) return undefined;
  const found = list.find((s) => s.id === sId);
  if (found === undefined) return undefined;
  return found.name === null ? undefined : found.name;
}

function filenameForExposure(list: Exposure[] | undefined,
                              eId: number | undefined): string | undefined
{
  if (eId === undefined || list === undefined) return undefined;
  const found = list.find((e) => e.id === eId);
  if (found === undefined) return undefined;
  return found.filename === null ? undefined : found.filename;
}

function buildUrl(
  page: PageId,
  _experiment: string | undefined,
  _sample: string | undefined,
  _exposure: string | undefined,
  current: string,
): string {
  // Every legacy AppShell surface is retired — Inspect (#163), Index (#181),
  // Compare (#177) — so `activePage` is the inert `"none"` sentinel and there
  // is no Zustand-derived URL to emit. Returning `current` keeps the equality
  // guard happy and prevents an accidental redirect, so this hook never emits.
  // (The hook + the whole Zustand→URL bridge retire with the dual-nav model in
  // I5.1.)
  void page; void _experiment; void _sample; void _exposure;
  return current;
}

export function useUrlFromState(): void {
  const navigate = useNavigate();
  const location = useLocation();
  const qc = useQueryClient();

  const activePage = useAppState((s) => s.activePage);
  const activeExperimentId = useAppState((s) => s.activeExperimentId);
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const resolving = useAppState((s) => s.resolving);
  const staleUrlContext = useAppState((s) => s.staleUrlContext);

  // Subscribe to experiments + (when an experiment is active) samples
  // queries so SSE-driven cache rewrites trigger a re-render. We use
  // `useQuery` directly here (rather than `useSamples`) so we can gate
  // on `enabled: activeExperimentId !== undefined` — `useSamples` doesn't
  // expose `enabled`, so it would fire a GET on every cold mount.
  const { data: experiments } = useExperiments();
  const samplesQuery = useQuery({
    queryKey: queryKeys.samples(activeExperimentId),
    queryFn: () => api.listSamples(activeExperimentId as number),
    enabled: activeExperimentId !== undefined,
  });
  const samples = activeExperimentId !== undefined ? samplesQuery.data : undefined;
  // I1.7 (#163): Inspect — the only surface whose URL carried `?exposure=` —
  // is retired. No live surface emits the exposure into the URL, so the
  // exposure cache is never needed here and the query stays disabled.
  const exposuresEnabled = false;
  const exposuresQuery = useQuery({
    queryKey: queryKeys.exposures(activeSampleId),
    queryFn: () => api.listExposures(activeSampleId as number),
    enabled: exposuresEnabled,
  });
  const exposures = exposuresEnabled ? exposuresQuery.data : undefined;

  // Track the previous resolved slug pair so we can detect SSE-driven
  // disappearance (a slug went from defined → undefined because the
  // entity was deleted from the cache). Per spec §4.3 + §7, that case
  // should emit replace, not push (otherwise back-button stops at the
  // broken URL).
  const prevSlugsRef = useRef<{ exp: string | undefined; sample: string | undefined }>(
    { exp: undefined, sample: undefined },
  );

  useEffect(() => {
    // While useStateFromUrl is resolving the address bar (cold-mount deep
    // URL, popstate to a new slug), the URL is the authoritative input.
    // Tugging it back here from a half-populated Zustand snapshot would
    // race the resolve and clobber the user's deep link before
    // useStateFromUrl can apply the resolved id pair. The `resolving` flag
    // is precisely for this — released when applySuccess / staleUrlContext
    // dispatch lands. We re-read it via getState() because useStateFromUrl
    // (declared first in App.tsx, so its effect runs first) sets resolving
    // synchronously inside its effect; that value won't propagate into this
    // effect's closure until the next render commit, but Zustand's
    // getState() sees it immediately.
    if (resolving || useAppState.getState().resolving) return;
    // Same logic for stale URLs: the user is parked on an unresolvable
    // address (404 from /api/resolve, or unknown_path). Letting the URL be
    // overwritten by a half-populated state would clear the stale page
    // before they can act on it. StaleUrlPage / NavModal are responsible
    // for clearing staleUrlContext once recovery commits.
    if (staleUrlContext !== null) return;

    const expName = nameForExperiment(experiments, activeExperimentId);
    const sampleName = nameForSample(samples, activeSampleId);
    const exposureName = filenameForExposure(exposures, activeExposureId);
    // Compare against `window.location` (the DOM truth), not `useLocation()`.
    // When Zustand's useSyncExternalStore-driven re-render fires before
    // BrowserRouter's React.useState commit (e.g., PlotCard auto-pick
    // landing in the same commit cycle as a `,`/`.` sample-step keypress),
    // `useLocation()` lags by one render — the equality guard would miss
    // and emit a duplicate navigate to the URL we just wrote (issue #118).
    // Effect deps still include location.pathname/search so we re-run on
    // popstate and other URL changes; we just compare against the truth.
    // Same precedent: useStateFromUrl's mid-fetch race-detector reads
    // window.location directly.
    const current = window.location.pathname + window.location.search;
    const target = buildUrl(activePage, expName, sampleName, exposureName, current);

    // SSE-driven invalidation detection: a previously-resolvable slug
    // is now undefined (the row vanished from cache). Force replace.
    const prev = prevSlugsRef.current;
    const slugDisappeared =
      (prev.exp !== undefined && expName === undefined && activeExperimentId !== undefined) ||
      (prev.sample !== undefined && sampleName === undefined && activeSampleId !== undefined);
    prevSlugsRef.current = { exp: expName, sample: sampleName };

    // Wait for cache hydration before emitting. Without this, a cold-mount
    // deep link can be destroyed if /api/resolve returns before
    // /api/experiments: applySuccess populates Zustand active ids; this
    // hook fires with experiments still undefined; nameForExperiment
    // returns undefined; buildUrl emits /index; useStateFromUrl re-parses
    // /index and wipes the just-populated active ids.
    //
    // (Exposures used to gate this on Inspect; Inspect is retired (#163) and
    // no surface emits the exposure into the URL, so there's nothing to wait
    // for on that axis.)
    const cacheReady =
      (activeExperimentId === undefined || experiments !== undefined) &&
      (activeSampleId === undefined || samples !== undefined);
    if (!cacheReady) return;

    // Consume the emit-mode flag BEFORE the equality guard. Otherwise a
    // settle-emit (e.g. applySuccess flips to /index/lipid/JC001 — same as
    // the resolved URL) skips navigate but leaves `nextEmitMode = replace`
    // armed; the next user-driven URL change (tab click → push) silently
    // replaces instead, breaking back-button history.
    const explicitMode = consumeEmitMode();
    if (target === current) return;        // equality guard
    const mode = slugDisappeared ? "replace" : explicitMode;
    navigate(target, { replace: mode === "replace" });
  }, [
    activePage, activeExperimentId, activeSampleId, activeExposureId,
    resolving, staleUrlContext,
    experiments, samples, exposures,
    location.pathname, location.search,
    navigate, qc,
  ]);
}
