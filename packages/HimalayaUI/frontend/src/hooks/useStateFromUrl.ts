import { useEffect } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import { useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
import { parseLocation } from "../lib/url/parseLocation";
import * as api from "../api";
import type { ResolveSuccess, Experiment, Sample } from "../api";
import { queryKeys } from "../queries";

// Spec §4.2 — URL → Zustand. Reads `useLocation()` so popstate AND
// useNavigate both flow through the hook. Origin-tagged fetches avoid
// the Zustand-mid-flight race (Zustand mutations don't change `location`,
// so AbortController alone is insufficient).

const VALID_PAGES = new Set(["index", "inspect", "compare"]);

export function useStateFromUrl(): void {
  const location = useLocation();
  const navigate = useNavigate();
  const qc = useQueryClient();

  useEffect(() => {
    let cancelled = false;
    const origin = location.pathname + location.search;
    const parsed = parseLocation(location.pathname, location.search);

    if (parsed.kind === "stale") {
      useAppState.getState().setStaleUnknownPath(parsed.raw);
      return;
    }

    if (parsed.kind === "compare") {
      // Compare uses numeric ids resolved by react-router useParams in the
      // ComparePage component itself; just set the active page + clear
      // resolving. `setActivePage` already clears `staleUrlContext` as a
      // side effect, so we don't need a separate call for that.
      const a = useAppState.getState();
      a.setActivePage("compare");
      a.setResolving(false);
      return;
    }

    if (parsed.kind === "root") {
      // §5 redirect: build slug URL from persisted Zustand and navigate
      // (replace). All emits go through useNavigate so react-router's
      // location subscription stays consistent.
      const s = useAppState.getState();
      const page = VALID_PAGES.has(s.activePage) ? s.activePage : "index";
      const expId = s.activeExperimentId;
      const sId = s.activeSampleId;
      if (page === "compare") {
        navigate(expId !== undefined
          ? `/experiments/${expId}/compare`
          : "/compare/all", { replace: true });
        return;
      }
      if (expId === undefined) {
        navigate(`/${page}`, { replace: true });
        return;
      }
      // Synchronous cache hit?
      const exps = qc.getQueryData<Experiment[]>(queryKeys.experiments) ?? [];
      const expName = exps.find((e) => e.id === expId)?.name;
      const samples = sId !== undefined
        ? qc.getQueryData<Sample[]>(queryKeys.samples(expId)) ?? []
        : [];
      const sName = sId !== undefined ? samples.find((x) => x.id === sId)?.name : null;
      if (expName !== null && expName !== undefined &&
          (sId === undefined || (sName !== null && sName !== undefined))) {
        const path = sName !== undefined && sName !== null
          ? `/${page}/${encodeURIComponent(expName)}/${encodeURIComponent(sName)}`
          : `/${page}/${encodeURIComponent(expName)}`;
        navigate(path, { replace: true });
        return;
      }
      // Cold-mount fallback: resolve-by-id.
      // Set `resolving: true` so useUrlFromState's guard suppresses its
      // settle emit while the async fetch is in flight. Without this the
      // Zustand→URL hook sees an empty experiments cache (queries haven't
      // resolved yet either), computes target=`/index`, and navigates away
      // from `/` BEFORE the resolve-by-id call lands — clobbering the
      // seeded slug pair we're trying to recover. Symmetric with the
      // resolving:true set on the recognized-kind branch below.
      useAppState.getState().setResolving(true);
      (async () => {
        const q: api.ResolveQuery = { experiment_id: expId };
        if (sId !== undefined) q.sample_id = sId;
        let body;
        try {
          body = await api.resolve(q);
        } catch {
          if (!cancelled) {
            useAppState.getState().setResolving(false);
            navigate("/index", { replace: true });
          }
          return;
        }
        if (cancelled) return;
        useAppState.getState().setResolving(false);
        if ("error" in body) {
          navigate("/index", { replace: true });
          return;
        }
        const path = body.sample_name !== undefined
          ? `/${page}/${encodeURIComponent(body.experiment_name)}/${encodeURIComponent(body.sample_name)}`
          : `/${page}/${encodeURIComponent(body.experiment_name)}`;
        navigate(path, { replace: true });
      })();
      // Cleanup also runs for this branch — see the cleanup at the end of
      // the effect for the rationale on why we always reset resolving.
      return () => {
        cancelled = true;
        useAppState.getState().setResolving(false);
      };
    }

    // index or inspect — fetch /api/resolve with whichever slugs are present.
    if (parsed.experiment === undefined) {
      // Bare /index or /inspect — clear active selection atomically.
      useAppState.getState().setResolveSuccess({
        page: parsed.kind,
        experimentId: undefined,
        sampleId: undefined,
        exposureId: undefined,
      });
      return;
    }

    // Pre-fetch clear of staleUrlContext + set resolving.
    {
      const a = useAppState.getState();
      a.setStaleUrlContext(null);
      a.setResolving(true);
    }

    // Capture window.location at fetch start. We re-read at fetch end
    // and bail if it changed, catching the raw `history.replaceState` /
    // `pushState` case where the URL was changed without going through
    // react-router (no popstate, no cleanup). For react-router-driven
    // changes (TabRocker, NavModal, useNavigate), `cancelled` is set by
    // the effect cleanup. The combination handles both production
    // (BrowserRouter) and tests (MemoryRouter — where window.location is
    // stable, so this check no-ops and we rely on `cancelled`).
    const startWindowUrl = window.location.pathname + window.location.search;

    const ctl = new AbortController();
    const q: api.ResolveQuery = { experiment: parsed.experiment };
    if (parsed.sample !== undefined) q.sample = parsed.sample;
    if (parsed.kind === "inspect" && parsed.exposure !== undefined) {
      q.exposure = parsed.exposure;
    }

    (async () => {
      let body;
      try {
        body = await api.resolve(q, ctl.signal);
      } catch (e) {
        if ((e as Error).name === "AbortError") return;
        if (!cancelled) useAppState.getState().setResolving(false);
        return;
      }
      if (cancelled) return;
      const currentWindowUrl = window.location.pathname + window.location.search;
      if (currentWindowUrl !== startWindowUrl) return;
      if ("error" in body && body.error === "not_found") {
        useAppState.getState().setStaleNotFound({
          kind: "not_found",
          missing: body.missing,
          missing_value: body.missing_value,
          experiment_resolved: body.experiment_resolved,
          sample_resolved: body.sample_resolved,
        });
        return;
      }
      if ("error" in body) {
        useAppState.getState().setStaleUnknownPath(origin);
        return;
      }
      const ok = body as ResolveSuccess;
      useAppState.getState().setResolveSuccess({
        page: parsed.kind,
        experimentId: ok.experiment_id,
        sampleId: ok.sample_id,
        exposureId: ok.exposure_id,
      });
    })();

    return () => {
      cancelled = true;
      ctl.abort();
      // Defensive reset: if the effect is unmounting (no remount immediately
      // following), `resolving:true` would otherwise leak past the
      // component's lifetime. If a remount IS following (location changed
      // mid-fetch), the new effect's `setResolving(true)` runs synchronously
      // after this cleanup within React's deferred-effects flush, so no
      // intermediate render observes the false→true transition.
      useAppState.getState().setResolving(false);
    };
  }, [location.pathname, location.search, navigate, qc]);
}
