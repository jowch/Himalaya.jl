import { useEffect } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import { useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
import { parseLocation } from "../lib/url/parseLocation";
import * as api from "../api";
import type { ResolveSuccess, Experiment, Sample } from "../api";
import { queryKeys } from "../queries";
import { emitReplaceNext } from "../lib/url/emitMode";

// Spec §4.2 — URL → Zustand. Reads `useLocation()` so popstate AND
// useNavigate both flow through the hook. Origin-tagged fetches avoid
// the Zustand-mid-flight race (Zustand mutations don't change `location`,
// so AbortController alone is insufficient).

/**
 * Atomic apply of a 200 resolve response. Single setState commit so
 * useUrlFromState recomputes once, no cascading partial URL emits.
 *
 * App-init URL sync (per spec §4.3 trigger table) → replace, so the
 * initial state→URL emit doesn't push a redundant history entry over
 * the URL the user just landed on.
 */
function applySuccess(body: ResolveSuccess, page: "index" | "inspect" | "compare") {
  emitReplaceNext();
  useAppState.setState({
    activePage: page,
    activeExperimentId: body.experiment_id,
    activeSampleId: body.sample_id,
    activeExposureId: body.exposure_id,
    staleUrlContext: null,
    resolving: false,
  });
}

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
      useAppState.setState({
        staleUrlContext: { kind: "unknown_path", raw: parsed.raw },
        resolving: false,
      });
      return;
    }

    if (parsed.kind === "compare") {
      // Compare uses numeric ids resolved by react-router useParams in the
      // ComparePage component itself; just set the active page.
      useAppState.setState({
        activePage: "compare",
        staleUrlContext: null,
        resolving: false,
      });
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
      useAppState.setState({ resolving: true });
      (async () => {
        const q: api.ResolveQuery = { experiment_id: expId };
        if (sId !== undefined) q.sample_id = sId;
        let body;
        try {
          body = await api.resolve(q);
        } catch {
          if (!cancelled) {
            useAppState.setState({ resolving: false });
            navigate("/index", { replace: true });
          }
          return;
        }
        if (cancelled) return;
        useAppState.setState({ resolving: false });
        if ("error" in body) {
          navigate("/index", { replace: true });
          return;
        }
        const path = body.sample_name !== undefined
          ? `/${page}/${encodeURIComponent(body.experiment_name)}/${encodeURIComponent(body.sample_name)}`
          : `/${page}/${encodeURIComponent(body.experiment_name)}`;
        navigate(path, { replace: true });
      })();
      return;
    }

    // index or inspect — fetch /api/resolve with whichever slugs are present.
    if (parsed.experiment === undefined) {
      useAppState.setState({
        activePage: parsed.kind,
        activeExperimentId: undefined,
        activeSampleId: undefined,
        activeExposureId: undefined,
        staleUrlContext: null,
        resolving: false,
      });
      return;
    }

    // Pre-fetch clear of staleUrlContext + set resolving.
    useAppState.setState({ staleUrlContext: null, resolving: true });

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
        if (!cancelled) useAppState.setState({ resolving: false });
        return;
      }
      if (cancelled) return;
      const currentWindowUrl = window.location.pathname + window.location.search;
      if (currentWindowUrl !== startWindowUrl) return;
      if ("error" in body && body.error === "not_found") {
        useAppState.setState({
          staleUrlContext: {
            kind: "not_found",
            missing: body.missing,
            missing_value: body.missing_value,
            experiment_resolved: body.experiment_resolved,
            sample_resolved: body.sample_resolved,
          },
          resolving: false,
        });
        return;
      }
      if ("error" in body) {
        useAppState.setState({
          staleUrlContext: { kind: "unknown_path", raw: origin },
          resolving: false,
        });
        return;
      }
      applySuccess(body as ResolveSuccess, parsed.kind);
    })();

    return () => {
      cancelled = true;
      ctl.abort();
    };
  }, [location.pathname, location.search, navigate, qc]);
}
