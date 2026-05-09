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
      (async () => {
        const q: api.ResolveQuery = { experiment_id: expId };
        if (sId !== undefined) q.sample_id = sId;
        let body;
        try {
          body = await api.resolve(q);
        } catch {
          if (!cancelled) navigate("/index", { replace: true });
          return;
        }
        if (cancelled) return;
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
      // Origin-tag check: did the URL change during the fetch?
      if (cancelled || (window.location.pathname + window.location.search) !== origin) {
        return;
      }
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
