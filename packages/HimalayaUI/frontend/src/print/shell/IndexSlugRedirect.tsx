import { useEffect, useState } from "react";
import { Navigate, useParams } from "react-router-dom";
import * as api from "../../api";
import { ResolvingFallback } from "./ResolvingFallback";

/**
 * IndexSlugRedirect — I4.4 (#181). Old `/index/:experiment/:sample` deep-links
 * carry slug names; the focus workspace route wants a numeric `:sampleId`.
 * Resolve the slugs via `/api/resolve`, then redirect to `/sample/:id`. On
 * 404 / error / a success body without a sample, land on the corpus contact
 * experiments home (`/experiments`). Index is retired; this only preserves old permalinks.
 */
export function IndexSlugRedirect(): JSX.Element {
  const { experiment, sample } = useParams();
  const [target, setTarget] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;
    (async () => {
      if (experiment === undefined || sample === undefined) {
        if (!cancelled) setTarget("/experiments");
        return;
      }
      try {
        const body = await api.resolve({ experiment, sample });
        if (cancelled) return;
        // Undefined-guard: ResolveSuccess.sample_id is `number | undefined`
        // (api.ts). A success body without a sample (e.g. an experiment-only
        // resolve) would otherwise build "/sample/undefined" — treat it like
        // an error and fall back to the corpus.
        setTarget(
          "error" in body || body.sample_id === undefined
            ? "/experiments"
            : `/sample/${body.sample_id}`,
        );
      } catch {
        if (!cancelled) setTarget("/experiments");
      }
    })();
    return () => { cancelled = true; };
  }, [experiment, sample]);

  if (target === null) return <ResolvingFallback />;
  return <Navigate to={target} replace />;
}
