import { useEffect, useState } from "react";
import { Navigate, useParams } from "react-router-dom";
import * as api from "../../api";
import { ResolvingFallback } from "./ResolvingFallback";
import { LINK_FAILED, STALE_LINK, workspaceHref } from "./SampleRedirect";
import { showToast } from "../../lib/toast";

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
        if ("error" in body || body.sample_id === undefined) {
          showToast(STALE_LINK, "warning");
          setTarget("/experiments");
        } else {
          setTarget(workspaceHref(body.experiment_id, body.sample_id, body.exposure_id));
        }
      } catch {
        if (!cancelled) {
          showToast(LINK_FAILED, "error");
          setTarget("/experiments");
        }
      }
    })();
    return () => { cancelled = true; };
  }, [experiment, sample]);

  if (target === null) return <ResolvingFallback />;
  return <Navigate to={target} replace />;
}
