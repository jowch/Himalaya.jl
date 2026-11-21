import { useMemo } from "react";
import { Navigate, useParams } from "react-router-dom";
import { useQuery } from "@tanstack/react-query";
import * as api from "../../api";
import { useExperiment, useTriggerScan } from "../../queries";
import { useScanState } from "../../hooks/useScanState";
import { ScanFailedPage } from "../pages/ScanFailedPage";
import { ResolvingFallback } from "./ResolvingFallback";

export function ExperimentEntry(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const exp = useExperiment(expId);
  const scan = useScanState(expId);
  const triggerScan = useTriggerScan(expId);

  const patterns = useMemo(() => ({
    ...(exp.data?.image_pattern != null && { image: exp.data.image_pattern }),
    ...(exp.data?.metadata_pattern != null && { metadata: exp.data.metadata_pattern }),
    ...(exp.data?.integration_pattern != null && { integration: exp.data.integration_pattern }),
  }), [exp.data?.image_pattern, exp.data?.metadata_pattern, exp.data?.integration_pattern]);

  const manifest = useQuery({
    queryKey: ["manifest", exp.data?.data_dir ?? "", patterns, exp.data?.analysis_dir ?? ""],
    queryFn: () => api.fetchManifest(exp.data!.data_dir, patterns, undefined, exp.data!.analysis_dir ?? undefined),
    enabled: scan.failed && !!exp.data?.data_dir,
  });

  if (!scan.loaded && !exp.isError) return <ResolvingFallback />;
  if (scan.scanning) return <Navigate to={`/experiments/${expId}/grouping`} replace />;
  if (scan.failed) {
    return (
      <ScanFailedPage
        experimentId={expId}
        unmatched={manifest.data?.unmatched ?? []}
        parsedCount={manifest.data?.matched.image ?? 0}
        onIngestParsed={() => triggerScan.mutate(true)}
        dataDir={exp.data?.data_dir ?? ""}
        {...(exp.data?.analysis_dir != null ? { analysisDir: exp.data.analysis_dir } : {})}
        patterns={patterns}
      />
    );
  }
  return <Navigate to={`/experiments/${expId}/workspace`} replace />;
}
