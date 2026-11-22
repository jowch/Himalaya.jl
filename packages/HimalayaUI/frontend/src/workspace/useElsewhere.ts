import { useMemo, useState } from "react";
import { useCorpusSamples, useExperiments } from "../queries";
import { sampleDisplayName } from "../lib/sample/displayName";
import { workspaceHref } from "../lib/url/workspaceHref";

export interface Elsewhere {
  experiments: { id: number; name: string }[];
  samples: { id: number; name: string; experiment: string; href: string }[];
}

export function useElsewhere(experimentId: number): { elsewhere: Elsewhere; onType: (query: string) => void } {
  const experimentsQ = useExperiments();
  const [searched, setSearched] = useState(false);
  const corpusQ = useCorpusSamples({ enabled: searched });
  const elsewhere = useMemo(() => {
    const experiments = (experimentsQ.data ?? []).filter((e) => e.id !== experimentId).map((e) => ({ id: e.id, name: e.name ?? `Experiment ${e.id}` }));
    const names = new Map(experiments.map((e) => [e.id, e.name]));
    const samples = (corpusQ.data ?? [])
      .filter((s) => s.experiment_id !== experimentId)
      .map((s) => ({ id: s.id, name: sampleDisplayName(s), experiment: names.get(s.experiment_id) ?? `Experiment ${s.experiment_id}`, href: workspaceHref(s.experiment_id, s.id) }));
    return { experiments, samples };
  }, [experimentsQ.data, corpusQ.data, experimentId]);
  const onType = (query: string): void => {
    if (query.trim()) setSearched(true);
  };
  return { elsewhere, onType };
}
