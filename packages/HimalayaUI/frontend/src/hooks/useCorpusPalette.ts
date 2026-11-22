import { useMemo, useState } from "react";
import { useCorpusSamples, useExperiments } from "../queries";
import { sampleDisplayName } from "../lib/sample/displayName";
import { workspaceHref } from "../lib/url/workspaceHref";

export interface CorpusExperiment {
  id: number;
  name: string;
}

export interface CorpusSampleRow {
  id: number;
  name: string;
  experimentId: number;
  experiment: string;
  href: string;
}

export function useCorpusPalette(): { experiments: CorpusExperiment[]; samples: CorpusSampleRow[]; onType: (query: string) => void } {
  const experimentsQ = useExperiments();
  const [searched, setSearched] = useState(false);
  const corpusQ = useCorpusSamples({ enabled: searched });
  const experiments = useMemo(
    () => (experimentsQ.data ?? []).map((e) => ({ id: e.id, name: e.name ?? `Experiment ${e.id}` })),
    [experimentsQ.data],
  );
  const samples = useMemo(() => {
    const names = new Map(experiments.map((e) => [e.id, e.name]));
    return (corpusQ.data ?? []).map((s) => ({
      id: s.id,
      name: sampleDisplayName(s),
      experimentId: s.experiment_id,
      experiment: names.get(s.experiment_id) ?? `Experiment ${s.experiment_id}`,
      href: workspaceHref(s.experiment_id, s.id),
    }));
  }, [corpusQ.data, experiments]);
  const onType = (query: string): void => {
    if (query.trim()) setSearched(true);
  };
  return { experiments, samples, onType };
}
