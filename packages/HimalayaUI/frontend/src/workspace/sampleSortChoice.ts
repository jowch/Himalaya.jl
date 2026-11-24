import { useCallback, useState } from "react";
import type { SampleSort } from "./sampleSort";

const keyFor = (experimentId: number): string => `himalaya-ui:sample-sort:${experimentId}`;

function read(experimentId: number): SampleSort {
  try {
    return window.localStorage.getItem(keyFor(experimentId)) === "name" ? "name" : "measured";
  } catch {
    return "measured";
  }
}

export function useSampleSort(experimentId: number): [SampleSort, (sort: SampleSort) => void] {
  const [chosen, setChosen] = useState<{ experimentId: number; sort: SampleSort }>(() => ({ experimentId, sort: read(experimentId) }));
  const sort = chosen.experimentId === experimentId ? chosen.sort : read(experimentId);
  const choose = useCallback((next: SampleSort) => {
    setChosen({ experimentId, sort: next });
    try {
      if (next === "measured") window.localStorage.removeItem(keyFor(experimentId));
      else window.localStorage.setItem(keyFor(experimentId), next);
    } catch {
      return;
    }
  }, [experimentId]);
  return [sort, choose];
}
