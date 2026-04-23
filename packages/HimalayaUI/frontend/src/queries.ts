import { useQuery } from "@tanstack/react-query";
import * as api from "./api";

export const queryKeys = {
  experiment: (id: number) => ["experiment", id] as const,
  samples:    (experimentId: number) => ["experiment", experimentId, "samples"] as const,
};

export function useExperiment(id: number) {
  return useQuery({
    queryKey: queryKeys.experiment(id),
    queryFn: () => api.getExperiment(id),
  });
}

export function useSamples(experimentId: number) {
  return useQuery({
    queryKey: queryKeys.samples(experimentId),
    queryFn: () => api.listSamples(experimentId),
  });
}
