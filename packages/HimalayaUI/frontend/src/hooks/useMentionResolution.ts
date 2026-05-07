import { useMemo } from "react";
import { useQueries } from "@tanstack/react-query";
import * as api from "../api";
import { queryKeys } from "../queries";
import type { MentionToken } from "../lib/renderMentions";

export type ResolvedMention =
  | { type: "peak";       data: api.Peak }
  | { type: "index";      data: api.IndexEntry }
  | { type: "exposure";   data: api.Exposure }
  | { type: "sample";     data: api.Sample }
  | { type: "experiment"; data: api.Experiment }
  | { type: "comparison"; data: api.Comparison };

type ResolutionEntry = ResolvedMention | "loading" | "dead";

function queryForToken(token: MentionToken) {
  const { type, id } = token;
  switch (type) {
    case "peak":
      return { queryKey: queryKeys.peak(id), queryFn: () => api.getPeak(id), retry: false };
    case "index":
      return { queryKey: queryKeys.index(id), queryFn: () => api.getIndex(id), retry: false };
    case "exposure":
      return { queryKey: queryKeys.exposure(id), queryFn: () => api.getExposure(id), retry: false };
    case "sample":
      return { queryKey: queryKeys.sample(id), queryFn: () => api.getSample(id), retry: false };
    case "experiment":
      return { queryKey: queryKeys.experiment(id), queryFn: () => api.getExperiment(id), retry: false };
    case "comparison":
      return { queryKey: queryKeys.comparison(id), queryFn: () => api.getComparison(id), retry: false };
  }
}

export function useMentionResolution(tokens: MentionToken[]): Map<string, ResolutionEntry> {
  const queries = useQueries({ queries: tokens.map(queryForToken) });

  return useMemo(() => {
    const map = new Map<string, ResolutionEntry>();
    tokens.forEach((token, i) => {
      const key = `${token.type}:${token.id}`;
      const q   = queries[i];
      if (q.isError) {
        const err = q.error;
        map.set(key, err instanceof api.ApiError && err.status === 404 ? "dead" : "loading");
      } else if (q.isSuccess && q.data !== undefined) {
        map.set(key, { type: token.type, data: q.data } as ResolvedMention);
      } else {
        map.set(key, "loading");
      }
    });
    return map;
  }, [tokens, queries]);
}
