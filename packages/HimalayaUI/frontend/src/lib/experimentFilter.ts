import type { CorpusSample, Experiment } from "../api";

/**
 * The one honest label both surfaces use for a ?experiment that names nothing.
 * Exported so the experiment filter UI and the contact-sheet EmptyState stay in
 * literal agreement (SA-F5): a single predicate AND a single copy string.
 */
export const UNKNOWN_EXPERIMENT_LABEL = "Unknown experiment";

export interface ExperimentFilter {
  /** Parsed ?experiment id. Absent or malformed params filter nothing. */
  id: number | undefined;
  /** True only once BOTH lists are loaded and the id names no experiment and
   *  no sample carries it. While either list is still loading we claim
   *  nothing — an unknown verdict must never be a loading artifact. */
  unknown: boolean;
  /** The experiment's REAL name when one exists. Undefined for nameless or
   *  sample-vouched experiments — fallback labels are presentation, so each
   *  consumer formats its own ("experiment N" h1 vs "Experiment N" option). */
  name: string | undefined;
}

/**
 * Shared resolver for the `?experiment=` filter (SA-F5). Both the
 * SamplesPage body and the experiment filter UI judge the filter through this
 * single predicate so the two surfaces can never disagree about whether the
 * URL names a real experiment.
 */
export function resolveExperimentFilter(
  raw: string | null,
  experiments: Experiment[] | undefined,
  samples: CorpusSample[] | undefined,
): ExperimentFilter {
  const id = raw !== null && /^\d+$/.test(raw) ? Number(raw) : undefined;
  if (id === undefined) return { id: undefined, unknown: false, name: undefined };

  const record = experiments?.find((e) => e.id === id);
  if (record) return { id, unknown: false, name: record.name ?? undefined };

  // No /experiments record, but the corpus can still vouch for the id: a
  // sample carrying this experiment_id makes the filter real (the two
  // endpoints can disagree, e.g. transiently after a reingest).
  const vouched = samples?.some((s) => s.experiment_id === id) ?? false;
  if (vouched) return { id, unknown: false, name: undefined };

  const loaded = experiments !== undefined && samples !== undefined;
  return { id, unknown: loaded, name: undefined };
}
