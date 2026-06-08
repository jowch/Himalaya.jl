import type { GroupingMode } from "./lib/comparison/coloring";

export interface User {
  id: number;
  username: string;
  first_name: string | null;
  last_name: string | null;
}

export interface Experiment {
  id: number;
  name: string | null;
  path: string;
  data_dir: string;
  analysis_dir: string;
  manifest_path: string | null;
  created_at: string;
  q_units: string | null;
}

export interface SampleTag {
  id: number;
  key: string;
  value: string;
  source: string;
}

export interface Sample {
  id: number;
  experiment_id: number;
  name: string | null;
  display_name: string | null;
  notes: string | null;
  tags: SampleTag[];
}

// Corpus samples carry q_units (resolved from the owning experiment's
// config) — the per-experiment Sample does not. Phase 3 normalization
// reads this field. Returned by the corpus-wide GET /api/samples route.
//
// `screened` and `phase` are optional, forward-looking seams for the contact
// sheet (R1 / #224):
//   - `screened`: the per-sample triage flag owned by #162's backend (not yet
//     wired). Until it lands, the contact sheet derives screened-ness from the
//     exposures (see lib/sample/screened.ts); this field, when present, wins.
//   - `phase`: the sample's resolved liquid-crystalline phase, surfaced by a
//     future indexing-rollup route. When present the status cell shows a phase
//     chip; when absent it shows the hollow-dot "Not indexed" affordance (M-6).
export interface CorpusSample extends Sample {
  q_units: string;
  screened?: boolean;
  phase?: string | null;
}

export class ApiError extends Error {
  constructor(public status: number, message: string, public body: unknown) {
    super(message);
    this.name = "ApiError";
  }
}

export interface AuthOpts {
  username?: string;
  clientId?: string;
  clientOpId?: string;
}

export async function request<T>(
  method: string,
  path: string,
  body?: unknown,
  opts?: AuthOpts,
): Promise<T> {
  const headers: Record<string, string> = {};
  if (body !== undefined) headers["Content-Type"] = "application/json";
  if (opts?.username) headers["X-Username"] = opts.username;
  if (opts?.clientId) headers["X-Client-Id"] = opts.clientId;
  if (opts?.clientOpId) headers["X-Client-Op-Id"] = opts.clientOpId;

  const init: RequestInit = { method, headers };
  if (body !== undefined) init.body = JSON.stringify(body);

  const res = await fetch(path, init);
  if (res.status === 204) return undefined as T;

  const text = await res.text();
  const parsed = text ? safeJson(text) : null;
  if (!res.ok) {
    const msg = parsed && typeof parsed === "object" && parsed !== null && "error" in parsed
      ? String((parsed as { error: unknown }).error)
      : `${method} ${path} failed with ${res.status}`;
    throw new ApiError(res.status, msg, parsed);
  }
  return parsed as T;
}

function safeJson(s: string): unknown {
  try { return JSON.parse(s); } catch { return s; }
}

// Users
export const listUsers  = () => request<User[]>("GET", "/api/users");
export const createUser = (
  username: string,
  fields?: { first_name?: string; last_name?: string },
  opts?: AuthOpts,
) => request<User>("POST", "/api/users", { username, ...fields }, opts);

// Experiments
export const listExperiments = () =>
  request<Experiment[]>("GET", "/api/experiments");
export const getExperiment = (id: number) =>
  request<Experiment>("GET", `/api/experiments/${id}`);
export const updateExperiment = (
  id: number,
  patch: Record<string, never>,
  opts?: AuthOpts,
) => request<Experiment>("PATCH", `/api/experiments/${id}`, patch, opts);

// Samples
export const listSamples    = (experiment_id: number) =>
  request<Sample[]>("GET", `/api/experiments/${experiment_id}/samples`);
export const listCorpusSamples = (): Promise<CorpusSample[]> =>
  request<CorpusSample[]>("GET", "/api/samples");
export const updateSample   = (id: number, patch: { display_name?: string; notes?: string }, opts?: AuthOpts) =>
  request<Sample>("PATCH", `/api/samples/${id}`, patch, opts);
export const addSampleTag   = (id: number, key: string, value: string, opts?: AuthOpts) =>
  request<SampleTag>("POST", `/api/samples/${id}/tags`, { key, value }, opts);
export const removeSampleTag = (id: number, tag_id: number, opts?: AuthOpts) =>
  request<void>("DELETE", `/api/samples/${id}/tags/${tag_id}`, undefined, opts);

// Exposures
export interface ExposureTag {
  id: number;
  key: string;
  value: string;
  source: string;
}

export interface Exposure {
  id: number;
  sample_id: number;
  filename: string | null;
  kind: "file" | "averaged" | "background_subtracted";
  selected: boolean;
  status: "accepted" | "rejected" | null;
  image_path: string | null;
  /**
   * Cache-busting token combining IMAGE_PROCESSING_VERSION (server-side
   * constant) with the source TIFF's mtime. Empty string when no image.
   * Append as `?v=<image_version>` to the image URL so the browser's
   * cache key changes when the image content actually changes.
   */
  image_version: string;
  tags: ExposureTag[];
  sources: unknown[];
  trace_hash: string | null;
  analysis_inputs_hash: string | null;
}

export const listExposures = (sample_id: number) =>
  request<Exposure[]>("GET", `/api/samples/${sample_id}/exposures`);

export const setExposureStatus = (
  id: number,
  status: "accepted" | "rejected" | null,
  opts?: AuthOpts,
) => request<{ id: number; status: string | null }>(
  "PATCH", `/api/exposures/${id}/status`, { status }, opts);

export const selectExposure = (id: number, opts?: AuthOpts) =>
  request<{ id: number; selected: boolean }>(
    "PATCH", `/api/exposures/${id}/select`, {}, opts);

export const addExposureTag = (
  id: number, key: string, value: string, opts?: AuthOpts,
) => request<ExposureTag>(
  "POST", `/api/exposures/${id}/tags`, { key, value }, opts);

export const removeExposureTag = (
  id: number, tag_id: number, opts?: AuthOpts,
) => request<void>(
  "DELETE", `/api/exposures/${id}/tags/${tag_id}`, undefined, opts);

// Trace
export interface Trace {
  q: number[];
  I: number[];
  sigma: number[];
}

export const getTrace = (exposure_id: number) =>
  request<Trace>("GET", `/api/exposures/${exposure_id}/trace`);

/** Batch member traces for a series, keyed by exposure_id. Matches the
 *  `toWaterfallRows(members, tracesById)` contract (a plain Record, number index).
 *  Unresolvable members (no exposure / derived / missing .dat) are absent from the map. */
export const getSeriesTraces = (series_id: number) =>
  request<Record<number, Trace>>("GET", `/api/series/${series_id}/traces`);

// Peaks
export interface Peak {
  id: number;
  exposure_id: number;
  q: number;
  intensity: number | null;
  prominence: number | null;
  sharpness: number | null;
  source: "auto" | "manual";
  /** When true (only meaningful for auto peaks), the peak is soft-disabled by the user. */
  excluded: boolean;
}

/**
 * Backend response for POST /api/exposures/:id/peaks. Carries the inserted
 * peak plus event metadata + the post-state hash so the client can mark the
 * exposure cache fresh without a refetch round-trip.
 */
/**
 * Backend response for POST /api/exposures/:id/peaks. The peak fields are
 * inlined (matches `routes_peaks.jl` which JSON3.writes a flat Dict). Earlier
 * draft typed this as `{peak: Peak, ...}` and mutators read `response.peak.id`
 * — that would throw in production; only the unit tests passed because their
 * mock fixture matched the (wrong) type. Now extends `Peak` like
 * `PeakUpdatedResponse` does.
 */
export interface PeakAddResponse extends Peak {
  event_id: number;
  view_row_id: number;
  analysis_inputs_hash: string;
}

/**
 * Backend response for PATCH /api/peaks/:id. The peak fields are inlined; the
 * event metadata fields are nullable for the no-op case where excluded did
 * not actually change.
 */
export interface PeakUpdatedResponse extends Peak {
  event_id: number | null;
  view_row_id: number | null;
  analysis_inputs_hash: string;
}

/**
 * Backend response for DELETE /api/peaks/:id. Carries the post-state hash so
 * the client can mark the exposure cache fresh inline with the optimistic
 * delete and avoid a transient StaleIndicesBanner flash before the SSE frame
 * arrives.
 */
export interface PeakRemoveResponse {
  event_id: number;
  view_row_id: number | null;
  analysis_inputs_hash: string;
}

export const listPeaks = (exposure_id: number) =>
  request<Peak[]>("GET", `/api/exposures/${exposure_id}/peaks`);
export const addPeak = (exposure_id: number, q: number, opts?: AuthOpts) =>
  request<PeakAddResponse>("POST", `/api/exposures/${exposure_id}/peaks`, { q }, opts);
export const removePeak = (peak_id: number, opts?: AuthOpts) =>
  request<PeakRemoveResponse>("DELETE", `/api/peaks/${peak_id}`, undefined, opts);
export const setPeakExcluded = (peak_id: number, excluded: boolean, opts?: AuthOpts) =>
  request<PeakUpdatedResponse>(
    "PATCH", `/api/peaks/${peak_id}`, { excluded }, opts);

// Indices
export interface IndexPeakRef {
  peak_id: number;
  ratio_position: number;
  residual: number;
  q_observed: number;
}

/** Gauss–Bonnet coexistence flag for an index candidate vs the current
 *  assignment. Display-and-ranking only — never folded into `score` (which
 *  stays coverage×consistency). Recomputed per request, never persisted. */
export interface BonnetFlag {
  predicted_a: number;
  consistent: boolean;
}

export interface IndexEntry {
  id: number;
  exposure_id: number;
  phase: string;
  basis: number;
  score: number | null;
  r_squared: number | null;
  lattice_d: number | null;
  ngc: number | null;
  status: "candidate";
  kind: "auto" | "speculative";
  inputs_hash: string | null;
  peaks: IndexPeakRef[];
  predicted_q: number[];
  /** Bonnet coexistence flag vs the assignment (null when N/A). Rendered as
   *  the ⭙ Bonnet badge in Plan D. */
  bonnet?: BonnetFlag | null;
}

export const listIndices = (exposure_id: number) =>
  request<IndexEntry[]>("GET", `/api/exposures/${exposure_id}/indices`);

export const deleteIndex = (id: number, opts?: AuthOpts) =>
  request<{ deleted: number }>("DELETE", `/api/indices/${id}`, undefined, opts);

// Speculative indices (user-built, sub-minpeaks). Anchor + snap UX.
export interface SpeculativeSnap {
  ratio_position: number;
  predicted_q: number;
  suggested_peak_id: number | null;
  suggested_q: number | null;
  suggested_residual: number | null;
  is_anchor: boolean;
}

export const getSpeculativeSnap = (
  exposure_id: number,
  phase: string,
  anchor_peak_id: number,
  anchor_ratio: number,
) => {
  const qs = `?phase=${encodeURIComponent(phase)}&anchor_peak_id=${anchor_peak_id}&anchor_ratio=${anchor_ratio}`;
  return request<SpeculativeSnap[]>("GET", `/api/exposures/${exposure_id}/speculative-snap${qs}`);
};

export interface SpeculativeAdditional { ratio_position: number; peak_id: number }

export const createSpeculative = (
  exposure_id: number,
  body: {
    phase: string;
    anchor_peak_id: number;
    anchor_ratio: number;
    additional: SpeculativeAdditional[];
    active?: boolean;
  },
  opts?: AuthOpts,
) => request<IndexEntry>("POST", `/api/exposures/${exposure_id}/speculative`, body, opts);

// Custom index (Plan D-9): a client-fitted lattice hypothesis. `basis` is the
// q₁ slope the modal computes via physics (2π/a × first(phaseratios(P))). The
// route persists a speculative index and adds it to the assignment. Response is
// the new IndexEntry (+ queue metadata).
export type CustomIndexResponse = IndexEntry & {
  event_id: number;
  view_row_id: number | null;
};
export const createCustomIndex = (
  exposure_id: number, phase: string, basis: number, opts?: AuthOpts,
) => request<CustomIndexResponse>(
  "POST", `/api/exposures/${exposure_id}/custom-index`, { phase, basis }, opts);

// Assignment (Plan D) — the durable per-exposure 3-state phase assignment cart.
// `state` is explicit, never inferred from members.length: an `indexed`
// assignment with 0 members is a "call in progress"; `form_factor`/`null`
// always carry 0 members. Replaces the retired single-active group model.
export type AssignmentState = "indexed" | "form_factor" | "null";

export interface Assignment {
  exposure_id: number;
  state: AssignmentState;
  members: number[]; // index ids, ascending
}

/**
 * Assignment mutation responses carry queue-framework metadata (event_id,
 * view_row_id) alongside the canonical Assignment body — the mutator's
 * onSuccess strips these before writing the row into the cache.
 */
export type AssignmentMutationResponse = Assignment & {
  event_id: number;
  view_row_id: number | null;
};

export const getAssignment = (exposure_id: number) =>
  request<Assignment>("GET", `/api/exposures/${exposure_id}/assignment`);

export const setAssignmentState = (
  exposure_id: number, state: AssignmentState, opts?: AuthOpts,
) => request<AssignmentMutationResponse>(
  "POST", `/api/exposures/${exposure_id}/assignment/state`, { state }, opts);

export const addAssignmentPhase = (
  exposure_id: number, index_id: number, opts?: AuthOpts,
) => request<AssignmentMutationResponse>(
  "POST", `/api/exposures/${exposure_id}/assignment/members`, { index_id }, opts);

export const removeAssignmentPhase = (
  exposure_id: number, index_id: number, opts?: AuthOpts,
) => request<AssignmentMutationResponse>(
  "DELETE", `/api/exposures/${exposure_id}/assignment/members/${index_id}`, undefined, opts);

// Sample messages (chat log)
export interface SampleMessage {
  id: number;
  sample_id: number;
  author_id: number | null;
  /** Null if the author's user row has been deleted. */
  author: string | null;
  body: string;
  created_at: string;
}

export const listSampleMessages = (sample_id: number) =>
  request<SampleMessage[]>("GET", `/api/samples/${sample_id}/messages`);

export const postSampleMessage = (sample_id: number, body: string, opts?: AuthOpts) =>
  request<SampleMessage>("POST", `/api/samples/${sample_id}/messages`, { body }, opts);

// Analysis
export interface ReanalyzeResponse {
  id: number;
  analyzed: boolean;
  /** New post-analyze hash. Used to clear StaleIndicesBanner inline with
   *  the HTTP response, before the SSE post_state arrives. */
  analysis_inputs_hash: string;
}
export const reanalyzeExposure = (exposure_id: number, opts?: AuthOpts) =>
  request<ReanalyzeResponse>("POST", `/api/exposures/${exposure_id}/analyze`, {}, opts);

// Single-entity fetchers for mention resolution
export const getPeak     = (id: number) => request<Peak>("GET", `/api/peaks/${id}`);
export const getIndex    = (id: number) => request<IndexEntry>("GET", `/api/indices/${id}`);
export const getExposure = (id: number) => request<Exposure>("GET", `/api/exposures/${id}`);
export const getSample   = (id: number) => request<Sample>("GET", `/api/samples/${id}`);

// ─── Comparisons (Plan §Phase 3) ────────────────────────────────────────────
//
// Shapes mirror the Julia route emit / `fetch_comparison_with_members` in
// `comparisons.jl`. The `MemberSnapshot` type lives here because it must be
// the single source of truth for both the HTTP response parser AND the SSE
// `applyRemoteToCache` handler — both paths must produce the same parsed
// shape to avoid cache divergence during reconciliation. (The Compare-era
// client-side `lib/comparison/snapshot.ts` deriver was removed in I5.3 #184.)

export interface MemberSnapshotPeak {
  id: number;
  q: number;
  intensity: number | null;
  sharpness: number;
  source: "auto" | "manual";
}

export interface MemberSnapshotConfirmedIndex {
  id: number;
  phase: string;
  lattice_d: number;
  r_squared: number;
  ngc: number;
  peak_ids: number[];
}

/** One assigned phase + its fitted lattice (Plan E E-4). `lattice_d` is the
 *  index lattice parameter (`a` for cubics, `d` for lamellar/hexagonal). */
export interface MemberSnapshotPhase {
  phase: string;
  lattice_d: number | null;
}

export interface MemberSnapshot {
  effective_peaks: MemberSnapshotPeak[];
  confirmed_index: MemberSnapshotConfirmedIndex | null;
  /** Durable 3-state assignment (Plan E E-7). Older snapshots predate this
   *  field; treat a missing value as "indexed". */
  assignment_state?: AssignmentState;
  /** Distinct phases the member's assignment carries (Plan E E-4), each with
   *  its fitted lattice. Drives the coexistence reading / member rows / strip
   *  cells (both lattices under coexistence). Empty for form-factor / null
   *  members. Missing on older snapshots → derive from confirmed_index. */
  confirmed_phases?: MemberSnapshotPhase[];
  analysis_inputs_hash: string;
}

/** Per-member input shape sent to `POST /api/comparisons` and `/submit`. */
export interface ComparisonMemberInput {
  /** Existing member id (UPDATE on submit) or null/undefined (INSERT/create). */
  id?: number | null;
  exposure_id: number | null;
  display_order: number;
  band_height?: number;
  y_offset?: number;
  normalization?: string;
  color_override?: string | null;
  label_override?: string | null;
  q_window_min?: number | null;
  q_window_max?: number | null;
  peak_display?: unknown;
  snapshot: MemberSnapshot;
}

/** Per-member shape returned by GET / POST endpoints. Mirrors `fetch_comparison_with_members`. */
export interface ComparisonMember {
  id: number;
  comparison_id: number;
  exposure_id: number | null;
  display_order: number;
  band_height: number;
  y_offset: number;
  normalization: string;
  color_override: string | null;
  label_override: string | null;
  q_window_min: number | null;
  q_window_max: number | null;
  peak_display: unknown;
  snapshot: MemberSnapshot | null;
  is_stale: boolean;
  created_by: number | null;
  created_at: string | null;
}

/**
 * Per-member shape returned by the series GET / POST endpoints.
 * Mirrors `fetch_series_with_plate` (packages/HimalayaUI/src/series.jl).
 * Field-for-field identical to `ComparisonMember` except `comparison_id` →
 * `series_id`. This is the render pipeline's going-forward input type.
 *
 * The two types are intentionally twinned until `ComparisonMember` is retired
 * (I3.6 / I5.3) — tighten a shared field (e.g. `peak_display: unknown`) in
 * both, or the render pipeline and Compare diverge.
 */
export interface SeriesMember {
  id: number;
  series_id: number;
  exposure_id: number | null;
  display_order: number;
  band_height: number;
  y_offset: number;
  normalization: string;
  color_override: string | null;
  label_override: string | null;
  q_window_min: number | null;
  q_window_max: number | null;
  peak_display: unknown;
  snapshot: MemberSnapshot | null;
  is_stale: boolean;
  created_by: number | null;
  created_at: string | null;
}

export interface Comparison {
  id: number;
  title: string;
  description: string | null;
  content_hash: string;
  created_by: number | null;
  created_at: string | null;
  updated_at: string | null;
  forked_from_id: number | null;
  forked_at_hash: string | null;
  forked_from_title: string | null;
  /** Author's persisted view choices (spec §6.4); NULL = author never picked. */
  view_grouping_mode: string | null;
  view_show_peak_ticks: boolean | null;
  view_show_peak_labels: boolean | null;
  last_event_at: string | null;
  members: ComparisonMember[];
}

/** Lightweight summary row used by the listing endpoints. */
export interface ComparisonSummary {
  id: number;
  title: string;
  description: string | null;
  content_hash: string;
  created_by: number | null;
  created_at: string | null;
  updated_at: string | null;
  forked_from_id: number | null;
  forked_at_hash: string | null;
  /**
   * Author's persisted view choices (spec §6.4); NULL = author never picked.
   * `view_grouping_mode` is intentionally loose (`string | null`) on read —
   * permissive-read / strict-write; `SaveComparisonBody` uses `GroupingMode`.
   * Populated by #137 (backend); reads `undefined` until that lands.
   */
  view_grouping_mode: string | null;
  view_show_peak_ticks: boolean | null;
  view_show_peak_labels: boolean | null;
  last_event_at: string | null;
  /**
   * Listing projection fields (see `_comparison_listing_rows`).
   * Populated by #137 (backend); reads `undefined` until that lands.
   */
  author_username: string | null;
  member_count: number;
  /** Backend-capped top-3 distinct phases (`_topk_phases`). */
  member_phases: string[];
  /** True distinct-phase total — drives the `+N more` overflow in the sidebar. */
  member_phase_count: number;
  has_stale_members: boolean;
}

export interface ComparisonMessage {
  id: number;
  comparison_id: number;
  author_id: number | null;
  author: string | null;
  body: string;
  created_at: string;
}

/**
 * Body shape posted to `POST /api/comparisons` (create) and
 * `POST /api/comparisons/:id/submit` (update). The mutator picks the route
 * based on the presence of `id` in the payload.
 */
export interface SaveComparisonBody {
  title: string;
  description?: string | null;
  members: ComparisonMemberInput[];
  /** Required on submit (existing comparison); absent on create. */
  expected_content_hash?: string;
  /** Set when forking; both fields ride together or not at all. */
  forked_from_id?: number | null;
  forked_at_hash?: string | null;
  /**
   * Author's view choices (spec §6.4); omitted = not changed by this save.
   * Strict `GroupingMode` on write vs. permissive `string | null` on the read
   * types is deliberate — see the read-type comments.
   */
  view_grouping_mode?: GroupingMode | null;
  view_show_peak_ticks?: boolean | null;
  view_show_peak_labels?: boolean | null;
}

/**
 * Thrown by a fetcher when the server returns 409 (content_hash drift).
 * Carries the server's `current_hash` and `current_state`. `status` is 409 so
 * the queue's failure-class router treats it as a validation error (no retry,
 * surfaces in `onError`); `useQueueMutation` suppresses the toast on it.
 *
 * There is no longer a conflict surface that consumes it. Compare is retired
 * (#177, replay-only) and the series-commit route is last-write-wins (Plan 6a,
 * no longer 409s). The type is kept because `saveComparison` / `commitSeriesPlate`
 * still defensively parse a 409 should one ever arrive.
 */
export class ConflictError extends Error {
  status = 409 as const;
  constructor(
    public current_hash: string | null,
    // A comparison-submit 409 carries a `Comparison`; a series-commit 409
    // (`commitSeriesPlate`) carries a `Series`. No discriminator.
    public current_state: Comparison | Series | null,
    message?: string,
  ) {
    super(message ?? "content_hash conflict");
    this.name = "ConflictError";
  }
}

/**
 * Save a comparison — create (no id) or submit (id present). Internally calls
 * `POST /api/comparisons` or `POST /api/comparisons/:id/submit`. On 409,
 * throws `ConflictError` with the server's current state attached.
 *
 * Why this branches on id (rather than two separate fetchers): the queue
 * mutator's `request` function takes a single payload — branching here keeps
 * the OpKind-to-route mapping in one place and the mutator's `request`
 * trivially testable.
 */
export async function saveComparison(
  body: SaveComparisonBody,
  comparisonId: number | undefined,
  opts?: AuthOpts,
): Promise<Comparison> {
  const path = comparisonId === undefined
    ? "/api/comparisons"
    : `/api/comparisons/${comparisonId}/submit`;
  try {
    return await request<Comparison>("POST", path, body, opts);
  } catch (err) {
    if (err instanceof ApiError && err.status === 409) {
      const b = err.body as
        | { current_hash?: string; current_state?: Comparison }
        | null;
      throw new ConflictError(
        b?.current_hash ?? null,
        b?.current_state ?? null,
        err.message,
      );
    }
    throw err;
  }
}

export const deleteComparison = (id: number, opts?: AuthOpts) =>
  request<{ id: number; deleted: boolean; event_id: number }>(
    "DELETE", `/api/comparisons/${id}`, undefined, opts);

export const getComparison = (id: number) =>
  request<Comparison>("GET", `/api/comparisons/${id}`);

export const listComparisons = () =>
  request<ComparisonSummary[]>("GET", "/api/comparisons");

export const listExperimentComparisons = (experiment_id: number) =>
  request<ComparisonSummary[]>("GET", `/api/experiments/${experiment_id}/comparisons`);

export const getComparisonForks = (id: number) =>
  request<ComparisonSummary[]>("GET", `/api/comparisons/${id}/forks`);

export const listComparisonMessages = (comparison_id: number) =>
  request<ComparisonMessage[]>("GET", `/api/comparisons/${comparison_id}/messages`);

export const postComparisonMessage = (
  comparison_id: number, body: string, opts?: AuthOpts,
) => request<ComparisonMessage>(
  "POST", `/api/comparisons/${comparison_id}/messages`, { body }, opts);

// ─── Picker / scoping shared types ──────────────────────────────────────────

/** Per-pair shape returned by `GET /api/sample-tags` (corpus) and
 *  the experiment-scoped `/api/experiments/:eid/sample-tags`. */
export interface SampleTagPair {
  key: string;
  value: string;
}

/** Per-row shape returned by `GET /api/experiments/:eid/picker-samples`
 *  and the corpus-wide `GET /api/picker-samples`. */
export interface PickerSampleRow {
  sample: Sample;
  indexing_exposure_id: number | null;
  all_exposures: PickerSampleExposure[];
}

export interface PickerSampleExposure {
  id: number;
  sample_id: number;
  filename: string | null;
  selected: boolean;
}

// ─── Corpus scoping reads (I0.2 / I0.3 corpus siblings; I3.4 consumer) ──────
/** Corpus-wide distinct (key,value) tag pairs — GET /api/sample-tags. Reuses
 *  the existing SampleTagPair shape (same {key,value} contract). */
export const getCorpusSampleTags = (): Promise<SampleTagPair[]> =>
  request<SampleTagPair[]>("GET", "/api/sample-tags");

/** Corpus-wide picker projection — GET /api/picker-samples. Same PickerSampleRow
 *  shape as the experiment-scoped getPickerSamples. */
export const getCorpusPickerSamples = (): Promise<PickerSampleRow[]> =>
  request<PickerSampleRow[]>("GET", "/api/picker-samples");

// ─── Series scoping batch write (I2.6 route, source='scoping') ──────────────
export interface BatchSampleTagEntry { sample_id: number; value: string }
export interface BatchSampleTagResult {
  id: number; sample_id: number; key: string; value: string; source: string;
}
/** Batch sample-tag insert — POST /api/samples/tags/batch. One `key`, N
 *  {sample_id,value}, one `source`. Used by scoping with source='scoping'. */
export const batchSampleTags = (
  key: string,
  tags: BatchSampleTagEntry[],
  source: string,
  opts?: AuthOpts,
): Promise<BatchSampleTagResult[]> =>
  request<BatchSampleTagResult[]>(
    "POST", "/api/samples/tags/batch", { key, tags, source }, opts);

// ─── Comparison pins (Plan §Phase 13, Task 13.2) ────────────────────────────
//
// Per-user pinned comparisons surface at the top of the sidebar. Pin/unpin
// are trivial idempotent state toggles — no `with_idempotency`, no SSE — so
// the API is a straightforward POST/DELETE pair. The list endpoint reads
// `X-Username` and returns a flat array of comparison ids in
// most-recently-pinned-first order.

export const listComparisonPins = (opts?: AuthOpts): Promise<number[]> =>
  request<number[]>("GET", "/api/users/me/comparison-pins", undefined, opts);

export const pinComparison = (id: number, opts?: AuthOpts) =>
  request<{ comparison_id: number; pinned: boolean }>(
    "POST", `/api/comparisons/${id}/pin`, undefined, opts);

export const unpinComparison = (id: number, opts?: AuthOpts) =>
  request<{ comparison_id: number; pinned: boolean }>(
    "DELETE", `/api/comparisons/${id}/pin`, undefined, opts);

// ─── Series (#166 / #167 / #168 — series event-kind cluster) ────────────────
//
// Queue-side scaffolding (#198) + the folio read layer (#173 / I3.3): the
// listing summary type `SeriesSummary` and the read fetchers
// (listSeries / getSeries / forksOfSeries) below back the read hooks
// useSeriesList / useSeries in queries.ts. See the Decision Record in
// docs/superpowers/plans/2026-05-18-series-event-kinds.md.
// Shapes mirror `fetch_series_with_plate` / `_series_listing_rows` in series.jl.

/**
 * Lightweight summary row returned by GET /api/series (the corpus listing).
 * Mirrors `_series_listing_rows` in packages/HimalayaUI/src/series.jl — keep
 * the field set in lockstep with that projection. NOT the full `Series` type
 * (that is the GET-detail shape with nested members/samples).
 *
 * `content_hash === ""` is the draft sentinel (an uncommitted series), the
 * same convention the full `Series` type documents.
 */
export interface SeriesSummary {
  id: number;
  title: string;
  description: string | null;
  content_hash: string;
  created_by: number | null;
  created_at: string | null;
  updated_at: string | null;
  forked_from_id: number | null;
  forked_at_hash: string | null;
  view_grouping_mode: string | null;
  view_show_peak_ticks: boolean | null;
  view_show_peak_labels: boolean | null;
  /** Normalised (SQLite `datetime()`) recency sort key; backend sorts DESC. */
  last_event_at: string | null;
  author_username: string | null;
  member_count: number;
  /** Backend-capped top-3 distinct phases (`_topk_phases`). */
  member_phases: string[];
  /** True distinct-phase total — for a `+N more` overflow if shown. */
  member_phase_count: number;
  has_stale_members: boolean;
  /** Recipe ordering variable (e.g. "LL37 : lipid ratio"); null until scoped. */
  ordering_variable: string | null;
  /** True when the members resolve to >1 distinct `samples.experiment_id`.
   *  Valid because q is absolute (Å⁻¹) — series may legitimately span beamtimes. */
  spans_experiments: boolean;
}

/** The recipe membership — one `series_samples` row. */
export interface SeriesSample {
  id: number;
  series_id: number;
  sample_id: number;
  position: number;
  pinned: boolean;
  excluded: boolean;
}

/** The plate — one `series_members` row. Mirrors `ComparisonMember`. */
export interface SeriesMember {
  id: number;
  series_id: number;
  exposure_id: number | null;
  display_order: number;
  band_height: number;
  y_offset: number;
  normalization: string;
  color_override: string | null;
  label_override: string | null;
  q_window_min: number | null;
  q_window_max: number | null;
  peak_display: unknown;
  snapshot: MemberSnapshot | null;
  is_stale: boolean;
  created_by: number | null;
  created_at: string | null;
}

/** How a series orders its samples. Shared by `Series`, `SaveSeriesBody`, and
 *  the `saveSeries` mutator input so the literal union has one source. */
export type OrderRule = "ascending" | "descending" | "manual";

/** Full nested response from `GET/POST/PATCH /api/series*`. */
export interface Series {
  id: number;
  title: string;
  description: string | null;
  /** Empty string `""` is the draft sentinel: `fetch_series_with_plate`
   *  projects a NULL `content_hash` (an uncommitted draft) to `""`. */
  content_hash: string;
  created_by: number | null;
  created_at: string | null;
  updated_at: string | null;
  forked_from_id: number | null;
  forked_at_hash: string | null;
  forked_from_title: string | null;
  view_grouping_mode: string | null;
  view_show_peak_ticks: boolean | null;
  view_show_peak_labels: boolean | null;
  ordering_variable: string | null;
  order_rule: string;
  state: string;
  members: SeriesMember[];
  samples: SeriesSample[];
}

/**
 * Structural type-guard: is `x` a full `Series` projection (both nested
 * arrays + `state`), as opposed to the partial shape `saveSeries`'s
 * `synthesizeFromSse` yields on the SSE-wins path? Used by the series
 * mutators' `onSuccess` to choose between a cache splice and an invalidate.
 */
export function isFullSeries(x: unknown): x is Series {
  if (typeof x !== "object" || x === null) return false;
  const o = x as Record<string, unknown>;
  return typeof o.id === "number"
    && typeof o.state === "string"
    && Array.isArray(o.members)
    && Array.isArray(o.samples);
}

/** Per-recipe-row input for `POST /api/series` and `PATCH /api/series/:id`. */
export interface SeriesSampleInput {
  sample_id: number;
  position?: number;
  pinned?: boolean;
  excluded?: boolean;
}

/** Per-plate-member input for `POST /api/series/:id/commit`. Members carry no
 *  id — the dispatcher mints them. `snapshot` is server-filled when omitted. */
export interface SeriesMemberInput {
  exposure_id: number | null;
  display_order: number;
  band_height?: number;
  y_offset?: number;
  normalization?: string;
  color_override?: string | null;
  label_override?: string | null;
  q_window_min?: number | null;
  q_window_max?: number | null;
  peak_display?: unknown;
  snapshot?: MemberSnapshot;
}

/** Body for `POST /api/series` (create) and `PATCH /api/series/:id` (recipe). */
export interface SaveSeriesBody {
  title: string;
  description?: string | null;
  samples: SeriesSampleInput[];
  ordering_variable?: string | null;
  order_rule?: OrderRule;
  forked_from_id?: number | null;
  forked_at_hash?: string | null;
  view_grouping_mode?: string | null;
  view_show_peak_ticks?: boolean | null;
  view_show_peak_labels?: boolean | null;
}

/** Body for `POST /api/series/:id/commit`. */
export interface CommitSeriesPlateBody {
  members: SeriesMemberInput[];
}

/**
 * Save a series — create (no id ⇒ `POST /api/series`) or recipe-edit
 * (id present ⇒ `PATCH /api/series/:id`). Branches on id so the queue
 * mutator's `request` stays a single-payload call (mirrors `saveComparison`).
 */
export async function saveSeries(
  body: SaveSeriesBody,
  seriesId: number | undefined,
  opts?: AuthOpts,
): Promise<Series> {
  return seriesId === undefined
    ? request<Series>("POST", "/api/series", body, opts)
    : request<Series>("PATCH", `/api/series/${seriesId}`, body, opts);
}

/**
 * Commit the plate (the old "submit"). On 409 (content_hash drift) throws the
 * typed `ConflictError` carrying the server's `current_hash` + `current_state`
 * (a `Series`), mirroring `saveComparison` (I3.5b). This is the ONLY series
 * fetcher that throws `ConflictError`: recipe-save (`PATCH /api/series/:id`)
 * never reads `expected_content_hash` and never 409s, so `saveSeries` is left
 * untouched.
 */
export async function commitSeriesPlate(
  id: number, body: CommitSeriesPlateBody, opts?: AuthOpts,
): Promise<Series> {
  try {
    return await request<Series>("POST", `/api/series/${id}/commit`, body, opts);
  } catch (err) {
    if (err instanceof ApiError && err.status === 409) {
      const b = err.body as
        | { current_hash?: string; current_state?: Series }
        | null;
      throw new ConflictError(
        b?.current_hash ?? null,
        b?.current_state ?? null,
        err.message,
      );
    }
    throw err;
  }
}

export const deleteSeries = (id: number, opts?: AuthOpts) =>
  request<{ id: number; deleted: boolean; event_id: number }>(
    "DELETE", `/api/series/${id}`, undefined, opts);

// ── Series read fetchers (#173 / I3.3 — the folio read layer) ──

/** Corpus-wide listing for the series folio (#173). GET /api/series. */
export const listSeries = () =>
  request<SeriesSummary[]>("GET", "/api/series");

/** Full nested detail for one series (#175 builder reuses this). */
export const getSeries = (id: number) =>
  request<Series>("GET", `/api/series/${id}`);

/** Forks of one series — same SeriesSummary[] shape as the listing. */
export const forksOfSeries = (id: number) =>
  request<SeriesSummary[]>("GET", `/api/series/${id}/forks`);

// ─── Permalink resolve (Plan §Task 8) ───────────────────────────────────────

export interface ResolveSuccess {
  experiment_id: number;
  experiment_name: string;
  sample_id: number | undefined;
  sample_name: string | undefined;
  exposure_id: number | undefined;
  exposure_filename: string | undefined;
}

export interface ResolveError404 {
  error: "not_found";
  missing: "experiment" | "sample" | "exposure";
  missing_value: string;
  experiment_resolved: { id: number; name: string } | undefined;
  sample_resolved:     { id: number; name: string } | undefined;
}

export interface ResolveError400 {
  error: "ambiguous_params" | "missing_experiment" | "missing_sample";
}

export interface ResolveQuery {
  experiment?: string;
  sample?: string;
  exposure?: string;
  experiment_id?: number;
  sample_id?: number;
  exposure_id?: number;
}

/**
 * Calls /api/resolve. Returns either a `ResolveSuccess` (200), a
 * `ResolveError404`, or a `ResolveError400`. Caller distinguishes on
 * `error` field. Read-only — no auth headers.
 *
 * `signal` exposed for AbortController; caller is responsible for
 * the origin-tag race check (§4.2).
 */
export async function resolve(
  q: ResolveQuery,
  signal: AbortSignal | undefined = undefined,
): Promise<ResolveSuccess | ResolveError404 | ResolveError400> {
  const params = new URLSearchParams();
  if (q.experiment !== undefined) params.set("experiment", q.experiment);
  if (q.sample !== undefined) params.set("sample", q.sample);
  if (q.exposure !== undefined) params.set("exposure", q.exposure);
  if (q.experiment_id !== undefined) params.set("experiment_id", String(q.experiment_id));
  if (q.sample_id !== undefined) params.set("sample_id", String(q.sample_id));
  if (q.exposure_id !== undefined) params.set("exposure_id", String(q.exposure_id));

  const init: RequestInit = {};
  if (signal !== undefined) init.signal = signal;
  const res = await fetch(`/api/resolve?${params.toString()}`, init);
  const body = await res.json();
  return body as ResolveSuccess | ResolveError404 | ResolveError400;
}
