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
export interface CorpusSample extends Sample {
  q_units: string;
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

// Groups
export interface GroupEntry {
  id: number;
  exposure_id: number;
  kind: "auto" | "custom";
  active: boolean;
  members: number[];
}

/**
 * Mutation responses on group routes carry queue-framework metadata
 * (event_id, view_row_id) alongside the row. The mutator's onSuccess MUST
 * destructure these out before writing the row into the cache — otherwise
 * `GroupEntry` rows get polluted with queue plumbing fields.
 */
export type GroupMutationResponse = GroupEntry & {
  event_id: number;
  view_row_id: number | null;
};

export const listGroups = (exposure_id: number) =>
  request<GroupEntry[]>("GET", `/api/exposures/${exposure_id}/groups`);
export const addIndexToGroup = (group_id: number, index_id: number, opts?: AuthOpts) =>
  request<GroupMutationResponse>("POST", `/api/groups/${group_id}/members`, { index_id }, opts);
export const removeIndexFromGroup = (group_id: number, index_id: number, opts?: AuthOpts) =>
  request<GroupMutationResponse>("DELETE", `/api/groups/${group_id}/members/${index_id}`, undefined, opts);

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
// `comparisons.jl`. The `MemberSnapshot` type lives here (not in
// `lib/comparison/snapshot.ts`) because it must be the single source of
// truth for both the HTTP response parser AND the SSE `applyRemoteToCache`
// handler — both paths must produce the same parsed shape to avoid cache
// divergence during reconciliation.

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

export interface MemberSnapshot {
  effective_peaks: MemberSnapshotPeak[];
  confirmed_index: MemberSnapshotConfirmedIndex | null;
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
 * Thrown by `saveComparison` when the server returns 409 (content_hash drift).
 * Carries the server's `current_hash` and `current_state` so the conflict
 * modal can render the diff. `status` is set to 409 so the queue's failure-
 * class router treats it as a validation error (no retry, surfaces in
 * `onError`); the modal opens off the typed throw, not the toast.
 */
export class ConflictError extends Error {
  status = 409 as const;
  constructor(
    public current_hash: string | null,
    public current_state: Comparison | null,
    message?: string,
  ) {
    super(message ?? "comparison content_hash conflict");
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

// ─── Picker support routes (Plan §Phase 5, Task 5.2) ───────────────────────
//
// Read-only GETs feeding the comparison picker. `recently-picked` returns
// a flat exposure-id list in most-recent-first order; `sample-tags` returns
// distinct (key, value) pairs scoped to one experiment.

/** Per-pair shape returned by `GET /api/experiments/:eid/sample-tags`. */
export interface SampleTagPair {
  key: string;
  value: string;
}

export const getRecentlyPickedExposures = (
  user_id: number, limit?: number,
): Promise<number[]> => {
  const qs = limit !== undefined ? `?limit=${limit}` : "";
  return request<number[]>("GET", `/api/users/${user_id}/recently-picked-exposures${qs}`);
};

export const getSampleTags = (experiment_id: number): Promise<SampleTagPair[]> =>
  request<SampleTagPair[]>("GET", `/api/experiments/${experiment_id}/sample-tags`);

/** Per-row shape returned by `GET /api/experiments/:eid/picker-samples`. */
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

export const getPickerSamples = (
  experiment_id: number,
): Promise<PickerSampleRow[]> =>
  request<PickerSampleRow[]>(
    "GET", `/api/experiments/${experiment_id}/picker-samples`);

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
