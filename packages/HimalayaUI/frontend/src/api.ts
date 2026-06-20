export interface User {
  id: number;
  username: string;
  first_name: string | null;
  last_name: string | null;
}

export type GeometrySource = "prp" | "setup" | "user" | "default" | "computed";
export type IngestStatus = "idle" | "scanning" | "analyzing" | "complete" | "failed";

export interface Experiment {
  id: number;
  name: string | null;
  description: string | null;     // Phase E1 additive column (Task 1c migration)
  path: string;
  data_dir: string;
  analysis_dir: string;
  manifest_path: string | null;
  created_at: string;
  q_units: string | null;
  beam_center_x: number | null;
  beam_center_y: number | null;
  pixel_size_um: number | null;
  energy_kev: number | null;
  flight_path_m: number | null;
  // Phase A typed-geometry per-field provenance + scan bookkeeping.
  energy_kev_source: GeometrySource;
  flight_path_m_source: GeometrySource;
  beam_center_x_source: GeometrySource;
  beam_center_y_source: GeometrySource;
  pixel_size_um_source: GeometrySource;
  q_units_source: GeometrySource;
  last_scanned_at: string | null;
  scan_signature: string | null;
  ingest_status: IngestStatus;
  // Phase E1 additive columns: editable file-pattern globs (Task 1c migration).
  // NULL = use legacy experiment.toml fallback.
  image_pattern: string | null;
  metadata_pattern: string | null;
  integration_pattern: string | null;
  // Roll-up counts returned by GET /api/experiments/:id (not the list endpoint).
  stats?: { loads: number; samples: number; exposures: number; sessions: number };
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
  name: string;            // non-null after the collapse (was `string | null` + `display_name`)
  notes: string | null;
  tags: SampleTag[];
}

/** A merge/split discrepancy flag the auto-grouper raised on a slot (spec §8.8).
 *  `null` when the slot is clean. The structural-edit dismissal arm
 *  (`grouping_flag_dismissed`, Phase D) clears it. */
export type GroupingFlag =
  | { kind: "merge"; merge_with_sample_id: number; merge_with_label: string }
  | { kind: "split"; split_at_index: number; jump_from: number; jump_to: number }
  | null;

/** One exposure leaf under a load's sample slot. */
export interface LoadExposure {
  id: number;
  filename: string;
  horizontal_position: number | null;
  timestamp: string | null;
}

/** One sample slot inside a load (a (load, slot) coordinate). `name_source`/
 *  `grouping_source` are provenance tags ("user" | "computed" | …);
 *  `merged_into_id` is non-null when this slot was merged into a sibling. */
export interface LoadSample {
  sample_id: number;
  name: string;
  slot_index: number;
  grouping_source: string;
  name_source: string;
  merged_into_id: number | null;
  flag: GroupingFlag;
  exposures: LoadExposure[];
}

/** One rack-load roll-up (Phase A `loads` table), returned NESTED by
 *  `GET /api/experiments/:id/loads` (see Task 2): Load ▸ Sample ▸ Exposures.
 *  Drives E2's LoadFold/SampleFold/ExposureLeaf + the grouping-review count
 *  (samples whose `flag` is non-null). */
export interface Load {
  load_id: number;
  load_index: number;
  session_id: number | null;
  start_time: string | null;
  end_time: string | null;
  frame_count: number;
  note: string | null;
  samples: LoadSample[];
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
  /** The representative exposure's durable assignment state (selected=1 else
   *  highest-id). Drives the contact-sheet status: `form_factor` shows a
   *  distinct "Form factor" status, `indexed` with a `phase` shows the chip,
   *  everything else reads "Not indexed". Absent on older payloads → treated as
   *  unindexed. */
  assignment_state?: AssignmentState;
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

/** Canonical PATCH body for `PATCH /api/experiments/:id`. **E1 DEFINES,
 *  E2 IMPORTS — never redefine.** All fields are optional; the backend
 *  writes what is present:
 *  - name/description: plain writes, NO *_source stamp, NO rescan.
 *  - Geometry ×6: each field written + *_source stamped 'user' server-side
 *    (already built in Phase C — this widens the same route).
 *  - File patterns ×3: plain write + scan_signature invalidated server-side
 *    so the next scan re-discovers with the new glob.
 *  - data_dir/analysis_dir/path are READ-ONLY (400 if sent). */
export interface ExperimentPatch {
  name?: string;
  description?: string | null;
  energy_kev?: number;
  flight_path_m?: number;
  beam_center_x?: number;
  beam_center_y?: number;
  pixel_size_um?: number;
  q_units?: string;
  image_pattern?: string;
  metadata_pattern?: string;
  integration_pattern?: string;
}

// Experiments
export const listExperiments = () =>
  request<Experiment[]>("GET", "/api/experiments");
export const getExperiment = (id: number) =>
  request<Experiment>("GET", `/api/experiments/${id}`);

export const updateExperiment = (
  id: number,
  patch: ExperimentPatch,
  opts?: AuthOpts,
) => request<Experiment>("PATCH", `/api/experiments/${id}`, patch, opts);

/** Create-from-directory (spec §9.2). Returns the new experiment id
 *  immediately; the first scan runs async with progress over SSE. */
export interface CreateExperimentBody {
  path: string;
  name?: string;
  patterns?: { image?: string; metadata?: string; integration?: string };
}
export const createExperiment = (body: CreateExperimentBody, opts?: AuthOpts) =>
  request<Experiment>("POST", "/api/experiments", body, opts);

export const deleteExperiment = (id: number, opts?: AuthOpts) =>
  request<void>("DELETE", `/api/experiments/${id}`, undefined, opts);

/** Rescan: cheap change-check then additive ingest of new files. Idempotent. */
export const triggerScan = (id: number, opts?: AuthOpts) =>
  request<Experiment>("POST", `/api/experiments/${id}/scan`, {}, opts);

/** The Load ▸ Sample ▸ Exposures roll-up for the grouping-review surface
 *  (spec §9.2 — a dedicated endpoint, distinct from the flat corpus samples). */
export const listLoads = (id: number) =>
  request<Load[]>("GET", `/api/experiments/${id}/loads`);

/** Directory-picker path autocomplete (spec §9.2, read-only). */
export interface PathSuggestResponse { suggestions: string[] }
export const suggestPaths = (prefix: string) =>
  request<PathSuggestResponse>(
    "GET", `/api/fs/suggest?prefix=${encodeURIComponent(prefix)}`);

/** Directory-picker validate-path probe (spec §9.2). `matched`/`scanned` drive
 *  the validation line; `ok=false` + `message` powers the failed-scan preview. */
export interface ValidatePathResponse {
  ok: boolean;
  matched: number;
  scanned: number;
  message: string | null;
}
export const validatePath = (path: string) =>
  request<ValidatePathResponse>(
    "GET", `/api/fs/validate?path=${encodeURIComponent(path)}`);

// Samples
export const listSamples    = (experiment_id: number) =>
  request<Sample[]>("GET", `/api/experiments/${experiment_id}/samples`);
export const listCorpusSamples = (): Promise<CorpusSample[]> =>
  request<CorpusSample[]>("GET", "/api/samples");
export const updateSample   = (id: number, patch: { name?: string; notes?: string }, opts?: AuthOpts) =>
  request<Sample>("PATCH", `/api/samples/${id}`, patch, opts);
export const addSampleTag   = (id: number, key: string, value: string, opts?: AuthOpts) =>
  request<SampleTag>("POST", `/api/samples/${id}/tags`, { key, value }, opts);
export const removeSampleTag = (id: number, tag_id: number, opts?: AuthOpts) =>
  request<void>("DELETE", `/api/samples/${id}/tags/${tag_id}`, undefined, opts);
export const editSampleTag = (id: number, tag_id: number, patch: { key?: string; value?: string }, opts?: AuthOpts) =>
  request<SampleTag>("PATCH", `/api/samples/${id}/tags/${tag_id}`, patch, opts);

// Structural sample/exposure edits (Phase E2 grouping-review surface)
export const renameSample = (id: number, name: string, opts?: AuthOpts): Promise<Sample> =>
  request<Sample>("PATCH", `/api/samples/${id}/name`, { name }, opts);

export const moveExposure = (exposureId: number, sampleId: number, opts?: AuthOpts): Promise<Exposure> =>
  request<Exposure>("POST", `/api/exposures/${exposureId}/move`, { sample_id: sampleId }, opts);

export interface MergeSamplesResponse { loser_id: number; survivor_id: number }
export const mergeSamples = (loserId: number, survivorId: number, opts?: AuthOpts): Promise<MergeSamplesResponse> =>
  request<MergeSamplesResponse>("POST", `/api/samples/${loserId}/merge`, { survivor_id: survivorId }, opts);

export interface SplitSampleResponse { new_sample_id: number }
export const splitSample = (sampleId: number, exposureIds: number[], name: string, opts?: AuthOpts): Promise<SplitSampleResponse> =>
  request<SplitSampleResponse>("POST", `/api/samples/${sampleId}/split`, { exposure_ids: exposureIds, name }, opts);

/** "Keep separate" — durable dismissal of a backend-produced grouping flag
 *  (spec §9.3: grouping_flag_dismissed; suppressed in get_loads_rollup, so it
 *  stays gone across rescans). */
export interface DismissGroupingFlagBody { flag_kind: "merge" | "split"; merge_with_sample_id?: number }
export const dismissGroupingFlag = (sampleId: number, body: DismissGroupingFlagBody, opts?: AuthOpts): Promise<void> =>
  request<void>("POST", `/api/samples/${sampleId}/dismiss-flag`, body, opts);

/** Undo a previous dismissal — re-shows the flag so the sample re-enters
 *  "Needs review". The backend route is POST /api/samples/:id/dismiss-flag/undo.
 *  Symmetric inverse of dismissGroupingFlag (re-show ↔ suppress). */
export const undoDismissGroupingFlag = (sampleId: number, opts?: AuthOpts): Promise<void> =>
  request<void>("POST", `/api/samples/${sampleId}/dismiss-flag/undo`, {}, opts);

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
 * delete and avoid a transient stale-indices alert flash before the SSE frame
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
  /** New post-analyze hash. Used to clear the stale-indices alert inline with
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

// ─── Member snapshots (shared by the Series plate) ──────────────────────────
//
// The `MemberSnapshot` type lives here because it must be the single source of
// truth for both the HTTP response parser AND the SSE `applyRemoteToCache`
// handler — both paths must produce the same parsed shape to avoid cache
// divergence during reconciliation. Consumed by `SeriesMember` (the going-
// forward render-pipeline input type) and the series reading/figure helpers.

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

/** Lightweight summary row used by the listing endpoints. */
/**
 * Thrown by a fetcher when the server returns 409 (content_hash drift).
 * Carries the server's `current_hash` and `current_state`. `status` is 409 so
 * the queue's failure-class router treats it as a validation error (no retry,
 * surfaces in `onError`); `useQueueMutation` suppresses the toast on it.
 *
 * There is no longer a conflict surface that consumes it. The series-commit
 * route is last-write-wins (Plan 6a, no longer 409s). The type is kept because
 * `commitSeriesPlate` still defensively parses a 409 should one ever arrive.
 */
export class ConflictError extends Error {
  status = 409 as const;
  constructor(
    public current_hash: string | null,
    // A series-commit 409 (`commitSeriesPlate`) carries a `Series`.
    public current_state: Series | null,
    message?: string,
  ) {
    super(message ?? "content_hash conflict");
    this.name = "ConflictError";
  }
}

// ─── Picker / scoping shared types ──────────────────────────────────────────

/** Per-pair shape returned by `GET /api/sample-tags` (corpus) and
 *  the experiment-scoped `/api/experiments/:eid/sample-tags`.
 *  `count` is the number of distinct samples carrying this (key, value) pair —
 *  used by the Manage-tags modal to rank suggestions by frequency.
 *  proposeOrdering (scoping) ignores this field and ranks by distinct-value
 *  count instead. */
export interface SampleTagPair {
  key: string;
  value: string;
  count?: number;
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

// ─── Series (#166 / #167 / #168 — series event-kind cluster) ────────────────
//
// Queue-side scaffolding (#198) + the folio read layer (#173 / I3.3): the
// listing summary type `SeriesSummary` and the read fetchers
// (listSeries / getSeries / forksOfSeries) below back the read hooks
// useSeriesList / useSeries in queries.ts. See the Decision Record in
// docs/event-log.md.
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
  /** Beamtime provenance: the members' single experiment's `name` when the
   *  series does NOT span experiments; null when spanning, memberless, or the single experiment has no name. */
  experiment_name: string | null;
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

/** The plate — one `series_members` row. */
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
 * mutator's `request` stays a single-payload call.
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
 * (a `Series`). This is the ONLY series
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
