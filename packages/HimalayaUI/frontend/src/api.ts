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
  ingest_status: IngestStatus;
  ingest_error?: string | null;
  image_pattern: string | null;
  metadata_pattern: string | null;
  integration_pattern: string | null;
  // Roll-up counts returned by both list and detail endpoints.
  stats?: { loads: number; samples: number; exposures: number; sessions: number; span_hours: number; started_at: string | null };
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

// Corpus samples carry q_units (resolved from the owning experiment's
// config) — the per-experiment Sample does not. Phase 3 normalization
// reads this field. Returned by the corpus-wide GET /api/samples route.
//
// `phase` is an optional, forward-looking seam for the contact sheet (R1 / #224):
// the sample's resolved liquid-crystalline phase, surfaced by a future
// indexing-rollup route. When present the status cell shows a phase chip; when
// absent it shows the "Not indexed" affordance (M-6).
export interface CorpusSample extends Sample {
  q_units: string;
  phase?: string | null;
  /** The representative exposure's durable assignment state (`representative_exposure_id`:
   *  the selected frame if kept, else the newest kept frame). Drives the contact-sheet status: `form_factor` shows a
   *  distinct "Form factor" status, `indexed` with a `phase` shows the chip,
   *  everything else reads "Not indexed". Absent on older payloads → treated as
   *  unindexed. */
  assignment_state?: AssignmentState;
  call_phases?: string[];
  representative_exposure_id?: number | null;
}

export class ApiError extends Error {
  constructor(public status: number, message: string, public body: unknown) {
    super(message);
    this.name = "ApiError";
  }
}

export class NetworkError extends Error {
  constructor(public method: string, public path: string, public override cause?: unknown) {
    super(`${method} ${path} could not reach the server`);
    this.name = "NetworkError";
  }
}

export interface AuthOpts {
  username?: string;
  clientId?: string;
  clientOpId?: string;
  signal?: AbortSignal;
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

  new Headers(headers);
  const init: RequestInit = { method, headers };
  if (opts?.signal) init.signal = opts.signal;
  if (body !== undefined) init.body = JSON.stringify(body);

  let res: Response;
  try {
    res = await fetch(path, init);
  } catch (e) {
    if (e instanceof DOMException && e.name === "AbortError") throw e;
    throw new NetworkError(method, path, e);
  }
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

/** Create-from-directory (spec §9.2). Returns the new experiment id
 *  immediately; the first scan runs async with progress over SSE. */
export interface CreateExperimentBody {
  /** Legacy single-path field (= data_dir). Optional now that the funnel sends
   *  explicit data_dir/analysis_dir/name confirmed by the user in Configuration. */
  path?: string;
  data_dir?: string;
  analysis_dir?: string;
  name?: string;
  patterns?: { image?: string; metadata?: string; integration?: string };
  /** First-run geometry, committed at create. The funnel preview already derived
   *  geometry (with the confirmed setup file), so Approve sends the WHOLE thing —
   *  each value with its honest source ('setup'/'prp'/'computed', or 'user' for a
   *  manual edit). Persisted verbatim; the scan only fills fields left unset, so
   *  geometry is derived once. A value with no source defaults to 'user'. */
  geometry?: {
    beam_center_x?: number; beam_center_x_source?: string;
    beam_center_y?: number; beam_center_y_source?: string;
    flight_path_m?: number; flight_path_m_source?: string;
    pixel_size_um?: number; pixel_size_um_source?: string;
    energy_kev?: number;    energy_kev_source?: string;
  };
}
export const createExperiment = (body: CreateExperimentBody, opts?: AuthOpts) =>
  request<Experiment>("POST", "/api/experiments", body, opts);

/** Rescan: cheap change-check then additive ingest of new files. Idempotent.
 *  Pass `force = true` to bypass the cheap change-check and always run a full
 *  scan — used when a pattern field edit changes what files are discovered. */
export const triggerScan = (id: number, opts?: AuthOpts, force = false) =>
  request<Experiment>("POST", `/api/experiments/${id}/scan`, { force }, opts);

export interface ProblemExposure {
  exposure_id: number;
  filename: string | null;
  sample_id: number | null;
  sample_name: string | null;
  problem: "no_trace" | "failed";
  reason: string;
}
export interface IngestProblems { ingest_error: string | null; exposures: ProblemExposure[] }
export const listProblems = (id: number) =>
  request<IngestProblems>("GET", `/api/experiments/${id}/problems`);

/** Directory-picker path autocomplete (spec §9.2, read-only). */
export interface PathSuggestResponse { suggestions: string[] }
export const suggestPaths = (prefix: string) =>
  request<PathSuggestResponse>(
    "GET", `/api/fs/suggest?prefix=${encodeURIComponent(prefix)}`);

export interface SetupPreview {
  root: string;
  name: string;
  data_dir: string;
  analysis_dir: string | null;
  setup_file: string | null;
  setup_ambiguous: boolean;
  layout_detected: boolean;
  patterns: { image: string; metadata: string; integration: string };
  exposures: number;
  without_metadata: number;
  without_trace: number;
  geometry: ManifestGeometry;
  discrepancies: ManifestDiscrepancy[];
  duplicate_of: { id: number; name: string } | null;
}
export interface FolderPreview { path: string; setups: SetupPreview[] }
export const previewFolder = (path: string) =>
  request<FolderPreview>("GET", `/api/fs/preview?path=${encodeURIComponent(path)}`);

/** Phase-1 manifest for the Configuration first-run step (spec §6.5).
 *  `GET /api/fs/manifest` with query params derived from the draft path +
 *  optional glob overrides. Returns per-type matched counts and the list of
 *  unmatched files so the user can tune patterns before creating the
 *  experiment. `patterns` keys are optional — omitting them lets the backend
 *  use its defaults. */
export interface ManifestUnmatched { file: string; miss: string; near?: string }
export interface ManifestGeometry {
  energy_kev: number | null; energy_kev_source: string;
  flight_path_m: number | null; flight_path_m_source: string;
  beam_center_x: number | null; beam_center_x_source: string;
  beam_center_y: number | null; beam_center_y_source: string;
  pixel_size_um: number | null; pixel_size_um_source: string;
}
export interface ManifestDiscrepancy { field: string; message: string }
export interface ManifestResponse {
  total: number;
  matched: { image: number; metadata: number; integration: number };
  unmatched: ManifestUnmatched[];
  geometry?: ManifestGeometry;
  discrepancies?: ManifestDiscrepancy[];
  matched_files?: string[];
}
export const fetchManifest = (
  path: string,
  patterns: { image?: string; metadata?: string; integration?: string } = {},
  setupFile?: string,
  analysisDir?: string,
): Promise<ManifestResponse> => {
  const params = new URLSearchParams({ path });
  if (setupFile)            params.set("setup_file",          setupFile);
  // Integration (.dat) is matched against the analysis subtree, mirroring the
  // real scan — so the preview's integration count reflects where .dat lives.
  if (analysisDir)          params.set("analysis_dir",        analysisDir);
  if (patterns.image)       params.set("image_pattern",       patterns.image);
  if (patterns.metadata)    params.set("metadata_pattern",    patterns.metadata);
  if (patterns.integration) params.set("integration_pattern", patterns.integration);
  return request<ManifestResponse>("GET", `/api/fs/manifest?${params.toString()}`);
};

// Samples
export const listSamples    = (experiment_id: number) =>
  request<Sample[]>("GET", `/api/experiments/${experiment_id}/samples`);
export const listCorpusSamples = (experimentId?: number): Promise<CorpusSample[]> =>
  request<CorpusSample[]>("GET", experimentId === undefined ? "/api/samples" : `/api/samples?experiment_id=${experimentId}`);
export const updateSample   = (id: number, patch: { name?: string; notes?: string }, opts?: AuthOpts) =>
  request<Sample>("PATCH", `/api/samples/${id}`, patch, opts);
export const addSampleTag   = (id: number, key: string, value: string, opts?: AuthOpts) =>
  request<SampleTag>("POST", `/api/samples/${id}/tags`, { key, value }, opts);
export const removeSampleTag = (id: number, tag_id: number, opts?: AuthOpts) =>
  request<void>("DELETE", `/api/samples/${id}/tags/${tag_id}`, undefined, opts);
export const editSampleTag = (id: number, tag_id: number, patch: { key?: string; value?: string }, opts?: AuthOpts) =>
  request<SampleTag>("PATCH", `/api/samples/${id}/tags/${tag_id}`, patch, opts);

export interface SampleBridge {
  sample_id: number;
  group_id: number;
}

export const listBridges = (experimentId: number): Promise<SampleBridge[]> =>
  request<SampleBridge[]>("GET", `/api/experiments/${experimentId}/bridges`);

export const bridgeSamples = (sampleId: number, toSampleId: number, opts?: AuthOpts): Promise<SampleBridge[]> =>
  request<SampleBridge[]>("POST", `/api/samples/${sampleId}/bridge`, { to_sample_id: toSampleId }, opts);

export const unbridgeSample = (sampleId: number, opts?: AuthOpts): Promise<SampleBridge[]> =>
  request<SampleBridge[]>("DELETE", `/api/samples/${sampleId}/bridge`, undefined, opts);

export type SampleLabelState = "suggested" | "confirmed";

export interface SampleLabel {
  id: number;
  sample_id: number;
  label: string;
  state: SampleLabelState;
  source: string;
}

export const listLabels = (experimentId: number): Promise<SampleLabel[]> =>
  request<SampleLabel[]>("GET", `/api/experiments/${experimentId}/labels`);

export const addLabel = (sampleId: number, label: string, opts?: AuthOpts): Promise<SampleLabel[]> =>
  request<SampleLabel[]>("POST", `/api/samples/${sampleId}/labels`, { label }, opts);

export const removeLabel = (sampleId: number, label: string, opts?: AuthOpts): Promise<SampleLabel[]> =>
  request<SampleLabel[]>("DELETE", `/api/samples/${sampleId}/labels?label=${encodeURIComponent(label)}`, undefined, opts);

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
  exposure_time: number | null;
  counting_time: number | null;
  transmitted_counts: number | null;
  incident_monitor: number | null;
  timestamp: string | null;
  beam_energy_ev: number | null;
  pipe_length_m: number | null;
  detector: string | null;
}

export const listExposures = (sample_id: number, opts?: { signal?: AbortSignal }) =>
  request<Exposure[]>("GET", `/api/samples/${sample_id}/exposures`, undefined, opts);

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
export interface TraceMonitors {
  transmitted_counts: number | null;
  exposure_time: number | null;
  counting_time: number | null;
}

export interface Trace {
  q: number[];
  I: number[];
  sigma: number[];
  monitors?: TraceMonitors;
}

export const getTrace = (exposure_id: number, opts?: { signal?: AbortSignal }) =>
  request<Trace>("GET", `/api/exposures/${exposure_id}/trace`, undefined, opts);

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
 * Backend response for PATCH /api/peaks/auto/:id. The peak fields are inlined; the
 * event metadata fields are nullable for the no-op case where excluded did
 * not actually change.
 */
export interface PeakUpdatedResponse extends Peak {
  event_id: number | null;
  view_row_id: number | null;
  analysis_inputs_hash: string;
}

export interface PeakRemoveResponse {
  event_id: number;
  view_row_id: number | null;
  analysis_inputs_hash: string;
}

export const listPeaks = (exposure_id: number, opts?: { signal?: AbortSignal }) =>
  request<Peak[]>("GET", `/api/exposures/${exposure_id}/peaks`, undefined, opts);
export const addPeak = (exposure_id: number, q: number, opts?: AuthOpts) =>
  request<PeakAddResponse>("POST", `/api/exposures/${exposure_id}/peaks`, { q }, opts);
export const removeManualPeak = (peak_id: number, opts?: AuthOpts) =>
  request<PeakRemoveResponse>("DELETE", `/api/peaks/manual/${peak_id}`, undefined, opts);
export const setAutoPeakExcluded = (peak_id: number, excluded: boolean, opts?: AuthOpts) =>
  request<PeakUpdatedResponse>(
    "PATCH", `/api/peaks/auto/${peak_id}`, { excluded }, opts);

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

export const listIndices = (exposure_id: number, opts?: { signal?: AbortSignal }) =>
  request<IndexEntry[]>("GET", `/api/exposures/${exposure_id}/indices`, undefined, opts);

export const deleteIndex = (id: number, opts?: AuthOpts) =>
  request<{ deleted: number }>("DELETE", `/api/indices/${id}`, undefined, opts);

// Speculative indices (user-built, sub-minpeaks).

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
  exposure_id: number, phase: string, basis: number, ratios: number[], opts?: AuthOpts,
) => request<CustomIndexResponse>(
  "POST", `/api/exposures/${exposure_id}/custom-index`, { phase, basis, ratios }, opts);

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

export const getAssignment = (exposure_id: number, opts?: { signal?: AbortSignal }) =>
  request<Assignment>("GET", `/api/exposures/${exposure_id}/assignment`, undefined, opts);

export interface SampleCallMember {
  index_id: number;
  exposure_id: number;
  phase: string;
}

export interface SampleCall {
  sample_id: number;
  state: AssignmentState;
  phase: string | null;
  call_phases: string[];
  members: SampleCallMember[];
}

export type SampleCallMutationResponse = SampleCall & {
  event_id: number;
  view_row_id: number | null;
};

export const getSampleCall = (sample_id: number) =>
  request<SampleCall>("GET", `/api/samples/${sample_id}/call`);

export const setSampleCallState = (
  sample_id: number, state: AssignmentState, opts?: AuthOpts,
) => request<SampleCallMutationResponse>(
  "POST", `/api/samples/${sample_id}/call/state`, { state }, opts);

export const addAssignmentPhase = (
  exposure_id: number, index_id: number, opts?: AuthOpts,
) => request<AssignmentMutationResponse>(
  "POST", `/api/exposures/${exposure_id}/assignment/members`, { index_id }, opts);

export const removeAssignmentPhase = (
  exposure_id: number, index_id: number, opts?: AuthOpts,
) => request<AssignmentMutationResponse>(
  "DELETE", `/api/exposures/${exposure_id}/assignment/members/${index_id}`, undefined, opts);

// Analysis
export interface ReanalyzeResponse {
  id: number;
  analyzed: boolean;
  analysis_inputs_hash: string;
}
export const reanalyzeExposure = (exposure_id: number, opts?: AuthOpts) =>
  request<ReanalyzeResponse>("POST", `/api/exposures/${exposure_id}/analyze`, {}, opts);

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


// ─── Corpus scoping reads (I0.2 / I0.3 corpus siblings; I3.4 consumer) ──────
/** Corpus-wide distinct (key,value) tag pairs — GET /api/sample-tags. Reuses
 *  the existing SampleTagPair shape (same {key,value} contract). */
export const getCorpusSampleTags = (): Promise<SampleTagPair[]> =>
  request<SampleTagPair[]>("GET", "/api/sample-tags");


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
