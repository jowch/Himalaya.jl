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
  label: string | null;
  name: string | null;
  notes: string | null;
  tags: SampleTag[];
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

async function request<T>(
  method: string,
  path: string,
  body?: unknown,
  opts?: AuthOpts,
): Promise<T> {
  const headers: Record<string, string> = {};
  if (body !== undefined) headers["Content-Type"] = "application/json";
  if (opts?.username && method !== "GET") headers["X-Username"] = opts.username;
  if (opts?.clientId && method !== "GET") headers["X-Client-Id"] = opts.clientId;
  if (opts?.clientOpId && method !== "GET") headers["X-Client-Op-Id"] = opts.clientOpId;

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
  patch: Partial<Pick<Experiment, "name" | "data_dir" | "analysis_dir" | "manifest_path">>,
  opts?: AuthOpts,
) => request<Experiment>("PATCH", `/api/experiments/${id}`, patch, opts);

// Samples
export const listSamples    = (experiment_id: number) =>
  request<Sample[]>("GET", `/api/experiments/${experiment_id}/samples`);
export const updateSample   = (id: number, patch: { name?: string; notes?: string }, opts?: AuthOpts) =>
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

export const listExposures = (
  sample_id: number,
  opts?: { excludeRejected?: boolean },
) => {
  const qs = opts?.excludeRejected ? "?exclude_rejected=true" : "";
  return request<Exposure[]>("GET", `/api/samples/${sample_id}/exposures${qs}`);
};

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
export interface PeakAddResponse {
  peak: Peak;
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
export const reanalyzeExposure = (exposure_id: number, opts?: AuthOpts) =>
  request<{ id: number; analyzed: boolean }>("POST", `/api/exposures/${exposure_id}/analyze`, {}, opts);

// Single-entity fetchers for mention resolution
export const getPeak     = (id: number) => request<Peak>("GET", `/api/peaks/${id}`);
export const getIndex    = (id: number) => request<IndexEntry>("GET", `/api/indices/${id}`);
export const getExposure = (id: number) => request<Exposure>("GET", `/api/exposures/${id}`);
export const getSample   = (id: number) => request<Sample>("GET", `/api/samples/${id}`);
