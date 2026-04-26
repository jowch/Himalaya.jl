export interface User { id: number; username: string; }

export interface Experiment {
  id: number;
  name: string | null;
  path: string;
  data_dir: string;
  analysis_dir: string;
  manifest_path: string | null;
  created_at: string;
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

export interface AuthOpts { username?: string }

async function request<T>(
  method: string,
  path: string,
  body?: unknown,
  opts?: AuthOpts,
): Promise<T> {
  const headers: Record<string, string> = {};
  if (body !== undefined) headers["Content-Type"] = "application/json";
  if (opts?.username && method !== "GET") headers["X-Username"] = opts.username;

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
export const createUser = (username: string, opts?: AuthOpts) =>
  request<User>("POST", "/api/users", { username }, opts);

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
  tags: ExposureTag[];
  sources: unknown[];
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

export interface PeakCreated extends Peak { stale_indices: number }

export const listPeaks = (exposure_id: number) =>
  request<Peak[]>("GET", `/api/exposures/${exposure_id}/peaks`);
export const addPeak = (exposure_id: number, q: number, opts?: AuthOpts) =>
  request<PeakCreated>("POST", `/api/exposures/${exposure_id}/peaks`, { q }, opts);
export const removePeak = (peak_id: number, opts?: AuthOpts) =>
  request<void>("DELETE", `/api/peaks/${peak_id}`, undefined, opts);
export const setPeakExcluded = (peak_id: number, excluded: boolean, opts?: AuthOpts) =>
  request<Peak & { stale_indices: number }>(
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
  status: "candidate" | "stale";
  peaks: IndexPeakRef[];
  predicted_q: number[];
}

export const listIndices = (exposure_id: number) =>
  request<IndexEntry[]>("GET", `/api/exposures/${exposure_id}/indices`);

// Groups
export interface GroupEntry {
  id: number;
  exposure_id: number;
  kind: "auto" | "custom";
  active: boolean;
  members: number[];
}

export const listGroups = (exposure_id: number) =>
  request<GroupEntry[]>("GET", `/api/exposures/${exposure_id}/groups`);
export const addIndexToGroup = (group_id: number, index_id: number, opts?: AuthOpts) =>
  request<GroupEntry>("POST", `/api/groups/${group_id}/members`, { index_id }, opts);
export const removeIndexFromGroup = (group_id: number, index_id: number, opts?: AuthOpts) =>
  request<GroupEntry>("DELETE", `/api/groups/${group_id}/members/${index_id}`, undefined, opts);

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
