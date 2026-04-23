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

let _username: string | undefined = undefined;
export function setUsername(name: string | undefined): void { _username = name; }
export function getUsername(): string | undefined { return _username; }

async function request<T>(method: string, path: string, body?: unknown): Promise<T> {
  const headers: Record<string, string> = {};
  if (body !== undefined) headers["Content-Type"] = "application/json";
  if (_username && method !== "GET") headers["X-Username"] = _username;

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
export const createUser = (username: string) =>
  request<User>("POST", "/api/users", { username });

// Experiments
export const getExperiment = (id: number) =>
  request<Experiment>("GET", `/api/experiments/${id}`);
export const updateExperiment = (id: number, patch: Partial<Pick<Experiment,
  "name" | "data_dir" | "analysis_dir" | "manifest_path">>) =>
  request<Experiment>("PATCH", `/api/experiments/${id}`, patch);

// Samples
export const listSamples    = (experiment_id: number) =>
  request<Sample[]>("GET", `/api/experiments/${experiment_id}/samples`);
export const updateSample   = (id: number, patch: { name?: string; notes?: string }) =>
  request<Sample>("PATCH", `/api/samples/${id}`, patch);
export const addSampleTag   = (id: number, key: string, value: string) =>
  request<SampleTag>("POST", `/api/samples/${id}/tags`, { key, value });
export const removeSampleTag = (id: number, tag_id: number) =>
  request<void>("DELETE", `/api/samples/${id}/tags/${tag_id}`);
