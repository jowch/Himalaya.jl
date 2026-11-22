import type { OpKind } from "./types";
import { NetworkError } from "../../api";

interface ErrorWithStatus {
  status?: number;
  message?: string;
}

function asErrorWithStatus(err: unknown): ErrorWithStatus | null {
  if (err === null || err === undefined) return null;
  if (typeof err !== "object") return null;
  return err as ErrorWithStatus;
}

/**
 * 4xx HTTP responses indicate the server rejected the request semantically.
 * The request will not succeed on retry without user input — surface as a
 * Validation toast so the user knows what to change.
 */
export function isValidationError(err: unknown): boolean {
  const e = asErrorWithStatus(err);
  if (!e || typeof e.status !== "number") return false;
  return e.status >= 400 && e.status < 500;
}

/**
 * 404 specifically. Mutators that opt into `treats404AsSuccess` use this to
 * recognise the "already gone on the server" case — typical when a remove
 * retries after a successful first attempt that the client missed (5xx /
 * network blip), where the desired end state already exists.
 */
export function is404Error(err: unknown): boolean {
  const e = asErrorWithStatus(err);
  return !!e && e.status === 404;
}

/**
 * 5xx HTTP, network failures (TypeError from fetch), and timeout errors are
 * infrastructure problems — likely transient. Retry behaviour applies; the
 * InfrastructureBanner (mounted at App.tsx by M2) surfaces persistent failures.
 */
export function isInfrastructureError(err: unknown): boolean {
  const e = asErrorWithStatus(err);
  if (!e) return false;
  if (typeof e.status === "number") {
    return e.status >= 500 && e.status < 600;
  }
  // No status → network or runtime failure. Treat as infrastructure.
  return true;
}

const UNREACHABLE_STATUSES = new Set([502, 503, 504]);

export function isUnreachableError(err: unknown): boolean {
  if (err instanceof NetworkError) return true;
  const e = asErrorWithStatus(err);
  return !!e && typeof e.status === "number" && UNREACHABLE_STATUSES.has(e.status);
}

export function isServerError(err: unknown): boolean {
  const e = asErrorWithStatus(err);
  return !!e && typeof e.status === "number" && e.status >= 500 && e.status < 600;
}

const VALIDATION_COPY: Partial<Record<OpKind, string>> = {
  peak_added: "Couldn't add peak",
  peak_excluded: "Couldn't exclude peak",
  peak_unexcluded: "Couldn't restore peak",
  peak_removed: "Couldn't remove peak",
  speculative_created: "Couldn't create speculative index",
  speculative_deleted: "Couldn't delete speculative index",
  set_exposure_status: "Couldn't update exposure status",
  select_exposure: "Couldn't switch exposure",
  update_sample: "Couldn't update sample",
  add_label: "Couldn't add label",
  remove_label: "Couldn't remove label",
  restore_label_suggestion: "Couldn't put the label back to suggested",
  reanalyze_exposure: "Couldn't re-analyze exposure",
  delete_index: "Couldn't delete index",
  assignment_add: "Couldn't call phase",
  assignment_remove: "Couldn't remove phase",
  patch_experiment: "Couldn't save the experiment",
  patch_setup: "Couldn't save the setup",
};

/**
 * Build a user-facing toast message for a validation failure. Per-kind copy
 * lives in VALIDATION_COPY; add new kinds when introducing new mutators.
 */
export function buildValidationMessage(kind: OpKind, err: unknown): string {
  const head = VALIDATION_COPY[kind] ?? "Action failed";
  const detail =
    asErrorWithStatus(err)?.message ?? "see browser console for details";
  return `${head}: ${detail}`;
}

export function buildServerErrorMessage(kind: OpKind): string {
  const head = VALIDATION_COPY[kind] ?? "Action failed";
  return `${head}: the server hit an error, so the change was undone`;
}
