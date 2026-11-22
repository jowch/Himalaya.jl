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

export function isValidationError(err: unknown): boolean {
  const e = asErrorWithStatus(err);
  if (!e || typeof e.status !== "number") return false;
  return e.status >= 400 && e.status < 500;
}
export function is404Error(err: unknown): boolean {
  const e = asErrorWithStatus(err);
  return !!e && e.status === 404;
}
export function isInfrastructureError(err: unknown): boolean {
  const e = asErrorWithStatus(err);
  if (!e) return false;
  if (typeof e.status === "number") {
    return e.status >= 500 && e.status < 600;
  }
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
  set_exposure_status: "Couldn't drop or restore frame",
  update_sample: "Couldn't update sample",
  add_label: "Couldn't add label",
  remove_label: "Couldn't remove label",
  restore_label_suggestion: "Couldn't put the label back to suggested",
  reanalyze_exposure: "Couldn't re-analyze frame",
  delete_index: "Couldn't delete index",
  assignment_add: "Couldn't call phase",
  assignment_remove: "Couldn't remove phase",
  patch_experiment: "Couldn't save the experiment",
  patch_setup: "Couldn't save the setup",
};
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
