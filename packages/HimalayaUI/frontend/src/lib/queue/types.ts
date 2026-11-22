import type { QueryClient } from "@tanstack/react-query";
import type { Assignment, SampleCall } from "../../api";

export type OpKind =
  | "peak_added" | "peak_excluded" | "peak_unexcluded" | "peak_removed"
  | "speculative_created" | "speculative_deleted"
  | "set_exposure_status"
  | "update_sample"
  | "reanalyze_exposure"
  | "delete_index"
  | "assignment_add" | "assignment_remove" | "sample_call_set_state"
  | "custom_index_commit"
  | "bridge_samples" | "unbridge_sample"
  | "patch_compares"
  | "add_label" | "remove_label" | "restore_label_suggestion"
  | "set_tag_color"
  | "patch_experiment"
  | "patch_setup";
export interface OpPayload<T = unknown> {
  kind: OpKind;
  payload: T;
  clientOpId: string;
  exposureId?: number;
}
export interface RollbackContext {
  restore: () => void;
}
export interface PendingDeferred<T> {
  promise: Promise<T>;
  resolve: (value: T) => void;
  reject: (reason: unknown) => void;
  controller?: AbortController;
}
export interface CurationPostState {
  analysis_inputs_hash: string;
  indices: unknown[];
  assignment?: Pick<Assignment, "state" | "members">;
  assignment_dropped?: string[];
}
export interface AssignmentPostState {
  assignment: Pick<Assignment, "state" | "members">;
}
export interface SampleCallPostState {
  call: Pick<SampleCall, "state" | "members">;
}
export interface SseEvent {
  id: number;
  kind: string;
  entity_type: string;
  entity_id: number;
  actor?: string | null;
  client_id?: string | null;
  client_op_id?: string | null;
  ts?: string;
  payload?: unknown;
  post_state?: CurationPostState | AssignmentPostState | SampleCallPostState;
}
export interface QueueResponseMeta {
  event_id: number;
  client_op_id: string | null | undefined;
  analysis_inputs_hash: string | undefined;
  view_row_id?: number;
}
export type FlatPayload<TInput, TScope> = Omit<OpPayload<TInput>, "payload"> & TScope & TInput;
export interface Mutator<TInput, TScope, TResponse> {
  kind: OpKind;
  onMutate: (payload: FlatPayload<TInput, TScope>, qc: QueryClient) => RollbackContext;
  request: (payload: FlatPayload<TInput, TScope>, signal: AbortSignal) => Promise<TResponse>;
  onSuccess: (payload: FlatPayload<TInput, TScope>, response: TResponse, qc: QueryClient) => void;
  synthesizeFromSse?: (remote: SseEvent, base: QueueResponseMeta) => TResponse | undefined;
  treats404AsSuccess?: boolean;
}