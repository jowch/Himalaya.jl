import type { QueueResponseMeta } from "./types";

export function stripQueueMetadata<T extends object>(
  response: T,
): { meta: QueueResponseMeta; payload: Omit<T, keyof QueueResponseMeta> } {
  const {
    event_id,
    view_row_id,
    analysis_inputs_hash,
    client_op_id,
    ...payload
  } = response as T & QueueResponseMeta;
  const meta: QueueResponseMeta = { event_id, analysis_inputs_hash, client_op_id };
  if (view_row_id !== undefined) {
    meta.view_row_id = view_row_id;
  }
  return {
    meta,
    payload: payload as Omit<T, keyof QueueResponseMeta>,
  };
}
