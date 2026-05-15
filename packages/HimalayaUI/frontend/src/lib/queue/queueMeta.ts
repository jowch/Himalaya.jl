import type { QueueResponseMeta } from "./types";

/**
 * Split a queue mutator response into its framework-owned metadata and
 * the mutator-specific payload. Use in `onSuccess` so the per-mutator
 * cache write doesn't have to spell out every framework field by hand —
 * and a new field added to `QueueResponseMeta` later flows through
 * automatically.
 *
 *   const { meta, payload } = stripQueueMetadata(response);
 *   qc.setQueryData(...payload row...)
 *   writeExposureHash(qc, exposureId, meta.analysis_inputs_hash);
 */
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
