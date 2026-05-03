export type {
  OpKind,
  OpPayload,
  RollbackContext,
  PendingDeferred,
  SseEvent,
  Mutator,
} from "./types";

export {
  pendingDeferreds,
  makeDeferred,
  getDeferred,
  clearDeferred,
} from "./deferred";

export { handleRemoteEvent } from "./replayCoordinator";

export {
  attachPersistence,
  rehydrate,
  STORAGE_KEY,
  SCHEMA_VERSION,
} from "./persistence";
export type { RehydrateResult } from "./persistence";

export {
  useExposureHasPendingPeakOps,
  useQueueOpStatus,
} from "./hooks";

export {
  useQueueMutation,
} from "./useQueueMutation";
export type { UseQueueMutationResult } from "./useQueueMutation";
export {
  isValidationError,
  isInfrastructureError,
  buildValidationMessage,
} from "./errors";
