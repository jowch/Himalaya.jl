// Spec §4.3 — push/replace mode for the next useUrlFromState emit.
// Hoisted into its own module to avoid a state.ts ↔ useUrlFromState cycle
// and to keep test isolation simpler (consume resets to push on read).

let nextEmitMode: "push" | "replace" = "push";

export function emitReplaceNext(): void {
  nextEmitMode = "replace";
}

export function consumeEmitMode(): "push" | "replace" {
  const mode = nextEmitMode;
  nextEmitMode = "push";
  return mode;
}

// For tests: reset between cases.
export function _resetEmitMode(): void {
  nextEmitMode = "push";
}
