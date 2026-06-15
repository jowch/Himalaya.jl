// SR-only status announcer — the quiet sibling of lib/toast.ts.
//
// `showToast` surfaces consequential events visibly. `announce` is for
// frequent, low-consequence, immediately-visible state changes (peak add/remove,
// frame flips) where a visible toast on every action would be spam, but a
// screen-reader user still needs WCAG 4.1.3 status feedback. The message is
// written into a visually-hidden aria-live region (see print/ui/LiveRegion.tsx).

export interface AnnounceOptions {
  /** true → assertive region (interrupts); default false → polite. */
  assertive?: boolean;
}

let activeImpl: (msg: string, opts?: AnnounceOptions) => void = (msg) => {
  // eslint-disable-next-line no-console
  console.debug(`[announce] ${msg}`);
};

/**
 * Announce a status message to screen readers via the mounted LiveRegion.
 * No-op (console.debug) when no region is mounted.
 */
export function announce(msg: string, opts?: AnnounceOptions): void {
  activeImpl(msg, opts);
}

/**
 * Replace the announce implementation. LiveRegion calls this on mount.
 * Passing `null` restores the console.debug fallback.
 */
export function setAnnounceImpl(
  impl: ((msg: string, opts?: AnnounceOptions) => void) | null,
): void {
  if (impl === null) {
    activeImpl = (msg) => {
      // eslint-disable-next-line no-console
      console.debug(`[announce] ${msg}`);
    };
  } else {
    activeImpl = impl;
  }
}
