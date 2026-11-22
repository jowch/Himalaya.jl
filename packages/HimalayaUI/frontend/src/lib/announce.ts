export interface AnnounceOptions {
  assertive?: boolean;
}

let activeImpl: (msg: string, opts?: AnnounceOptions) => void = (msg) => {
  // eslint-disable-next-line no-console -- dev fallback when no live region is mounted
  console.debug(`[announce] ${msg}`);
};

export function announce(msg: string, opts?: AnnounceOptions): void {
  activeImpl(msg, opts);
}
export function setAnnounceImpl(
  impl: ((msg: string, opts?: AnnounceOptions) => void) | null,
): void {
  if (impl === null) {
    activeImpl = (msg) => {
      // eslint-disable-next-line no-console -- dev fallback when no live region is mounted
      console.debug(`[announce] ${msg}`);
    };
  } else {
    activeImpl = impl;
  }
}