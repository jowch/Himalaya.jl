import { useEffect, useRef, useState } from "react";
import { setAnnounceImpl, type AnnounceOptions } from "../../lib/announce";

/**
 * LiveRegion — the global SR-only status surface.
 *
 * Mounts two visually-hidden aria-live nodes (polite + assertive) and wires
 * `setAnnounceImpl` into them on mount, restoring the console.debug fallback on
 * unmount (mirrors ToastContainer's setToastImpl lifecycle). Consumers call
 * `announce(msg, { assertive })` from lib/announce; the message lands in the
 * matching node so screen readers speak it (WCAG 4.1.3).
 *
 * Render once at the root of the app (App.tsx), next to ToastContainer.
 *
 * Appearance note: the only style here is the `sr-only` utility — the project's
 * visually-hidden mechanism. There is no visible chrome by design.
 */
export function LiveRegion(): JSX.Element {
  const [polite, setPolite] = useState("");
  const [assertive, setAssertive] = useState("");
  // A toggled trailing space so two identical consecutive announcements still
  // change the node's textContent — screen readers only re-speak on a real
  // text mutation. The space is invisible and reads as nothing.
  //
  // PER-REGION flips: the polite and assertive nodes each own their own toggle.
  // A single shared toggle would let an interleaved write to the OTHER region
  // flip it back, so an identical repeat could land on the exact prior text and
  // never mutate, and the SR would silently skip it.
  const flipPolite = useRef(false);
  const flipAssertive = useRef(false);

  useEffect(() => {
    setAnnounceImpl((msg: string, opts?: AnnounceOptions) => {
      if (opts?.assertive) {
        flipAssertive.current = !flipAssertive.current;
        setAssertive(flipAssertive.current ? msg : `${msg} `);
      } else {
        flipPolite.current = !flipPolite.current;
        setPolite(flipPolite.current ? msg : `${msg} `);
      }
    });
    return () => setAnnounceImpl(null);
  }, []);

  return (
    <>
      <div
        data-testid="live-region-polite"
        aria-live="polite"
        aria-atomic="true"
        className="sr-only"
      >
        {polite}
      </div>
      <div
        data-testid="live-region-assertive"
        aria-live="assertive"
        aria-atomic="true"
        className="sr-only"
      >
        {assertive}
      </div>
    </>
  );
}
