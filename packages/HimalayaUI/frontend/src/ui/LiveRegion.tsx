import { useEffect, useRef, useState } from "react";
import { setAnnounceImpl, type AnnounceOptions } from "../lib/announce";

export function LiveRegion(): JSX.Element {
  const [polite, setPolite] = useState("");
  const [assertive, setAssertive] = useState("");
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