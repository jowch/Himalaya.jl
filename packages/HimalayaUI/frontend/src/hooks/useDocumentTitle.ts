import { useEffect } from "react";

/** App name suffix; the tab reads "<lead> · Himalaya". */
const BASE = "Himalaya";

/**
 * Set `document.title` for the lifetime of the calling component and restore the
 * previous title on unmount. Pass the page-specific lead (sample name, series
 * title, page label); it is suffixed with the app name. A nullish/blank lead
 * falls back to the bare app name (e.g. while data is still loading) so the tab
 * never shows a stale, empty, or "undefined" title.
 *
 * Resolves the "static document.title" re-score nit shared across the page
 * surfaces (BU-RESCORE3 N3, FO-RESCORE2 F14, the folio/builder flags): the tab
 * label now names what you are looking at, so multiple open tabs are
 * distinguishable and browser history reads meaningfully.
 */
export function useDocumentTitle(lead: string | null | undefined): void {
  useEffect(() => {
    const prev = document.title;
    document.title = lead && lead.trim() !== "" ? `${lead} · ${BASE}` : BASE;
    return () => {
      document.title = prev;
    };
  }, [lead]);
}
