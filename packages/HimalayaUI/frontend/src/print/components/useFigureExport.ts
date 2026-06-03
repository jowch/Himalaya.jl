import { useCallback, useRef, useState } from "react";
import type { ExportSpec } from "../../lib/figure-export/types";
import { buildExportPng, buildExportSvg, canExportPng } from "../../lib/figure-export/renderer";
import { canCopyPngToClipboard, copyPngToClipboard } from "../../lib/figure-export/clipboard";
import { downloadBlob } from "../../lib/figure-export/download";
import { buildFilename } from "../../lib/figure-export/filename";
import { showToast } from "../../lib/toast";

export interface UseFigureExport {
  onCopy: () => void;
  onDownloadPng: () => void;
  onDownloadSvg: () => void;
  copyDisabled: boolean;
  pngDisabled: boolean;
  pending: boolean;
}

/**
 * The figure-export logic hook — the ONLY consumer of `lib/figure-export/*`.
 * `spec` is a thunk evaluated at click time so it captures fresh state. The
 * clean scientific styling lives in the engine adapter that builds the
 * `ExportSpec`; this hook is style-agnostic.
 */
export function useFigureExport(
  spec: () => ExportSpec,
  filenameStem: string,
  ariaContext: string,
): UseFigureExport {
  const [pending, setPending] = useState(false);
  // A synchronous in-flight ref: `pending` state flips async, so two clicks in
  // one tick would both pass a state-only guard. The ref closes that window;
  // the state still drives the UI.
  const inFlight = useRef(false);

  const run = useCallback(
    async (work: () => Promise<void> | void, failMsg: string): Promise<void> => {
      if (inFlight.current) return;
      inFlight.current = true;
      setPending(true);
      try {
        await work();
      } catch (err) {
        showToast(failMsg, "error");
        // eslint-disable-next-line no-console
        console.warn(err);
      } finally {
        inFlight.current = false;
        setPending(false);
      }
    },
    [],
  );

  const onCopy = useCallback(() => {
    void run(async () => {
      const blob = await buildExportPng(spec());
      await copyPngToClipboard(blob);
      showToast(`Copied ${ariaContext} to clipboard`, "success");
    }, "Couldn't copy figure. Try Download instead.");
  }, [run, spec, ariaContext]);

  const onDownloadPng = useCallback(() => {
    void run(async () => {
      const blob = await buildExportPng(spec());
      downloadBlob(blob, buildFilename(filenameStem, "png"));
    }, "Couldn't render figure for download.");
  }, [run, spec, filenameStem]);

  const onDownloadSvg = useCallback(() => {
    void run(() => {
      const svg = buildExportSvg(spec());
      const xml = new XMLSerializer().serializeToString(svg);
      const blob = new Blob([xml], { type: "image/svg+xml" });
      downloadBlob(blob, buildFilename(filenameStem, "svg"));
    }, "Couldn't render figure for download.");
  }, [run, spec, filenameStem]);

  const canCopy = canCopyPngToClipboard();
  const canPng = canExportPng();
  return {
    onCopy,
    onDownloadPng,
    onDownloadSvg,
    copyDisabled: !canCopy || !canPng || pending,
    pngDisabled: !canPng || pending,
    pending,
  };
}
