import { useCallback, useRef, useState } from "react";
import {
  canRasterizePng,
  svgStringToBlob,
  svgStringToPng,
} from "../../lib/figure-export/raster";
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
 * The figure-export logic hook. `renderSvg` is a thunk evaluated at click time
 * so it captures fresh state; it returns the export figure as a standalone SVG
 * markup string (the greenfield CleanFigure builder). This hook is style- and
 * engine-agnostic: it blobs the string for SVG, and rasterizes it for PNG.
 */
export function useFigureExport(
  renderSvg: () => string,
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
      const blob = await svgStringToPng(renderSvg());
      await copyPngToClipboard(blob);
      showToast(`Copied ${ariaContext} to clipboard`, "success");
    }, "Couldn't copy figure. Try Download instead.");
  }, [run, renderSvg, ariaContext]);

  const onDownloadPng = useCallback(() => {
    void run(async () => {
      const blob = await svgStringToPng(renderSvg());
      downloadBlob(blob, buildFilename(filenameStem, "png"));
    }, "Couldn't render figure for download.");
  }, [run, renderSvg, filenameStem]);

  const onDownloadSvg = useCallback(() => {
    void run(() => {
      const blob = svgStringToBlob(renderSvg());
      downloadBlob(blob, buildFilename(filenameStem, "svg"));
    }, "Couldn't render figure for download.");
  }, [run, renderSvg, filenameStem]);

  const canCopy = canCopyPngToClipboard();
  const canPng = canRasterizePng();
  return {
    onCopy,
    onDownloadPng,
    onDownloadSvg,
    copyDisabled: !canCopy || !canPng || pending,
    pngDisabled: !canPng || pending,
    pending,
  };
}
