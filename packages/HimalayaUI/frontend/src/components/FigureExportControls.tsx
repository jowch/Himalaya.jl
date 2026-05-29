// FigureExportControls.tsx — Copy + split-Download (PNG/SVG) buttons.
import { useEffect, useRef, useState } from "react";
import type { ExportSpec } from "../lib/figure-export/types";
import { buildExportPng, buildExportSvg, canExportPng } from "../lib/figure-export/renderer";
import { canCopyPngToClipboard, copyPngToClipboard } from "../lib/figure-export/clipboard";
import { downloadBlob } from "../lib/figure-export/download";
import { buildFilename } from "../lib/figure-export/filename";
import { showToast } from "../lib/toast";

export interface FigureExportControlsProps {
  /** Thunk evaluated at click time so it captures fresh state (xDomain,
   *  current peaks, etc.) rather than stale state at component mount. */
  spec: () => ExportSpec;
  /** Already-resolved, slugified filename stem. Component appends the date
   *  + extension internally. */
  filenameStem: string;
  /** Used in aria-labels: "Copy {ariaContext} to clipboard". */
  ariaContext: string;
  /** Parent disables both buttons when the underlying data is incomplete
   *  (e.g. `!traceQ.data || !peaksQ.data`). */
  disabled?: boolean;
}

export function FigureExportControls({
  spec, filenameStem, ariaContext, disabled = false,
}: FigureExportControlsProps): JSX.Element {
  const [pending, setPending] = useState(false);
  const [menuOpen, setMenuOpen] = useState(false);
  const triggerRef = useRef<HTMLButtonElement>(null);
  const panelRef = useRef<HTMLDivElement>(null);
  const canCopy = canCopyPngToClipboard();
  const canPng = canExportPng();

  // Esc + outside-click dismiss (mirrors ForksPopover pattern; spec §Popover).
  useEffect(() => {
    if (!menuOpen) return;
    const onKey = (e: KeyboardEvent): void => {
      if (e.key === "Escape") setMenuOpen(false);
    };
    document.addEventListener("keydown", onKey);
    return () => document.removeEventListener("keydown", onKey);
  }, [menuOpen]);

  useEffect(() => {
    if (!menuOpen) return;
    const onDown = (e: MouseEvent): void => {
      const t = e.target as Node | null;
      if (t === null) return;
      if (panelRef.current?.contains(t)) return;
      if (triggerRef.current?.contains(t)) return;
      setMenuOpen(false);
    };
    document.addEventListener("mousedown", onDown);
    return () => document.removeEventListener("mousedown", onDown);
  }, [menuOpen]);

  const onCopy = async (): Promise<void> => {
    if (disabled || !canCopy || !canPng || pending) return;
    setPending(true);
    try {
      const blob = await buildExportPng(spec());
      await copyPngToClipboard(blob);
      showToast("Copied figure to clipboard", "success");
    } catch (err) {
      showToast("Couldn't copy figure. Try Download instead.", "error");
      // eslint-disable-next-line no-console
      console.warn(err);
    } finally {
      setPending(false);
    }
  };

  const onDownloadPng = async (): Promise<void> => {
    if (disabled || !canPng || pending) return;
    setMenuOpen(false);
    setPending(true);
    try {
      const blob = await buildExportPng(spec());
      downloadBlob(blob, buildFilename(filenameStem, "png"));
    } catch (err) {
      showToast("Couldn't render figure for download.", "error");
      // eslint-disable-next-line no-console
      console.warn(err);
    } finally {
      setPending(false);
    }
  };

  const onDownloadSvg = (): void => {
    if (disabled || pending) return;
    setMenuOpen(false);
    setPending(true);
    try {
      const svg = buildExportSvg(spec());
      const xml = new XMLSerializer().serializeToString(svg);
      const blob = new Blob([xml], { type: "image/svg+xml" });
      downloadBlob(blob, buildFilename(filenameStem, "svg"));
    } catch (err) {
      showToast("Couldn't render figure for download.", "error");
      // eslint-disable-next-line no-console
      console.warn(err);
    } finally {
      setPending(false);
    }
  };

  // Copy needs both clipboard support AND the PNG renderer (which produces the
  // blob). Save-PNG (and the chevron's PNG menu row) needs the renderer; SVG
  // download stays accessible regardless because it doesn't go through
  // OffscreenCanvas.
  const copyDisabled = disabled || !canCopy || !canPng || pending;
  const pngDisabled = disabled || !canPng || pending;
  const downloadAnyDisabled = disabled || pending;

  const copyTitle = !canCopy
    ? "Clipboard requires HTTPS"
    : !canPng
      ? "Browser doesn't support PNG export"
      : `Copy ${ariaContext} to clipboard`;
  const pngTitle = !canPng
    ? "Browser doesn't support PNG export"
    : `Download ${ariaContext} as PNG`;

  return (
    <span
      data-testid="figure-export-controls"
      className="relative inline-flex items-center gap-1"
    >
      <button
        type="button"
        data-testid="figure-export-copy"
        aria-label={`Copy ${ariaContext} to clipboard`}
        title={copyTitle}
        disabled={copyDisabled}
        onClick={() => { void onCopy(); }}
        className="px-1.5 py-0.5 rounded text-xs text-ink-faint hover:text-ink
                   hover:bg-paper-sunk border border-transparent hover:border-hair-strong
                   disabled:opacity-40 disabled:cursor-default"
      >
        Copy
      </button>
      <span className="inline-flex items-stretch border border-hair-strong rounded overflow-hidden">
        <button
          type="button"
          data-testid="figure-export-download-png"
          aria-label={`Download ${ariaContext} as PNG`}
          title={pngTitle}
          disabled={pngDisabled}
          onClick={() => { void onDownloadPng(); }}
          className="px-1.5 py-0.5 text-xs text-ink-faint hover:text-ink hover:bg-paper-sunk
                     disabled:opacity-40 disabled:cursor-default"
        >
          Save
        </button>
        <span className="w-px bg-hair-strong" />
        <button
          ref={triggerRef}
          type="button"
          data-testid="figure-export-download-menu-trigger"
          aria-label="Other download formats"
          aria-haspopup="true"
          aria-expanded={menuOpen}
          disabled={downloadAnyDisabled}
          onClick={() => setMenuOpen((v) => !v)}
          className="px-1 text-xs text-ink-faint hover:text-ink hover:bg-paper-sunk
                     disabled:opacity-40 disabled:cursor-default"
        >
          ▾
        </button>
      </span>
      {menuOpen && (
        <div
          ref={panelRef}
          data-testid="figure-export-download-menu"
          className="absolute z-50 top-full mt-1 right-0 min-w-[140px] card border border-hair-strong
                     bg-plate shadow-lg p-1"
        >
          <button
            type="button"
            data-testid="figure-export-download-menu-png"
            disabled={!canPng}
            title={pngTitle}
            onClick={() => { void onDownloadPng(); }}
            className="block w-full text-left px-2 py-1 text-xs text-ink hover:bg-paper-sunk rounded
                       disabled:opacity-40 disabled:cursor-default"
          >
            Download as PNG
          </button>
          <button
            type="button"
            data-testid="figure-export-download-menu-svg"
            onClick={onDownloadSvg}
            className="block w-full text-left px-2 py-1 text-xs text-ink hover:bg-paper-sunk rounded"
          >
            Download as SVG
          </button>
        </div>
      )}
    </span>
  );
}
