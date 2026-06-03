import { useEffect, useRef, useState } from "react";
import { Button, IconButton, Menu } from "../ui";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ExportButtonProps {
  onCopy: () => void;
  onDownloadPng: () => void;
  onDownloadSvg: () => void;
  /** Copy unavailable (no clipboard / no PNG renderer / a render in flight). */
  copyDisabled?: boolean;
  /** PNG download unavailable (no PNG renderer / a render in flight). */
  pngDisabled?: boolean;
  /** A render is in flight (blocks SVG too). */
  pending?: boolean;
  /** Page-level gate (data not ready). ORs into every action. */
  disabled?: boolean;
  /** Fills aria-labels: "Copy {ariaContext} to clipboard". */
  ariaContext: string;
  /** PLACEMENT ONLY. */
  className?: string;
}

/**
 * The figure-export split button (mockup: series-builder rail-foot "Copy as
 * PNG"): a bordered group with a primary **Copy** action and a `▾` chevron that
 * opens a two-item download menu (PNG / SVG). Presentational — every side
 * effect is a prop (wire it with `useFigureExport`). Owns ONLY its menu-open
 * state + outside-pointerdown dismissal (the `Field` precedent; `Menu` owns
 * Escape + arrow-nav).
 */
export function ExportButton({
  onCopy,
  onDownloadPng,
  onDownloadSvg,
  copyDisabled = false,
  pngDisabled = false,
  pending = false,
  disabled = false,
  ariaContext,
  className,
}: ExportButtonProps): JSX.Element {
  const [open, setOpen] = useState(false);
  const wrapRef = useRef<HTMLSpanElement | null>(null);

  // Outside-pointerdown closes the menu (mirrors Field/Popover). Bound only
  // while open.
  useEffect(() => {
    if (!open) return;
    const onPointerDown = (e: PointerEvent): void => {
      if (wrapRef.current && !wrapRef.current.contains(e.target as Node)) {
        setOpen(false);
      }
    };
    document.addEventListener("pointerdown", onPointerDown);
    return () => document.removeEventListener("pointerdown", onPointerDown);
  }, [open]);

  const copyOff = disabled || copyDisabled;
  const pngOff = disabled || pngDisabled;
  const svgOff = disabled || pending;

  return (
    <span
      ref={wrapRef}
      data-testid="export-button"
      className={cx("relative inline-flex", className)}
    >
      <span className="inline-flex items-stretch border border-hair-strong rounded overflow-hidden">
        <Button
          variant="ghost"
          data-testid="export-copy"
          aria-label={`Copy ${ariaContext} to clipboard`}
          disabled={copyOff}
          onClick={onCopy}
        >
          Copy
        </Button>
        <span className="w-px bg-hair-strong" aria-hidden="true" />
        <IconButton
          label="Download formats"
          data-testid="export-menu-trigger"
          aria-haspopup="menu"
          aria-expanded={open}
          disabled={disabled || pending}
          onClick={() => setOpen((o) => !o)}
        >
          ▾
        </IconButton>
      </span>
      <Menu<"png" | "svg">
        open={open}
        options={[
          { value: "png", label: "Download as PNG", disabled: pngOff },
          { value: "svg", label: "Download as SVG", disabled: svgOff },
        ]}
        onSelect={(v) => {
          if (v === "png") onDownloadPng();
          else onDownloadSvg();
        }}
        onClose={() => setOpen(false)}
        aria-label="Download formats"
        className="right-0 top-full"
      />
    </span>
  );
}
