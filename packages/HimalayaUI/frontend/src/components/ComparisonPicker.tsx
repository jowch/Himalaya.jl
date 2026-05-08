import { useEffect, useRef, useState } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
import { useFocusTrap } from "../hooks/useFocusTrap";
import { ComparisonPickerBody, type Pick } from "./ComparisonPickerBody";

interface Props {
  isOpen: boolean;
  onClose: () => void;
  experimentId: number | undefined;
}

export function ComparisonPicker({
  isOpen, onClose, experimentId,
}: Props): JSX.Element | null {
  const dialogRef = useRef<HTMLDivElement>(null);
  const inputRef  = useRef<HTMLInputElement>(null);
  useFocusTrap(dialogRef, isOpen);

  const qc = useQueryClient();
  const draft = useAppState((s) => s.activeDraft);
  const addMember = useAppState((s) => s.addMember);

  const [picks, setPicks] = useState<Pick[]>([]);
  useEffect(() => {
    if (isOpen) {
      setPicks([]);
      inputRef.current?.focus();
    }
  }, [isOpen]);

  if (!isOpen) return null;

  const alreadyAddedExposureIds = new Set(
    (draft?.members ?? [])
      .map((m) => m.exposure_id)
      .filter((id): id is number => id !== null),
  );

  const onAddSelected = (): void => {
    for (const p of picks) addMember(p.exposure_id, qc);
    onClose();
  };

  return (
    <div
      data-testid="comparison-picker-overlay"
      className="fixed inset-0 z-50 flex items-start justify-center pt-[10vh]
                 bg-[oklch(0.05_0_0/0.65)] backdrop-blur-sm anim-pal-in"
      role="presentation"
      onClick={(e) => { if (e.target === e.currentTarget) onClose(); }}
    >
      <div
        ref={dialogRef}
        data-testid="comparison-picker"
        role="dialog"
        aria-modal="true"
        aria-labelledby="comparison-picker-title"
        onKeyDown={(e) => { if (e.key === "Escape") { e.preventDefault(); onClose(); } }}
        className="w-[min(720px,calc(100vw-48px))] max-h-[78vh]
                   bg-bg-elevated border border-border rounded-xl shadow-2xl
                   flex flex-col overflow-hidden anim-pal-scale"
      >
        <div className="flex items-center justify-between px-4 py-3 border-b border-border">
          <h2 id="comparison-picker-title" className="text-base font-medium text-fg">
            Add traces
          </h2>
          <button
            type="button"
            data-testid="comparison-picker-close"
            onClick={onClose}
            className="text-fg-muted hover:text-fg text-sm px-2 py-1"
            aria-label="Close picker"
          >
            esc
          </button>
        </div>

        <ComparisonPickerBody
          experimentId={experimentId}
          picks={picks}
          onPicksChange={setPicks}
          alreadyAddedExposureIds={alreadyAddedExposureIds}
          searchInputRef={inputRef}
        />

        <div className="flex items-center gap-2 px-4 py-3 border-t border-border">
          <span className="text-xs text-fg-dim flex-1">{picks.length} selected</span>
          <button type="button"
            data-testid="comparison-picker-cancel"
            onClick={onClose}
            className="px-3 py-1 rounded border border-border text-sm text-fg-muted">
            Cancel
          </button>
          <button type="button"
            data-testid="comparison-picker-add"
            onClick={onAddSelected}
            disabled={picks.length === 0}
            className="px-3 py-1 rounded bg-accent text-bg text-sm font-medium disabled:opacity-50">
            {picks.length === 0 ? "Add selected" : `Add ${picks.length} selected`}
          </button>
        </div>
      </div>
    </div>
  );
}
