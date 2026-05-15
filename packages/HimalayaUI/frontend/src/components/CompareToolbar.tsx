import { useEffect, useRef, useState } from "react";
import { useFocusTrap } from "../hooks/useFocusTrap";

interface Props {
  groupingControl: React.ReactNode;
  annotationControl: React.ReactNode;
  forksCount: number;
  onCopyLink: () => void;
  onDelete: () => void;
  onDiscardChanges: (() => void) | null;
  onFork: () => void;
  exportControl: React.ReactNode;
  saveControl: React.ReactNode;
}

export function CompareToolbar(p: Props): JSX.Element {
  const [open, setOpen] = useState(false);
  const menuRef = useRef<HTMLDivElement>(null);
  // Trap Tab/Shift+Tab inside the open menu — canonical pattern from
  // NavModal / OnboardingFlow. Restores focus to the trigger on close.
  useFocusTrap(menuRef, open);

  useEffect(() => {
    if (!open) return;
    const onDoc = (e: MouseEvent) => {
      if (!menuRef.current?.contains(e.target as Node)) setOpen(false);
    };
    const onEsc = (e: KeyboardEvent) => { if (e.key === "Escape") setOpen(false); };
    document.addEventListener("mousedown", onDoc);
    document.addEventListener("keydown", onEsc);
    return () => {
      document.removeEventListener("mousedown", onDoc);
      document.removeEventListener("keydown", onEsc);
    };
  }, [open]);

  return (
    <div data-testid="compare-toolbar" className="flex items-center gap-2 flex-wrap">
      {p.groupingControl}
      {p.annotationControl}
      <div className="relative" ref={menuRef}>
        <button
          type="button"
          data-testid="compare-toolbar-more"
          data-interactable="button"
          aria-expanded={open}
          aria-haspopup="menu"
          onClick={() => setOpen((v) => !v)}
          className="px-2 py-0.5 rounded text-xs text-fg-dim hover:text-fg hover:bg-bg-hover border border-transparent hover:border-border"
        >
          ⋯ more
        </button>
        {open && (
          <div role="menu" className="absolute z-50 mt-1 right-0 min-w-[200px] card border border-border bg-bg-elevated shadow-lg p-1">
            <button
              type="button"
              data-testid="compare-toolbar-forks"
              data-interactable="button"
              onClick={() => { setOpen(false); }}
              className="w-full text-left px-2 py-1 hover:bg-bg-hover text-sm"
            >
              Forks ({p.forksCount})
            </button>
            <button type="button" onClick={() => { p.onCopyLink(); setOpen(false); }} className="w-full text-left px-2 py-1 hover:bg-bg-hover text-sm">Copy link</button>
            <button type="button" onClick={() => { p.onFork(); setOpen(false); }} className="w-full text-left px-2 py-1 hover:bg-bg-hover text-sm">Fork</button>
            <button type="button" onClick={() => { p.onDelete(); setOpen(false); }} className="w-full text-left px-2 py-1 hover:bg-bg-hover text-sm text-danger">Delete</button>
            {p.onDiscardChanges !== null && (
              <button type="button" onClick={() => { p.onDiscardChanges!(); setOpen(false); }} className="w-full text-left px-2 py-1 hover:bg-bg-hover text-sm">Discard changes</button>
            )}
          </div>
        )}
      </div>
      <span className="ml-auto inline-flex items-center gap-2">
        {p.exportControl}
        {p.saveControl}
      </span>
    </div>
  );
}
