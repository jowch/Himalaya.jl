import { useEffect, useRef, useState } from "react";
import { Link } from "react-router-dom";
import { useFocusTrap } from "../hooks/useFocusTrap";

/** A child fork, pre-resolved to its review-page path by the caller. */
export interface ForkLink {
  id: number;
  title: string;
  href: string;
}

interface Props {
  groupingControl: React.ReactNode;
  annotationControl: React.ReactNode;
  /** Child forks of this comparison, each with a route-ready `href`. */
  forksList: ForkLink[];
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
          <div role="menu" className="absolute z-50 mt-1 right-0 min-w-[220px] card border border-border bg-bg-elevated shadow-lg p-1">
            {/* Forks sub-section (Compare UX C-17): the child forks of this
                comparison, folded in from the old standalone ForksPopover.
                Always shown inside the open menu; lists each fork as a link
                to its review page, or a "No forks yet" empty state. */}
            <div data-testid="compare-toolbar-forks" className="px-2 py-1">
              <div className="text-xs text-fg-muted">
                Forks ({p.forksList.length})
              </div>
              {p.forksList.length === 0 ? (
                <div className="text-xs text-fg-muted italic mt-0.5">
                  No forks yet
                </div>
              ) : (
                <ul className="flex flex-col mt-0.5 max-h-[200px] overflow-y-auto">
                  {p.forksList.map((f) => (
                    <li key={f.id}>
                      <Link
                        to={f.href}
                        data-testid="compare-toolbar-fork-link"
                        onClick={() => setOpen(false)}
                        className="block px-1 py-0.5 rounded text-sm text-fg
                                   hover:bg-bg-hover hover:text-accent truncate"
                      >
                        {f.title || `Comparison #${f.id}`}
                      </Link>
                    </li>
                  ))}
                </ul>
              )}
            </div>
            <div className="my-1 border-t border-border" />
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
