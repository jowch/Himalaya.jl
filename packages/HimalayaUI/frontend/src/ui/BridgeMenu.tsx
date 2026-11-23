import { useEffect, useRef, useState } from "react";
import { Tooltip } from "./Tooltip";

export interface BridgeMate {
  id: number;
  title: string;
}

interface Props {
  mates: readonly BridgeMate[];
  editable: boolean;
  onGo: (id: number) => void;
  onRemove: (id: number) => void;
}

export function BridgeMenu({ mates, editable, onGo, onRemove }: Props): JSX.Element {
  const [open, setOpen] = useState(false);
  const root = useRef<HTMLSpanElement | null>(null);
  const opener = useRef<HTMLButtonElement | null>(null);
  const name = `Same as ${mates.map((m) => m.title).join(", ")}`;

  useEffect(() => {
    if (!open) return;
    root.current?.querySelector<HTMLElement>("[role='dialog'] button")?.focus();
    const onDown = (ev: MouseEvent): void => {
      if (root.current && ev.target instanceof Node && !root.current.contains(ev.target)) setOpen(false);
    };
    document.addEventListener("mousedown", onDown, true);
    return () => document.removeEventListener("mousedown", onDown, true);
  }, [open]);

  const close = (): void => {
    setOpen(false);
    opener.current?.focus();
  };

  return (
    <span ref={root} className="relative inline-flex shrink-0 self-center">
      <Tooltip content={name}>
        <button
          ref={opener}
          type="button"
          aria-label={name}
          aria-expanded={open}
          data-same-sample
          onClick={() => setOpen((o) => !o)}
          className="inline-flex items-center rounded px-1 text-base leading-6 text-ink-faint tabular-nums hover:text-ink focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"
        >
          <span aria-hidden="true">⇄{mates.length}</span>
        </button>
      </Tooltip>
      {open && (
        <div
          role="dialog"
          aria-label="Same sample as"
          data-keys-trap
          data-testid="workspace-same-samples"
          className="absolute left-0 top-full z-20 mt-1 flex min-w-56 max-w-80 flex-col rounded-sm border border-hair bg-plate py-1 shadow-lg"
          onKeyDown={(e) => {
            if (e.key !== "Escape") return;
            e.preventDefault();
            e.stopPropagation();
            close();
          }}
        >
          <p className="px-2.5 pb-1 pt-0.5 text-caption text-ink-faint">Same sample as</p>
          {mates.map((m) => (
            <div key={m.id} className="flex items-center gap-1 px-1">
              <button
                type="button"
                aria-label={`Go to ${m.title}`}
                onClick={() => {
                  setOpen(false);
                  onGo(m.id);
                }}
                className="min-w-0 flex-1 truncate rounded px-1.5 py-1 text-left text-base text-ink hover:bg-hair/40 focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent"
              >
                {m.title}
              </button>
              {editable && (
                <Tooltip content={`Not the same sample as ${m.title}`}>
                  <button
                    type="button"
                    aria-label={`Not the same sample as ${m.title}`}
                    onClick={() => onRemove(m.id)}
                    className="shrink-0 rounded px-1.5 py-1 text-ink-faint hover:text-ink focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent"
                  >
                    ×
                  </button>
                </Tooltip>
              )}
            </div>
          ))}
        </div>
      )}
    </span>
  );
}
