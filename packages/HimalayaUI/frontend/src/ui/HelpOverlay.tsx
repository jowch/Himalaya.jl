import { Fragment, useEffect, useRef } from "react";
import { IconButton } from "./IconButton";
import { Kbd } from "./Kbd";
import { Overlay } from "./Overlay";

export type HelpRow = { keys: string[]; gesture?: never; what: string } | { gesture: string; keys?: never; what: string };
export type HelpSection = { title: string; rows: HelpRow[] };

const TITLE = "Keyboard shortcuts";

export function HelpOverlay({ sections, onClose }: { sections: readonly HelpSection[]; onClose: () => void }): JSX.Element {
  const body = useRef<HTMLDivElement | null>(null);
  useEffect(() => { body.current?.focus(); }, []);
  return (
    <Overlay label={TITLE} onClose={onClose} testId="workspace-help">
      <div className="flex max-h-[calc(88vh-1.5rem)] flex-col">
        <header data-testid="help-header" className="flex shrink-0 items-center justify-between gap-4 border-b border-hair py-2 pl-5 pr-3">
          <h2 className="text-label text-ink-soft">{TITLE}</h2>
          <IconButton label="Close" dismiss tabIndex={-1} onClick={onClose} />
        </header>
        <div
          ref={body}
          tabIndex={0}
          role="region"
          aria-label={TITLE}
          data-testid="help-body"
          className="grid min-h-0 grid-cols-[auto_1fr] items-center gap-x-4 gap-y-2 overflow-y-auto overscroll-contain px-5 py-4 text-body outline-none"
        >
          {sections.map((s) => (
            <Fragment key={s.title}>
              <h3 className="col-span-2 mt-2 text-label text-ink-soft first:mt-0">{s.title}</h3>
              {s.rows.map((r) => (
                <Fragment key={r.what}>
                  <span className="inline-flex items-center gap-1">
                    {r.keys ? r.keys.map((k) => <Kbd key={k}>{k}</Kbd>) : <span className="text-ink-soft">{r.gesture}</span>}
                  </span>
                  <span>{r.what}</span>
                </Fragment>
              ))}
            </Fragment>
          ))}
        </div>
      </div>
    </Overlay>
  );
}
