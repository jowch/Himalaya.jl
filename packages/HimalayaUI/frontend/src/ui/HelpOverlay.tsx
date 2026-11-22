import { Fragment } from "react";
import { Kbd } from "./Kbd";
import { Overlay } from "./Overlay";

export type HelpRow = { keys: string[]; gesture?: never; what: string } | { gesture: string; keys?: never; what: string };
export type HelpSection = { title: string; rows: HelpRow[] };

export function HelpOverlay({ sections, onClose }: { sections: readonly HelpSection[]; onClose: () => void }): JSX.Element {
  return (
    <Overlay label="Keyboard shortcuts" onClose={onClose} testId="workspace-help">
      <div className="grid grid-cols-[auto_1fr] items-center gap-x-4 gap-y-2 px-5 py-4 text-body">
        {sections.map((s) => (
          <Fragment key={s.title}>
            <h2 className="col-span-2 mt-2 text-label text-ink-soft first:mt-0">{s.title}</h2>
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
    </Overlay>
  );
}
