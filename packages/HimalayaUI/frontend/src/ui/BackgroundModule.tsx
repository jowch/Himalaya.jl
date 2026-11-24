import type { ReactNode } from "react";
import type { CompareScale } from "../workspace/compare";
import { SCALE_NOTES } from "./CompareBar";
import { Kbd } from "./Kbd";
import { Section } from "./Section";
import { Tooltip } from "./Tooltip";

export type BackgroundView =
  | { kind: "loading" }
  | { kind: "buffer" }
  | { kind: "off"; reason: "chosen" | "no-buffer"; canUseRunBuffer: boolean }
  | {
      kind: "set";
      name: string;
      source: "default" | "user";
      trace: "loading" | "ready" | "missing" | "hidden";
      scale: CompareScale | undefined;
      canUseRunBuffer: boolean;
    };

interface Props {
  view: BackgroundView;
  editable: boolean;
  picker?: ReactNode;
  onChange: () => void;
  onNone: () => void;
  onDefault: () => void;
}

const NOTE = "h-5 truncate whitespace-nowrap text-caption text-ink-soft";
const ACTION = "shrink-0 whitespace-nowrap rounded-sm border border-hair-strong bg-plate px-2 py-0.5 text-caption text-ink hover:bg-paper-sunk disabled:text-ink-faint disabled:hover:bg-plate";

const MATCHED_AT_HIGH_Q = "matched at high q";

const SCALE_LABELS: Record<CompareScale["source"], string | undefined> = {
  normalised: undefined,
  monitor: "scaled",
  matched: MATCHED_AT_HIGH_Q,
  unknown_normalisation: MATCHED_AT_HIGH_Q,
  unknown_mode: MATCHED_AT_HIGH_Q,
  capillary: MATCHED_AT_HIGH_Q,
  other_setup: "other setup, matched",
  beam_off: "beam off, matched",
  i2_collapsed: "I₂ dropped, matched",
};

function offLabel(view: BackgroundView): string | undefined {
  if (view.kind === "buffer") return "Buffer";
  if (view.kind === "off") return view.reason === "chosen" ? "None" : "No buffer in this run";
  return undefined;
}

function nameLine(view: BackgroundView): ReactNode {
  if (view.kind === "set") {
    return (
      <>
        <Tooltip content={view.name} overflowOnly><span className="min-w-0 truncate text-body text-ink">{view.name}</span></Tooltip>
        <span className="shrink-0 text-caption text-ink-faint">{view.source === "default" ? "run buffer" : "chosen"}</span>
      </>
    );
  }
  if (view.kind === "loading") return <span className="text-body text-ink-faint">Loading</span>;
  return <span className="min-w-0 truncate text-body text-ink-soft">{offLabel(view)}</span>;
}

function traceNote(view: BackgroundView): { text: string; description?: string } | null {
  if (view.kind !== "set") return null;
  if (view.trace === "loading") return { text: "loading" };
  if (view.trace === "missing") return { text: "no trace" };
  if (view.trace === "hidden") return { text: "setup hidden" };
  if (!view.scale) return null;
  const label = SCALE_LABELS[view.scale.source];
  if (!label) return null;
  return { text: `${label} ×${view.scale.scale.toPrecision(3)}`, description: SCALE_NOTES[view.scale.source] };
}

export function BackgroundModule({ view, editable, picker, onChange, onNone, onDefault }: Props): JSX.Element {
  const loading = view.kind === "loading";
  const idle = loading || !editable;
  const hasControls = view.kind !== "buffer";
  const offerDefault = (view.kind === "set" || view.kind === "off") && view.canUseRunBuffer;
  const canTurnOff = view.kind !== "off";
  const note = traceNote(view);
  return (
    <Section title="Background" className="relative shrink-0">
      <div className="flex h-6 items-baseline justify-between gap-2" data-testid="background-name">{nameLine(view)}</div>
      <p className={NOTE} data-testid="background-note" {...(note?.description ? { "aria-description": note.description } : {})}>{note?.text}</p>
      <div className="mt-1 flex h-7 items-center gap-1.5 overflow-hidden" data-testid="background-controls">
        {hasControls && (
          <>
            <button type="button" className={ACTION} disabled={idle} onClick={onChange}>
              Change <Kbd>B</Kbd>
            </button>
            {canTurnOff && (
              <button type="button" className={ACTION} disabled={idle} onClick={onNone}>None</button>
            )}
            {offerDefault && (
              <Tooltip content="Use the run's buffer as the background"><button type="button" className={ACTION} disabled={idle} onClick={onDefault}>Run buffer</button></Tooltip>
            )}
          </>
        )}
      </div>
      {picker}
    </Section>
  );
}
