import type { ReactNode } from "react";
import type { CompareScale } from "../workspace/compare";
import { SCALE_NOTES } from "./CompareBar";
import { Kbd } from "./Kbd";
import { Section } from "./Section";

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

const LINE = "min-h-5 text-caption";
const ACTION = "rounded-sm border border-hair-strong bg-plate px-2 py-0.5 text-caption text-ink hover:bg-paper-sunk disabled:text-ink-faint disabled:hover:bg-plate";

function nameLine(view: BackgroundView): ReactNode {
  if (view.kind === "set") {
    return (
      <>
        <span className="min-w-0 truncate text-body text-ink">{view.name}</span>
        <span className="shrink-0 text-caption text-ink-faint">{view.source === "default" ? "The run's first buffer" : "Your choice"}</span>
      </>
    );
  }
  if (view.kind === "buffer") return <span className="text-body text-ink-soft">This is a buffer, so nothing is subtracted from it.</span>;
  if (view.kind === "off") {
    return <span className="text-body text-ink-soft">{view.reason === "chosen" ? "No background. Nothing is subtracted." : "No buffer in this run. Choose one to subtract."}</span>;
  }
  return <span className="text-body text-ink-faint">Loading</span>;
}

function traceLine(view: BackgroundView): ReactNode {
  if (view.kind !== "set") return null;
  if (view.trace === "loading") return "Loading the trace";
  if (view.trace === "missing") return "No trace for this background";
  if (view.trace === "hidden") return "The setup of this background is hidden";
  if (!view.scale) return null;
  return (
    <>
      <span data-scale-source={view.scale.source}>{SCALE_NOTES[view.scale.source]}</span>
      {" "}
      <span data-scale-value className="text-ink-faint">{`×${view.scale.scale.toPrecision(3)}`}</span>
    </>
  );
}

export function BackgroundModule({ view, editable, picker, onChange, onNone, onDefault }: Props): JSX.Element {
  const loading = view.kind === "loading";
  const idle = loading || !editable;
  const hasControls = view.kind !== "buffer";
  const offerDefault = (view.kind === "set" || view.kind === "off") && view.canUseRunBuffer;
  const canTurnOff = view.kind !== "off" || view.reason !== "chosen";
  return (
    <Section title="Background" className="relative shrink-0">
      <div className="flex min-h-6 items-baseline justify-between gap-2" data-testid="background-name">{nameLine(view)}</div>
      <p className={LINE + " text-ink-soft"}>{traceLine(view)}</p>
      {!hasControls && <div className="mt-1 min-h-6" aria-hidden="true" />}
      {hasControls && (
        <div className="mt-1 flex min-h-6 flex-wrap items-center gap-1.5">
          <button type="button" className={ACTION} disabled={idle} onClick={onChange}>
            Change <Kbd>B</Kbd>
          </button>
          {canTurnOff && (
            <button type="button" className={ACTION} disabled={idle} onClick={onNone}>None</button>
          )}
          {offerDefault && (
            <button type="button" className={ACTION} disabled={idle} onClick={onDefault}>Use the run's buffer</button>
          )}
        </div>
      )}
      {picker}
    </Section>
  );
}
