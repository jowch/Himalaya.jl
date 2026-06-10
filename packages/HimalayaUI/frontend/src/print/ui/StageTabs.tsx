import { Dot } from "./Dot";

export type StageKey = "samples" | "index" | "series";

interface StageTabsProps {
  active: StageKey;
  onChange: (s: StageKey) => void;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

const STAGES: ReadonlyArray<{ key: StageKey; label: string }> = [
  { key: "samples", label: "Samples" },
  { key: "index", label: "Index" },
  { key: "series", label: "Series" },
];

export function StageTabs({ active, onChange, className }: StageTabsProps): JSX.Element {
  return (
    <div
      role="tablist"
      aria-label="Workspace stage"
      data-testid="stage-tabs"
      className={cx("flex items-center gap-1", className)}
    >
      {STAGES.map((s) => (
        <button
          key={s.key}
          type="button"
          role="tab"
          aria-selected={active === s.key}
          data-stage={s.key}
          data-active={active === s.key ? "true" : "false"}
          onClick={() => onChange(s.key)}
          // text-sm (~11.5px) stands in for the mockup's 11px label (-0.5px, within
          // the named scale; no off-scale token added). The accessible name is the
          // proper-case label; `uppercase` only restyles it visually (textContent
          // stays "Samples").
          className={cx(
            "inline-flex items-center text-sm font-semibold uppercase tracking-wide px-3 py-1.5 rounded-sm transition-colors",
            "focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent",
            active === s.key ? "text-ink bg-paper-sunk" : "text-ink-soft hover:text-ink",
          )}
        >
          <Dot tone={active === s.key ? "accent" : "neutral"} size="xs" aria-hidden className="mr-1.5" />
          {s.label}
        </button>
      ))}
    </div>
  );
}
