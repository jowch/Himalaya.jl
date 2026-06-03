import type { ReactNode } from "react";
import { Swatch, PhaseLabel } from "../ui";

export interface SeriesMemberRowProps {
  /** Dominant + any coexisting phase names, display order. `[]` = form factor. */
  phases: string[];
  /** The variable value (e.g. "1:0.5"), right-aligned, mono. */
  variableValue: ReactNode;
  /** Page-derived lattice + first-peak line (e.g. "a = 195 Å · q₁ 0.057 Å⁻¹"). */
  dataLine: ReactNode;
  /** Page-owned highlight (synced with the waterfall hot row). */
  hot?: boolean;
  onHover?: () => void;
  onLeave?: () => void;
  /** PLACEMENT-ONLY. Appended last. */
  className?: string;
}

// The series-plot reading "member" (mockup `.member`) — a self-decoding row:
// swatch shows phase(s), names are colored per phase, the variable sits right,
// and a mono data line carries lattice + q₁. Presentational: `hot` + handlers
// are page-owned so hovering here / the waterfall / nothing stays in sync.
export function SeriesMemberRow({
  phases,
  variableValue,
  dataLine,
  hot = false,
  onHover,
  onLeave,
  className,
}: SeriesMemberRowProps): JSX.Element {
  const formFactor = phases.length === 0;
  return (
    <div
      data-testid="series-member-row"
      {...(hot ? { "data-hot": "true" } : {})}
      onMouseEnter={onHover}
      onMouseLeave={onLeave}
      className={`flex items-center gap-[9px] px-[9px] py-[7px] rounded cursor-pointer border ${hot ? "bg-plate border-hair" : "border-transparent hover:bg-plate"}${className ? ` ${className}` : ""}`}
    >
      {formFactor ? (
        <Swatch phase="" empty size="md" shape="circle" />
      ) : phases.length >= 2 ? (
        <Swatch phase={phases[0]!} coexistWith={phases[1]!} size="md" shape="circle" />
      ) : (
        <Swatch phase={phases[0]!} size="md" shape="circle" />
      )}
      <div className="flex-1 min-w-0 flex flex-col gap-px">
        <div className="flex items-baseline justify-between gap-2">
          <span className="inline-flex items-baseline min-w-0">
            {formFactor ? (
              <span className="text-meta font-semibold text-ink-faint">Form factor</span>
            ) : (
              phases.map((p, i) => (
                <span key={p} className="inline-flex items-baseline">
                  {i > 0 && <span className="text-meta text-ink-faint mx-0.5"> + </span>}
                  <PhaseLabel phase={p} className="text-meta font-semibold">{p}</PhaseLabel>
                </span>
              ))
            )}
          </span>
          <span className="text-data text-ink-soft shrink-0">{variableValue}</span>
        </div>
        <div className="text-data text-ink-faint truncate">{dataLine}</div>
      </div>
    </div>
  );
}
