import type { KnownPhase } from "../phases";

interface ShortName {
  base: string;
  sub?: string;
}

export const LAMELLAR_SHORT: ShortName = { base: "L", sub: "α" };
export const HEXAGONAL_SHORT: ShortName = { base: "H" };

const SHORT_NAMES: Partial<Record<KnownPhase, ShortName>> = {
  Lamellar: LAMELLAR_SHORT,
  Hexagonal: HEXAGONAL_SHORT,
};

export function isShortened(phase: string): boolean {
  return phase in SHORT_NAMES;
}

export function PhaseShortName({ phase }: { phase: string }): JSX.Element {
  const short = SHORT_NAMES[phase as KnownPhase];
  if (!short) return <span>{phase}</span>;
  return (
    <span>
      {short.base}
      {short.sub && <sub>{short.sub}</sub>}
    </span>
  );
}
