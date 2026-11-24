import type { ReactNode } from "react";
import { Tooltip } from "./Tooltip";

interface TermProps {
  explain: string;
  children: ReactNode;
  className?: string;
}

const LOOK =
  "cursor-help rounded-sm underline decoration-dotted decoration-hair-strong decoration-1 underline-offset-3 hover:decoration-ink-faint focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent";

export function Term({ explain, children, className }: TermProps): JSX.Element {
  return (
    <Tooltip content={explain}>
      <span tabIndex={0} aria-description={explain} className={className ? `${LOOK} ${className}` : LOOK}>
        {children}
      </span>
    </Tooltip>
  );
}
