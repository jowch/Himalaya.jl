import { useId } from "react";
import type { ReactNode } from "react";
import { Tooltip } from "./Tooltip";

interface TermProps {
  explain: string;
  children: ReactNode;
  className?: string;
}

const LOOK =
  "cursor-help rounded-sm underline decoration-dotted decoration-hair-strong decoration-1 underline-offset-3 hover:decoration-ink-faint focus-visible:outline-2 focus-visible:outline-offset-0 focus-visible:outline-accent";

export function Term({ explain, children, className }: TermProps): JSX.Element {
  const id = useId();
  return (
    <>
      <Tooltip content={explain} describes={false}>
        <span tabIndex={0} aria-describedby={id} className={className ? `${LOOK} ${className}` : LOOK}>
          {children}
        </span>
      </Tooltip>
      <span id={id} hidden>
        {explain}
      </span>
    </>
  );
}
