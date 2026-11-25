import { useId } from "react";
import type { ReactNode } from "react";
import { createPortal } from "react-dom";
import { Tooltip } from "./Tooltip";

interface TermProps {
  explain: string;
  children: ReactNode;
  className?: string;
  nested?: boolean;
}

const MARK = "cursor-help underline decoration-dotted decoration-hair-strong decoration-1 underline-offset-3 hover:decoration-ink-faint";
const FOCUSABLE = `${MARK} rounded-sm focus-visible:outline-2 focus-visible:outline-offset-0 focus-visible:outline-accent`;

const withClass = (look: string, className: string | undefined): string => (className ? `${look} ${className}` : look);

export function Explanation({ id, texts }: { id: string; texts: readonly string[] }): JSX.Element {
  return createPortal(
    <span id={id} hidden>
      {texts.join(" ")}
    </span>,
    document.body,
  );
}

export function Term({ explain, children, className, nested = false }: TermProps): JSX.Element {
  const id = useId();
  if (nested) {
    return (
      <Tooltip content={explain} describes={false} focusWithin>
        <span className={withClass(MARK, className)}>{children}</span>
      </Tooltip>
    );
  }
  return (
    <>
      <Tooltip content={explain} describes={false}>
        <span tabIndex={0} aria-describedby={id} className={withClass(FOCUSABLE, className)}>
          {children}
        </span>
      </Tooltip>
      <Explanation id={id} texts={[explain]} />
    </>
  );
}
