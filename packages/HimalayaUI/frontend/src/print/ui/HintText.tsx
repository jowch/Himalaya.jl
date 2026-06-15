interface HintTextProps {
  children: React.ReactNode;
  /** Placement-only className (spacing, max-width, etc). Appearance is fixed. */
  className?: string;
}

export function HintText({ children, className = "" }: HintTextProps): JSX.Element {
  return (
    <p className={`text-ink-soft text-base italic ${className}`.trim()}>{children}</p>
  );
}
