interface HintTextProps {
  children: React.ReactNode;
}

export function HintText({ children }: HintTextProps): JSX.Element {
  return (
    <p className="text-fg-dim text-base italic">{children}</p>
  );
}
