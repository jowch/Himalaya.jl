interface SectionLabelProps {
  children: React.ReactNode;
}

export function SectionLabel({ children }: SectionLabelProps): JSX.Element {
  return (
    <h3 className="text-fg-muted text-[11px] font-semibold uppercase tracking-widest mb-2">
      {children}
    </h3>
  );
}
