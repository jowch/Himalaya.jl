interface SectionLabelProps {
  children: React.ReactNode;
}

export function SectionLabel({ children }: SectionLabelProps): JSX.Element {
  return (
    <h3 className="text-label mb-2">
      {children}
    </h3>
  );
}
