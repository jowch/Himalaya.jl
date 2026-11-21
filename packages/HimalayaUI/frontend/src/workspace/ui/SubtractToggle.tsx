interface Props {
  name: string;
  on: boolean;
  onChange: (on: boolean) => void;
}

export function SubtractToggle({ name, on, onChange }: Props): JSX.Element {
  const option = (pressed: boolean, label: string, value: boolean): JSX.Element => (
    <button
      type="button"
      title={label}
      aria-pressed={pressed}
      onClick={() => onChange(value)}
      className={"max-w-40 truncate rounded-sm px-2 py-0.5 " + (pressed ? "bg-plate text-ink shadow-sm" : "text-ink-soft hover:text-ink")}
    >
      {label}
    </button>
  );
  return (
    <div
      role="group"
      aria-label="Subtraction"
      data-testid="workspace-subtraction-toggle"
      className="absolute right-3 top-1 z-10 flex items-center gap-0.5 rounded-sm border border-hair bg-paper p-0.5 text-caption"
    >
      {option(!on, "Raw", false)}
      {option(on, `Minus ${name}`, true)}
    </div>
  );
}
