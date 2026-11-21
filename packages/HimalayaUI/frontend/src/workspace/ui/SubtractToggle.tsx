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
      className={"max-w-48 truncate rounded-sm px-2 " + (pressed ? "bg-plate text-ink shadow-sm" : "text-ink-soft hover:text-ink")}
    >
      {label}
    </button>
  );
  return (
    <div
      role="group"
      aria-label="Subtraction"
      data-testid="workspace-subtraction-toggle"
      className="flex items-center rounded-sm border border-hair bg-paper text-caption"
    >
      {option(!on, "Raw", false)}
      {option(on, `Minus ${name}`, true)}
    </div>
  );
}
