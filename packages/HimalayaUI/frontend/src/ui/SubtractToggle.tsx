import { Tooltip } from "./Tooltip";

interface Props {
  name: string;
  on: boolean;
  onChange: (on: boolean) => void;
}

export function SubtractToggle({ name, on, onChange }: Props): JSX.Element {
  const option = (pressed: boolean, label: string, value: boolean, full = label): JSX.Element => (
    <Tooltip content={full === label ? undefined : full}>
      <button
        type="button"
        aria-label={full}
        aria-pressed={pressed}
        onClick={() => onChange(value)}
        className={"rounded-sm px-2 " + (pressed ? "bg-plate text-ink shadow-sm" : "text-ink-soft hover:text-ink")}
      >
        {label}
      </button>
    </Tooltip>
  );
  return (
    <div
      role="group"
      aria-label="Subtraction"
      data-testid="workspace-subtraction-toggle"
      className="flex items-center rounded-sm border border-hair bg-paper text-caption"
    >
      {option(!on, "Raw", false)}
      {option(on, "Minus", true, `Minus ${name}`)}
    </div>
  );
}
