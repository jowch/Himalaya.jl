import { useId } from "react";
import { minusExplanation, RAW } from "./explanations";
import { Explanation, Term } from "./Term";

interface Props {
  name: string;
  on: boolean;
  onChange: (on: boolean) => void;
  reserved?: boolean;
}

export function SubtractToggle({ name, on, onChange, reserved = false }: Props): JSX.Element {
  const id = useId();
  const option = (pressed: boolean, label: string, value: boolean, explain: string, full = label): JSX.Element => (
    <button
      type="button"
      aria-label={full}
      aria-pressed={pressed}
      aria-describedby={`${id}-${label}`}
      onClick={() => onChange(value)}
      className={"rounded-sm px-2 " + (pressed ? "bg-plate text-ink shadow-sm" : "text-ink-soft hover:text-ink")}
    >
      <Term explain={explain} nested>{label}</Term>
      <Explanation id={`${id}-${label}`} texts={[explain]} />
    </button>
  );
  return (
    <div
      role="group"
      aria-label="Subtraction"
      data-testid="workspace-subtraction-toggle"
      {...(reserved ? { "aria-hidden": true } : {})}
      className={"flex items-center rounded-sm border border-hair bg-paper text-caption" + (reserved ? " invisible" : "")}
    >
      {option(!on, "Raw", false, RAW)}
      {option(on, "Minus", true, minusExplanation(name), `Minus ${name}`)}
    </div>
  );
}
