import { buttonSizeClass } from "./Button";
import { Kbd } from "./Kbd";
import { Tooltip } from "./Tooltip";

interface Props {
  open: boolean;
  disabled: boolean;
  onClick: () => void;
}

export function CustomIndexButton({ open, disabled, onClick }: Props): JSX.Element {
  return (
    <Tooltip content={open ? "Close the custom index" : undefined}>
      <button
        type="button"
        aria-pressed={open}
        disabled={disabled && !open}
        onClick={onClick}
        className={
          buttonSizeClass.sm + " border disabled:opacity-50 " +
          (open ? "border-ink-soft bg-paper-sunk text-ink" : "border-dashed border-hair-strong text-ink-soft hover:text-ink")
        }
      >
        Custom index <Kbd>{open ? "Esc" : "I"}</Kbd>
      </button>
    </Tooltip>
  );
}
