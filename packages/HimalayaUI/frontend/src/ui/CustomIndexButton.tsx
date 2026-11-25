import { Button } from "./Button";
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
      <Button type="button" variant={open ? "sunk" : "dashed"} size="sm" aria-pressed={open} disabled={disabled && !open} onClick={onClick}>
        Custom index <Kbd>{open ? "Esc" : "I"}</Kbd>
      </Button>
    </Tooltip>
  );
}
