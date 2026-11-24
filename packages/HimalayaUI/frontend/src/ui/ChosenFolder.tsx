import { cx } from "../lib/cx";
import { Button } from "./Button";

export interface ChosenFolderProps {
  label: string;
  path: string;
  onChange: () => void;
  testId?: string;
  className?: string;
}

export function ChosenFolder({ label, path, onChange, testId = "chosen-folder", className = "" }: ChosenFolderProps): JSX.Element {
  return (
    <div data-testid={testId} className={cx("flex flex-wrap items-center gap-x-3 gap-y-1 rounded-sm border border-hair-strong bg-plate px-3 py-1.5", className)}>
      <span className="text-meta text-ink-soft">{label}</span>
      <span data-testid={`${testId}-path`} className="min-w-0 break-all font-semibold text-ink">{path}</span>
      <Button variant="ghost" className="ml-auto" onClick={onChange}>Change</Button>
    </div>
  );
}
