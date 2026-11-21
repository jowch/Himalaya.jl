import { Kbd } from "../Kbd";

export function PaletteFooter({ keys }: { keys: readonly (readonly [string, string])[] }): JSX.Element {
  return (
    <div className="flex min-h-8 flex-wrap items-center gap-4 border-t border-hair px-3 text-xs text-ink-faint">
      {keys.map(([key, what]) => (
        <span key={key} className="inline-flex items-center gap-1.5">
          <Kbd>{key}</Kbd> {what}
        </span>
      ))}
    </div>
  );
}
