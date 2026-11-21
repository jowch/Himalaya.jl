import { useEffect, useRef } from "react";
import { TAG_COLORS, type TagColor } from "../../api";
import { TAG_SWATCH } from "./tagSwatches";

interface Props {
  name: string;
  current: TagColor;
  isDefault: boolean;
  left: number;
  onPick: (color: TagColor | null) => void;
  onClose: () => void;
}

export function TagColorMenu({ name, current, isDefault, left, onPick, onClose }: Props): JSX.Element {
  const root = useRef<HTMLDivElement | null>(null);
  useEffect(() => {
    root.current?.querySelector<HTMLButtonElement>("[aria-checked='true']")?.focus();
    const onDown = (ev: MouseEvent): void => {
      if (root.current && ev.target instanceof Node && !root.current.contains(ev.target)) onClose();
    };
    document.addEventListener("mousedown", onDown, true);
    return () => document.removeEventListener("mousedown", onDown, true);
  }, [onClose]);
  const pick = (color: TagColor | null): void => {
    onPick(color);
    onClose();
  };
  return (
    <div
      ref={root}
      role="dialog"
      aria-label={`Colour of tag ${name}`}
      data-keys-trap
      data-testid="workspace-tag-colors"
      style={{ left }}
      className="absolute top-full z-20 mt-1 rounded-sm border border-hair bg-plate p-2 shadow-lg"
      onKeyDown={(e) => {
        if (e.key === "Escape") {
          e.preventDefault();
          onClose();
        }
      }}
    >
      <div role="radiogroup" aria-label={`Colour of tag ${name}`} className="flex items-center gap-1.5">
        {TAG_COLORS.map((c) => (
          <button
            key={c}
            type="button"
            role="radio"
            aria-checked={c === current}
            aria-label={TAG_SWATCH[c].name}
            title={TAG_SWATCH[c].name}
            data-tag-color={c}
            onClick={() => pick(c)}
            className={`h-5 w-5 rounded-full ${TAG_SWATCH[c].solid} ${c === current ? "ring-2 ring-ink ring-offset-2 ring-offset-plate" : ""}`}
          />
        ))}
      </div>
      {!isDefault && (
        <button type="button" onClick={() => pick(null)} className="mt-2 w-full text-left text-caption text-ink-soft hover:text-ink">
          Default colour
        </button>
      )}
    </div>
  );
}
