import { useEffect, useMemo, useRef, useState } from "react";
import { labelOptions } from "../../lib/labels";
import { useListHighlight } from "./useListHighlight";

interface Props {
  vocabulary: readonly string[];
  onSample: readonly string[];
  onPick: (label: string) => void;
  onRemoveLast: () => void;
  onClose: () => void;
}

type Row = { label: string; create: boolean };

export function LabelPicker({ vocabulary, onSample, onPick, onRemoveLast, onClose }: Props): JSX.Element {
  const [q, setQ] = useState("");
  const input = useRef<HTMLInputElement | null>(null);
  const root = useRef<HTMLDivElement | null>(null);
  useEffect(() => input.current?.focus(), []);
  useEffect(() => {
    const onDown = (ev: MouseEvent): void => {
      if (root.current && ev.target instanceof Node && !root.current.contains(ev.target)) onClose();
    };
    document.addEventListener("mousedown", onDown, true);
    return () => document.removeEventListener("mousedown", onDown, true);
  }, [onClose]);
  const rows: Row[] = useMemo(() => {
    const o = labelOptions(vocabulary, q, onSample);
    return [...o.matches.map((label) => ({ label, create: false })), ...(o.create ? [{ label: o.create, create: true }] : [])];
  }, [vocabulary, q, onSample]);
  const pick = (r: Row | undefined): void => {
    if (!r) return;
    onPick(r.label);
    onClose();
  };
  const highlight = useListHighlight({ count: rows.length, onCommit: (i) => pick(rows[i]), onEscape: onClose, resetKey: q });
  const at = highlight.at;
  return (
    <div
      ref={root}
      role="dialog"
      aria-label="Add a label"
      data-keys-trap
      data-testid="workspace-label-picker"
      className="absolute right-0 top-full z-20 mt-1 w-64 rounded-sm border border-hair bg-plate shadow-lg"
      onKeyDown={(e) => {
        if (highlight.onKeyDown(e)) return;
        if (e.key === "Backspace" && q === "" && !e.repeat) {
          e.preventDefault();
          onRemoveLast();
        }
      }}
    >
      <input
        ref={input}
        type="search"
        value={q}
        onChange={(e) => setQ(e.target.value)}
        placeholder="Label"
        aria-label="Label name"
        className="w-full border-b border-hair bg-plate px-2 py-1 text-body placeholder:text-ink-faint focus:outline-none"
      />
      <ul className="max-h-64 overflow-y-auto py-1">
        {rows.map((r, i) => (
          <li key={(r.create ? "+" : "") + r.label}>
            <button
              type="button"
              data-label-option={r.create ? "create" : "existing"}
              tabIndex={-1}
              aria-current={i === at ? "true" : undefined}
              onMouseEnter={() => highlight.setAt(i)}
              onMouseDown={(e) => e.preventDefault()}
              onClick={() => pick(r)}
              className={"w-full px-3 py-1 text-left text-body " + (i === at ? "bg-row-cursor" : "") + (r.create ? " text-ink-soft" : "")}
            >
              {r.create ? `Create “${r.label}”` : r.label}
            </button>
          </li>
        ))}
        {rows.length === 0 && <li className="px-3 py-2 text-caption">No more labels to add.</li>}
      </ul>
    </div>
  );
}
