import { useCallback, useRef, useState } from "react";
import { prepare, layout } from "@chenglou/pretext";
import { MentionPicker } from "./MentionPicker";

interface MentionComposeProps {
  disabled: boolean;
  onSubmit: (body: string) => void;
}

function computeTextareaHeight(text: string, width: number, font: string): number {
  if (width <= 0) return 40;
  const prepared = prepare(text || " ", font, { whiteSpace: "pre-wrap" });
  const { height } = layout(prepared, width, 20);
  return Math.max(40, Math.min(height + 20, 200));
}

export function MentionCompose({ disabled, onSubmit }: MentionComposeProps): JSX.Element {
  const [text, setText] = useState("");
  const [pickerQuery, setPickerQuery] = useState<string | null>(null);
  const [atStart, setAtStart] = useState(0);
  const ref = useRef<HTMLTextAreaElement>(null);
  const wrapRef = useRef<HTMLDivElement>(null);

  const updateMeta = useCallback((val: string, cursor: number) => {
    const before = val.slice(0, cursor);
    const atIdx = before.lastIndexOf("@");
    if (atIdx !== -1 && !before.slice(atIdx + 1).includes(" ")) {
      setAtStart(atIdx);
      setPickerQuery(before.slice(atIdx + 1));
    } else {
      setPickerQuery(null);
    }

    if (ref.current && wrapRef.current) {
      const w = wrapRef.current.clientWidth - 20;
      const font = getComputedStyle(ref.current).font;
      try {
        ref.current.style.height = `${computeTextareaHeight(val, w, font)}px`;
      } catch {
        // pretext can throw in jsdom (no canvas measurement) — fall back to fixed height
        ref.current.style.height = "40px";
      }
    }
  }, []);

  const handleChange = (e: React.ChangeEvent<HTMLTextAreaElement>): void => {
    const val = e.target.value;
    setText(val);
    updateMeta(val, e.target.selectionStart ?? val.length);
  };

  const handleSelect = useCallback((token: string) => {
    const cursor = ref.current?.selectionStart ?? text.length;
    const before = text.slice(0, atStart);
    const after  = text.slice(cursor);
    const next   = before + token + after;
    setText(next);
    setPickerQuery(null);
    if (ref.current) {
      ref.current.focus();
      const pos = before.length + token.length;
      ref.current.setSelectionRange(pos, pos);
    }
  }, [text, atStart]);

  const trySubmit = useCallback(() => {
    const trimmed = text.trim();
    if (!trimmed || disabled) return;
    onSubmit(trimmed);
    setText("");
    setPickerQuery(null);
    if (ref.current) ref.current.style.height = "40px";
  }, [text, disabled, onSubmit]);

  const onKeyDown = (e: React.KeyboardEvent<HTMLTextAreaElement>): void => {
    if (e.key === "Enter" && !e.shiftKey && pickerQuery === null) {
      e.preventDefault();
      trySubmit();
    }
    // Esc and arrow keys handled inside MentionPicker via window listener
  };

  return (
    <div ref={wrapRef} className="flex-shrink-0 border-t border-border bg-bg px-2.5 py-2 relative">
      {pickerQuery !== null && (
        <MentionPicker
          query={pickerQuery}
          onSelect={handleSelect}
          onDismiss={() => setPickerQuery(null)}
        />
      )}
      <textarea
        ref={ref}
        value={text}
        onChange={handleChange}
        onKeyDown={onKeyDown}
        placeholder={disabled ? "Sign in to post…" : "Write a note… (@ to mention)"}
        data-testid="chat-compose"
        className="w-full resize-none bg-transparent text-fg text-base font-sans
                   placeholder:text-fg-dim outline-0 border-0"
        style={{ height: "40px" }}
      />
      <div className="flex items-center justify-between text-xs text-fg-dim">
        <span>
          <kbd className="border border-border rounded px-1">⏎</kbd> send
          {" · "}
          <kbd className="border border-border rounded px-1">⇧⏎</kbd> newline
          {" · "}
          <kbd className="border border-border rounded px-1">@</kbd> mention
        </span>
      </div>
    </div>
  );
}
