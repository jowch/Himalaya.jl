import { useEffect, useId, useRef, useState, type KeyboardEvent } from "react";
import { cycleOption, fieldErrors, type PaletteForm } from "./model";
import { PaletteFooter } from "./PaletteFooter";

interface Props {
  form: PaletteForm;
  onBack: () => void;
  onDone: () => void;
}

export function PaletteFormView({ form, onBack, onDone }: Props): JSX.Element {
  const [values, setValues] = useState<Record<string, string>>(() => Object.fromEntries(form.fields.map((f) => [f.id, f.value])));
  const [focus, setFocus] = useState(0);
  const [tried, setTried] = useState(false);
  const listId = useId();
  const errorId = useId();
  const refs = useRef<(HTMLElement | null)[]>([]);
  const errors = fieldErrors(form.fields, values);
  const field = form.fields[focus]!;
  const wide = form.fields.length === 1 && field.kind === "text";
  const shownError = tried ? errors[field.id] : undefined;

  const opening = useRef(!form.fields[0]?.typed);
  useEffect(() => {
    const el = refs.current[focus];
    el?.focus();
    if (!(el instanceof HTMLInputElement)) return;
    if (opening.current && focus === 0) el.select();
    else el.setSelectionRange(el.value.length, el.value.length);
  }, [focus]);
  useEffect(() => {
    if (focus !== 0) opening.current = false;
  }, [focus]);

  const set = (id: string, v: string): void => setValues((cur) => ({ ...cur, [id]: v }));
  const move = (to: number): void => setFocus((to + form.fields.length) % form.fields.length);
  const submit = (): void => {
    const bad = form.fields.findIndex((f) => errors[f.id]);
    if (bad >= 0) {
      setTried(true);
      setFocus(bad);
      return;
    }
    onDone();
    form.submit(values);
  };
  const onKeyDown = (e: KeyboardEvent<HTMLElement>): void => {
    if (e.metaKey || e.ctrlKey || e.altKey || e.nativeEvent.isComposing) return;
    const value = values[field.id] ?? "";
    if (e.key === "Tab") {
      e.preventDefault();
      move(focus + (e.shiftKey ? -1 : 1));
    } else if (e.key === "Enter" || e.key === "Escape") {
      e.preventDefault();
      if (e.repeat) return;
      if (e.key === "Enter") submit();
      else onBack();
    } else if (field.kind === "choice" && (e.key === "ArrowDown" || e.key === "ArrowUp")) {
      e.preventDefault();
      set(field.id, cycleOption(field.options ?? [], value, e.key === "ArrowDown" ? 1 : -1));
    } else if (e.key === "Backspace" && (field.kind === "choice" || value === "")) {
      e.preventDefault();
      if (e.repeat) return;
      if (focus > 0) move(focus - 1);
      else onBack();
    }
  };

  const options = field.kind === "choice" ? (field.options ?? []) : [];
  const chosen = options.findIndex((o) => o.value === values[field.id]);
  return (
    <>
      <div className="flex h-12 items-center gap-1.5 overflow-hidden border-b border-hair px-3" data-testid="workspace-palette-bar">
        <span className="mr-0.5 whitespace-nowrap text-lg font-semibold">{form.title}</span>
        {form.fields.map((f, i) => {
          const on = i === focus;
          const bad = tried && errors[f.id] !== undefined;
          const ring = bad ? "shadow-[inset_0_0_0_1px_var(--color-error)]" : on ? "shadow-[inset_0_0_0_1px_var(--color-accent)]" : "";
          const opt = f.options?.find((o) => o.value === values[f.id]);
          return (
            <span
              key={f.id}
              data-field={f.id}
              onMouseDown={(e) => {
                if (e.target instanceof HTMLInputElement) return setFocus(i);
                e.preventDefault();
                setFocus(i);
                refs.current[i]?.focus();
              }}
              className={"inline-flex h-7 min-w-0 items-center gap-1.5 whitespace-nowrap rounded-sm px-2 " + (wide ? "flex-1 " : "") + (on ? "bg-plate " : "bg-paper-sunk ") + ring}
            >
              <span className="text-xs text-ink-faint">{f.label}</span>
              {f.kind === "text" ? (
                <input
                  ref={(el) => {
                    refs.current[i] = el;
                  }}
                  type="text"
                  autoComplete="off"
                  spellCheck={false}
                  aria-label={f.label}
                  aria-invalid={bad || undefined}
                  aria-describedby={bad && on ? errorId : undefined}
                  value={values[f.id] ?? ""}
                  placeholder={f.placeholder}
                  size={wide ? 1 : Math.max(3, Math.min(28, (values[f.id] || f.placeholder || "").length + 1))}
                  onChange={(e) => set(f.id, e.target.value)}
                  onFocus={() => setFocus(i)}
                  onKeyDown={onKeyDown}
                  className={"min-w-[2ch] bg-transparent text-lg text-ink placeholder:text-ink-faint focus:outline-none" + (wide ? " flex-1" : "")}
                />
              ) : (
                <span
                  ref={(el) => {
                    refs.current[i] = el;
                  }}
                  role="combobox"
                  tabIndex={-1}
                  aria-label={f.label}
                  aria-expanded={on}
                  aria-controls={listId}
                  aria-activedescendant={on && chosen >= 0 ? `${listId}-${chosen}` : undefined}
                  onFocus={() => setFocus(i)}
                  onKeyDown={onKeyDown}
                  className="inline-flex cursor-pointer items-center gap-1.5 text-lg focus:outline-none"
                >
                  {opt?.swatch && <span aria-hidden="true" className={"h-2 w-2 shrink-0 rounded-full " + opt.swatch} />}
                  {opt?.label ?? ""}
                  <span aria-hidden="true" className="text-xs text-ink-faint">▾</span>
                </span>
              )}
            </span>
          );
        })}
      </div>
      <ul id={listId} role="listbox" aria-label={field.label} className="max-h-[50vh] overflow-y-auto p-1">
        {options.map((o, i) => (
          <li
            key={o.value}
            id={`${listId}-${i}`}
            role="option"
            aria-selected={i === chosen}
            aria-label={o.label}
            onMouseDown={(e) => e.preventDefault()}
            onClick={() => set(field.id, o.value)}
            className={"grid h-8 cursor-pointer grid-cols-[minmax(0,1fr)_auto] items-center gap-x-3 rounded-sm px-2 text-base " + (i === chosen ? "bg-row-cursor" : "")}
          >
            <span className="flex min-w-0 items-center gap-2">
              {o.swatch && <span aria-hidden="true" className={"h-2 w-2 shrink-0 rounded-full " + o.swatch} />}
              <span className="truncate">{o.label}</span>
            </span>
            <span aria-hidden="true" className="font-semibold text-accent">{i === chosen ? "✓" : ""}</span>
          </li>
        ))}
        {shownError && (
          <li id={errorId} role="presentation" className="px-2 py-2 text-xs text-error">
            {shownError}
          </li>
        )}
      </ul>
      <PaletteFooter
        keys={[
          ...(form.fields.length > 1 ? [["Tab", "next field"] as const] : []),
          ...(field.kind === "choice" ? [["↑ ↓", "choose"] as const] : []),
          ["↵", form.verb] as const,
          ["Esc", "back"] as const,
        ]}
      />
    </>
  );
}
