import { KbLegend } from "../ui";
import { SHORTCUTS, shortcutLabel, type ShortcutId } from "./shortcuts";

/**
 * KbdLegend — a registry-driven keyboard legend. Each entry's key cap and its
 * description are read from the shared shortcut registry by `id`, so the printed
 * hint can never disagree with the live binding (the whole point of the unified
 * library: one source of truth feeds handlers, `aria-keyshortcuts`, AND the
 * on-screen legend). Appearance lives in the `KbLegend` ui primitive; this is
 * the registry adapter. Pass placement-only `className`.
 */
export function KbdLegend({
  ids,
  className,
  testId,
}: {
  ids: ShortcutId[];
  /** PLACEMENT-ONLY. Forwarded to the KbLegend row. */
  className?: string;
  /** Override the row's test id when more than one legend coexists on a surface. */
  testId?: string;
}): JSX.Element {
  const shortcuts = ids.map((id) => ({
    keyLabel: shortcutLabel(id),
    description: SHORTCUTS[id].label,
  }));
  return (
    <KbLegend
      shortcuts={shortcuts}
      {...(className ? { className } : {})}
      {...(testId ? { testId } : {})}
    />
  );
}
