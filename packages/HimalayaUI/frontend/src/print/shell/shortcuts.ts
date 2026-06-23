// Unified keyboard shortcut library — the single source of truth for the app's
// gesture vocabulary. Handlers (useShortcuts), on-screen hints, the legend
// (<KbdLegend>) and aria-keyshortcuts all derive from this registry, so they can
// never drift apart. Design: print/shell/AGENTS.md §"Keyboard shortcut registry"
//
// Rev-2 axes (2026-06-21): ↑/↓ = sample nav, ←/→ = detail nav (frame on
// Corpus/Loupe, candidate on Focus — same id, page-interpreted handler).
// SUPERSEDES the "Locked decisions (Jonathan 2026-06-13)" comment ([ ] sample,
// ←→ exposure, ↑↓ candidate) — that model has been retired by the ingestion
// redesign shell unification (M2, T2.1).

export type ShortcutId =
  | "prevSample"
  | "nextSample"
  | "prevDetail"
  | "nextDetail"
  | "openFocus"
  | "openLoupe"
  | "toggleSelect"
  | "extendPrev"
  | "extendNext"
  | "selectAll"
  | "restore"
  | "helpOverlay"
  | "drop"
  | "keep"
  | "representative"
  | "dismiss"
  | "undo"
  | "redo"
  | "reorderUp"
  | "reorderDown"
  | "addSample"
  | "confirm"
  | "addPeak"
  | "find";

export type ShortcutGroup = "Navigate" | "Screen" | "Edit" | "General";

export interface ShortcutDef {
  id: ShortcutId;
  /** Normalized combos (see `eventCombo`): `"x"`, `"["`, `"ArrowLeft"`,
   *  `"Mod+z"`, `"Mod+Shift+z"`, `"Alt+ArrowUp"`, `"/"`. `Mod` = ⌘ (mac) / Ctrl. */
  keys: string[];
  /** Human label for legends/tooltips (the action, not the key). */
  label: string;
  group: ShortcutGroup;
}

export const SHORTCUTS: Record<ShortcutId, ShortcutDef> = {
  // Navigate — rev-2 axes: ↑/↓ = sample, ←/→ = detail (page-interpreted).
  prevSample: { id: "prevSample", keys: ["ArrowUp"], label: "Previous sample", group: "Navigate" },
  nextSample: { id: "nextSample", keys: ["ArrowDown"], label: "Next sample", group: "Navigate" },
  prevDetail: { id: "prevDetail", keys: ["ArrowLeft"], label: "Previous detail", group: "Navigate" },
  nextDetail: { id: "nextDetail", keys: ["ArrowRight"], label: "Next detail", group: "Navigate" },
  openFocus: { id: "openFocus", keys: ["Enter"], label: "Open Focus workspace", group: "Navigate" },
  openLoupe: { id: "openLoupe", keys: ["l"], label: "Open Loupe", group: "Navigate" },
  // Select — multi-select verbs.
  toggleSelect: { id: "toggleSelect", keys: [" "], label: "Toggle selection", group: "Navigate" },
  extendPrev: { id: "extendPrev", keys: ["Shift+ArrowLeft"], label: "Extend selection left", group: "Navigate" },
  extendNext: { id: "extendNext", keys: ["Shift+ArrowRight"], label: "Extend selection right", group: "Navigate" },
  selectAll: { id: "selectAll", keys: ["Mod+a"], label: "Select all", group: "Navigate" },
  restore: { id: "restore", keys: ["Backspace"], label: "Restore / remove from selection", group: "Navigate" },
  // General — help + dismiss + find.
  helpOverlay: { id: "helpOverlay", keys: ["?"], label: "Keyboard shortcuts", group: "General" },
  dismiss: { id: "dismiss", keys: ["Escape"], label: "Back / dismiss", group: "General" },
  find: { id: "find", keys: ["Mod+k", "/"], label: "Find a sample", group: "General" },
  // Screen — exposure status verbs.
  drop: { id: "drop", keys: ["x"], label: "Drop", group: "Screen" },
  keep: { id: "keep", keys: ["k"], label: "Keep", group: "Screen" },
  representative: { id: "representative", keys: ["r"], label: "Set representative", group: "Screen" },
  // Edit — undoable mutations + list reorder.
  undo: { id: "undo", keys: ["Mod+z"], label: "Undo", group: "Edit" },
  redo: { id: "redo", keys: ["Mod+Shift+z"], label: "Redo", group: "Edit" },
  reorderUp: { id: "reorderUp", keys: ["Alt+ArrowUp"], label: "Move up", group: "Edit" },
  reorderDown: { id: "reorderDown", keys: ["Alt+ArrowDown"], label: "Move down", group: "Edit" },
  // Surface-local in-surface verbs (Series A / ⌘Enter, Focus P) — kept out of
  // OVERLAY_IDS (the overlay lists only the global/navigation keymap).
  addSample: { id: "addSample", keys: ["a"], label: "Add a sample", group: "Edit" },
  confirm: { id: "confirm", keys: ["Mod+Enter"], label: "Confirm / build", group: "Edit" },
  addPeak: { id: "addPeak", keys: ["p"], label: "Toggle add-peak mode", group: "Edit" },
};

/**
 * Normalize a keyboard event to a combo string in the registry's grammar:
 * `(Mod+)(Alt+)(Shift+)key`, where `Mod` collapses ⌘/Ctrl (cross-platform) and a
 * single character key is lowercased (so CapsLock-X and x both read `x`, while
 * Shift+X reads `Shift+x`). Named keys (ArrowLeft, Escape) pass through as-is.
 */
export function eventCombo(e: KeyboardEvent): string {
  // '?' is Shift+/ on US layouts but layout-variable; emit a stable token so the
  // help binding is layout-robust (mirrors the single-char lowercasing below).
  if (e.key === '?') return '?'
  const parts: string[] = [];
  if (e.metaKey || e.ctrlKey) parts.push("Mod");
  if (e.altKey) parts.push("Alt");
  if (e.shiftKey) parts.push("Shift");
  const k = e.key.length === 1 ? e.key.toLowerCase() : e.key;
  parts.push(k);
  return parts.join("+");
}

/** Resolve a keyboard event to the action it triggers, or null. */
export function matchShortcut(e: KeyboardEvent): ShortcutId | null {
  const combo = eventCombo(e);
  for (const def of Object.values(SHORTCUTS)) {
    if (def.keys.includes(combo)) return def.id;
  }
  return null;
}

const GLYPH: Record<string, string> = {
  ArrowLeft: "←",
  ArrowRight: "→",
  ArrowUp: "↑",
  ArrowDown: "↓",
  Escape: "Esc",
};

/** Render one normalized combo for display (mac glyphs vs spelled-out words). */
export function keyComboLabel(combo: string, isMac: boolean = defaultIsMac()): string {
  const parts = combo.split("+");
  const base = parts[parts.length - 1]!;
  const mods = parts.slice(0, -1);
  const baseLabel = GLYPH[base] ?? (base.length === 1 ? base.toUpperCase() : base);
  if (isMac) {
    const pre = mods
      .map((m) => (m === "Mod" ? "⌘" : m === "Shift" ? "⇧" : m === "Alt" ? "⌥" : m))
      .join("");
    return pre + baseLabel;
  }
  const pre = mods.map((m) => (m === "Mod" ? "Ctrl" : m)).join("+");
  return pre ? `${pre}+${baseLabel}` : baseLabel;
}

/** The display label for an action's primary (first) binding. */
export function shortcutLabel(id: ShortcutId, isMac: boolean = defaultIsMac()): string {
  return keyComboLabel(SHORTCUTS[id].keys[0]!, isMac);
}

function defaultIsMac(): boolean {
  if (typeof navigator === "undefined") return false;
  return /mac|iphone|ipad/i.test(navigator.userAgent || "");
}
