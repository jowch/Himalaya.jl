import type { Action, ActionGroup, CoreId } from "./types";

/** The fixed cross-page gestures. Keys are app-wide constants — a page supplies
 *  only the handler (and optionally enabled/dock/label). Evolved from the old
 *  shortcuts.ts global vocabulary, minus the page-interpreted ids. */
export const CORE: Record<CoreId, { label: string; keys: string[]; group: ActionGroup }> = {
  back:      { label: "Back",  keys: ["Escape"],      group: "Navigate" },
  openFocus: { label: "Focus", keys: ["Enter"],       group: "Navigate" },
  openLoupe: { label: "Loupe", keys: ["l"],           group: "Navigate" },
  undo:      { label: "Undo",  keys: ["Mod+z"],       group: "Edit" },
  redo:      { label: "Redo",  keys: ["Mod+Shift+z"], group: "Edit" },
  help:      { label: "Shortcuts", keys: ["?"],       group: "Screen" },
  find:      { label: "Find",  keys: ["/", "Mod+k"],  group: "Screen" },
};

type CoreOverride = {
  run: Action["run"];
  enabled?: () => boolean;
  dock?: Action["dock"];
  label?: string;
  mode?: string;
};

export function core(id: CoreId, over: CoreOverride): Action {
  const base = CORE[id];
  const a: Action = {
    id,
    label: over.label ?? base.label,
    keys: base.keys,
    group: base.group,
    run: over.run,
  };
  if (over.enabled) a.enabled = over.enabled;
  if (over.dock !== undefined) a.dock = over.dock;
  if (over.mode !== undefined) a.mode = over.mode;
  return a;
}

type PageDef = {
  label: string;
  keys?: string[];
  group: ActionGroup;
  run: Action["run"];
  enabled?: () => boolean;
  dock?: Action["dock"];
  mode?: string;
};

export function page(id: string, def: PageDef): Action {
  const a: Action = { id, label: def.label, group: def.group, run: def.run };
  if (def.keys) a.keys = def.keys;
  if (def.enabled) a.enabled = def.enabled;
  if (def.dock !== undefined) a.dock = def.dock;
  if (def.mode !== undefined) a.mode = def.mode;
  return a;
}

const CORE_KEYS = new Set(Object.values(CORE).flatMap((c) => c.keys));

/** Build-time guard: a page verb must not reuse a core key. Called inside
 *  usePageActions so a colliding declaration throws in dev/test immediately. */
export function assertNoCoreCollision(actions: Action[]): void {
  const coreIds = new Set(Object.keys(CORE));
  for (const a of actions) {
    if (coreIds.has(a.id)) continue;
    for (const k of a.keys ?? []) {
      if (CORE_KEYS.has(k)) {
        throw new Error(`Action "${a.id}" reuses core key "${k}"; use core("${k}") or pick another key.`);
      }
    }
  }
}
