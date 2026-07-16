import type { ReactNode, MouseEvent as ReactMouseEvent } from "react";

export type ActionGroup = "Navigate" | "Act" | "Screen" | "Edit";

/** Fixed cross-page gestures. Keys are app-wide constants (see core.ts). */
export type CoreId = "back" | "openFocus" | "openLoupe" | "undo" | "redo" | "help" | "find";

/** A page-local verb id is any other string (e.g. "cull", "merge", "addPeak"). */
export type ActionId = CoreId | (string & {});

/** A page's current interaction mode is its own discriminated union; the layer
 *  only needs the `kind` tag to gate actions, so it sees the bare string. */
export type ModeKind = string;

export interface Action {
  id: ActionId;
  /** dock button · legend · aria-keyshortcuts · palette */
  label: string;
  /** normalized combos: "x" · "Mod+z" · "Shift+ArrowUp" · "Enter" · "?" */
  keys?: string[];
  group: ActionGroup;
  /** pure synchronous closure; the layer reads it before firing (never stale). */
  enabled?: () => boolean;
  run: (e?: KeyboardEvent | ReactMouseEvent) => void;
  /** show in dock; "primary" = Enter / double-click target, rendered prominently. */
  dock?: boolean | "primary";
  /** live only in this page-mode; omit = always live. */
  mode?: ModeKind;
  glyph?: ReactNode;
}

/** Spread onto every cursorable row. tabIndex is the roving substrate. */
export interface RowProps {
  ref: (el: HTMLElement | null) => void;
  tabIndex: 0 | -1;
  onClick: (e: ReactMouseEvent) => void;
  onDoubleClick: (e: ReactMouseEvent) => void;
  role: "row";
  "aria-current": "true" | undefined;
  "data-cursored": "true" | "false";
}

/** Feeds the shell's DockStepper. Mirrors DockStepperProps minus appearance. */
export interface CursorStepperProps {
  label: string;
  axis: "vertical" | "horizontal";
  testIdBase: string;
  count: ReactNode;
  onPrev: () => void;
  onNext: () => void;
  prevDisabled: boolean;
  nextDisabled: boolean;
}

export interface ListCursor {
  cursorId: number | null;
  selected: Set<number>;
  setCursor: (id: number) => void;
  moveBy: (delta: number) => void;
  activate: () => void;
  toggleSelect: (id?: number) => void;
  rowProps: (id: number) => RowProps;
  stepperProps: () => CursorStepperProps;
}
