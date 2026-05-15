import type { CompareMode } from "../hooks/useCompareMode";

interface Props {
  dirty: boolean;
  mode: CompareMode;
  onSave: () => void;
  isSaving: boolean;
}

const COPY: Record<CompareMode["kind"], string> = {
  "viewing":             "",  // hidden anyway
  "viewing-stale":       "",  // hidden anyway
  "editing-mine":        "Save changes",
  "editing-as-fork-of":  "Save as fork…",
  "creating-blank":      "Save",
  "creating-from-fork":  "Save fork",  // post-morph; user already confirmed the title
};

export function SavePill({ dirty, mode, onSave, isSaving }: Props): JSX.Element | null {
  if (!dirty) return null;
  // Viewing modes never own a dirty draft; guard regardless of the caller so
  // the pill can't render an empty-label button if `dirty` is ever wrong.
  if (mode.kind === "viewing" || mode.kind === "viewing-stale") return null;
  return (
    <button
      type="button"
      data-testid="save-pill"
      data-interactable="button"
      data-mode={mode.kind}
      disabled={isSaving}
      onClick={onSave}
      title="Save (Cmd+Enter)"
      className="px-3 py-1 rounded bg-accent text-bg disabled:opacity-50 text-sm font-medium"
    >
      {isSaving ? "Saving…" : COPY[mode.kind]}
    </button>
  );
}
