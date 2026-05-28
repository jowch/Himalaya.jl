interface Props {
  flagCount: number;
  memberCount: number;
  keyLabel: string;
  canBuild: boolean;
  onBuild: () => void;
}

/**
 * The narrative gate footer (series-scoping.html `.scope-foot`): a confirmation
 * state line (amber "N to check" / sage "ready"), the metadata-as-byproduct
 * note, and the single durable action. Build button = Print ink (S-B), greyed
 * ink when gated (S-C) — never the ice-blue accent.
 */
export function ScopingFoot({ flagCount, memberCount, keyLabel, canBuild, onBuild }: Props): JSX.Element {
  const ready = flagCount === 0;
  const stateText = ready
    ? `All ${memberCount} values confirmed — ready to build`
    : `${flagCount} value${flagCount === 1 ? "" : "s"} to check before you can build`;
  return (
    <div className="mt-6 flex items-center justify-between gap-5 border-t border-hair pt-4">
      <div className="flex flex-col gap-1">
        <div
          data-testid="scoping-foot-state"
          className={`flex items-center gap-2 text-[12.5px] font-semibold ${ready ? "text-ink" : "text-print-accent"}`}
        >
          <span
            className="h-2 w-2 shrink-0 rounded-full"
            style={{ background: ready ? "var(--color-success)" : "var(--color-print-accent)" }}
          />
          {stateText}
        </div>
        <div className="max-w-[42ch] text-[10.5px] text-ink-faint">
          Confirming records the {keyLabel} on every sample — the next series that needs it already knows.
        </div>
      </div>
      <button
        type="button"
        data-testid="scoping-open-confirm"
        disabled={!canBuild}
        onClick={onBuild}
        title={canBuild ? undefined : "Check the flagged values above before building"}
        className="shrink-0 rounded-md border border-ink bg-ink px-[18px] py-[11px] text-[13px] font-semibold text-paper transition-colors hover:bg-ink/85 disabled:cursor-not-allowed disabled:border-hair-strong disabled:bg-paper-sunk disabled:text-ink-faint"
      >
        Confirm &amp; build →
      </button>
    </div>
  );
}
