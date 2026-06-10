// src/print/components/ColdAssignPanel.tsx
import type { ReactNode } from "react";
import { Input, Kicker, HintText } from "../ui";
import type { ColdAssignRow } from "../pages/scopingDerive";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ColdAssignPanelProps {
  rows: ColdAssignRow[];
  /** The ordering-variable key the user is naming. */
  variableKey: string;
  onKeyChange: (key: string) => void;
  onValueChange: (sampleId: number, value: string) => void;
  /** The intro paragraph under the kicker. Undefined renders the default
   *  cold-corpus copy ("These samples share no tag key yet..."); an explicit
   *  `null` suppresses the paragraph entirely (custom mode, where that copy
   *  would be false); any other node renders in the same slot and typography. */
  intro?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * ColdAssignPanel — the cold-corpus path inside SeriesScopingPage.
 *
 * When no tag key can be proposed (cold / unstructured selection), this panel
 * lets the user name the ordering variable once and fill in each sample's value.
 * Composed from `src/print/ui/` primitives only; no legacy imports.
 *
 * Layout: a labelled key-name Input at the top, then one row per sample showing
 * the sample name (text-body font-semibold text-ink) and a value Input.
 * controls-don't-lie: values are plain inputs — no fake "type-and-it-saves"
 * affordance; the foot-gate "Confirm & build" is the single commit action.
 */
export function ColdAssignPanel({
  rows,
  variableKey,
  onKeyChange,
  onValueChange,
  intro,
  className,
}: ColdAssignPanelProps): JSX.Element {
  // Undefined keeps the default cold-corpus paragraph (call sites unchanged);
  // null is a deliberate suppression; any other node fills the same slot.
  const introNode =
    intro === undefined ? (
      <>
        These samples share no tag key yet. Name the variable (e.g. "lipid ratio") and assign each
        sample's value, then Confirm &amp; build.
      </>
    ) : (
      intro
    );
  return (
    <div data-testid="cold-assign-panel" className={cx("space-y-4", className)}>
      <div>
        <Kicker tone="accent" className="mb-1">
          Name the ordering variable
        </Kicker>
        {introNode !== null ? <p className="text-body text-ink-soft mb-2">{introNode}</p> : null}
        <Input
          value={variableKey}
          onValueChange={onKeyChange}
          placeholder="e.g. lipid ratio"
          aria-label="Ordering variable name"
          data-testid="cold-key-input"
        />
      </div>
      <div className="space-y-1">
        {rows.map((r) => (
          <div
            key={r.sampleId}
            data-testid="cold-assign-row"
            data-sample-id={r.sampleId}
            className="flex items-center gap-3 rounded border border-hair px-3 py-2"
          >
            <div className="flex-1 min-w-0">
              <div className="text-body font-semibold text-ink truncate">{r.sampleName}</div>
              <HintText className="font-mono">smp_{r.sampleId}</HintText>
            </div>
            <Input
              value={r.value}
              onValueChange={(v) => onValueChange(r.sampleId, v)}
              placeholder="value"
              aria-label={`Value for ${r.sampleName}`}
              inputSize="sm"
              className="w-28 flex-shrink-0"
            />
          </div>
        ))}
      </div>
    </div>
  );
}
