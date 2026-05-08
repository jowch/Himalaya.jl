import { useMemo } from "react";
import type { RefObject } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
import { usePickerSamples } from "../queries";
import { ComparisonPickerBody, type Pick } from "./ComparisonPickerBody";

interface Props {
  experimentId: number | undefined;
  /** External ref the parent uses to focus the search input from a sibling CTA. */
  searchInputRef?: RefObject<HTMLInputElement>;
}

export function ComparisonPickerPanel({
  experimentId, searchInputRef,
}: Props): JSX.Element {
  const qc = useQueryClient();
  const draft = useAppState((s) => s.activeDraft);
  const addMember = useAppState((s) => s.addMember);
  const removeMember = useAppState((s) => s.removeMember);

  // Same query the body fetches — TanStack dedupes, no extra HTTP. Used here
  // to map exposure_id → sample_id so we can swap an existing same-sample
  // member when the user picks a different exposure for the same sample.
  const pickerQ = usePickerSamples(experimentId);

  const exposureToSample = useMemo(() => {
    const m = new Map<number, number>();
    for (const r of pickerQ.data ?? []) {
      for (const e of r.all_exposures) m.set(e.id, r.sample.id);
    }
    return m;
  }, [pickerQ.data]);

  const alreadyAddedExposureIds = useMemo(
    () => new Set(
      (draft?.members ?? [])
        .map((m) => m.exposure_id)
        .filter((id): id is number => id !== null),
    ),
    [draft],
  );

  const handlePick = (p: Pick): void => {
    if (alreadyAddedExposureIds.has(p.exposure_id)) return;

    // Swap semantics: if a different exposure of the same sample is already a
    // member, remove it before adding the new one. Iterates back-to-front so
    // the index stays valid for the splice in `removeMember(index)`.
    const members = draft?.members ?? [];
    for (let i = members.length - 1; i >= 0; i--) {
      const eid = members[i]?.exposure_id;
      if (eid === null || eid === undefined || eid === p.exposure_id) continue;
      if (exposureToSample.get(eid) === p.sample_id) {
        removeMember(i);
        break;   // only one same-sample member at a time
      }
    }

    addMember(p.exposure_id, qc);
  };

  return (
    <div
      data-testid="comparison-picker-panel"
      className="flex flex-col h-full overflow-hidden"
    >
      <div className="px-4 py-3 border-b border-border">
        <h3 className="text-sm font-medium text-fg">Add traces</h3>
      </div>
      <ComparisonPickerBody
        experimentId={experimentId}
        picks={[]}                        // immediate mode: ignored
        onPicksChange={() => {}}          // immediate mode: ignored
        onPick={handlePick}
        alreadyAddedExposureIds={alreadyAddedExposureIds}
        {...(searchInputRef ? { searchInputRef } : {})}
      />
    </div>
  );
}
