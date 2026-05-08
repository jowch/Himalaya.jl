import type { RefObject } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
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

  const alreadyAddedExposureIds = new Set(
    (draft?.members ?? [])
      .map((m) => m.exposure_id)
      .filter((id): id is number => id !== null),
  );

  const handlePick = (p: Pick): void => {
    if (alreadyAddedExposureIds.has(p.exposure_id)) return;
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
