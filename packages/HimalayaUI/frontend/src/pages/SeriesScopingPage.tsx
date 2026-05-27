import { useEffect, useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import {
  useCorpusSampleTags,
  useCorpusPickerSamples,
  useScopeSeries,
} from "../queries";
import { proposeOrdering, type OrderingRow } from "../lib/scoping/proposeOrdering";
import { ScopingConfirmModal } from "../components/ScopingConfirmModal";

// Static skeleton shape for boneyard's headless capture (docs/boneyard.md).
const SCOPING_FIXTURE = (
  <div className="space-y-2">
    {[0, 1, 2, 3].map((i) => (
      <div key={i} className="flex items-center gap-3 rounded border border-hair p-3">
        <div className="h-4 w-4 rounded bg-paper-sunk" />
        <div className="h-4 w-1/3 rounded bg-paper-sunk" />
        <div className="h-4 w-1/4 rounded bg-paper-sunk" />
      </div>
    ))}
  </div>
);

/**
 * SeriesScopingPage — the series scoping step at /series/new (#174 / I3.4).
 *
 * Per master plan §6.4 this is the confirm-and-build GATE that *writes* the
 * structured (key,value) sample_tags — NOT series creation. It reads corpus
 * tags to machine-propose an ordering variable and corpus picker samples as
 * member candidates, lets the user resolve flagged values + toggle which
 * members are included, and on confirm writes one batch (source='scoping')
 * through the queue, then navigates to the folio. I3.6 later upgrades this to
 * create + open the new series builder (the folio→scoping→builder stitch).
 */
export function SeriesScopingPage(): JSX.Element {
  const navigate = useNavigate();
  const tagsQ = useCorpusSampleTags();
  const pickerQ = useCorpusPickerSamples();
  const scopeSeries = useScopeSeries();

  const isLoading = tagsQ.isLoading || pickerQ.isLoading;
  const isError = tagsQ.isError || pickerQ.isError;

  const proposal = useMemo(
    () => proposeOrdering(tagsQ.data ?? [], pickerQ.data ?? []),
    [tagsQ.data, pickerQ.data],
  );

  // Local editable copy of the proposed rows — seeded from the proposal once
  // the queries resolve, then user-owned (value edits clear `flagged`, include
  // checkboxes toggle membership).
  const [rows, setRows] = useState<OrderingRow[]>([]);
  useEffect(() => {
    setRows(proposal.rows);
  }, [proposal.rows]);

  const [confirmOpen, setConfirmOpen] = useState(false);

  const included = rows.filter((r) => r.include && !r.flagged);
  // Build gate: every INCLUDED row must be non-flagged, and there must be at
  // least one included row to write.
  const canBuild =
    proposal.orderingKey !== undefined &&
    included.length > 0 &&
    rows.every((r) => !r.include || !r.flagged);

  function setValue(sampleId: number, value: string): void {
    setRows((prev) =>
      prev.map((r) =>
        r.sampleId === sampleId ? { ...r, value, flagged: value === "" } : r,
      ),
    );
  }
  function toggleInclude(sampleId: number): void {
    setRows((prev) =>
      prev.map((r) => (r.sampleId === sampleId ? { ...r, include: !r.include } : r)),
    );
  }

  function handleConfirm(): void {
    if (proposal.orderingKey === undefined) return;
    scopeSeries.mutate({
      key: proposal.orderingKey,
      tags: included.map((r) => ({ sampleId: r.sampleId, value: r.value })),
    });
    setConfirmOpen(false);
    // D1: land on the folio (exists since I3.3). I3.6 upgrades this to create
    // + open the new series builder (folio→scoping→builder stitch).
    navigate("/series");
  }

  return (
    <div data-testid="scoping-page" className="mx-auto flex max-w-[900px] flex-col gap-5 px-8 py-7">
      <header className="flex flex-col gap-1">
        <div className="text-xs font-semibold uppercase tracking-wide text-print-accent">
          New series — scoping
        </div>
        <p className="text-sm text-ink-faint">
          Ordering variable:{" "}
          <span data-testid="ordering-key" className="font-mono text-ink">
            {proposal.orderingKey ?? "—"}
          </span>
        </p>
      </header>

      {isError ? (
        <div data-testid="scoping-error" className="px-4 py-8 text-sm text-ink-soft">
          Could not load corpus tags. Try reloading the page.
        </div>
      ) : (
        <Skeleton
          name="scoping"
          className="w-full"
          loading={isLoading}
          stagger={50}
          transition={200}
          fixture={SCOPING_FIXTURE}
        >
          {rows.length === 0 ? (
            <div data-testid="scoping-empty" className="px-4 py-12 text-center text-sm text-ink-faint">
              No samples in the corpus to scope.
            </div>
          ) : (
            <div className="flex flex-col gap-2">
              {rows.map((r) => (
                <div
                  key={r.sampleId}
                  data-testid={`scoping-row-${r.sampleId}`}
                  data-flagged={r.flagged ? "true" : undefined}
                  className="flex items-center gap-3 rounded border border-hair px-3 py-2"
                >
                  <input
                    type="checkbox"
                    data-testid={`scoping-include-${r.sampleId}`}
                    aria-label={`Include ${r.sampleName}`}
                    checked={r.include}
                    onChange={() => toggleInclude(r.sampleId)}
                  />
                  <span className="flex-1 text-sm text-ink">{r.sampleName}</span>
                  <input
                    type="text"
                    data-testid={`scoping-value-${r.sampleId}`}
                    aria-label={`Value for ${r.sampleName}`}
                    value={r.value}
                    onChange={(e) => setValue(r.sampleId, e.target.value)}
                    placeholder={r.flagged ? "set a value" : undefined}
                    className={`w-28 rounded border px-2 py-1 text-sm ${
                      r.flagged ? "border-print-accent" : "border-hair"
                    }`}
                  />
                </div>
              ))}
            </div>
          )}
        </Skeleton>
      )}

      <div className="flex items-center justify-end gap-3">
        <span className="text-xs text-ink-faint">
          {included.length} of {rows.length} included
        </span>
        <button
          type="button"
          data-testid="scoping-open-confirm"
          disabled={!canBuild}
          onClick={() => setConfirmOpen(true)}
          className="rounded bg-accent px-3 py-1.5 text-sm font-semibold text-paper disabled:opacity-40"
          title={canBuild ? undefined : "Resolve every included flagged value to build"}
        >
          Confirm &amp; build →
        </button>
      </div>

      <ScopingConfirmModal
        open={confirmOpen}
        orderingKey={proposal.orderingKey}
        count={included.length}
        onConfirm={handleConfirm}
        onClose={() => setConfirmOpen(false)}
      />
    </div>
  );
}
