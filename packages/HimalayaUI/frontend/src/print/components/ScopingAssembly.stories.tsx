import type { Meta, StoryObj } from "@storybook/react-vite";
import { useEffect, useState } from "react";
import { ScopePlate } from "./ScopePlate";
import { ScopeSampleRow } from "./ScopeSampleRow";
import { ScopeCandidateRow } from "./ScopeCandidateRow";
import { useDragReorder, reorder } from "./useDragReorder";
import type { PhaseSegment } from "../ui";
import { realTraces } from "../fixtures/realTraces";
import { buildFootState } from "../pages/scopingDerive";

/**
 * Page simulation (NOT a component): assembles `ScopePlate` with mapped
 * `ScopeSampleRow` / `ScopeCandidateRow` children into the series-scoping
 * worksheet as the mockup lays it out — the machine's grouping, the parsed
 * ordering values, the candidate it found-but-did-not-assume, the phase preview,
 * and the gated build action.
 *
 * This story owns the cross-component state ScopePlate refuses to hold: the
 * per-sample flagged map, the candidate→series fold, and an in-session undo
 * stack (⌘Z / the undo button). It DERIVES everything ScopePlate renders —
 * preview segments, count, foot state, the build gate — exactly as the future
 * Layer-4 page will. ScopePlate stays presentational.
 */

interface Member {
  id: string;
  name: string;
  value: string;
  /** Numeric sort key behind the ratio (low → high). */
  key: number;
  phase: string;
  coexistWith?: string[];
  flagged: boolean;
}

// Cycle the 5 real fixture traces across the members.
const TRACE_IDS = [37, 65, 66, 67, 93];
const traceFor = (i: number) => realTraces[TRACE_IDS[i % TRACE_IDS.length]!]!;

const INITIAL_SERIES: Member[] = [
  { id: "smp_04", name: "Lipid 1-2, no LL37", value: "1 : 0", key: 0, phase: "Pn3m", flagged: true },
  { id: "smp_07", name: "Lipid 1-2 + LL37 1:0.25", value: "1 : 0.25", key: 0.25, phase: "Pn3m", flagged: false },
  { id: "smp_09", name: "Lipid 1-2 + LL37 1:0.5", value: "1 : 0.5", key: 0.5, phase: "Pn3m", coexistWith: ["Lamellar"], flagged: false },
  { id: "smp_12", name: "Lipid 1-2 + LL37 1:1", value: "1 : 1", key: 1, phase: "Lamellar", coexistWith: ["Pn3m"], flagged: false },
  { id: "smp_15", name: "Lipid 1-2 + LL37 1:2", value: "1 : 2", key: 2, phase: "Lamellar", flagged: false },
  { id: "smp_18", name: "Lipid 1-2 + LL37 1:4", value: "1 : 4", key: 4, phase: "Lamellar", flagged: false },
];

const INITIAL_CANDIDATES: Member[] = [
  { id: "smp_21", name: "Lipid 1-1 + LL37 1:1", value: "1 : 1", key: 1, phase: "Pn3m", flagged: false },
];

const ORDER_OPTIONS = [
  { value: "LL37 : lipid ratio", label: "LL37 : lipid ratio" },
  { value: "Time", label: "Time" },
  { value: "Dose", label: "Dose" },
  { value: "Temperature", label: "Temperature" },
  { value: "Define your own…", label: "Define your own…" },
];

type HistoryEntry =
  | { type: "flag"; id: string; prev: boolean; label: string }
  | { type: "add"; id: string; label: string };

function ScopingView(): JSX.Element {
  const [series, setSeries] = useState<Member[]>(INITIAL_SERIES);
  const [candidates, setCandidates] = useState<Member[]>(INITIAL_CANDIDATES);
  const [history, setHistory] = useState<HistoryEntry[]>([]);
  const [orderedBy, setOrderedBy] = useState("LL37 : lipid ratio");
  // The displayed order is seeded value-sorted (low → high), but the grip is a
  // real MANUAL OVERRIDE: dragging a row rewrites this id list, and manual order
  // then wins (candidates append; we never re-sort).
  const [order, setOrder] = useState<string[]>(() =>
    [...INITIAL_SERIES].sort((a, b) => a.key - b.key).map((s) => s.id),
  );

  const { dragItemProps, dropEdge } = useDragReorder((from, to) =>
    setOrder((o) => reorder(o, from, to)),
  );

  // Keyboard reorder (SC-KBD) converges on the same setOrder(reorder) path as
  // drag; `reorder` no-ops on out-of-range boundaries. The real page also
  // announces moves to an SR live region — invisible in Storybook, so the
  // story wires the move only.
  const moveRow = (i: number, delta: -1 | 1): void => setOrder((o) => reorder(o, i, i + delta));

  // The shown members follow the page-owned manual `order` against the lookup.
  const byId = new Map(series.map((s) => [s.id, s]));
  const sorted = order.map((id) => byId.get(id)).filter((s): s is Member => s != null);

  // Clicking a value toggles whether this read is skipped from the batch
  // write. Recorded so it steps back with Undo / ⌘Z.
  const toggleFlag = (id: string): void => {
    const m = series.find((s) => s.id === id);
    if (!m) return;
    setHistory((h) => [
      ...h,
      { type: "flag", id, prev: m.flagged, label: (m.flagged ? "restored " : "skipped ") + id },
    ]);
    setSeries((cur) => cur.map((s) => (s.id === id ? { ...s, flagged: !s.flagged } : s)));
  };

  // Adding a candidate folds it into the series. Manual order wins, so it is
  // APPENDED to the displayed order rather than re-sorted in by key.
  const addCandidate = (id: string): void => {
    const c = candidates.find((x) => x.id === id);
    if (!c) return;
    setCandidates((cur) => cur.filter((x) => x.id !== id));
    setSeries((cur) => [...cur, c]);
    setOrder((o) => [...o, id]);
    setHistory((h) => [...h, { type: "add", id, label: "added " + id }]);
  };

  const undo = (): void => {
    setHistory((h) => {
      const e = h[h.length - 1];
      if (!e) return h;
      if (e.type === "flag") {
        setSeries((cur) => cur.map((s) => (s.id === e.id ? { ...s, flagged: e.prev } : s)));
      } else {
        setSeries((cur) => {
          const m = cur.find((s) => s.id === e.id);
          if (m) setCandidates((cs) => [...cs, m]);
          return cur.filter((s) => s.id !== e.id);
        });
        setOrder((o) => o.filter((id) => id !== e.id));
      }
      return h.slice(0, -1);
    });
  };

  // ⌘Z / Ctrl-Z steps back, mirroring the mockup.
  useEffect(() => {
    const onKey = (e: KeyboardEvent): void => {
      if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "z") {
        e.preventDefault();
        undo();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  });

  // Derivations the page owns — ScopePlate just renders them. The foot + gate
  // mirror the REAL page contract (scopingDerive.ts): a skipped (flagged) read
  // is excluded from the write, never blocks; warn only when nothing is kept.
  const n = sorted.length;
  const skipped = sorted.filter((s) => s.flagged).length;
  const kept = n - skipped;
  const preview: PhaseSegment[] = sorted.map((s) =>
    s.coexistWith ? { phase: s.phase, coexistWith: s.coexistWith } : { phase: s.phase },
  );
  const footState = buildFootState(kept, skipped);
  const lastLabel = history.length ? history[history.length - 1]!.label : undefined;

  return (
    <div className="bg-paper p-10 flex justify-center min-h-screen">
      <ScopePlate
        seriesName="LL37 titration of lipid 1-2"
        grouping={
          <>
            You selected <strong>{n} samples</strong> on the contact sheet. Himalaya grouped them from
            their names and read the order from the <strong>LL37 : lipid ratio</strong>. Click a value
            to skip that read from the write — one starts skipped here to show the exclusion look.
          </>
        }
        orderedBy={orderedBy}
        orderOptions={ORDER_OPTIONS}
        onOrderSelect={setOrderedBy}
        orderNote="Read from the sample names. Change it to time, dose, temperature, or define your own."
        count={`${n} samples · low to high`}
        {...(history.length ? { onUndo: undo, ...(lastLabel ? { undoLabel: `Step back: ${lastLabel}` } : {}) } : {})}
        rows={sorted.map((s, i) => {
          const dprops = dragItemProps(i);
          const edge = dropEdge(i);
          return (
            <div
              key={s.id}
              {...dprops}
              className={`relative cursor-grab${dprops["data-dragging"] ? " opacity-50" : ""}`}
            >
              {edge && (
                <span
                  aria-hidden="true"
                  className={`pointer-events-none absolute left-0 right-0 z-10 h-0.5 rounded-full bg-accent ${edge === "top" ? "-top-px" : "-bottom-px"}`}
                />
              )}
              <ScopeSampleRow
                name={s.name}
                sampleId={s.id}
                trace={traceFor(i)}
                phase={s.phase}
                value={s.value}
                {...(s.flagged ? { flagged: true } : {})}
                onToggleFlag={() => toggleFlag(s.id)}
                onMoveBy={(delta) => moveRow(i, delta)}
              />
            </div>
          );
        })}
        candidates={
          candidates.length ? (
            candidates.map((c, i) => (
              <ScopeCandidateRow
                key={c.id}
                name={c.name}
                why={
                  <>
                    has LL37, but the{" "}
                    <strong className="text-accent font-semibold">1-1 lipid line</strong> — its own series?
                  </>
                }
                trace={traceFor(i)}
                phase={c.phase}
                onAdd={() => addCandidate(c.id)}
              />
            ))
          ) : (
            <div className="text-meta text-ink-faint italic">
              Nothing else in the corpus matches this grouping.
            </div>
          )
        }
        preview={preview}
        footState={footState}
        footNote={
          <>
            Confirming records the LL37 : lipid ratio on every kept sample — the next series that
            needs it already knows.
          </>
        }
        {...(kept === 0 ? { buildDisabled: true } : {})}
        onBuild={() => undefined}
      />
    </div>
  );
}

const meta: Meta<typeof ScopePlate> = {
  title: "components/ScopingAssembly",
  component: ScopePlate,
};

export default meta;
type Story = StoryObj<typeof meta>;

export const Page: Story = { render: () => <ScopingView /> };
