import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { BuilderRail } from "./BuilderRail";
import { SeriesPlate, type SeriesScale } from "./SeriesPlate";
import { MemberRow } from "./MemberRow";
import { RailBack } from "./RailBack";
import { Dock } from "./Dock";
import { useDragReorder, reorder } from "./useDragReorder";
import { TRANSITION } from "../waterfall/waterfall.fixtures";

/**
 * Page simulation (NOT a component): assembles the figure (`SeriesPlate`) on the
 * left and the `BuilderRail` editing rail on the right into the series-builder
 * view. The page owns the cross-component state the components are forbidden to
 * hold — offset, scale (the figure↔rail link is the SAME state), and the
 * collapsed/full-bleed mode.
 *
 * When collapsed the rail is UNMOUNTED (the figure goes full-bleed) and the
 * floating `RailBack` tab + `Dock` offset control take over — mirroring the
 * mockup's full-bleed reading mode.
 *
 * The Layer-4 builder page (plate shell, nav, top bar) is deferred; this story
 * only simulates the page's state ownership. Rows reuse the shared
 * `TRANSITION` WaterfallRow fixture — the exact construction SeriesPlate.stories
 * uses (`toWaterfallRows(transitionSeries, realTraces)`).
 */

interface TraceDatum {
  id: string;
  name: string;
  dose: string;
  phase: string;
  coexistWith?: string[];
}

// Mirrors the mockup SERIES: six members ordered by LL37 : lipid ratio, with the
// two middle members in coexistence.
const TRACES: TraceDatum[] = [
  { id: "smp_04", name: "Pn3m", dose: "1 : 0", phase: "Pn3m" },
  { id: "smp_08", name: "Pn3m", dose: "1 : 1", phase: "Pn3m" },
  { id: "smp_11", name: "Pn3m", dose: "1 : 2", phase: "Pn3m", coexistWith: ["Im3m"] },
  { id: "smp_14", name: "Im3m", dose: "1 : 3", phase: "Im3m", coexistWith: ["Pn3m"] },
  { id: "smp_16", name: "Lamellar", dose: "1 : 3.5", phase: "Lamellar" },
  { id: "smp_18", name: "Lamellar", dose: "1 : 4", phase: "Lamellar" },
];

const ORDER_OPTIONS = [
  { value: "LL37 : lipid ratio", label: "LL37 : lipid ratio" },
  { value: "Time", label: "Time" },
  { value: "Dose", label: "Dose" },
  { value: "Temperature", label: "Temperature" },
  { value: "Define your own…", label: "Define your own…" },
];

function SeriesBuilderView(): JSX.Element {
  const [offset, setOffset] = useState(1.2);
  const [scale, setScale] = useState<SeriesScale>("log");
  const [collapsed, setCollapsed] = useState(false);
  const [orderedBy, setOrderedBy] = useState("LL37 : lipid ratio");
  const [traceOrder, setTraceOrder] = useState<TraceDatum[]>(TRACES);

  // The rail's "Traces — drag to reorder" label is honest: dragging a row
  // rewrites the page-owned `traceOrder`.
  const { dragItemProps } = useDragReorder((from, to) =>
    setTraceOrder((o) => reorder(o, from, to)),
  );

  const traces = traceOrder.map((t, i) => {
    const props = dragItemProps(i);
    return (
      <div key={t.id} {...props} className={`cursor-grab${props["data-dragging"] ? " opacity-50" : ""}`}>
        <MemberRow
          name={t.name}
          sub={`${t.id} · ${t.dose}`}
          phase={t.phase}
          {...(t.coexistWith ? { coexistWith: t.coexistWith } : {})}
        />
      </div>
    );
  });

  return (
    <div className="bg-paper min-h-screen">
      <div className="flex">
        <div className="flex-1 min-w-0 p-6">
          <SeriesPlate
            kicker="Series"
            title="LL37 titration of lipid 1-2"
            subtitle="6 exposures · variable: LL37 : lipid · SSRL Apr 2026"
            rows={TRANSITION}
            offsetScale={offset}
            scale={scale}
            onScaleChange={setScale}
            legendPhases={["Pn3m", "Im3m", "Lamellar"]}
            footHint="peaks are light anchors — hover a trace to read its indexing"
            footNote={`offset ×${offset.toFixed(2)} · ${scale === "log" ? "log" : "linear"} q`}
          />
        </div>

        {!collapsed && (
          <div className="w-[336px] shrink-0">
            <BuilderRail
              grouping={
                <>
                  Himalaya read <strong>6 samples</strong> as one series from their
                  names, ordered by <strong>LL37 : lipid ratio</strong>.
                </>
              }
              orderedBy={orderedBy}
              orderOptions={ORDER_OPTIONS}
              onOrderSelect={setOrderedBy}
              orderNote="Confirmed once, then reused by every series that needs it."
              offset={offset}
              onOffsetChange={setOffset}
              scale={scale}
              onScaleChange={setScale}
              traces={traces}
              onCollapse={() => setCollapsed(true)}
            />
          </div>
        )}
      </div>

      {collapsed && (
        <>
          <RailBack onClick={() => setCollapsed(false)} />
          <Dock offset={offset} onOffsetChange={setOffset} />
        </>
      )}
    </div>
  );
}

const meta: Meta<typeof BuilderRail> = {
  title: "components/SeriesBuilderAssembly",
  component: BuilderRail,
};

export default meta;
type Story = StoryObj<typeof meta>;

export const Page: Story = { render: () => <SeriesBuilderView /> };
