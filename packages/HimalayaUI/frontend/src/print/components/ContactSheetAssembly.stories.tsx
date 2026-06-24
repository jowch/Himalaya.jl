import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { SheetTable } from "./SheetTable";
import { SampleTableRow } from "./SampleTableRow";
import type { GalleryExposure } from "./ThumbnailGallery";
import { Kicker, ProgressBar } from "../ui";
import type { Tag } from "../ui";
import thumb37 from "../fixtures/thumbs/37.png?url";
import thumb65 from "../fixtures/thumbs/65.png?url";
import thumb66 from "../fixtures/thumbs/66.png?url";
import thumb67 from "../fixtures/thumbs/67.png?url";
import thumb93 from "../fixtures/thumbs/93.png?url";

/**
 * Page simulation (NOT a component): assembles `SheetTable` + `SampleTableRow`
 * children + `CullBar` into the contact-sheet view as the mockup lays it out —
 * a flat sample list with batch-cull selection, a screened-progress bar in the
 * page head, and a floating dark CullBar that appears on first selection.
 *
 * The Layer-4 contact-sheet page (plate shell, nav, back button) is deferred;
 * this story owns the cross-component state (selected exposure ids, screened
 * count) that the page will own.
 */

interface SampleDatum {
  sampleId: string;
  name: string;
  screened: boolean;
  phase: string | null;
  kept: number;
  total: number;
  dropped?: number;
  exposures: GalleryExposure[];
  tags: Tag[];
}

const TAGS_A: Tag[] = [{ key: "LL37" }, { key: "temperature", value: "37C" }];
const TAGS_B: Tag[] = [{ key: "DOPC" }, { key: "chol", value: "20%" }];
const TAGS_C: Tag[] = [{ key: "MO" }, { key: "buffer", value: "PBS" }];
const TAGS_D: Tag[] = [{ key: "DPPE" }];
const TAGS_E: Tag[] = [];

const SAMPLES: SampleDatum[] = [
  {
    sampleId: "s-001",
    name: "POPC + 20% cholesterol",
    screened: true,
    phase: "Pn3m",
    kept: 4,
    total: 5,
    dropped: 1,
    exposures: [
      { id: 101, src: thumb37, frameNo: "37", representative: true, kept: true },
      { id: 102, src: thumb65, frameNo: "65", kept: true },
      { id: 103, src: thumb66, frameNo: "66", rejected: true },
      { id: 104, src: thumb67, frameNo: "67", kept: true },
      { id: 105, src: thumb93, frameNo: "93", kept: true },
    ],
    tags: TAGS_A,
  },
  {
    sampleId: "s-002",
    name: "POPC + 40% cholesterol",
    screened: true,
    phase: "Im3m",
    kept: 3,
    total: 3,
    exposures: [
      { id: 201, src: thumb65, frameNo: "65", representative: true, kept: true },
      { id: 202, src: thumb66, frameNo: "66", kept: true },
      { id: 203, src: thumb67, frameNo: "67", kept: true },
    ],
    tags: TAGS_B,
  },
  {
    sampleId: "s-003",
    name: "MO + buffer (unscreened)",
    screened: false,
    phase: null,
    kept: 0,
    total: 4,
    dropped: 0,
    exposures: [
      { id: 301, src: thumb93, frameNo: "93" },
      { id: 302, src: thumb37, frameNo: "37" },
      { id: 303, src: thumb66, frameNo: "66" },
      { id: 304, src: thumb67, frameNo: "67" },
    ],
    tags: TAGS_C,
  },
  {
    sampleId: "s-004",
    name: "DPPE monolayer",
    screened: false,
    phase: null,
    kept: 0,
    total: 2,
    exposures: [
      { id: 401, src: thumb65, frameNo: "65" },
      { id: 402, src: thumb93, frameNo: "93" },
    ],
    tags: TAGS_D,
  },
  {
    sampleId: "s-005",
    name: "Ia3d reference",
    screened: true,
    phase: "Ia3d",
    kept: 5,
    total: 5,
    exposures: [
      { id: 501, src: thumb37, frameNo: "37", representative: true, kept: true },
      { id: 502, src: thumb65, frameNo: "65", kept: true },
      { id: 503, src: thumb66, frameNo: "66", kept: true },
      { id: 504, src: thumb67, frameNo: "67", kept: true },
      { id: 505, src: thumb93, frameNo: "93", kept: true },
    ],
    tags: TAGS_E,
  },
];

function ContactSheetView(): JSX.Element {
  // selected = set of selected exposure ids across ALL rows (drives CullBar).
  const [selected, setSelected] = useState<Set<number>>(new Set());

  const screenedCount = SAMPLES.filter((s) => s.screened).length;

  // Immutably toggle an exposure id in/out of the selection set.
  const toggle = (_sample: SampleDatum, id: number): void => {
    setSelected((prev) => {
      const next = new Set(prev);
      if (next.has(id)) {
        next.delete(id);
      } else {
        next.add(id);
      }
      return next;
    });
  };

  return (
    <div className="bg-paper p-6">
      <div className="max-w-[1240px] mx-auto">
        {/* Page head: title/sub on the left, screened-count + progress on the right */}
        <div className="flex items-end justify-between mb-5">
          <div>
            <Kicker tone="accent">Contact sheet</Kicker>
            <p className="text-body text-ink-soft mt-1">
              Dark frames on light paper — screen, tag, and cull before indexing.
            </p>
          </div>
          <div className="text-right">
            <p className="text-body text-ink-soft">
              {screenedCount} / {SAMPLES.length} screened
            </p>
            <ProgressBar
              value={screenedCount}
              total={SAMPLES.length}
              label="samples screened"
              className="w-40 mt-1"
            />
          </div>
        </div>

        {/* Sample table */}
        <SheetTable
          empty={
            <div className="p-10 text-center text-ink-faint">
              No samples screened yet
            </div>
          }
        >
          {SAMPLES.map((s) => {
            // Cull model is MULTI-select: hand the whole selection set down. The
            // gallery only renders THIS row's exposures, so only its own members
            // light up — every selected frame, not just the first.
            return (
              <SampleTableRow
                key={s.sampleId}
                name={s.name}
                sampleId={s.sampleId}
                {...(s.screened !== undefined ? { screened: s.screened } : {})}
                exposures={s.exposures}
                selectedExposureIds={selected}
                onSelectExposure={(id) => toggle(s, id)}
                kept={s.kept}
                total={s.total}
                {...(s.dropped != null ? { dropped: s.dropped } : {})}
                tags={s.tags}
                // phase: pass only when defined — StatusCell reads undefined as
                // "not yet indexed" whereas null = "screened but no phase"
                {...(s.phase !== undefined ? { phase: s.phase } : {})}
              />
            );
          })}
        </SheetTable>

      </div>
    </div>
  );
}

const meta: Meta<typeof SheetTable> = {
  title: "components/ContactSheetAssembly",
  component: SheetTable,
};

export default meta;
type Story = StoryObj<typeof meta>;

export const Page: Story = { render: () => <ContactSheetView /> };
