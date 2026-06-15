import type { Meta, StoryObj } from "@storybook/react-vite";
import { SampleTableRow, SAMPLE_TABLE_COLS } from "./SampleTableRow";
import type { GalleryExposure } from "./ThumbnailGallery";
import type { Tag } from "../ui/tag";
import thumb37 from "../fixtures/thumbs/37.png?url";
import thumb65 from "../fixtures/thumbs/65.png?url";
import thumb66 from "../fixtures/thumbs/66.png?url";
import thumb67 from "../fixtures/thumbs/67.png?url";
import thumb93 from "../fixtures/thumbs/93.png?url";

// Status-shaped fixtures: the strip's kept/rejected flags and the row's count
// cells must tell the same story the real pages derive together from exposure
// status (toGalleryExposures / toSampleRowModel).
type FrameStatus = "accepted" | "rejected" | null;

function expo(
  id: number,
  src: string,
  status: FrameStatus,
  representative = false,
): GalleryExposure {
  return {
    id,
    src,
    frameNo: String(id),
    representative,
    rejected: status === "rejected",
    kept: status === "accepted",
  };
}

// kept 4 / total 5 / dropped 1 — matches the meta count cells.
const EXPOSURES: GalleryExposure[] = [
  expo(37, thumb37, "accepted", true),
  expo(65, thumb65, "accepted"),
  expo(66, thumb66, "rejected"),
  expo(67, thumb67, "accepted"),
  expo(93, thumb93, "accepted"),
];

const EXPOSURES_ALL_KEPT: GalleryExposure[] = EXPOSURES.map((e) => ({
  ...e,
  rejected: false,
  kept: true,
}));

const EXPOSURES_UNSCREENED: GalleryExposure[] = EXPOSURES.map((e) => ({
  ...e,
  rejected: false,
  kept: false,
}));

const TAGS: Tag[] = [{ key: "LL37" }, { key: "temperature", value: "37C" }];

const meta = {
  title: "components/SampleTableRow",
  component: SampleTableRow,
  args: {
    name: "Sample A",
    sampleId: "s-001",
    exposures: EXPOSURES,
    kept: 4,
    total: 5,
    tags: TAGS,
  },
} satisfies Meta<typeof SampleTableRow>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Full row: indexed phase, screened, one dropped, representative + rejected thumbs. */
export const Indexed: Story = {
  args: {
    screened: true,
    phase: "Pn3m",
    dropped: 1,
  },
};

/** Not indexed and NOT screened — resting tint + "Not indexed"; every frame
 *  unscreened, so the strip shows no kept dots and the count cell says 0. */
export const Unindexed: Story = {
  args: {
    screened: false,
    phase: null,
    exposures: EXPOSURES_UNSCREENED,
    kept: 0,
    dropped: 0,
  },
};

/** dropped=0 → no dropped callout. */
export const AllKept: Story = {
  args: {
    screened: true,
    phase: "Im3m",
    exposures: EXPOSURES_ALL_KEPT,
    kept: 5,
    total: 5,
    dropped: 0,
  },
};

const STACK_ROWS = [
  {
    name: "POPC + 20% chol",
    sampleId: "s-001",
    screened: true,
    phase: "Pn3m" as const,
    exposures: EXPOSURES,
    kept: 4,
    total: 5,
    dropped: 1,
  },
  {
    name: "POPC + 40% chol",
    sampleId: "s-002",
    screened: true,
    phase: "Im3m" as const,
    exposures: EXPOSURES_ALL_KEPT,
    kept: 5,
    total: 5,
    dropped: 0,
  },
  {
    name: "MO + buffer",
    sampleId: "s-003",
    screened: false,
    phase: null,
    exposures: EXPOSURES_UNSCREENED.slice(0, 3),
    kept: 0,
    total: 3,
    dropped: 0,
  },
  {
    name: "MO + PEG",
    sampleId: "s-004",
    screened: false,
    phase: "Ia3d" as const,
    exposures: EXPOSURES_UNSCREENED.slice(0, 3),
    kept: 0,
    total: 3,
    dropped: 0,
  },
];

const SRC_CYCLE = [thumb37, thumb65, thumb66, thumb67, thumb93];
// kept 10 / total 12 / dropped 2 — matches the Overflowing count cells.
const MANY_EXPOSURES: GalleryExposure[] = Array.from(
  { length: 12 },
  (_, i): GalleryExposure =>
    expo(
      200 + i,
      SRC_CYCLE[i % SRC_CYCLE.length]!,
      i === 4 || i === 7 ? "rejected" : "accepted",
      i === 0,
    ),
);

const MANY_TAGS: Tag[] = [
  { key: "LL37" },
  { key: "temperature", value: "37C" },
  { key: "buffer", value: "PBS" },
  { key: "lipid", value: "DOPC" },
  { key: "pH", value: "7.4" },
  { key: "chol", value: "20%" },
];

/** Real row context with ~12 exposures and ~6 tags — the exposures cell scrolls
 *  horizontally within its column and the tags cell shows the "+N" overflow. */
export const Overflowing: Story = {
  parameters: { layout: "fullscreen" },
  render: () => (
    <div className="bg-plate min-w-[1000px]">
      <SampleTableRow
        name="POPC + 20% chol (long batch)"
        sampleId="s-999"
        screened
        phase="Pn3m"
        exposures={MANY_EXPOSURES}
        kept={10}
        total={12}
        dropped={2}
        tags={MANY_TAGS}
      />
    </div>
  ),
};

const HEADER_LABELS = ["SAMPLE", "EXPOSURES", "KEPT", "TAGS", "STATUS"];

/** Header row (story-local scaffolding) + 4 rows on the SAME grid track —
 *  proves column alignment. */
export const Stack: Story = {
  parameters: { layout: "fullscreen" },
  render: () => (
    <div className="bg-plate min-w-[1000px]">
      <div
        className="grid border-b border-hair-strong"
        style={{ gridTemplateColumns: SAMPLE_TABLE_COLS }}
      >
        {HEADER_LABELS.map((label) => (
          <div key={label} className="px-4 py-2">
            <span className="text-kicker text-kicker-faint">{label}</span>
          </div>
        ))}
      </div>
      {STACK_ROWS.map((row) => (
        <SampleTableRow
          key={row.sampleId}
          name={row.name}
          sampleId={row.sampleId}
          screened={row.screened}
          phase={row.phase}
          exposures={row.exposures}
          kept={row.kept}
          total={row.total}
          dropped={row.dropped}
          tags={TAGS}
        />
      ))}
    </div>
  ),
};
