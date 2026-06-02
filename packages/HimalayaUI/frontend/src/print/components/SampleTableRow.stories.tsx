import type { Meta, StoryObj } from "@storybook/react-vite";
import { SampleTableRow, SAMPLE_TABLE_COLS } from "./SampleTableRow";
import type { GalleryExposure } from "./ThumbnailGallery";
import type { Tag } from "../ui/tag";
import thumb37 from "../fixtures/thumbs/37.png?url";
import thumb65 from "../fixtures/thumbs/65.png?url";
import thumb66 from "../fixtures/thumbs/66.png?url";
import thumb67 from "../fixtures/thumbs/67.png?url";
import thumb93 from "../fixtures/thumbs/93.png?url";

const EXPOSURES: GalleryExposure[] = [
  { id: 37, src: thumb37, frameNo: "37", representative: true },
  { id: 65, src: thumb65, frameNo: "65" },
  { id: 66, src: thumb66, frameNo: "66", rejected: true },
  { id: 67, src: thumb67, frameNo: "67" },
  { id: 93, src: thumb93, frameNo: "93" },
];

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

/** Not indexed and NOT screened — shows the resting tint + "Not indexed". */
export const Unindexed: Story = {
  args: {
    screened: false,
    phase: null,
    dropped: 1,
  },
};

/** dropped=0 → no dropped callout. */
export const AllKept: Story = {
  args: {
    screened: true,
    phase: "Im3m",
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
    kept: 4,
    total: 5,
    dropped: 1,
  },
  {
    name: "POPC + 40% chol",
    sampleId: "s-002",
    screened: true,
    phase: "Im3m" as const,
    kept: 5,
    total: 5,
    dropped: 0,
  },
  {
    name: "MO + buffer",
    sampleId: "s-003",
    screened: false,
    phase: null,
    kept: 2,
    total: 3,
    dropped: 1,
  },
  {
    name: "MO + PEG",
    sampleId: "s-004",
    screened: false,
    phase: "Ia3d" as const,
    kept: 3,
    total: 3,
    dropped: 0,
  },
];

const SRC_CYCLE = [thumb37, thumb65, thumb66, thumb67, thumb93];
const MANY_EXPOSURES: GalleryExposure[] = Array.from(
  { length: 12 },
  (_, i): GalleryExposure => ({
    id: 200 + i,
    src: SRC_CYCLE[i % SRC_CYCLE.length]!,
    frameNo: String(200 + i),
    ...(i === 0 ? { representative: true } : {}),
    ...(i === 4 ? { rejected: true } : {}),
  }),
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
          exposures={EXPOSURES}
          kept={row.kept}
          total={row.total}
          dropped={row.dropped}
          tags={TAGS}
        />
      ))}
    </div>
  ),
};
