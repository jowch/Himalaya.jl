import type { Meta, StoryObj } from "@storybook/react-vite";
import { SheetTable } from "./SheetTable";
import { SampleTableRow } from "./SampleTableRow";
import type { GalleryExposure } from "./ThumbnailGallery";
import type { Tag } from "../ui/tag";
import thumb37 from "../fixtures/thumbs/37.png?url";
import thumb65 from "../fixtures/thumbs/65.png?url";
import thumb66 from "../fixtures/thumbs/66.png?url";
import thumb67 from "../fixtures/thumbs/67.png?url";
import thumb93 from "../fixtures/thumbs/93.png?url";

// Status-shaped fixtures: strip flags and count cells derive from one truth
// (the pages derive both from exposure status via toGalleryExposures).
// kept 4 / total 5 / dropped 1.
const EXPOSURES: GalleryExposure[] = [
  { id: 37, src: thumb37, frameNo: "37", representative: true },
  { id: 65, src: thumb65, frameNo: "65" },
  { id: 66, src: thumb66, frameNo: "66", rejected: true },
  { id: 67, src: thumb67, frameNo: "67" },
  { id: 93, src: thumb93, frameNo: "93" },
];

const EXPOSURES_ALL_KEPT: GalleryExposure[] = EXPOSURES.map((e) => ({
  ...e,
  rejected: false,
}));

const EXPOSURES_THREE_KEPT: GalleryExposure[] = EXPOSURES.slice(0, 3).map(
  (e) => ({ ...e, rejected: false }),
);

const TAGS: Tag[] = [{ key: "LL37" }, { key: "temperature", value: "37C" }];

const meta = {
  title: "components/SheetTable",
  component: SheetTable,
  args: {},
} satisfies Meta<typeof SheetTable>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Full table: 4 rows — indexed with one dropped frame (Pn3m), all-kept indexed
 *  (Im3m), not-indexed (all kept), indexed (Ia3d, all kept). */
export const Sheet: Story = {
  parameters: { layout: "fullscreen" },
  render: () => (
    <div className="bg-paper p-6">
      <SheetTable>
        <SampleTableRow
          name="POPC + 20% chol"
          sampleId="s-001"
          phase="Pn3m"
          exposures={EXPOSURES}
          kept={4}
          total={5}
          tags={TAGS}
        />
        <SampleTableRow
          name="POPC + 40% chol"
          sampleId="s-002"
          phase="Im3m"
          exposures={EXPOSURES_ALL_KEPT}
          kept={5}
          total={5}
          tags={TAGS}
        />
        <SampleTableRow
          name="MO + buffer"
          sampleId="s-003"
          phase={null}
          exposures={EXPOSURES_THREE_KEPT}
          kept={3}
          total={3}
          tags={[]}
        />
        <SampleTableRow
          name="MO + PEG"
          sampleId="s-004"
          phase="Ia3d"
          exposures={EXPOSURES_THREE_KEPT}
          kept={3}
          total={3}
          tags={[{ key: "lipid", value: "DOPC" }]}
        />
      </SheetTable>
    </div>
  ),
};

/** Empty state: no children + an empty node. */
export const Empty: Story = {
  parameters: { layout: "fullscreen" },
  render: () => (
    <div className="bg-paper p-6">
      <SheetTable
        empty={
          <div className="flex items-center justify-center py-12">
            <span className="text-body-faint">No samples screened yet</span>
          </div>
        }
      />
    </div>
  ),
};
