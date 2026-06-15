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
  { id: 37, src: thumb37, frameNo: "37", representative: true, kept: true },
  { id: 65, src: thumb65, frameNo: "65", kept: true },
  { id: 66, src: thumb66, frameNo: "66", rejected: true },
  { id: 67, src: thumb67, frameNo: "67", kept: true },
  { id: 93, src: thumb93, frameNo: "93", kept: true },
];

const EXPOSURES_ALL_KEPT: GalleryExposure[] = EXPOSURES.map((e) => ({
  ...e,
  rejected: false,
  kept: true,
}));

const EXPOSURES_UNSCREENED: GalleryExposure[] = EXPOSURES.slice(0, 3).map(
  (e) => ({ ...e, rejected: false, kept: false }),
);

const TAGS: Tag[] = [{ key: "LL37" }, { key: "temperature", value: "37C" }];

const meta = {
  title: "components/SheetTable",
  component: SheetTable,
  args: {},
} satisfies Meta<typeof SheetTable>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Full table: 4 rows — screened+indexed (Pn3m, 1 dropped), screened+indexed
 *  all-kept (Im3m), unscreened+not-indexed, unscreened+indexed (Ia3d). */
export const Sheet: Story = {
  parameters: { layout: "fullscreen" },
  render: () => (
    <div className="bg-paper p-6">
      <SheetTable>
        <SampleTableRow
          name="POPC + 20% chol"
          sampleId="s-001"
          screened
          phase="Pn3m"
          exposures={EXPOSURES}
          kept={4}
          total={5}
          dropped={1}
          tags={TAGS}
        />
        <SampleTableRow
          name="POPC + 40% chol"
          sampleId="s-002"
          screened
          phase="Im3m"
          exposures={EXPOSURES_ALL_KEPT}
          kept={5}
          total={5}
          dropped={0}
          tags={TAGS}
        />
        <SampleTableRow
          name="MO + buffer"
          sampleId="s-003"
          screened={false}
          phase={null}
          exposures={EXPOSURES_UNSCREENED}
          kept={0}
          total={3}
          dropped={0}
          tags={[]}
        />
        <SampleTableRow
          name="MO + PEG"
          sampleId="s-004"
          screened={false}
          phase="Ia3d"
          exposures={EXPOSURES_UNSCREENED}
          kept={0}
          total={3}
          dropped={0}
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
