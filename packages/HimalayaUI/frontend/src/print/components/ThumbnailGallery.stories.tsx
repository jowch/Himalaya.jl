import type { Meta, StoryObj } from "@storybook/react-vite";
import { ThumbnailGallery, type GalleryExposure } from "./ThumbnailGallery";
import thumb37 from "../fixtures/thumbs/37.png?url";
import thumb65 from "../fixtures/thumbs/65.png?url";
import thumb66 from "../fixtures/thumbs/66.png?url";
import thumb67 from "../fixtures/thumbs/67.png?url";
import thumb93 from "../fixtures/thumbs/93.png?url";

// Status-shaped fixtures, derived the same way the pages derive (binary:
// rejected = dropped, null = kept/default — only the dropped channel is marked).
const FIXTURES = [
  { id: 37, src: thumb37, frameNo: "37", representative: true, status: null },
  { id: 65, src: thumb65, frameNo: "65", representative: false, status: null },
  { id: 66, src: thumb66, frameNo: "66", representative: false, status: "rejected" },
  { id: 67, src: thumb67, frameNo: "67", representative: false, status: null },
  { id: 93, src: thumb93, frameNo: "93", representative: false, status: null },
] as const;

const EXPOSURES: GalleryExposure[] = FIXTURES.map((f) => ({
  id: f.id,
  src: f.src,
  frameNo: f.frameNo,
  representative: f.representative,
  rejected: f.status === "rejected",
}));

const meta = {
  title: "components/ThumbnailGallery",
  component: ThumbnailGallery,
  args: {
    exposures: EXPOSURES,
  },
} satisfies Meta<typeof ThumbnailGallery>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Default sheet strip — horizontal scroll: four kept (one of them the
 *  representative; kept frames carry no marker) and one dropped. */
export const Sheet: Story = {};

/** Loupe strip — larger thumbnails, centered row. */
export const Loupe: Story = {
  args: { size: "lg", align: "center" },
};

/** One exposure selected via selectedId. */
export const Selected: Story = {
  args: { selectedId: 65 },
};

const SRC_CYCLE = [thumb37, thumb65, thumb66, thumb67, thumb93];
const MANY_EXPOSURES: GalleryExposure[] = Array.from(
  { length: 12 },
  (_, i): GalleryExposure => ({
    id: 100 + i,
    src: SRC_CYCLE[i % SRC_CYCLE.length]!,
    frameNo: String(100 + i),
    ...(i === 0 ? { representative: true } : {}),
    ...(i === 3 ? { rejected: true } : {}),
  }),
);

/** Many exposures inside a fixed-width wrapper — the strip stays one row and
 *  scrolls horizontally instead of widening its container. */
export const Overflowing: Story = {
  args: { exposures: MANY_EXPOSURES },
  render: (args) => (
    <div style={{ width: 360 }} className="border border-hair">
      <ThumbnailGallery {...args} className="w-full" />
    </div>
  ),
};
