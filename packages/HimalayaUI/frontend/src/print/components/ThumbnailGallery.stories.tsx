import type { Meta, StoryObj } from "@storybook/react-vite";
import { ThumbnailGallery, type GalleryExposure } from "./ThumbnailGallery";
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

const meta = {
  title: "components/ThumbnailGallery",
  component: ThumbnailGallery,
  args: {
    exposures: EXPOSURES,
  },
} satisfies Meta<typeof ThumbnailGallery>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Default sheet strip — horizontal scroll, one representative, one rejected. */
export const Sheet: Story = {};

/** Loupe strip — centered row. */
export const Loupe: Story = {
  args: { variant: "loupe" },
};

/** One exposure selected via selectedId. */
export const Selected: Story = {
  args: { selectedId: 65 },
};
