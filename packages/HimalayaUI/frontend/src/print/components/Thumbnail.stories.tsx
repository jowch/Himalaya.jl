import type { Meta, StoryObj } from "@storybook/react-vite";
import { Thumbnail } from "./Thumbnail";
import thumb65 from "../fixtures/thumbs/65.png?url";

const meta = {
  title: "components/Thumbnail",
  component: Thumbnail,
  args: {
    src: thumb65,
    frameNo: "65",
  },
} satisfies Meta<typeof Thumbnail>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Normal: Story = {};

export const Representative: Story = {
  args: { representative: true },
};

export const Rejected: Story = {
  args: { rejected: true },
};

export const Selected: Story = {
  args: { selected: true },
};

export const RepresentativeAndSelected: Story = {
  args: { representative: true, selected: true },
};

export const Loupe: Story = {
  args: { size: "loupe" },
};

export const NoImage: Story = {
  args: { src: null, frameNo: "66" },
};

export const AllStates: Story = {
  args: { representative: true, rejected: true, selected: true },
};
