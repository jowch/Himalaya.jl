import type { Meta, StoryObj } from "@storybook/react-vite";
import type { ReactNode } from "react";
import { Thumbnail } from "./Thumbnail";
import thumb37 from "../fixtures/thumbs/37.png?url";
import thumb65 from "../fixtures/thumbs/65.png?url";
import thumb66 from "../fixtures/thumbs/66.png?url";
import thumb67 from "../fixtures/thumbs/67.png?url";

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

export const Large: Story = {
  args: { size: "lg" },
};

export const NoImage: Story = {
  args: { src: null, frameNo: "66" },
};

/** The state language at a glance: each state in its own labelled column
 *  (all `size="sm"`), then one `size="lg"` example for scale. */
export const AllStates: Story = {
  render: () => {
    const Cell = ({
      label,
      children,
    }: {
      label: string;
      children: ReactNode;
    }) => (
      <div className="flex flex-col items-center gap-1.5">
        {children}
        <span className="text-kicker text-kicker-faint">{label}</span>
      </div>
    );
    return (
      <div className="flex flex-wrap items-start gap-5">
        <Cell label="Default">
          <Thumbnail src={thumb65} frameNo="65" size="sm" />
        </Cell>
        <Cell label="Representative">
          <Thumbnail src={thumb37} frameNo="37" representative size="sm" />
        </Cell>
        <Cell label="Selected">
          <Thumbnail src={thumb67} frameNo="67" selected size="sm" />
        </Cell>
        <Cell label="Rejected">
          <Thumbnail src={thumb66} frameNo="66" rejected size="sm" />
        </Cell>
        <Cell label="No image">
          <Thumbnail src={null} frameNo="—" size="sm" />
        </Cell>
        <Cell label="Large (lg)">
          <Thumbnail src={thumb65} frameNo="65" representative size="lg" />
        </Cell>
      </div>
    );
  },
};
