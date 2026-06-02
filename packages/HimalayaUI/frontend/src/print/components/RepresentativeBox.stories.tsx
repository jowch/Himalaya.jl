import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { RepresentativeBox } from "./RepresentativeBox";

const meta = {
  title: "components/RepresentativeBox",
  component: RepresentativeBox,
  args: {
    isRepresentative: false,
  },
} satisfies Meta<typeof RepresentativeBox>;

export default meta;
type Story = StoryObj<typeof meta>;

export const NotRepresentative: Story = {
  args: { isRepresentative: false, onSetRepresentative: () => {} },
  render: (args) => (
    <div className="bg-paper p-4 w-80">
      <RepresentativeBox {...args} />
    </div>
  ),
};

export const IsRepresentative: Story = {
  args: { isRepresentative: true, onSetRepresentative: () => {} },
  render: (args) => (
    <div className="bg-paper p-4 w-80">
      <RepresentativeBox {...args} />
    </div>
  ),
};

// Stateful wrapper so toggling works interactively in Storybook.
function RepresentativeBoxToggle(): JSX.Element {
  const [isRep, setIsRep] = useState(false);
  return (
    <div className="bg-paper p-4 w-80">
      <RepresentativeBox
        isRepresentative={isRep}
        onSetRepresentative={() => setIsRep(true)}
      />
    </div>
  );
}

export const Interactive: Story = {
  render: () => <RepresentativeBoxToggle />,
};
