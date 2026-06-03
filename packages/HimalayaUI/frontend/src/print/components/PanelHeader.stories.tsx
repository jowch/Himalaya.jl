// src/print/components/PanelHeader.stories.tsx
import type { Meta, StoryObj } from "@storybook/react-vite";
import { PanelHeader } from "./PanelHeader";
import { SegmentedControl } from "../ui";

const meta: Meta<typeof PanelHeader> = {
  title: "components/PanelHeader",
  component: PanelHeader,
};
export default meta;
type Story = StoryObj<typeof meta>;

export const LabelOnly: Story = {
  args: { label: "Detector image" },
};

export const WithTools: Story = {
  render: () => (
    <PanelHeader label="Reflections — comb">
      <SegmentedControl
        size="xs"
        options={[
          { value: "comb", label: "comb" },
          { value: "resid", label: "indexing space" },
        ]}
        value="comb"
        onChange={() => {}}
        aria-label="comb view"
      />
    </PanelHeader>
  ),
};
