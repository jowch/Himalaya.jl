import type { Meta, StoryObj } from "@storybook/react-vite";
import { Chip } from "./Chip";

const meta = {
  title: "ui/Chip",
  component: Chip,
  args: { children: "Pn3m" },
} satisfies Meta<typeof Chip>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Static: Story = { args: { variant: "static" } };
export const Removable: Story = {
  args: { variant: "removable", onRemove: () => {} },
};
export const Add: Story = { args: { variant: "add", children: "+ add" } };
export const ToggleOff: Story = {
  args: { variant: "toggle", active: false, children: "Kept" },
};
export const ToggleOn: Story = {
  args: { variant: "toggle", active: true, children: "Kept" },
};
export const Trigger: Story = {
  args: { variant: "trigger", children: "Beamtime" },
};

/** Size is orthogonal to variant: the SAME variant rendered at `sm` and `md`. */
export const Sizes: Story = {
  render: () => (
    <div style={{ display: "flex", alignItems: "center", gap: 16 }}>
      <Chip variant="static" size="sm">
        static sm
      </Chip>
      <Chip variant="static" size="md">
        static md
      </Chip>
      <Chip variant="toggle" size="sm" active>
        toggle sm
      </Chip>
      <Chip variant="toggle" size="md" active>
        toggle md
      </Chip>
    </div>
  ),
};
