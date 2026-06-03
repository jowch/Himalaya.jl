import type { Meta, StoryObj } from "@storybook/react-vite";
import { ModalFieldRow } from "./ModalFieldRow";

const meta = {
  title: "components/ModalFieldRow",
  component: ModalFieldRow,
} satisfies Meta<typeof ModalFieldRow>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {
  args: {
    label: "Symmetry",
    children: <input type="text" placeholder="e.g. Pn3m" />,
  },
};

export const WithSuffix: Story = {
  args: {
    label: "Lattice",
    labelSuffix: "a",
    children: <input type="text" placeholder="e.g. 81.2" />,
  },
};
