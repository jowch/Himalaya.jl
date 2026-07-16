import type { Meta, StoryObj } from "@storybook/react-vite";
import { ModalHead } from "./ModalHead";

const noop = () => {};

const meta = {
  title: "components/ModalHead",
  component: ModalHead,
} satisfies Meta<typeof ModalHead>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {
  args: {
    kicker: "Speculative",
    title: "Custom index",
    onClose: noop,
  },
};
