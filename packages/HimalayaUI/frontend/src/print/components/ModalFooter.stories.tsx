import type { Meta, StoryObj } from "@storybook/react-vite";
import { Button } from "../ui";
import { ModalFooter } from "./ModalFooter";

const meta = {
  title: "components/ModalFooter",
  component: ModalFooter,
} satisfies Meta<typeof ModalFooter>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {
  args: {
    note: "Changes apply immediately.",
    actions: (
      <>
        <Button variant="ghost">Cancel</Button>
        <Button variant="accent">Add</Button>
      </>
    ),
  },
};
