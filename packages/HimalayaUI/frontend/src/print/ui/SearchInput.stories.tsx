import type { Meta, StoryObj } from "@storybook/react-vite";
import { SearchInput } from "./SearchInput";

const meta = {
  title: "ui/SearchInput",
  component: SearchInput,
  args: { value: "", placeholder: "Search series", onChange: () => {} },
} satisfies Meta<typeof SearchInput>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Empty: Story = {
  args: { value: "", placeholder: "Search series", onChange: () => {} },
};
export const WithValue: Story = {
  args: { value: "lipid-A", placeholder: "Search series", onChange: () => {} },
};
