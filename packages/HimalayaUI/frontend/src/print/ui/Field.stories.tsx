import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { Field } from "./Field";

const meta = {
  title: "ui/Field",
  component: Field,
  args: { value: "LL37 : lipid ratio" },
} satisfies Meta<typeof Field>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const Placeholder: Story = {
  args: { value: "", placeholder: "Choose a variable" },
};

const ORDER_OPTIONS = [
  { value: "LL37 : lipid ratio", label: "LL37 : lipid ratio" },
  { value: "Time", label: "Time" },
  { value: "Dose", label: "Dose" },
  { value: "Temperature", label: "Temperature" },
  { value: "Define your own…", label: "Define your own…" },
] as const;

function DropdownDemo(): JSX.Element {
  const [value, setValue] = useState("LL37 : lipid ratio");
  return (
    <div className="w-[320px]">
      <Field value={value} options={ORDER_OPTIONS} onSelect={setValue} menuLabel="Ordering variable" />
    </div>
  );
}

export const Dropdown: Story = { render: () => <DropdownDemo /> };
