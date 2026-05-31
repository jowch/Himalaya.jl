import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { Input } from "./Input";

const meta = {
  title: "ui/Input",
  component: Input,
  args: { value: "", onValueChange: () => {}, placeholder: "Lipid name" },
  render: (args) => {
    const [v, setV] = useState(args.value);
    return <Input {...args} value={v} onValueChange={setV} />;
  },
} satisfies Meta<typeof Input>;

export default meta;
type Story = StoryObj<typeof meta>;

function Magnifier(): JSX.Element {
  return (
    <svg
      width={14}
      height={14}
      viewBox="0 0 14 14"
      fill="none"
      aria-hidden="true"
      className="flex-shrink-0"
    >
      <circle cx="6" cy="6" r="4.3" stroke="var(--color-ink-faint)" strokeWidth={1.5} />
      <line
        x1="9.2"
        y1="9.2"
        x2="12.5"
        y2="12.5"
        stroke="var(--color-ink-faint)"
        strokeWidth={1.5}
        strokeLinecap="round"
      />
    </svg>
  );
}

export const Default: Story = {};

export const WithLeadingIcon: Story = {
  args: { leading: <Magnifier />, placeholder: "Search series" },
};

export const Invalid: Story = {
  args: { invalid: true, value: "??", placeholder: "q value" },
};

export const Small: Story = {
  args: { inputSize: "sm", placeholder: "order" },
};
