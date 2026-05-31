import { useEffect } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { ToastContainer } from "./Toast";
import { showToast } from "../../lib/toast";

// ToastContainer renders from an internal queue that is fed via the global
// `showToast` impl (wired in ToastContainer's mount effect). With no toasts it
// renders an empty fixed surface, so the stories below seed one toast per kind
// after mount via showToast — the public API — so the surface is visible.
const meta = {
  title: "ui/Toast",
  component: ToastContainer,
} satisfies Meta<typeof ToastContainer>;

export default meta;
type Story = StoryObj<typeof meta>;

function Seeded({
  msg,
  kind,
}: {
  msg: string;
  kind: "info" | "success" | "warning" | "error";
}): JSX.Element {
  // Seed after mount so ToastContainer's setToastImpl effect has run.
  useEffect(() => {
    const t = setTimeout(() => showToast(msg, kind), 0);
    return () => clearTimeout(t);
  }, [msg, kind]);
  return <ToastContainer />;
}

export const Info: Story = {
  render: () => <Seeded msg="Re-analysis queued." kind="info" />,
};

export const Success: Story = {
  render: () => <Seeded msg="Phase call saved." kind="success" />,
};

export const Warning: Story = {
  render: () => <Seeded msg="Some peaks fell outside the fit window." kind="warning" />,
};

export const Error: Story = {
  render: () => <Seeded msg="Could not reach the server." kind="error" />,
};
