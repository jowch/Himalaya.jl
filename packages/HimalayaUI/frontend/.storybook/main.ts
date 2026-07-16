import { defineMain } from "@storybook/react-vite/node";

// No `core.builder` — the Vite builder is implied by the react-vite framework.
// builder-vite auto-merges this project's vite.config.ts, so @tailwindcss/vite
// and @vitejs/plugin-react apply in stories with no extra wiring.
export default defineMain({
  framework: "@storybook/react-vite",
  stories: ["../src/print/**/*.stories.@(ts|tsx)"],
});
