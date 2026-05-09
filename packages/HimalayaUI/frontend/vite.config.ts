import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tailwindcss from "@tailwindcss/vite";
import { boneyardPlugin } from "boneyard-js/vite";

// Backend port defaults to 8080 (matching prod); override with VITE_API_PORT
// when running a second dev pair on a worktree (the prod serve already
// holds :8080). Frontend port is also overridable via the standard --port
// CLI flag.
const apiPort = process.env.VITE_API_PORT ?? "8080";

export default defineConfig({
  plugins: [react(), tailwindcss(), boneyardPlugin()],
  server: {
    port: 5173,
    proxy: {
      "/api": { target: `http://127.0.0.1:${apiPort}`, changeOrigin: false },
    },
    // Don't HMR-track e2e specs or Playwright artifacts. Without this, a
    // spec edit during a live Playwright run triggers a page reload that
    // crashes the Vite dev server (no further log; child exits silently).
    // e2e specs aren't part of the served bundle anyway.
    //
    // src/bones/** holds boneyard skeleton captures + the auto-generated
    // registry. The boneyard plugin writes these on cold loads of any
    // captured Skeleton; without ignoring them, every cold load forces a
    // full page reload mid-interaction.
    watch: {
      ignored: [
        "**/e2e/**",
        "**/test-results/**",
        "**/playwright-report/**",
        "**/src/bones/**",
      ],
    },
  },
  build: {
    outDir: "dist",
    emptyOutDir: true,
    sourcemap: true,
  },
});
