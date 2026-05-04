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
  },
  build: {
    outDir: "dist",
    emptyOutDir: true,
    sourcemap: true,
  },
});
