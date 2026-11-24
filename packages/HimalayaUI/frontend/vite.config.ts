import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tailwindcss from "@tailwindcss/vite";

// Backend port defaults to 8080 (matching prod); override with VITE_API_PORT
// when running a second dev pair on a worktree (the prod serve already
// holds :8080). Frontend port is also overridable via the standard --port
// CLI flag.
const apiPort = process.env.VITE_API_PORT ?? "8080";

export default defineConfig({
  plugins: [react(), tailwindcss()],
  server: {
    port: 5173,
    proxy: {
      "/api": {
        target: `http://127.0.0.1:${apiPort}`,
        changeOrigin: false,
        configure: (proxy) => {
          proxy.on("error", (_err, _req, res) => {
            if (res && "writeHead" in res && !res.headersSent) {
              res.writeHead(502, { "Content-Type": "text/plain" });
              res.end("backend unreachable");
            }
          });
          proxy.on("proxyRes", (upstream, _req, res) => {
            upstream.on("close", () => {
              if (!res.writableEnded) res.destroy();
            });
          });
        },
      },
    },
    // Don't HMR-track e2e specs or Playwright artifacts. Without this, a
    // spec edit during a live Playwright run triggers a page reload that
    // crashes the Vite dev server (no further log; child exits silently).
    // e2e specs aren't part of the served bundle anyway.
    watch: {
      ignored: [
        "**/e2e/**",
        "**/test-results/**",
        "**/playwright-report/**",
      ],
    },
  },
  build: {
    outDir: "dist",
    emptyOutDir: true,
    sourcemap: true,
    rolldownOptions: {
      output: {
        codeSplitting: {
          groups: [{ name: "vendor", test: /node_modules/ }],
        },
      },
    },
  },
});
