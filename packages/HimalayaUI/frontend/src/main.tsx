import "./bones/registry";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { configureBoneyard } from "boneyard-js/react";
import { App } from "./App";
import { ErrorBoundary } from "./ErrorBoundary";

// Global Skeleton defaults — kept here (not in registry.ts) because the Vite
// HMR plugin rewrites registry.ts without preserving configureBoneyard calls.
configureBoneyard({
  color: "rgba(30, 31, 38, 1)",
  darkColor: "rgba(42, 44, 52, 1)",
  animate: "pulse",
  boneClass: "bone",
});

const queryClient = new QueryClient({
  defaultOptions: {
    queries: { staleTime: 30_000, retry: 1, refetchOnWindowFocus: false },
  },
});

const root = document.getElementById("app");
if (!root) throw new Error("#app root missing");
createRoot(root).render(
  <StrictMode>
    <ErrorBoundary>
      <QueryClientProvider client={queryClient}>
        <App />
      </QueryClientProvider>
    </ErrorBoundary>
  </StrictMode>,
);
