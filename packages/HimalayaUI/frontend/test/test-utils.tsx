import type { ReactElement, ReactNode } from "react";
import { render, type RenderOptions } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";

export function makeClient(): QueryClient {
  return new QueryClient({
    // gcTime: Infinity so observer-less seeded data (setQueryData without a
    // useQuery) survives across mutate→success transitions. Tests create a
    // fresh client each run, so cross-test memory pressure is irrelevant.
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

export function renderWithProviders(
  ui: ReactElement,
  options: RenderOptions & { client?: QueryClient } = {},
) {
  const { client = makeClient(), ...rest } = options;
  const Wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, ...render(ui, { wrapper: Wrapper, ...rest }) };
}
