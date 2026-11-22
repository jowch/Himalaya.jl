import type { ReactNode } from "react";

export function Note({ children }: { children: ReactNode }): JSX.Element {
  return <p className="text-caption">{children}</p>;
}
