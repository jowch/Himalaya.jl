import type { ReactNode } from "react";

export function FieldError({ children, testId }: { children: ReactNode; testId?: string }): JSX.Element {
  return <span role="alert" data-testid={testId} className="text-caption text-error">{children}</span>;
}
