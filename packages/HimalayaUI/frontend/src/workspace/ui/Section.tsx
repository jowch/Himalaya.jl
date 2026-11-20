import type { ReactNode } from "react";

export function Section({ title, aside, children }: { title: string; aside?: ReactNode; children: ReactNode }): JSX.Element {
  return (
    <section className="flex flex-col gap-2">
      <header className="flex items-baseline justify-between gap-2">
        <h2 className="text-label">{title}</h2>
        {aside && <span className="text-caption">{aside}</span>}
      </header>
      {children}
    </section>
  );
}
