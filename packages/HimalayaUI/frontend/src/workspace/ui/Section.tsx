import type { ReactNode } from "react";

export function Section({ title, aside, className, children }: { title: string; aside?: ReactNode; className?: string; children: ReactNode }): JSX.Element {
  return (
    <section className={"flex flex-col px-3.5 py-3" + (className ? ` ${className}` : "")}>
      <header className="mb-2 flex items-baseline justify-between gap-2">
        <h2 className="text-sm font-semibold text-ink-soft">{title}</h2>
        {aside && <span className="text-xs text-ink-faint">{aside}</span>}
      </header>
      {children}
    </section>
  );
}

export function MainSection({ title, hint, children }: { title: string; hint?: ReactNode; children: ReactNode }): JSX.Element {
  return (
    <section className="flex flex-col">
      <header className="flex items-baseline gap-2.5 px-3.5 pb-1.5 pt-3">
        <h2 className="text-sm font-semibold text-ink-soft">{title}</h2>
        {hint && <span className="text-xs text-ink-faint">{hint}</span>}
      </header>
      <div className="px-3.5 pb-4">{children}</div>
    </section>
  );
}

export function SampleHeading({ name, meta }: { name: string; meta?: ReactNode }): JSX.Element {
  return (
    <div className="flex min-w-0 items-baseline gap-3 px-3.5 pb-1 pt-2.5">
      <h2 className="m-0 min-w-0 truncate text-xl font-normal leading-tight text-ink" title={name}>{name}</h2>
      {meta && <span className="shrink-0 whitespace-nowrap text-base text-ink-faint">{meta}</span>}
    </div>
  );
}
