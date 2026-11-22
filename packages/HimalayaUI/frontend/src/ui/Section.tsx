import type { ReactNode } from "react";
import { Tooltip } from "./Tooltip";

interface Collapsible {
  open: boolean;
  onOpen?: (() => void) | undefined;
}

export function Section({ title, aside, folded, className, collapsible, children }: { title: string; aside?: ReactNode; folded?: string | undefined; className?: string; collapsible?: Collapsible | undefined; children: ReactNode }): JSX.Element {
  const open = collapsible?.open ?? true;
  return (
    <section className={"flex flex-col px-3.5 py-3" + (className ? ` ${className}` : "")} data-section={title} data-collapsed={open ? undefined : "true"}>
      <header className={(open ? "mb-2 " : "") + "flex items-baseline justify-between gap-2"}>
        {collapsible?.onOpen ? (
          <h2 className="text-sm font-semibold text-ink-soft">
            <button type="button" aria-expanded={open} onClick={collapsible.onOpen} className="inline-flex items-baseline gap-1 hover:text-ink">
              <span aria-hidden="true" className="inline-block w-2.5 text-ink-faint">{open ? "▾" : "▸"}</span>
              {title}
            </button>
          </h2>
        ) : collapsible ? (
          <h2 className="inline-flex items-baseline gap-1 text-sm font-semibold text-ink-faint">
            <span aria-hidden="true" className="inline-block w-2.5">{open ? "▾" : "▸"}</span>
            {title}
          </h2>
        ) : (
          <h2 className="text-sm font-semibold text-ink-soft">{title}</h2>
        )}
        {aside && open && <span className="text-xs text-ink-faint">{aside}</span>}
        {folded && !open && (
          <Tooltip content={folded} overflowOnly>
            <span className="min-w-0 flex-1 truncate text-right text-xs text-ink-faint" data-section-folded="">{folded}</span>
          </Tooltip>
        )}
      </header>
      {open && children}
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

export function SampleHeading({ name, meta, aside }: { name: string; meta?: ReactNode; aside?: ReactNode }): JSX.Element {
  return (
    <div className="@container flex min-w-0 items-center gap-3 px-3.5 pb-1 pt-2.5">
      <div className="flex min-w-0 items-baseline gap-3">
        <Tooltip content={name} overflowOnly>
          <h2 className="m-0 min-w-0 truncate text-xl font-normal leading-tight text-ink">{name}</h2>
        </Tooltip>
        {meta && <span className="shrink-0 whitespace-nowrap text-base text-ink-faint @max-[34rem]:hidden">{meta}</span>}
      </div>
      {aside && <div className="flex min-w-min flex-1 basis-0 justify-end">{aside}</div>}
    </div>
  );
}
