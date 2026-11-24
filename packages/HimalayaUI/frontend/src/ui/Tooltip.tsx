import { cloneElement, useCallback, useEffect, useId, useLayoutEffect, useRef, useState } from "react";
import type { FocusEvent, MouseEvent, MutableRefObject, PointerEvent, ReactElement, ReactNode, Ref } from "react";
import { createPortal } from "react-dom";
import { Kbd } from "./Kbd";

export const TOOLTIP_OPEN_DELAY_MS = 400;
export const TOOLTIP_WARM_WINDOW_MS = 300;
const SCROLL_GRACE_MS = 100;
const GAP = 6;
const MARGIN = 8;

type Side = "top" | "bottom";

interface TooltipProps {
  content: ReactNode;
  children: ReactElement;
  shortcut?: string;
  side?: Side;
  overflowOnly?: boolean;
}

type TriggerProps = {
  onMouseEnter?: (e: MouseEvent<HTMLElement>) => void;
  onMouseLeave?: (e: MouseEvent<HTMLElement>) => void;
  onPointerDown?: (e: PointerEvent<HTMLElement>) => void;
  onFocus?: (e: FocusEvent<HTMLElement>) => void;
  onBlur?: (e: FocusEvent<HTMLElement>) => void;
  "aria-describedby"?: string;
  "aria-label"?: string;
  "aria-description"?: string;
};

type Suppression = "pointer" | "key" | null;

let lastClosedAt = Number.NEGATIVE_INFINITY;
let active: { close: () => void; onEscape: () => void } | null = null;
let escapeListening = false;

function onWindowKeyDown(e: KeyboardEvent): void {
  if (e.key === "Escape") active?.onEscape();
}

function listenForEscape(): void {
  if (escapeListening || typeof window === "undefined") return;
  window.addEventListener("keydown", onWindowKeyDown, true);
  escapeListening = true;
}

export function resetTooltipWarmth(): void {
  lastClosedAt = Number.NEGATIVE_INFINITY;
  active = null;
}

function isWarm(): boolean {
  return Date.now() - lastClosedAt < TOOLTIP_WARM_WINDOW_MS;
}

function focusVisible(el: HTMLElement): boolean {
  try {
    return el.matches(":focus-visible");
  } catch {
    return true;
  }
}

function isTruncated(el: HTMLElement): boolean {
  return el.scrollWidth > el.clientWidth || el.scrollHeight > el.clientHeight;
}

function assignRef<T>(ref: Ref<T> | undefined, value: T): void {
  if (typeof ref === "function") ref(value);
  else if (ref) (ref as MutableRefObject<T>).current = value;
}

export function placeTooltip(
  trigger: { top: number; bottom: number; left: number; width: number },
  tip: { width: number; height: number },
  viewport: { width: number; height: number },
  preferred: Side,
): { top: number; left: number; side: Side } {
  const fitsTop = trigger.top - GAP - tip.height >= MARGIN;
  const fitsBottom = trigger.bottom + GAP + tip.height <= viewport.height - MARGIN;
  const side: Side = preferred === "top" ? (fitsTop || !fitsBottom ? "top" : "bottom") : (fitsBottom || !fitsTop ? "bottom" : "top");
  const top = side === "top" ? trigger.top - GAP - tip.height : trigger.bottom + GAP;
  const centred = trigger.left + trigger.width / 2 - tip.width / 2;
  const left = Math.max(MARGIN, Math.min(centred, viewport.width - MARGIN - tip.width));
  return { top: Math.max(MARGIN, top), left, side };
}

export function Tooltip({ content, children, shortcut, side = "top", overflowOnly = false }: TooltipProps): JSX.Element {
  const id = useId();
  const [open, setOpen] = useState(false);
  const [pos, setPos] = useState<{ top: number; left: number; side: Side } | null>(null);
  const openRef = useRef(false);
  const openedAt = useRef(0);
  const triggerRef = useRef<HTMLElement | null>(null);
  const tipRef = useRef<HTMLSpanElement | null>(null);
  const timer = useRef<ReturnType<typeof setTimeout> | null>(null);
  const suppressed = useRef<Suppression>(null);
  const empty = content === null || content === undefined || content === "";

  const clearTimer = (): void => {
    if (timer.current !== null) clearTimeout(timer.current);
    timer.current = null;
  };

  const handle = useRef<{ close: () => void; onEscape: () => void } | null>(null);

  const close = useCallback((): void => {
    clearTimer();
    if (openRef.current) lastClosedAt = Date.now();
    openRef.current = false;
    if (active === handle.current) active = null;
    setOpen(false);
    setPos(null);
  }, []);

  if (handle.current === null) {
    handle.current = {
      close: () => close(),
      onEscape: () => {
        suppressed.current = "key";
        close();
      },
    };
  }

  const claim = (): void => {
    if (active && active !== handle.current) active.close();
    active = handle.current;
    listenForEscape();
  };

  const show = useCallback((): void => {
    timer.current = null;
    const el = triggerRef.current;
    if (!el || suppressed.current || empty) return;
    if (overflowOnly && !isTruncated(el)) return;
    claim();
    openRef.current = true;
    openedAt.current = Date.now();
    setOpen(true);
  }, [empty, overflowOnly]);

  const schedule = (): void => {
    clearTimer();
    if (isWarm()) {
      show();
      return;
    }
    claim();
    timer.current = setTimeout(show, TOOLTIP_OPEN_DELAY_MS);
  };

  useEffect(() => {
    if (empty) close();
  }, [empty, close]);

  useEffect(() => () => {
    clearTimer();
    if (active === handle.current) active = null;
  }, []);

  useEffect(() => {
    if (import.meta.env.DEV && !empty && triggerRef.current === null) {
      console.error("Tooltip: the trigger must be a DOM element or a component that forwards its ref to one.");
    }
  }, [empty]);

  useEffect(() => {
    if (!open) return;
    const onScroll = (): void => {
      if (Date.now() - openedAt.current >= SCROLL_GRACE_MS) close();
    };
    window.addEventListener("scroll", onScroll, true);
    window.addEventListener("resize", close);
    return () => {
      window.removeEventListener("scroll", onScroll, true);
      window.removeEventListener("resize", close);
    };
  }, [open, close]);

  useLayoutEffect(() => {
    if (!open) return;
    const el = triggerRef.current;
    const tip = tipRef.current;
    if (!el || !tip) return;
    const r = el.getBoundingClientRect();
    const t = tip.getBoundingClientRect();
    setPos(placeTooltip(r, t, { width: window.innerWidth, height: window.innerHeight }, side));
  }, [open, side, content, shortcut]);

  const own = children.props as TriggerProps;
  const childRef = (children as unknown as { ref?: Ref<HTMLElement> }).ref;
  const namesTrigger = overflowOnly || (typeof content === "string" && own["aria-label"] === content);
  const repeatsText = typeof content === "string" && (triggerRef.current?.textContent ?? "").includes(content);
  const describedAlready = typeof content === "string" && own["aria-description"] === content;
  const describes = open && !empty && !namesTrigger && !repeatsText && !describedAlready;
  const trigger = cloneElement(children, {
    ref: (node: HTMLElement | null) => {
      triggerRef.current = node;
      assignRef(childRef, node);
    },
    onMouseEnter: (e: MouseEvent<HTMLElement>) => {
      own.onMouseEnter?.(e);
      if (suppressed.current === "key") suppressed.current = null;
      schedule();
    },
    onMouseLeave: (e: MouseEvent<HTMLElement>) => {
      own.onMouseLeave?.(e);
      suppressed.current = null;
      close();
    },
    onPointerDown: (e: PointerEvent<HTMLElement>) => {
      own.onPointerDown?.(e);
      suppressed.current = "pointer";
      close();
    },
    onFocus: (e: FocusEvent<HTMLElement>) => {
      own.onFocus?.(e);
      if (focusVisible(e.currentTarget)) {
        clearTimer();
        show();
      }
    },
    onBlur: (e: FocusEvent<HTMLElement>) => {
      own.onBlur?.(e);
      suppressed.current = null;
      close();
    },
    ...(describes ? { "aria-describedby": [own["aria-describedby"], id].filter(Boolean).join(" ") } : {}),
  } as TriggerProps);

  return (
    <>
      {trigger}
      {open &&
        !empty &&
        createPortal(
          <span
            ref={tipRef}
            role="tooltip"
            id={id}
            data-testid="tooltip"
            data-side={pos?.side ?? side}
            className="pointer-events-none fixed z-[60] flex w-max max-w-[18rem] items-baseline gap-2 whitespace-normal break-words rounded-sm border border-hair bg-plate px-2 py-1 text-left text-xs font-normal text-ink shadow-lg"
            style={pos ? { top: pos.top, left: pos.left } : { top: 0, left: 0, visibility: "hidden" }}
          >
            <span className="min-w-0">{content}</span>
            {shortcut && <Kbd>{shortcut}</Kbd>}
          </span>,
          document.body,
        )}
    </>
  );
}
