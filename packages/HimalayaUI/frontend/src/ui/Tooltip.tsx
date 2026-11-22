import { cloneElement, useCallback, useEffect, useId, useLayoutEffect, useRef, useState } from "react";
import type { FocusEvent, MouseEvent, MutableRefObject, ReactElement, ReactNode, Ref } from "react";
import { createPortal } from "react-dom";
import { Kbd } from "./Kbd";

export const TOOLTIP_OPEN_DELAY_MS = 400;
export const TOOLTIP_WARM_WINDOW_MS = 300;
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
  onPointerDown?: (e: MouseEvent<HTMLElement>) => void;
  onFocus?: (e: FocusEvent<HTMLElement>) => void;
  onBlur?: (e: FocusEvent<HTMLElement>) => void;
  "aria-describedby"?: string;
};

let lastClosedAt = Number.NEGATIVE_INFINITY;
let closeOpen: (() => void) | null = null;

export function resetTooltipWarmth(): void {
  lastClosedAt = Number.NEGATIVE_INFINITY;
  closeOpen = null;
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
  const triggerRef = useRef<HTMLElement | null>(null);
  const tipRef = useRef<HTMLSpanElement | null>(null);
  const timer = useRef<ReturnType<typeof setTimeout> | null>(null);
  const suppressed = useRef(false);

  const clearTimer = (): void => {
    if (timer.current !== null) clearTimeout(timer.current);
    timer.current = null;
  };

  const close = useCallback((): void => {
    clearTimer();
    setOpen((was) => {
      if (was) lastClosedAt = Date.now();
      return false;
    });
    setPos(null);
  }, []);

  const show = useCallback((): void => {
    const el = triggerRef.current;
    if (!el || suppressed.current) return;
    if (overflowOnly && !isTruncated(el)) return;
    if (closeOpen && closeOpen !== close) closeOpen();
    closeOpen = close;
    setOpen(true);
  }, [close, overflowOnly]);

  const schedule = (): void => {
    clearTimer();
    if (isWarm()) show();
    else timer.current = setTimeout(show, TOOLTIP_OPEN_DELAY_MS);
  };

  useEffect(() => () => {
    clearTimer();
    if (closeOpen === close) closeOpen = null;
  }, [close]);

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent): void => {
      if (e.key === "Escape") {
        suppressed.current = true;
        close();
      }
    };
    window.addEventListener("keydown", onKey, true);
    window.addEventListener("scroll", close, true);
    window.addEventListener("resize", close);
    return () => {
      window.removeEventListener("keydown", onKey, true);
      window.removeEventListener("scroll", close, true);
      window.removeEventListener("resize", close);
      if (closeOpen === close) closeOpen = null;
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

  const empty = content === null || content === undefined || content === "";
  if (empty) return children;

  const own = children.props as TriggerProps;
  const childRef = (children as unknown as { ref?: Ref<HTMLElement> }).ref;
  const trigger = cloneElement(children, {
    ref: (node: HTMLElement | null) => {
      triggerRef.current = node;
      assignRef(childRef, node);
    },
    onMouseEnter: (e: MouseEvent<HTMLElement>) => {
      own.onMouseEnter?.(e);
      schedule();
    },
    onMouseLeave: (e: MouseEvent<HTMLElement>) => {
      own.onMouseLeave?.(e);
      suppressed.current = false;
      close();
    },
    onPointerDown: (e: MouseEvent<HTMLElement>) => {
      own.onPointerDown?.(e);
      suppressed.current = true;
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
      suppressed.current = false;
      close();
    },
    ...(open ? { "aria-describedby": [own["aria-describedby"], id].filter(Boolean).join(" ") } : {}),
  } as TriggerProps);

  return (
    <>
      {trigger}
      {open &&
        createPortal(
          <span
            ref={tipRef}
            role="tooltip"
            id={id}
            data-testid="tooltip"
            data-side={pos?.side ?? side}
            className="pointer-events-none fixed z-[60] flex w-max max-w-[18rem] items-baseline gap-2 whitespace-normal break-words rounded-sm bg-ink px-2 py-1 text-left text-xs font-normal text-paper shadow-lg"
            style={pos ? { top: pos.top, left: pos.left } : { top: 0, left: 0, visibility: "hidden" }}
          >
            <span className="min-w-0">{content}</span>
            {shortcut && <Kbd onFill>{shortcut}</Kbd>}
          </span>,
          document.body,
        )}
    </>
  );
}
