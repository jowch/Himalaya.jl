import { useEffect, useRef } from "react";
import type { ReactNode } from "react";
import { useFocusTrap } from "../../hooks/useFocusTrap";

export type ModalSize = "sm" | "md" | "lg";
export type ModalAlign = "center" | "top";
export type ModalVariant = "dialog" | "drawer";

export interface ModalShellProps {
  /** Render gate. Returns null when false (matches every call site's `if (!open) return null`). */
  open: boolean;
  onClose: () => void;
  size?: ModalSize;            // default "md"; ignored for variant="drawer"
  align?: ModalAlign;          // default "center"; ignored for variant="drawer"
  /** Esc dismiss. Default true. Set false where the parent owns Escape (e.g. NavModal). */
  closeOnEsc?: boolean;
  /** Backdrop-click dismiss. Default true. */
  closeOnOutsideClick?: boolean;
  /** "drawer" → right-edge sheet with a lower-z scrim. Default "dialog". */
  variant?: ModalVariant;
  /** a11y: the caller must supply at least one of these. */
  "aria-label"?: string;
  "aria-labelledby"?: string;
  "aria-describedby"?: string;
  /** Forwarded to the FRAME (role=dialog). Scrim gets `${testId}-scrim`. */
  testId?: string;
  /** PLACEMENT-ONLY on the frame (max-height, flex/grid, one-off width). No appearance utilities. */
  className?: string;
  children: ReactNode;
}

const sizeClass: Record<ModalSize, string> = {
  sm: "w-[min(440px,calc(100vw-48px))]",
  md: "w-[min(640px,calc(100vw-48px))]",
  lg: "w-[min(820px,calc(100vw-48px))]",
};

const alignClass: Record<ModalAlign, string> = {
  center: "items-center justify-center",
  top: "items-start justify-center pt-[12vh]",
};

/** Tiny placement-only class joiner (brief-sanctioned; no cva/clsx/tailwind-merge). */
function cx(...parts: (string | false | undefined)[]): string {
  return parts.filter(Boolean).join(" ");
}

export function ModalShell({
  open,
  onClose,
  size = "md",
  align = "center",
  closeOnEsc = true,
  closeOnOutsideClick = true,
  variant = "dialog",
  testId,
  className = "",
  children,
  ...aria
}: ModalShellProps): JSX.Element | null {
  const frameRef = useRef<HTMLDivElement>(null);
  useFocusTrap(frameRef, open);

  useEffect(() => {
    if (!open || !closeOnEsc) return;
    const onKey = (e: KeyboardEvent): void => {
      if (e.key === "Escape") {
        e.preventDefault();
        onClose();
      }
    };
    document.addEventListener("keydown", onKey);
    return () => document.removeEventListener("keydown", onKey);
  }, [open, closeOnEsc, onClose]);

  if (!open) return null;

  const onScrimClick = (e: React.MouseEvent): void => {
    if (closeOnOutsideClick && e.target === e.currentTarget) onClose();
  };

  const ariaProps = {
    "aria-label": aria["aria-label"],
    "aria-labelledby": aria["aria-labelledby"],
    "aria-describedby": aria["aria-describedby"],
  };

  if (variant === "drawer") {
    return (
      <div data-testid={testId} className="contents">
        <div
          data-testid={testId ? `${testId}-scrim` : undefined}
          role="presentation"
          onClick={onScrimClick}
          className="fixed inset-0 z-40 bg-scrim anim-pal-in"
        />
        <div
          ref={frameRef}
          role="dialog"
          aria-modal="true"
          data-variant="drawer"
          {...ariaProps}
          className={cx(
            "fixed right-0 top-14 bottom-0 z-50 w-[300px] max-w-[85vw] overflow-y-auto",
            "bg-plate border-l border-hair-strong shadow-2xl",
            className,
          )}
        >
          {children}
        </div>
      </div>
    );
  }

  return (
    <div
      data-testid={testId ? `${testId}-scrim` : undefined}
      role="presentation"
      onClick={onScrimClick}
      className={cx(
        "fixed inset-0 z-50 flex bg-scrim backdrop-blur-sm anim-pal-in",
        alignClass[align],
      )}
    >
      <div
        ref={frameRef}
        data-testid={testId}
        role="dialog"
        aria-modal="true"
        data-variant="dialog"
        data-size={size}
        data-align={align}
        {...ariaProps}
        className={cx(
          sizeClass[size],
          "bg-plate border border-hair-strong rounded-md shadow-2xl anim-pal-scale",
          "flex flex-col overflow-hidden",
          className,
        )}
      >
        {children}
      </div>
    </div>
  );
}
