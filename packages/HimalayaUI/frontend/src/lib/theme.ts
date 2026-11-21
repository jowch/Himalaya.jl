import { useCallback, useLayoutEffect, useSyncExternalStore } from "react";

export type WorkspaceTheme = "system" | "light" | "dark";

const KEY = "himalaya-ui:workspace-theme";
const CHANGED = "himalaya-ui:theme-changed";
const DARK_QUERY = "(prefers-color-scheme: dark)";

let unsaved: WorkspaceTheme = "system";

function readTheme(): WorkspaceTheme {
  try {
    const v = window.localStorage.getItem(KEY);
    return v === "light" || v === "dark" ? v : "system";
  } catch {
    return unsaved;
  }
}

function writeTheme(next: WorkspaceTheme): void {
  try {
    if (next === "system") window.localStorage.removeItem(KEY);
    else window.localStorage.setItem(KEY, next);
  } catch {
    unsaved = next;
  }
  window.dispatchEvent(new Event(CHANGED));
}

function darkQuery(): MediaQueryList | undefined {
  return typeof window.matchMedia === "function" ? window.matchMedia(DARK_QUERY) : undefined;
}

function subscribe(onChange: () => void): () => void {
  const mq = darkQuery();
  window.addEventListener(CHANGED, onChange);
  window.addEventListener("storage", onChange);
  mq?.addEventListener?.("change", onChange);
  return () => {
    window.removeEventListener(CHANGED, onChange);
    window.removeEventListener("storage", onChange);
    mq?.removeEventListener?.("change", onChange);
  };
}

function systemDark(): boolean {
  return darkQuery()?.matches ?? false;
}

export function effectiveTheme(theme: WorkspaceTheme, prefersDark: boolean): "light" | "dark" {
  if (theme === "system") return prefersDark ? "dark" : "light";
  return theme;
}

export function useWorkspaceTheme(): { theme: WorkspaceTheme; isDark: boolean; toggle: () => void; followSystem: () => void } {
  const theme = useSyncExternalStore(subscribe, readTheme);
  const prefersDark = useSyncExternalStore(subscribe, systemDark);
  const isDark = effectiveTheme(theme, prefersDark) === "dark";
  const toggle = useCallback(() => writeTheme(isDark ? "light" : "dark"), [isDark]);
  const followSystem = useCallback(() => writeTheme("system"), []);
  return { theme, isDark, toggle, followSystem };
}

export function useDocumentTheme(): void {
  const { theme, isDark } = useWorkspaceTheme();
  useLayoutEffect(() => {
    const root = document.documentElement;
    root.setAttribute("data-workspace-theme", theme);
    root.classList.toggle("dark", isDark);
  }, [theme, isDark]);
}
