import { useLayoutEffect, useSyncExternalStore } from "react";

export type AppTheme = "system" | "light" | "dark";

const KEY = "himalaya-ui:theme";
const OLD_KEY = "himalaya-ui:workspace-theme";
const CHANGED = "himalaya-ui:theme-changed";
const DARK_QUERY = "(prefers-color-scheme: dark)";

let unsaved: AppTheme = "system";

function readTheme(): AppTheme {
  try {
    const store = window.localStorage;
    const old = store.getItem(OLD_KEY);
    if (old !== null) {
      if (store.getItem(KEY) === null && (old === "light" || old === "dark")) store.setItem(KEY, old);
      store.removeItem(OLD_KEY);
    }
    const v = store.getItem(KEY);
    return v === "light" || v === "dark" ? v : "system";
  } catch {
    return unsaved;
  }
}

function writeTheme(next: AppTheme): void {
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

export function effectiveTheme(theme: AppTheme, prefersDark: boolean): "light" | "dark" {
  if (theme === "system") return prefersDark ? "dark" : "light";
  return theme;
}

export function useAppTheme(): { theme: AppTheme; isDark: boolean; choose: (next: AppTheme) => void } {
  const theme = useSyncExternalStore(subscribe, readTheme);
  const prefersDark = useSyncExternalStore(subscribe, systemDark);
  const isDark = effectiveTheme(theme, prefersDark) === "dark";
  return { theme, isDark, choose: writeTheme };
}

export function useDocumentTheme(): void {
  const { theme } = useAppTheme();
  useLayoutEffect(() => {
    document.documentElement.setAttribute("data-theme", theme);
  }, [theme]);
}
