import { useCallback, useEffect, useState } from "react";

export type WorkspaceTheme = "system" | "light" | "dark";

const KEY = "himalaya-ui:workspace-theme";
const CHANGED = "himalaya-ui:theme-changed";

function read(): WorkspaceTheme {
  try {
    const v = window.localStorage.getItem(KEY);
    return v === "light" || v === "dark" ? v : "system";
  } catch {
    return "system";
  }
}

function systemDark(): boolean {
  return typeof window.matchMedia === "function" && window.matchMedia("(prefers-color-scheme: dark)").matches;
}

export function effectiveTheme(theme: WorkspaceTheme, prefersDark: boolean): "light" | "dark" {
  if (theme === "system") return prefersDark ? "dark" : "light";
  return theme;
}

export function useWorkspaceTheme(): { theme: WorkspaceTheme; isDark: boolean; toggle: () => void; followSystem: () => void } {
  const [theme, setTheme] = useState<WorkspaceTheme>(read);
  const [prefersDark, setPrefersDark] = useState(systemDark);
  useEffect(() => {
    if (typeof window.matchMedia !== "function") return;
    const mq = window.matchMedia("(prefers-color-scheme: dark)");
    const on = (): void => setPrefersDark(mq.matches);
    mq.addEventListener?.("change", on);
    return () => mq.removeEventListener?.("change", on);
  }, []);
  const isDark = effectiveTheme(theme, prefersDark) === "dark";
  const choose = useCallback((next: WorkspaceTheme) => {
    setTheme(next);
    try {
      if (next === "system") window.localStorage.removeItem(KEY);
      else window.localStorage.setItem(KEY, next);
    } catch {
      return;
    } finally {
      window.dispatchEvent(new Event(CHANGED));
    }
  }, []);
  const toggle = useCallback(() => choose(isDark ? "light" : "dark"), [choose, isDark]);
  const followSystem = useCallback(() => choose("system"), [choose]);
  return { theme, isDark, toggle, followSystem };
}

export function useDocumentTheme(): void {
  useEffect(() => {
    const apply = (): void => document.documentElement.setAttribute("data-workspace-theme", read());
    apply();
    window.addEventListener(CHANGED, apply);
    window.addEventListener("storage", apply);
    return () => {
      window.removeEventListener(CHANGED, apply);
      window.removeEventListener("storage", apply);
    };
  }, []);
}
