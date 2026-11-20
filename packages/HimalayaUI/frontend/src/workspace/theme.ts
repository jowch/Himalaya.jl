import { useCallback, useEffect, useState } from "react";

export type WorkspaceTheme = "system" | "light" | "dark";

const KEY = "himalaya-ui:workspace-theme";

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

export function useWorkspaceTheme(): { theme: WorkspaceTheme; isDark: boolean; toggle: () => void } {
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
  const toggle = useCallback(() => {
    const next: WorkspaceTheme = isDark ? "light" : "dark";
    setTheme(next);
    try {
      window.localStorage.setItem(KEY, next);
    } catch {
      return;
    }
  }, [isDark]);
  return { theme, isDark, toggle };
}
