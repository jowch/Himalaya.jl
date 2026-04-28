import { useAppState } from "../state";

function userInitials(username: string | undefined): string {
  if (!username) return "?";
  const parts = username.split(/[\s_.-]+/).filter(Boolean);
  if (parts.length === 0) return "?";
  if (parts.length === 1) return parts[0]!.slice(0, 2).toUpperCase();
  return ((parts[0]![0] ?? "") + (parts[1]![0] ?? "")).toUpperCase();
}

/**
 * UtilityCluster — small cluster of app-wide utility controls in the top-right.
 * Currently: theme toggle + user-initials avatar (clicking opens the user picker).
 */
export function UtilityCluster(): JSX.Element {
  const theme         = useAppState((s) => s.theme);
  const setTheme      = useAppState((s) => s.setTheme);
  const username      = useAppState((s) => s.username);
  const clearUsername = useAppState((s) => s.clearUsername);

  const toggleTheme = (): void => setTheme(theme === "dark" ? "light" : "dark");
  const switchUser  = (): void => { clearUsername(); };

  const initials = userInitials(username);

  return (
    <div className="flex items-center gap-1" data-testid="utility-cluster">
      <button
        type="button"
        data-testid="theme-toggle"
        onClick={toggleTheme}
        aria-label="Toggle theme"
        title={`Theme: ${theme}`}
        className="w-7 h-7 grid place-items-center rounded-full
                   bg-transparent border border-transparent
                   text-fg-muted hover:text-fg hover:border-border hover:bg-bg-elevated"
      >
        {theme === "dark" ? (
          <svg viewBox="0 0 16 16" width="14" height="14" aria-hidden="true">
            <path d="M12 9 A5 5 0 1 1 7 4 A4 4 0 0 0 12 9 Z" fill="currentColor" />
          </svg>
        ) : (
          <svg viewBox="0 0 16 16" width="14" height="14" aria-hidden="true">
            <circle cx="8" cy="8" r="3" fill="currentColor" />
            <g stroke="currentColor" strokeWidth="1.2" strokeLinecap="round">
              <line x1="8" y1="1" x2="8" y2="3" />
              <line x1="8" y1="13" x2="8" y2="15" />
              <line x1="1" y1="8" x2="3" y2="8" />
              <line x1="13" y1="8" x2="15" y2="8" />
            </g>
          </svg>
        )}
      </button>
      <button
        type="button"
        data-testid="user-avatar"
        onClick={switchUser}
        aria-label="Switch user"
        title={username ?? "Sign in"}
        className="w-7 h-7 rounded-full font-sans text-xs font-semibold
                   text-white grid place-items-center bg-gradient-to-br
                   from-[oklch(0.68_0.05_240)] to-[oklch(0.56_0.06_200)]
                   hover:brightness-110"
      >
        {initials}
      </button>
    </div>
  );
}
