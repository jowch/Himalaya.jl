export interface NavbarProps {
  breadcrumb: string;
  username: string | undefined;
  onUserClick: () => void;
}

export function Navbar({ breadcrumb, username, onUserClick }: NavbarProps): JSX.Element {
  return (
    <header className="grid grid-cols-[auto_1fr_auto] items-center h-11 px-4 border-b border-border bg-bg">
      <div className="flex items-center gap-6">
        <span className="font-semibold tracking-wide" data-testid="nav-logo">Himalaya</span>
        <nav className="flex gap-4">
          <a className="text-fg border-b border-accent px-1 py-0.5" href="#">Analysis</a>
        </nav>
      </div>
      <div className="text-center text-fg-muted">
        <span className="font-mono text-[13px]" data-testid="nav-breadcrumb">{breadcrumb}</span>
      </div>
      <div className="flex">
        <button
          className="min-w-20 border border-border rounded-md px-2.5 py-1 hover:bg-bg-hover focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent"
          data-testid="nav-user"
          onClick={onUserClick}
        >
          {username ?? "Sign in"}
        </button>
      </div>
    </header>
  );
}
