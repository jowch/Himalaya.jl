export interface NavbarProps {
  breadcrumb: string;
  username: string | undefined;
  onUserClick: () => void;
}

export function Navbar({ breadcrumb, username, onUserClick }: NavbarProps): JSX.Element {
  return (
    <header className="nav">
      <div className="nav-left">
        <span className="nav-logo">Himalaya</span>
        <nav className="nav-links">
          <a className="nav-link active" href="#">Analysis</a>
        </nav>
      </div>
      <div className="nav-center">
        <span className="nav-breadcrumb">{breadcrumb}</span>
      </div>
      <div className="nav-right">
        <button className="nav-user" onClick={onUserClick}>
          {username ?? "Sign in"}
        </button>
      </div>
    </header>
  );
}
