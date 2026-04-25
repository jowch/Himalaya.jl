interface Pn3mIconProps {
  size?: number;
  className?: string;
}

/**
 * Single Im3m (Schwarz P / "plumber's nightmare") unit cell, viewed from [111]
 * in isometric projection — a cube with one elliptical channel opening per face.
 *
 * Geometry — 24×24 viewBox, hexagon circumradius = 9:
 *   Outer vertices: top (12,3), tr (19.79,7.5), br (19.79,16.5),
 *                   bot (12,21), bl (4.21,16.5), tl (4.21,7.5)
 *   Near corner (closest to viewer): (12,12)
 *   Face centers: top-face (12,7.5), right-face (15.9,14.25), left-face (8.1,14.25)
 *
 *   Three face fills at stepped opacity encode depth (top lightest, right darkest).
 *   Each face carries a filled ellipse — the channel opening — at higher opacity
 *   so it reads as a darker hole looking into the bicontinuous tunnel.
 */
export function Pn3mIcon({ size = 20, className = "" }: Pn3mIconProps): JSX.Element {
  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      viewBox="0 0 24 24"
      width={size}
      height={size}
      aria-label="Im3m bicontinuous cubic unit cell"
      role="img"
      className={className}
    >
      <title>Im3m bicontinuous cubic unit cell</title>

      {/* ── Face fills — stepped opacity gives isometric depth ── */}
      <polygon points="12,3 19.79,7.5 12,12 4.21,7.5"
               fill="currentColor" opacity="0.22"/>
      <polygon points="19.79,7.5 19.79,16.5 12,21 12,12"
               fill="currentColor" opacity="0.48"/>
      <polygon points="4.21,7.5 12,12 12,21 4.21,16.5"
               fill="currentColor" opacity="0.35"/>

      {/* ── Cube edges ── */}
      <polygon
        points="12,3 19.79,7.5 19.79,16.5 12,21 4.21,16.5 4.21,7.5"
        fill="none" stroke="currentColor" strokeWidth="1.25" strokeLinejoin="round"
      />
      <line x1="12"    y1="3"    x2="12"    y2="12"   stroke="currentColor" strokeWidth="1.25"/>
      <line x1="19.79" y1="16.5" x2="12"    y2="12"   stroke="currentColor" strokeWidth="1.25"/>
      <line x1="4.21"  y1="16.5" x2="12"    y2="12"   stroke="currentColor" strokeWidth="1.25"/>

      {/* ── Channel openings — filled ellipses (one per visible face)        ──
           Higher opacity than face fill → darker dot = looking into a tunnel. ── */}
      <ellipse cx="12"   cy="7.5"  rx="2.6" ry="1.3"
               fill="currentColor" opacity="0.72"/>
      <ellipse cx="15.9" cy="14.25" rx="2.6" ry="1.3"
               fill="currentColor" opacity="0.72"
               transform="rotate(-30 15.9 14.25)"/>
      <ellipse cx="8.1"  cy="14.25" rx="2.6" ry="1.3"
               fill="currentColor" opacity="0.72"
               transform="rotate(30 8.1 14.25)"/>
    </svg>
  );
}
