export interface GridCoord {
  row: number;
  col: number;
}

export interface GridDims {
  rows: number; // total rows INCLUDING the header row at index 0
  cols: number; // total columns
}

export interface NextCoordOpts {
  /** Rows moved by PageUp/PageDown. Default 10. */
  pageSize?: number;
  /** Ctrl/⌘ held — Home/End jump to grid corners instead of row ends. */
  ctrl?: boolean;
}

const clamp = (v: number, max: number): number => Math.max(0, Math.min(max, v));

/**
 * Compute the next focused cell coordinate for a keyboard navigation key in a
 * roving-tabindex grid. Clamps at all edges (NO wrap-around). Returns `null`
 * when `key` is not a navigation key the grid consumes.
 */
export function nextGridCoord(
  coord: GridCoord,
  key: string,
  dims: GridDims,
  opts: NextCoordOpts = {},
): GridCoord | null {
  const pageSize = opts.pageSize ?? 10;
  const ctrl = opts.ctrl ?? false;
  const maxRow = Math.max(0, dims.rows - 1);
  const maxCol = Math.max(0, dims.cols - 1);

  let { row, col } = coord;
  switch (key) {
    case "ArrowUp":
      row -= 1;
      break;
    case "ArrowDown":
      row += 1;
      break;
    case "ArrowLeft":
      col -= 1;
      break;
    case "ArrowRight":
      col += 1;
      break;
    case "Home":
      col = 0;
      if (ctrl) row = 0;
      break;
    case "End":
      col = maxCol;
      if (ctrl) row = maxRow;
      break;
    case "PageUp":
      row -= pageSize;
      break;
    case "PageDown":
      row += pageSize;
      break;
    default:
      return null;
  }

  return { row: clamp(row, maxRow), col: clamp(col, maxCol) };
}
