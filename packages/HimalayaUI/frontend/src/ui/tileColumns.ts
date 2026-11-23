export const TILE_GAP = 8;
export const TILE_ROW_CHROME = 36;
export const FULL_IMAGE_TILE = 160;
const FEWER_COLUMNS = [2, 3, 4] as const;
export const MOST_COLUMNS = 6;

export function tileWidth(width: number, columns: number): number {
  return (width - TILE_GAP * (columns - 1)) / columns;
}

export function tileColumns(width: number, freeHeight: number): number {
  if (width <= 0) return MOST_COLUMNS;
  for (const columns of FEWER_COLUMNS) {
    const twoRows = 2 * (tileWidth(width, columns) + TILE_ROW_CHROME) + TILE_GAP;
    if (twoRows <= freeHeight) return columns;
  }
  return MOST_COLUMNS;
}
