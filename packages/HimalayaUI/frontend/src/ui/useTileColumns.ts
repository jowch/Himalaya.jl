import { useLayoutEffect, useState, type RefObject } from "react";
import { MOST_COLUMNS, tileColumns, tileWidth } from "./tileColumns";

const BOTTOM_PADDING = 16;

export interface TileFit {
  columns: number;
  tile: number;
}

export function useTileColumns(root: RefObject<HTMLElement | null>): TileFit {
  const [fit, setFit] = useState<TileFit>({ columns: MOST_COLUMNS, tile: 0 });
  useLayoutEffect(() => {
    const el = root.current;
    const scroller = el?.closest("main");
    if (!el || !scroller || typeof ResizeObserver === "undefined") return;
    const measure = (): void => {
      const grid = el.querySelector("ul") ?? el;
      const top = grid.getBoundingClientRect().top - scroller.getBoundingClientRect().top + scroller.scrollTop;
      const width = grid.clientWidth;
      const columns = tileColumns(width, scroller.clientHeight - top - BOTTOM_PADDING);
      const tile = Math.floor(tileWidth(width, columns));
      setFit((f) => (f.columns === columns && f.tile === tile ? f : { columns, tile }));
    };
    measure();
    const observer = new ResizeObserver(measure);
    observer.observe(el);
    observer.observe(scroller);
    return () => observer.disconnect();
  }, [root]);
  return fit;
}
