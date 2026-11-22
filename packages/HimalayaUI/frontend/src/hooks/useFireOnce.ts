import { useEffect, useRef } from "react";

export function useFireOnce<T extends object>(token: T | undefined, ready: boolean, fire: (token: T) => void): void {
  const firedFor = useRef<T | undefined>(undefined);
  useEffect(() => {
    if (token === undefined || !ready || firedFor.current === token) return;
    firedFor.current = token;
    fire(token);
  }, [token, ready, fire]);
}
