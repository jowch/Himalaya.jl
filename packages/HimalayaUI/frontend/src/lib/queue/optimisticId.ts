import { getClientId } from "../clientId";

let cachedTabFingerprint: number | null = null;
function tabFingerprint(): number {
  if (cachedTabFingerprint !== null) return cachedTabFingerprint;
  const id = getClientId();
  let h = 0x811c9dc5;
  for (let i = 0; i < id.length; i++) {
    h ^= id.charCodeAt(i);
    h = Math.imul(h, 0x01000193);
  }
  const masked = (h >>> 8) & 0xffffff;
  cachedTabFingerprint = masked === 0 ? 1 : masked;
  return cachedTabFingerprint;
}
let optimisticCounter = 0;
export function nextOptimisticId(): number {
  optimisticCounter += 1;
  return -(tabFingerprint() * 0x10000000 + optimisticCounter);
}
export function __resetOptimisticIdForTest(): void {
  optimisticCounter = 0;
  cachedTabFingerprint = null;
}