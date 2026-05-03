import type { AuthOpts } from "../api";

/**
 * Build an `AuthOpts` object without leaking `undefined` keys — required by
 * TypeScript's `exactOptionalPropertyTypes: true`. Hoisted to its own module
 * so `lib/queue/mutators/*` can use it without importing from `queries.ts`
 * (which would create a circular import: queries imports queue mutators,
 * mutators would otherwise import authOpts back from queries).
 */
export function authOpts(
  username: string | undefined,
  clientId: string | undefined,
  clientOpId?: string,
): AuthOpts {
  const out: AuthOpts = {};
  if (username !== undefined) out.username = username;
  if (clientId !== undefined) out.clientId = clientId;
  if (clientOpId !== undefined) out.clientOpId = clientOpId;
  return out;
}
