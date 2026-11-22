import type { AuthOpts } from "../api";

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