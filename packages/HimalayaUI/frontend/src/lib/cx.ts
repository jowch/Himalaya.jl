// Join class-name parts, dropping falsy ones. Placement-only className helper
// shared across The Print's primitives and components.
export function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}
