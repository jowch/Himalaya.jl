const VALID_TYPES = ["peak", "index", "exposure", "sample", "experiment", "comparison"] as const;
type MentionType = typeof VALID_TYPES[number];

export type TextSegment  = { kind: "text"; text: string };
/**
 * `hash` is the 8-char truncated SHA-256 of a comparison's `content_hash` at
 * the moment the mention was inserted (eager-hash policy, see Phase 10).
 * Only present on `comparison` mentions; `parseMentions` populates it when
 * the input token has the `[[comparison:N@hhhhhhhh]]` form.
 */
export type MentionToken = {
  kind: "mention";
  type: MentionType;
  id: number;
  hash?: string;
};
export type BodySegment  = TextSegment | MentionToken;

// Token grammar:
//   [[<type>:<id>]]                — legacy + non-comparison types
//   [[<type>:<id>@<hash8>]]        — comparison eager-hash form
//
// `hash8` must be 8 lowercase hex chars; anything else falls through to the
// invalid-token branch (rendered as literal text). The strict regex here is
// load-bearing: a too-loose match would silently accept `@abc` or `@FOO12345`
// as a hash and the drift indicator would compare against garbage.
const TOKEN_RE = /\[\[(\w+):(\d+)(?:@([0-9a-f]{8}))?\]\]/g;

export function parseMentions(body: string): BodySegment[] {
  if (!body) return [];
  const segments: BodySegment[] = [];
  let last = 0;

  for (const match of body.matchAll(TOKEN_RE)) {
    const type = match[1];
    const id   = parseInt(match[2], 10);
    const hash = match[3]; // undefined when no `@hash` suffix
    if (!(VALID_TYPES as readonly string[]).includes(type) || Number.isNaN(id)) {
      // Treat invalid token as literal text — include the matched text in output
      const textChunk = body.slice(last, match.index! + match[0].length);
      if (textChunk) segments.push({ kind: "text", text: textChunk });
      last = match.index! + match[0].length;
      continue;
    }
    if (match.index! > last) {
      segments.push({ kind: "text", text: body.slice(last, match.index) });
    }
    const token: MentionToken = { kind: "mention", type: type as MentionType, id };
    if (hash !== undefined) token.hash = hash;
    segments.push(token);
    last = match.index! + match[0].length;
  }

  if (last < body.length) {
    segments.push({ kind: "text", text: body.slice(last) });
  }
  return segments;
}
