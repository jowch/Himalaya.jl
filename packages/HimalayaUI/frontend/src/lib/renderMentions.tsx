const VALID_TYPES = ["peak", "index", "exposure", "sample", "experiment"] as const;
type MentionType = typeof VALID_TYPES[number];

export type TextSegment  = { kind: "text"; text: string };
export type MentionToken = { kind: "mention"; type: MentionType; id: number };
export type BodySegment  = TextSegment | MentionToken;

const TOKEN_RE = /\[\[(\w+):(\d+)\]\]/g;

export function parseMentions(body: string): BodySegment[] {
  if (!body) return [];
  const segments: BodySegment[] = [];
  let last = 0;

  for (const match of body.matchAll(TOKEN_RE)) {
    const type = match[1];
    const id   = parseInt(match[2], 10);
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
    segments.push({ kind: "mention", type: type as MentionType, id });
    last = match.index! + match[0].length;
  }

  if (last < body.length) {
    segments.push({ kind: "text", text: body.slice(last) });
  }
  return segments;
}
