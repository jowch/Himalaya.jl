/** A tag: a key, optionally with a value. Backend model — some tags are
 *  key-only ("LL37"), some key+value ("temperature" = "37C"). */
export interface Tag {
  key: string;
  value?: string;
}
