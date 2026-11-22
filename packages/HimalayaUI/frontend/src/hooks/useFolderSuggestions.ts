import { useEffect, useState } from "react";
import * as api from "../api";

export const NO_FOLDERS: api.FolderSuggestions = { folders: [], more: false };

export const isFullPath = (p: string): boolean => p.startsWith("/") || p === "~" || p.startsWith("~/");

export function useFolderSuggestions(path: string, enabled = true): api.FolderSuggestions {
  const [suggested, setSuggested] = useState<api.FolderSuggestions>(NO_FOLDERS);
  const live = enabled && isFullPath(path);
  useEffect(() => {
    if (!live) { setSuggested(NO_FOLDERS); return; }
    let current = true;
    void api.suggestFolders(path)
      .then((r) => { if (current) setSuggested(r); })
      .catch(() => { if (current) setSuggested(NO_FOLDERS); });
    return () => { current = false; };
  }, [path, live]);
  return suggested;
}
