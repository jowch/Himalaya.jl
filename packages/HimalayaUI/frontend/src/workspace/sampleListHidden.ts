import { useCallback, useState } from "react";

const KEY = "himalaya-ui:workspace-sample-list-hidden";

function read(): boolean {
  try {
    return window.localStorage.getItem(KEY) === "1";
  } catch {
    return false;
  }
}

export function useSampleListHidden(): [boolean, (hidden: boolean) => void] {
  const [hidden, setHidden] = useState<boolean>(read);
  const choose = useCallback((next: boolean) => {
    setHidden(next);
    try {
      if (next) window.localStorage.setItem(KEY, "1");
      else window.localStorage.removeItem(KEY);
    } catch {
      return;
    }
  }, []);
  return [hidden, choose];
}
