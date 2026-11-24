import { useEffect } from "react";
import { Navigate } from "react-router";
import { showToast } from "../lib/toast";

export function LinkGone({ message }: { message: string }): JSX.Element {
  useEffect(() => { showToast(message, "warning"); }, [message]);
  return <Navigate to="/experiments" replace />;
}
