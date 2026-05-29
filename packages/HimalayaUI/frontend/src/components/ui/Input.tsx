import type { InputHTMLAttributes } from "react";

interface InputProps extends InputHTMLAttributes<HTMLInputElement> {
  className?: string;
}

export function Input({ className = "", ...props }: InputProps): JSX.Element {
  return (
    <input
      className={
        "bg-paper border border-hair-strong rounded-md px-2 py-1 " +
        "focus:outline focus:outline-1 focus:outline-accent focus:border-accent " +
        "transition-colors " +
        className
      }
      {...props}
    />
  );
}
