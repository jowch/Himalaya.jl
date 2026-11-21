import { Overlay } from "./Overlay";

const CLOSE_KEYS = ["Escape", " "] as const;

interface Props {
  src: string;
  caption: string;
  onClose: () => void;
}

export function DetectorOverlay({ src, caption, onClose }: Props): JSX.Element {
  return (
    <Overlay label="Detector image" onClose={onClose} testId="workspace-detector-big" closeKeys={CLOSE_KEYS}>
      <div className="flex flex-col items-center gap-2 p-3">
        <img src={src} alt="" className="aspect-square w-[min(420px,80vw)] rounded-[4px] bg-frame-edge object-contain" />
        <p className="m-0 text-sm text-ink-faint">{caption}</p>
      </div>
    </Overlay>
  );
}
