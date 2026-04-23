import { useState } from "react";
import "./styles.css";
import { useAppState } from "./state";
import { useExperiment, useSamples } from "./queries";
import { Navbar } from "./components/Navbar";
import { Layout } from "./components/Layout";
import { SampleList } from "./components/SampleList";
import { UserModal } from "./components/UserModal";

const EXPERIMENT_ID = 1;

export function App(): JSX.Element {
  const username        = useAppState((s) => s.username);
  const activeSampleId  = useAppState((s) => s.activeSampleId);
  const setUsername     = useAppState((s) => s.setUsername);
  const setActiveSample = useAppState((s) => s.setActiveSample);

  const [modalOpen, setModalOpen] = useState<boolean>(!username);

  const experimentQ = useExperiment(EXPERIMENT_ID);
  const samplesQ    = useSamples(EXPERIMENT_ID);

  const samples      = samplesQ.data ?? [];
  const activeSample = samples.find((s) => s.id === activeSampleId);
  const bootError    = experimentQ.error ?? samplesQ.error;

  const breadcrumb = bootError
    ? `Error: ${(bootError as Error).message}`
    : (experimentQ.data?.name ?? "experiment")
        + (activeSample
          ? ` › ${activeSample.label ?? ""} ${activeSample.name ?? ""}`.trimEnd()
          : "");

  return (
    <>
      <Navbar
        breadcrumb={breadcrumb}
        username={username}
        onUserClick={() => setModalOpen(true)}
      />
      <Layout
        left={
          <SampleList
            samples={samples}
            activeId={activeSampleId}
            onSelect={setActiveSample}
          />
        }
      />
      <UserModal
        open={modalOpen}
        onSelect={(name) => {
          setUsername(name);
          setModalOpen(false);
        }}
        onClose={() => setModalOpen(false)}
      />
    </>
  );
}
