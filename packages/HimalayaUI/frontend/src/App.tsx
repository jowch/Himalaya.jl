import { useEffect, useState } from "react";
import "./styles.css";
import * as api from "./api";
import type { Experiment, Sample } from "./api";
import { useAppState } from "./state";
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

  const [experiment, setExperiment] = useState<Experiment | null>(null);
  const [samples, setSamples]       = useState<Sample[]>([]);
  const [bootError, setBootError]   = useState<string | null>(null);
  const [modalOpen, setModalOpen]   = useState<boolean>(!username);

  useEffect(() => {
    void (async () => {
      try {
        const [exp, list] = await Promise.all([
          api.getExperiment(EXPERIMENT_ID),
          api.listSamples(EXPERIMENT_ID),
        ]);
        setExperiment(exp);
        setSamples(list);
      } catch (e) {
        setBootError((e as Error).message);
      }
    })();
  }, []);

  const activeSample = samples.find((s) => s.id === activeSampleId);
  const breadcrumb   = bootError
    ? `Error: ${bootError}`
    : (experiment?.name ?? "experiment")
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
