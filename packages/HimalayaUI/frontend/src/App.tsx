import { useEffect, useState } from "react";
import "./styles.css";
import { useAppState } from "./state";
import {
  useExperiment, useSamples, useExposures, useTrace, usePeaks,
  useAddPeak, useRemovePeak,
} from "./queries";
import { Navbar } from "./components/Navbar";
import { Layout } from "./components/Layout";
import { SampleList } from "./components/SampleList";
import { UserModal } from "./components/UserModal";
import { ExposureList } from "./components/ExposureList";
import { TraceViewer } from "./components/TraceViewer";
import { StaleIndicesBanner } from "./components/StaleIndicesBanner";

const EXPERIMENT_ID = 1;

export function App(): JSX.Element {
  const username          = useAppState((s) => s.username);
  const activeSampleId    = useAppState((s) => s.activeSampleId);
  const activeExposureId  = useAppState((s) => s.activeExposureId);
  const setUsername       = useAppState((s) => s.setUsername);
  const setActiveSample   = useAppState((s) => s.setActiveSample);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);

  const [modalOpen, setModalOpen] = useState<boolean>(!username);

  const experimentQ = useExperiment(EXPERIMENT_ID);
  const samplesQ    = useSamples(EXPERIMENT_ID);
  const exposuresQ  = useExposures(activeSampleId);
  const traceQ      = useTrace(activeExposureId);
  const peaksQ      = usePeaks(activeExposureId);

  const addPeak    = useAddPeak(activeExposureId ?? 0);
  const removePeak = useRemovePeak(activeExposureId ?? 0);

  // Auto-select the first exposure when samples or the active sample changes
  useEffect(() => {
    const exposures = exposuresQ.data ?? [];
    if (exposures.length === 0) return;
    const stillValid = exposures.some((e) => e.id === activeExposureId);
    if (!stillValid) setActiveExposure(exposures[0]!.id);
  }, [exposuresQ.data, activeExposureId, setActiveExposure]);

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
        centerTop={
          <div className="flex flex-col flex-1 min-h-0">
            <StaleIndicesBanner exposureId={activeExposureId} />
            {traceQ.data && peaksQ.data && activeExposureId !== undefined ? (
              <TraceViewer
                trace={traceQ.data}
                peaks={peaksQ.data}
                activeGroupIndices={[]}
                hoveredIndex={undefined}
                onAddPeak={(q) => addPeak.mutate(q)}
                onRemovePeak={(peakId) => removePeak.mutate(peakId)}
              />
            ) : (
              <p className="text-fg-muted italic flex-1 flex items-center justify-center">
                Select an exposure to view its trace.
              </p>
            )}
          </div>
        }
        centerBottom={<ExposureList />}
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
