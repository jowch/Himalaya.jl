import { test } from "node:test";
import assert from "node:assert/strict";
import { frameFromLabels, keptFromLabels, openSampleBy, statusFromText, workspaceStart } from "../views.mjs";
import { tracked } from "../settle.mjs";

test("the pressed thumbnail's frame number maps to the exposure list order", () => {
  const order = ["a_S1", "a_S2", "a_S3"];
  assert.equal(frameFromLabels(["Frame 2, representative"], order), "a_S2");
  assert.equal(frameFromLabels(["Frame 3, dropped"], order), "a_S3");
  assert.equal(frameFromLabels([], order), null);
  assert.equal(frameFromLabels(["Frame 1", "Frame 2"], order), null);
  assert.equal(frameFromLabels(["Frame 9"], order), null);
});

test("the workspace's frame cards give kept / total", () => {
  assert.equal(keptFromLabels(["Frame 1", "Frame 2, dropped", "Frame 3"]), "2/3");
  assert.equal(keptFromLabels(["Frame 1, dropped"]), "0/1");
  assert.equal(keptFromLabels(["Drop", "Restore"]), null);
  assert.equal(keptFromLabels([]), null);
});

test("status text maps onto the sample states", () => {
  assert.equal(statusFromText("Not indexed"), "unreviewed");
  assert.equal(statusFromText("Unreviewed"), "unreviewed");
  assert.equal(statusFromText("Index →"), "unreviewed");
  assert.equal(statusFromText("Form factor"), "form_factor");
  assert.equal(statusFromText("No exposures"), "no_exposures");
  assert.equal(statusFromText(" Pn3m "), "phase");
  assert.equal(statusFromText(null), null);
});

test("settle ignores SSE, fonts, the favicon and in-page blob or data URLs", () => {
  const req = (type, url) => ({ resourceType: () => type, url: () => url });
  assert.equal(tracked(req("fetch", "http://h/api/samples")), true);
  assert.equal(tracked(req("image", "http://h/api/exposures/1/image")), true);
  assert.equal(tracked(req("eventsource", "http://h/api/events")), false);
  assert.equal(tracked(req("font", "http://h/x.woff2")), false);
  assert.equal(tracked(req("image", "blob:http://h/1234")), false);
  assert.equal(tracked(req("image", "data:image/png;base64,AA")), false);
  assert.equal(tracked(req("other", "http://h/favicon.svg")), false);
  assert.equal(tracked(req("image", "http://h/favicon.ico?v=2")), false);
});

test("the workspace is read in place when the shared tab is already on it", () => {
  const ui = "http://localhost:5173";
  assert.deepEqual(workspaceStart(`${ui}/experiments/1/workspace/7`, ui, 1, 3), { reuse: true, route: null });
  assert.deepEqual(workspaceStart(`${ui}/experiments/1/workspace`, ui, 1, 3), { reuse: true, route: null });
  assert.deepEqual(workspaceStart(`${ui}/experiments/2/workspace/7`, ui, 1, 3), { reuse: false, route: "/experiments/1/workspace/3" });
  assert.deepEqual(workspaceStart(`${ui}/experiments/1/corpus`, ui, 1, 3), { reuse: false, route: "/experiments/1/workspace/3" });
  assert.deepEqual(workspaceStart("about:blank", ui, 1, 3), { reuse: false, route: "/experiments/1/workspace/3" });
  assert.deepEqual(workspaceStart("http://elsewhere:9000/experiments/1/workspace/7", ui, 1, 3), { reuse: false, route: "/experiments/1/workspace/3" });
});

test("a sample that is already open is reopened by switching away and back, so the frame read is the one it opens on, not the one the user stepped to", () => {
  const ui = "http://127.0.0.1:5173";
  assert.deepEqual(openSampleBy(`${ui}/experiments/1/workspace/7?exposure=71`, 1, 7), { reopen: true });
  assert.deepEqual(openSampleBy(`${ui}/experiments/1/workspace/7`, 1, 7), { reopen: true });
  assert.deepEqual(openSampleBy(`${ui}/experiments/1/workspace/7`, 1, 8), { click: true });
  assert.deepEqual(openSampleBy(`${ui}/experiments/2/workspace/7`, 1, 7), { click: true });
  assert.deepEqual(openSampleBy("about:blank", 1, 7), { click: true });
});
