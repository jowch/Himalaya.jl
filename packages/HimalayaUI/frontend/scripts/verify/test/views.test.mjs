import { test } from "node:test";
import assert from "node:assert/strict";
import { frameFromLabels, statusFromText } from "../views.mjs";
import { tracked } from "../settle.mjs";

test("the pressed thumbnail's frame number maps to the exposure list order", () => {
  const order = ["a_S1", "a_S2", "a_S3"];
  assert.equal(frameFromLabels(["Frame 2, representative"], order), "a_S2");
  assert.equal(frameFromLabels(["Frame 3, dropped"], order), "a_S3");
  assert.equal(frameFromLabels([], order), null);
  assert.equal(frameFromLabels(["Frame 1", "Frame 2"], order), null);
  assert.equal(frameFromLabels(["Frame 9"], order), null);
});

test("corpus status text maps onto the StatusCell states", () => {
  assert.equal(statusFromText("Not indexed"), "unreviewed");
  assert.equal(statusFromText("Unreviewed"), "unreviewed");
  assert.equal(statusFromText("Index →"), "unreviewed");
  assert.equal(statusFromText("Form factor"), "form_factor");
  assert.equal(statusFromText("No exposures"), "no_exposures");
  assert.equal(statusFromText(" Pn3m "), "phase");
  assert.equal(statusFromText(null), null);
});

test("settle ignores SSE, fonts and in-page blob or data URLs", () => {
  const req = (type, url) => ({ resourceType: () => type, url: () => url });
  assert.equal(tracked(req("fetch", "http://h/api/samples")), true);
  assert.equal(tracked(req("image", "http://h/api/exposures/1/image")), true);
  assert.equal(tracked(req("eventsource", "http://h/api/events")), false);
  assert.equal(tracked(req("font", "http://h/x.woff2")), false);
  assert.equal(tracked(req("image", "blob:http://h/1234")), false);
  assert.equal(tracked(req("image", "data:image/png;base64,AA")), false);
});
