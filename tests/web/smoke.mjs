/**
 * Headless smoke test for the browser demo.
 *
 * It runs the real `app.js` against a real DOM, with only the browser APIs jsdom lacks
 * stubbed out: canvas drawing and `fetch`. That is enough to catch what unit tests
 * cannot -- a mistyped element id, a property that does not exist, a listener wired to
 * the wrong event -- any of which would leave Elias looking at a blank page.
 *
 *   cd tests/web && npm install && npm test
 */

import { readFileSync } from "node:fs";
import { fileURLToPath, pathToFileURL } from "node:url";
import { dirname, resolve } from "node:path";
import { JSDOM } from "jsdom";

const here = dirname(fileURLToPath(import.meta.url));
const webDir = resolve(here, "../../web");

let failures = 0;
function check(label, condition, detail = "") {
  const ok = Boolean(condition);
  if (!ok) failures++;
  console.log(`${ok ? "  ok  " : "FAIL  "}${label}${detail ? ` — ${detail}` : ""}`);
}

/** A 2D context that records what was asked of it, so we can prove drawing happened. */
function recordingContext() {
  const calls = { stroke: 0, fill: 0, arc: 0, moveTo: 0, lineTo: 0, fillText: 0 };
  const noop = () => {};
  return new Proxy(
    {
      calls,
      canvas: null,
      setTransform: noop,
      clearRect: noop,
      beginPath: noop,
      closePath: noop,
      save: noop,
      restore: noop,
      setLineDash: noop,
      strokeRect: noop,
      measureText: () => ({ width: 10 }),
      stroke: () => calls.stroke++,
      fill: () => calls.fill++,
      arc: () => calls.arc++,
      moveTo: () => calls.moveTo++,
      lineTo: () => calls.lineTo++,
      fillText: () => calls.fillText++,
    },
    {
      get: (target, key) => (key in target ? target[key] : noop),
      set: (target, key, value) => ((target[key] = value), true),
    },
  );
}

const dom = new JSDOM(readFileSync(resolve(webDir, "index.html"), "utf8"), {
  url: "https://example.invalid/optic/",
  pretendToBeVisual: true,
});
const { window } = dom;

// Browser APIs jsdom does not provide.
const contexts = new Map();
window.HTMLCanvasElement.prototype.getContext = function getContext() {
  if (!contexts.has(this)) contexts.set(this, recordingContext());
  return contexts.get(this);
};
window.HTMLCanvasElement.prototype.getBoundingClientRect = () => ({
  width: 820,
  height: 300,
  top: 0,
  left: 0,
});

globalThis.window = window;
globalThis.document = window.document;
globalThis.getComputedStyle = window.getComputedStyle.bind(window);
globalThis.requestAnimationFrame = (cb) => setTimeout(cb, 0);
globalThis.fetch = async (url) => {
  const body = readFileSync(resolve(webDir, String(url)));
  return { ok: true, status: 200, statusText: "OK", arrayBuffer: async () => body };
};

const settle = () => new Promise((r) => setTimeout(r, 12));

const app = await import(pathToFileURL(resolve(webDir, "app.js")).href);
await app.ready;
await settle();

const doc = window.document;
const text = (id) => doc.getElementById(id).textContent;
const efl = () => {
  const match = text("first-order").match(/Focal length([\-\d.]+) mm/);
  return match ? Number(match[1]) : NaN;
};

console.log("\nbrowser demo smoke test\n");

check("the kernel loaded without error", doc.getElementById("error").hidden);
check("the preset list is populated", doc.getElementById("preset").options.length >= 2);
check(
  "the prescription table has a row per surface",
  doc.getElementById("prescription").rows.length === 8,
  `${doc.getElementById("prescription").rows.length} rows`,
);
check("first-order data is shown", text("first-order").includes("Focal length"));
check("the Cooke triplet focal length is right", Math.abs(efl() - 50.021) < 0.01, `EFL ${efl()}`);
check(
  "a spot diagram is shown per field",
  doc.getElementById("spots").querySelectorAll("figure").length === 3,
);
check(
  "the layout canvas was actually drawn",
  contexts.get(doc.getElementById("layout"))?.calls.stroke > 20,
);
check(
  "the spot canvases were actually drawn",
  contexts.get(doc.getElementById("spot-0"))?.calls.arc > 100,
);
check(
  "spot metrics are reported",
  /RMS [\d.]+/.test(doc.querySelector("#spots .metrics").innerHTML),
  doc.querySelector("#spots .metrics").textContent.trim(),
);

// Editing a radius must move the focal length.
const radius = doc.querySelector('#prescription input[data-field="radius"][data-row="0"]');
const before = efl();
radius.value = "25";
radius.dispatchEvent(new window.Event("input", { bubbles: true }));
await settle();
check("editing a radius re-analyses the system", Math.abs(efl() - before) > 1, `${before} -> ${efl()}`);

radius.value = "22.01359";
radius.dispatchEvent(new window.Event("input", { bubbles: true }));
await settle();
check("restoring the radius restores the focal length", Math.abs(efl() - before) < 1e-6);

// Defocus must blur the spot WITHOUT editing the design: it is a question asked of the
// prescription, not a change to it, and it must not fight the autofocus solve.
const thicknessCell = () =>
  Number(doc.querySelector('#prescription input[data-field="thickness"][data-row="6"]').value);
const rms = () => {
  const m = doc.querySelector("#spots .metrics").textContent.match(/RMS ([\d.]+)/);
  return m ? Number(m[1]) : NaN;
};

const focusBefore = thicknessCell();
const sharp = rms();
const slider = doc.getElementById("defocus");
slider.value = "1.5";
slider.dispatchEvent(new window.Event("input", { bubbles: true }));
await settle();
check("defocus blurs the spot", rms() > sharp * 2, `${sharp} -> ${rms()}`);
check("defocus leaves the prescription alone", thicknessCell() === focusBefore, `${thicknessCell()}`);
check("the defocus readout updates", text("defocus-value").includes("1.50"));

doc.getElementById("refocus").click();
await settle();
check("refocus restores the sharp spot", Math.abs(rms() - sharp) < 0.1, `${rms()}`);

// Solves must be visible and editable — Elias asked for exactly this.
const solveSelect = doc.querySelector('#prescription select[data-field="solve"][data-row="6"]');
check("the image distance shows its solve", solveSelect?.value === "marginal_ray_height", solveSelect?.value);
check(
  "a solved thickness is read-only",
  doc.querySelector('#prescription input[data-field="thickness"][data-row="6"]').readOnly,
);

solveSelect.value = "fixed";
solveSelect.dispatchEvent(new window.Event("change", { bubbles: true }));
await settle();
check(
  "releasing the solve makes the thickness editable",
  !doc.querySelector('#prescription input[data-field="thickness"][data-row="6"]').readOnly,
);
solveSelect.value = "marginal_ray_height";
solveSelect.dispatchEvent(new window.Event("change", { bubbles: true }));
await settle();
check("the solve can be put back", Math.abs(thicknessCell() - focusBefore) < 1e-6, `${thicknessCell()}`);

// The primary wavelength must be shown and changeable.
const primary = doc.getElementById("primary");
check("the primary wavelength is offered", primary.options.length === 3, `${primary.options.length}`);
check("it defaults to the middle line", text("first-order").includes("0.5876"));
primary.value = "0";
primary.dispatchEvent(new window.Event("change", { bubbles: true }));
await settle();
check("changing the primary wavelength takes effect", text("first-order").includes("0.4861"));
primary.value = "1";
primary.dispatchEvent(new window.Event("change", { bubbles: true }));
await settle();

// Fields must say which way they point. The on-axis field has no direction, so the
// label to inspect is one of the off-axis ones.
const offAxisLabel = () => [...doc.querySelectorAll(".field-label")].pop().textContent;
check("field labels name their axis", /Y\s*2?0\.0/.test(offAxisLabel()), offAxisLabel());
const axis = doc.getElementById("field-axis");
axis.value = "x";
axis.dispatchEvent(new window.Event("change", { bubbles: true }));
await settle();
check("switching the field axis to X is reflected", /X\s*2?0\.0/.test(offAxisLabel()), offAxisLabel());
axis.value = "y";
axis.dispatchEvent(new window.Event("change", { bubbles: true }));
await settle();

// Pupil sampling patterns.
const pattern = doc.getElementById("pupil-pattern");
const rmsSquare = rms();
for (const p of ["hexapolar", "dithered"]) {
  pattern.value = p;
  pattern.dispatchEvent(new window.Event("change", { bubbles: true }));
  await settle();
  check(`${p} sampling produces a spot`, Number.isFinite(rms()) && rms() > 1, `RMS ${rms()}`);
  check(`${p} agrees with the square grid to within 10%`,
    Math.abs(rms() - rmsSquare) / rmsSquare < 0.1, `${rmsSquare} vs ${rms()}`);
}
pattern.value = "square";
pattern.dispatchEvent(new window.Event("change", { bubbles: true }));
await settle();

// The Airy disc toggle.
const airy = doc.getElementById("show-airy");
check("the Airy disc is on by default", airy.checked);
const arcsWith = contexts.get(doc.getElementById("spot-0")).calls.arc;
airy.checked = false;
airy.dispatchEvent(new window.Event("change", { bubbles: true }));
await settle();
check("toggling the Airy disc redraws", contexts.get(doc.getElementById("spot-0")).calls.arc > arcsWith);

// Distortion, per field and wavelength.
const distortion = doc.getElementById("distortion");
check("distortion is tabulated", !distortion.hidden && distortion.rows.length >= 3,
  `${distortion.rows.length} rows`);
check("distortion covers every wavelength", distortion.rows[0].cells.length === 4,
  `${distortion.rows[0].cells.length} columns`);

// Switching preset must rebuild everything.
const preset = doc.getElementById("preset");
preset.value = "1";
preset.dispatchEvent(new window.Event("change", { bubbles: true }));
await settle();
check("switching preset reloads the system", Math.abs(efl() - 97.58) < 0.01, `EFL ${efl()}`);
check(
  "the table shrinks to the singlet's three surfaces",
  doc.getElementById("prescription").rows.length === 3,
);

// A broken prescription must surface a readable error, not a blank page.
const glass = doc.querySelector('#prescription select[data-field="glass"][data-row="0"]');
glass.value = "MIRROR";
glass.dispatchEvent(new window.Event("change", { bubbles: true }));
await settle();
check("a mirror in place of glass still analyses", doc.getElementById("error").hidden);

console.log(`\n${failures === 0 ? "all checks passed" : `${failures} check(s) failed`}\n`);
process.exit(failures === 0 ? 0 : 1);
