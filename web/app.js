/**
 * The demo application.
 *
 * It owns a prescription, hands it to the engine on every edit, and draws whatever comes
 * back. There is no incremental update path and none is wanted: a full analysis of a
 * triplet takes a few milliseconds, so recomputing everything keeps what you see and
 * what the kernel believes in permanent agreement.
 */

import { OpticEngine } from "./engine.js";
import { drawLayout, drawSpot, fieldColour } from "./draw.js";

const state = {
  engine: null,
  presets: [],
  glasses: [],
  spec: null,
  analysis: null,
  /** Thickness before the image plane with the defocus slider centred. */
  focusBase: 0,
  pending: false,
};

const el = (id) => document.getElementById(id);

async function boot() {
  try {
    state.engine = await OpticEngine.fromUrl("optic_wasm.wasm");
  } catch (error) {
    el("error").textContent = `The optical kernel failed to load: ${error.message}`;
    el("error").hidden = false;
    return;
  }

  const { presets, glasses } = state.engine.presets();
  state.presets = presets;
  state.glasses = glasses;

  el("preset").innerHTML = presets
    .map((p, i) => `<option value="${i}">${escapeHtml(p.title)}</option>`)
    .join("");
  el("preset").addEventListener("change", (e) => loadPreset(Number(e.target.value)));

  el("epd").addEventListener("input", (e) => {
    state.spec.entrance_pupil_diameter = Math.max(0.01, Number(e.target.value) || 0.01);
    schedule();
  });
  el("fields").addEventListener("change", (e) => {
    const parsed = parseNumbers(e.target.value);
    state.spec.fields = parsed.length ? parsed : [0];
    e.target.value = state.spec.fields.join(", ");
    schedule();
  });
  el("wavelengths").addEventListener("change", (e) => {
    const parsed = parseNumbers(e.target.value).filter((w) => w > 0.1 && w < 20);
    state.spec.wavelengths = parsed.length ? parsed : [0.5875618];
    e.target.value = state.spec.wavelengths.map((w) => w.toFixed(4)).join(", ");
    schedule();
  });

  el("defocus").addEventListener("input", () => {
    const offset = Number(el("defocus").value);
    el("defocus-value").textContent = `${offset >= 0 ? "+" : ""}${offset.toFixed(2)} mm`;
    lastSurface().thickness = state.focusBase + offset;
    schedule({ redrawTable: true });
  });

  el("refocus").addEventListener("click", () => update({ refocus: true }));

  // Bound once: renderTable replaces the rows, not the table body, so attaching there
  // would stack a new listener on every redraw.
  el("prescription").addEventListener("input", onCellEdit);
  el("prescription").addEventListener("change", onCellEdit);

  window.addEventListener("resize", () => render());

  loadPreset(0);
}

function lastSurface() {
  return state.spec.surfaces[state.spec.surfaces.length - 2];
}

function loadPreset(index) {
  state.spec = structuredClone(state.presets[index]);
  el("epd").value = state.spec.entrance_pupil_diameter;
  el("fields").value = state.spec.fields.join(", ");
  el("wavelengths").value = state.spec.wavelengths.map((w) => w.toFixed(4)).join(", ");
  recentreDefocus();
  renderTable();
  update();
}

function recentreDefocus() {
  state.focusBase = lastSurface().thickness;
  el("defocus").value = 0;
  el("defocus-value").textContent = "+0.00 mm";
}

/** Coalesce rapid edits into one analysis per frame. */
function schedule(options = {}) {
  if (state.pending) return;
  state.pending = true;
  requestAnimationFrame(() => {
    state.pending = false;
    update(options);
  });
}

function update({ refocus = false, redrawTable = false } = {}) {
  const result = state.engine.analyze({
    system: state.spec,
    rays_per_fan: 15,
    spot_grid: 23,
    refocus,
  });

  if (!result.ok) {
    el("error").textContent = result.error;
    el("error").hidden = false;
    return;
  }
  el("error").hidden = true;

  // The engine echoes the prescription it actually used, so a refocus is visible.
  if (refocus) {
    state.spec = result.system;
    recentreDefocus();
    renderTable();
  } else if (redrawTable) {
    syncTableValues();
  }

  state.analysis = result;
  render();
  renderReadout();
  renderWarnings();
}

function render() {
  if (!state.analysis) return;
  drawLayout(el("layout"), state.analysis);

  const spots = state.analysis.spots;
  const halfWidth = Math.max(
    5,
    1.12 * spots.reduce((m, s) => Math.max(m, Number.isFinite(s.geometric) ? s.geometric : 0), 0),
  );

  const container = el("spots");
  if (container.children.length !== spots.length) {
    container.innerHTML = spots
      .map(
        (s, i) => `
        <figure class="spot">
          <canvas id="spot-${i}"></canvas>
          <figcaption>
            <span class="swatch" style="background:${fieldColour(i)}"></span>
            <strong>${s.field.toFixed(1)}&deg;</strong>
            <span class="metrics"></span>
          </figcaption>
        </figure>`,
      )
      .join("");
  }

  spots.forEach((spot, i) => {
    drawSpot(el(`spot-${i}`), spot, halfWidth, state.spec.wavelengths);
    const caption = container.children[i].querySelector(".metrics");
    const limited = spot.rms < spot.airy;
    caption.innerHTML = Number.isFinite(spot.rms)
      ? `RMS ${spot.rms.toFixed(1)} &micro;m${limited ? " &middot; diffraction limited" : ""}`
      : "no rays";
    caption.classList.toggle("good", limited);
  });
}

function renderReadout() {
  const f = state.analysis.first_order;
  const rows = [
    ["Focal length", `${f.efl.toFixed(3)} mm`],
    ["Working f/#", `f/${f.fno.toFixed(2)}`],
    ["Entrance pupil", `${f.epd.toFixed(2)} mm`],
    ["Back focal distance", `${f.bfd.toFixed(3)} mm`],
    ["Total track", `${f.total_track.toFixed(2)} mm`],
  ];
  el("first-order").innerHTML = rows
    .map(([k, v]) => `<div><dt>${k}</dt><dd>${v}</dd></div>`)
    .join("");
}

function renderWarnings() {
  const list = state.analysis.warnings;
  el("warnings").hidden = list.length === 0;
  el("warnings").innerHTML = list.map((w) => `<li>${escapeHtml(w)}</li>`).join("");
}

function renderTable() {
  const rows = state.spec.surfaces.map((s, i) => {
    const isImage = i === state.spec.surfaces.length - 1;
    const name = isImage ? "IMG" : String(i + 1);
    const options = state.glasses
      .map(
        (g) =>
          `<option value="${g}"${g.toUpperCase() === s.glass.toUpperCase() ? " selected" : ""}>${g}</option>`,
      )
      .join("");

    return `
      <tr${isImage ? ' class="image-row"' : ""}>
        <th scope="row">${name}</th>
        <td><input type="number" step="0.1" data-field="radius" data-row="${i}"
                   value="${s.radius ?? ""}" placeholder="&#8734;" aria-label="Radius of surface ${name}"></td>
        <td>${
          isImage
            ? "&mdash;"
            : `<input type="number" step="0.1" data-field="thickness" data-row="${i}" value="${round(s.thickness)}" aria-label="Thickness after surface ${name}">`
        }</td>
        <td>${isImage ? "&mdash;" : `<select data-field="glass" data-row="${i}" aria-label="Glass after surface ${name}">${options}</select>`}</td>
        <td><input type="number" step="0.1" min="0" data-field="semi_diameter" data-row="${i}"
                   value="${s.semi_diameter ?? ""}" placeholder="auto" aria-label="Semi-diameter of surface ${name}"></td>
        <td>${isImage ? "" : `<input type="radio" name="stop" data-row="${i}"${s.stop ? " checked" : ""} aria-label="Aperture stop at surface ${name}">`}</td>
      </tr>`;
  });

  el("prescription").innerHTML = rows.join("");
}

function onCellEdit(event) {
  const target = event.target;
  const row = Number(target.dataset.row);
  if (Number.isNaN(row)) return;

  if (target.type === "radio") {
    state.spec.surfaces.forEach((s, i) => (s.stop = i === row));
  } else {
    const field = target.dataset.field;
    const raw = target.value.trim();
    if (field === "glass") {
      state.spec.surfaces[row].glass = raw;
    } else if (raw === "") {
      // Empty means "flat" for a radius and "work it out from the rays" for an aperture.
      state.spec.surfaces[row][field] = field === "thickness" ? 0 : null;
    } else {
      state.spec.surfaces[row][field] = Number(raw);
    }
    if (field === "thickness" && row === state.spec.surfaces.length - 2) {
      recentreDefocus();
    }
  }
  schedule();
}

/** Push engine-side changes back into the inputs without rebuilding the table. */
function syncTableValues() {
  const input = document.querySelector(
    `#prescription input[data-field="thickness"][data-row="${state.spec.surfaces.length - 2}"]`,
  );
  if (input) input.value = round(lastSurface().thickness);
}

function parseNumbers(text) {
  return text
    .split(/[,;\s]+/)
    .map((t) => Number(t))
    .filter((n) => Number.isFinite(n));
}

const round = (v) => Math.round(v * 1e6) / 1e6;

function escapeHtml(text) {
  return String(text).replace(
    /[&<>"']/g,
    (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c],
  );
}

/**
 * Kicked off on import. Exported so a headless test can await the first analysis
 * instead of guessing how long it takes.
 */
export const ready = boot();
