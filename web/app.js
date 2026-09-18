/**
 * The demo application.
 *
 * It owns a prescription, hands it to the engine on every edit, and draws whatever comes
 * back. There is no incremental update path and none is wanted: a full analysis of a
 * triplet takes a few milliseconds, so recomputing everything keeps what you see and
 * what the kernel believes in permanent agreement.
 *
 * Two things deliberately live outside the prescription, because they are questions
 * asked *of* a design rather than parts of it: the defocus shift and the pupil sampling
 * pattern. Neither edits the lens.
 */

import { OpticEngine } from "./engine.js";
import { drawLayout, drawSpot, fieldColour } from "./draw.js";

const state = {
  engine: null,
  presets: [],
  glasses: [],
  spec: null,
  analysis: null,
  /** Which transverse axis the field list refers to. */
  fieldAxis: "y",
  /** Image-plane shift in mm, applied after solves. Not part of the design. */
  defocus: 0,
  pupilPattern: "square",
  showAiry: true,
  pending: false,
  /** Warnings from the last .zmx import, shown alongside the analysis ones. */
  importWarnings: [],
};

const SOLVE_LABELS = {
  fixed: "—",
  marginal_ray_height: "Marginal ray",
  chief_ray_height: "Chief ray",
  pickup: "Pickup",
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
    setFields(parsed.length ? parsed : [0]);
    e.target.value = fieldMagnitudes().join(", ");
    schedule();
  });

  el("field-axis").addEventListener("change", (e) => {
    state.fieldAxis = e.target.value;
    setFields(fieldMagnitudes());
    schedule();
  });

  el("wavelengths").addEventListener("change", (e) => {
    const parsed = parseNumbers(e.target.value).filter((w) => w > 0.1 && w < 20);
    state.spec.wavelengths = parsed.length ? parsed : [0.5875618];
    if (state.spec.primary_wavelength >= state.spec.wavelengths.length) {
      state.spec.primary_wavelength = state.spec.wavelengths.length >> 1;
    }
    e.target.value = state.spec.wavelengths.map((w) => w.toFixed(4)).join(", ");
    renderPrimaryChoices();
    schedule();
  });

  el("primary").addEventListener("change", (e) => {
    state.spec.primary_wavelength = Number(e.target.value);
    schedule();
  });

  el("pupil-pattern").addEventListener("change", (e) => {
    state.pupilPattern = e.target.value;
    schedule();
  });

  el("show-airy").addEventListener("change", (e) => {
    state.showAiry = e.target.checked;
    render();
  });

  el("defocus").addEventListener("input", () => {
    state.defocus = Number(el("defocus").value);
    el("defocus-value").textContent = `${state.defocus >= 0 ? "+" : ""}${state.defocus.toFixed(2)} mm`;
    schedule();
  });

  el("refocus").addEventListener("click", () => {
    state.defocus = 0;
    el("defocus").value = 0;
    el("defocus-value").textContent = "+0.00 mm";
    // Only meaningful when nothing is solving the image distance for us.
    const solved = lastSurface().solve?.type !== "fixed";
    update({ refocus: !solved });
  });

  el("zmx-file").addEventListener("change", async (e) => {
    const file = e.target.files?.[0];
    if (file) await openZmx(file);
    e.target.value = "";
  });

  el("zmx-export").addEventListener("click", saveZmx);

  // Dropping a file anywhere on the page is the obvious gesture, so support it.
  for (const type of ["dragenter", "dragover"]) {
    window.addEventListener(type, (e) => {
      e.preventDefault();
      document.body.classList.add("dragging");
    });
  }
  for (const type of ["dragleave", "drop"]) {
    window.addEventListener(type, (e) => {
      e.preventDefault();
      if (type === "dragleave" && e.relatedTarget) return;
      document.body.classList.remove("dragging");
    });
  }
  window.addEventListener("drop", async (e) => {
    const file = e.dataTransfer?.files?.[0];
    if (file) await openZmx(file);
  });

  // Bound once: renderTable replaces the rows, not the table body, so attaching there
  // would stack a new listener on every redraw.
  el("prescription").addEventListener("input", onCellEdit);
  el("prescription").addEventListener("change", onCellEdit);

  window.addEventListener("resize", () => render());

  loadPreset(0);
}

/**
 * Load a Zemax prescription.
 *
 * Whatever the importer could not honour is kept and shown with the analysis warnings:
 * a design that silently loses a coordinate break or a vignetting factor is worse than
 * one that arrives with a list of caveats.
 */
async function openZmx(file) {
  let result;
  try {
    result = state.engine.importZmx(await file.arrayBuffer());
  } catch (error) {
    showError(`Could not read ${file.name}: ${error.message}`);
    return;
  }
  if (!result.ok) {
    showError(`${file.name}: ${result.error}`);
    return;
  }

  state.spec = result.system;
  state.importWarnings = result.warnings ?? [];
  state.defocus = 0;
  state.fieldAxis = state.spec.fields.some((f) => Math.abs(f.x) > Math.abs(f.y)) ? "x" : "y";

  el("preset").selectedIndex = -1;
  el("epd").value = state.spec.entrance_pupil_diameter;
  el("fields").value = fieldMagnitudes().join(", ");
  el("field-axis").value = state.fieldAxis;
  el("wavelengths").value = state.spec.wavelengths.map((w) => w.toFixed(4)).join(", ");
  el("defocus").value = 0;
  el("defocus-value").textContent = "+0.00 mm";

  renderPrimaryChoices();
  renderTable();
  update();
}

/** Hand the current prescription back as a .zmx file. */
function saveZmx() {
  const result = state.engine.exportZmx(state.spec);
  if (!result.ok) {
    showError(result.error);
    return;
  }
  const name = (state.spec.title || "optic").replace(/[^\w.-]+/g, "_");
  const url = URL.createObjectURL(new Blob([result.text], { type: "text/plain" }));
  const link = document.createElement("a");
  link.href = url;
  link.download = `${name}.zmx`;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function showError(message) {
  el("error").textContent = message;
  el("error").hidden = false;
}

function lastSurface() {
  return state.spec.surfaces[state.spec.surfaces.length - 2];
}

/** Field magnitudes along the currently selected axis. */
function fieldMagnitudes() {
  return state.spec.fields.map((f) => (Math.abs(f.x) > Math.abs(f.y) ? f.x : f.y));
}

function setFields(values) {
  state.spec.fields = values.map((v) =>
    state.fieldAxis === "x" ? { x: v, y: 0 } : { x: 0, y: v },
  );
}

function loadPreset(index) {
  state.spec = structuredClone(state.presets[index]);
  state.importWarnings = [];
  state.defocus = 0;
  state.fieldAxis = state.spec.fields.some((f) => Math.abs(f.x) > Math.abs(f.y)) ? "x" : "y";

  el("epd").value = state.spec.entrance_pupil_diameter;
  el("fields").value = fieldMagnitudes().join(", ");
  el("field-axis").value = state.fieldAxis;
  el("wavelengths").value = state.spec.wavelengths.map((w) => w.toFixed(4)).join(", ");
  el("defocus").value = 0;
  el("defocus-value").textContent = "+0.00 mm";

  renderPrimaryChoices();
  renderTable();
  update();
}

function renderPrimaryChoices() {
  const primary = state.spec.primary_wavelength ?? state.spec.wavelengths.length >> 1;
  state.spec.primary_wavelength = Math.min(primary, state.spec.wavelengths.length - 1);
  el("primary").innerHTML = state.spec.wavelengths
    .map(
      (w, i) =>
        `<option value="${i}"${i === state.spec.primary_wavelength ? " selected" : ""}>${w.toFixed(4)} &micro;m</option>`,
    )
    .join("");
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

function update({ refocus = false } = {}) {
  const result = state.engine.analyze({
    system: state.spec,
    rays_per_fan: 15,
    spot_grid: 23,
    refocus,
    pupil_pattern: state.pupilPattern,
    defocus: state.defocus,
  });

  if (!result.ok) {
    el("error").textContent = result.error;
    el("error").hidden = false;
    return;
  }
  el("error").hidden = true;

  // The engine echoes the prescription it actually used, so solved thicknesses and a
  // refocus show up in the editor rather than being applied invisibly.
  state.spec = result.system;
  state.analysis = result;

  syncSolvedThicknesses();
  render();
  renderReadout();
  renderDistortion();
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
            <strong class="field-label"></strong>
            <span class="metrics"></span>
          </figcaption>
        </figure>`,
      )
      .join("");
  }

  spots.forEach((spot, i) => {
    drawSpot(el(`spot-${i}`), spot, halfWidth, state.spec.wavelengths, state.showAiry);
    const figure = container.children[i];
    const axis = Math.abs(spot.field.x) > Math.abs(spot.field.y) ? "X" : "Y";
    const magnitude = Math.abs(spot.field.x) > Math.abs(spot.field.y) ? spot.field.x : spot.field.y;
    figure.querySelector(".field-label").textContent = `${axis} ${magnitude.toFixed(1)}°`;

    const caption = figure.querySelector(".metrics");
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
    ["Primary", `${state.analysis.primary_wavelength.toFixed(4)} µm`],
  ];
  el("first-order").innerHTML = rows
    .map(([k, v]) => `<div><dt>${k}</dt><dd>${v}</dd></div>`)
    .join("");
}

/**
 * Distortion per field and wavelength.
 *
 * It earns a table of its own because it is the one number that depends on nothing but
 * the chief ray — no aperture, no vignetting, no sampling — which makes it the sharpest
 * instrument available for comparing this kernel against another tool.
 */
function renderDistortion() {
  const rows = state.analysis.distortion;
  const table = el("distortion");
  if (!rows.length) {
    table.innerHTML = "";
    table.hidden = true;
    return;
  }
  table.hidden = false;

  const wavelengths = state.spec.wavelengths;
  const fields = [...new Set(rows.map((r) => r.field_index))].sort((a, b) => a - b);
  const find = (fi, wi) =>
    rows.find((r) => r.field_index === fi && Math.abs(r.wavelength - wavelengths[wi]) < 1e-12);

  const head = `<tr><th>Distortion</th>${wavelengths
    .map((w) => `<th>${w.toFixed(4)} &micro;m</th>`)
    .join("")}</tr>`;
  const body = fields
    .map((fi) => {
      const spot = state.analysis.spots[fi];
      const axis = Math.abs(spot.field.x) > Math.abs(spot.field.y) ? "X" : "Y";
      const magnitude =
        Math.abs(spot.field.x) > Math.abs(spot.field.y) ? spot.field.x : spot.field.y;
      const cells = wavelengths
        .map((_, wi) => {
          const cell = find(fi, wi);
          return `<td>${cell ? `${cell.percent >= 0 ? "+" : ""}${cell.percent.toFixed(4)}%` : "&mdash;"}</td>`;
        })
        .join("");
      return `<tr><th scope="row">${axis} ${magnitude.toFixed(1)}&deg;</th>${cells}</tr>`;
    })
    .join("");
  table.innerHTML = `<thead>${head}</thead><tbody>${body}</tbody>`;
}

function renderWarnings() {
  const list = [...state.importWarnings, ...state.analysis.warnings];
  el("warnings").hidden = list.length === 0;
  el("warnings").innerHTML = list.map((w) => `<li>${escapeHtml(w)}</li>`).join("");
}

function renderTable() {
  const rows = state.spec.surfaces.map((s, i) => {
    const isImage = i === state.spec.surfaces.length - 1;
    const name = isImage ? "IMG" : String(i + 1);
    const solveType = s.solve?.type ?? "fixed";
    const solved = solveType !== "fixed";

    const glassOptions = state.glasses
      .map(
        (g) =>
          `<option value="${g}"${g.toUpperCase() === s.glass.toUpperCase() ? " selected" : ""}>${g}</option>`,
      )
      .join("");

    const solveOptions = ["fixed", "marginal_ray_height", "chief_ray_height"]
      .map(
        (v) =>
          `<option value="${v}"${v === solveType ? " selected" : ""}>${SOLVE_LABELS[v]}</option>`,
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
            : `<input type="number" step="0.1" data-field="thickness" data-row="${i}"
                      value="${round(s.thickness)}"${solved ? " readonly class=\"solved\"" : ""}
                      aria-label="Thickness after surface ${name}">`
        }</td>
        <td>${isImage ? "&mdash;" : `<select data-field="glass" data-row="${i}" aria-label="Glass after surface ${name}">${glassOptions}</select>`}</td>
        <td>${isImage ? "&mdash;" : `<select data-field="solve" data-row="${i}" aria-label="Thickness solve for surface ${name}">${solveOptions}</select>`}</td>
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
    schedule();
    return;
  }

  const field = target.dataset.field;
  const raw = target.value.trim();

  if (field === "solve") {
    state.spec.surfaces[row].solve =
      raw === "fixed" ? { type: "fixed" } : { type: raw, height: 0 };
    renderTable();
    schedule();
    return;
  }

  if (field === "glass") {
    state.spec.surfaces[row].glass = raw;
  } else if (raw === "") {
    // Empty means "flat" for a radius and "work it out from the rays" for an aperture.
    state.spec.surfaces[row][field] = field === "thickness" ? 0 : null;
  } else {
    state.spec.surfaces[row][field] = Number(raw);
  }
  schedule();
}

/** Write solved thicknesses back into their inputs without rebuilding the table. */
function syncSolvedThicknesses() {
  state.spec.surfaces.forEach((s, i) => {
    if ((s.solve?.type ?? "fixed") === "fixed") return;
    const input = document.querySelector(
      `#prescription input[data-field="thickness"][data-row="${i}"]`,
    );
    if (input && document.activeElement !== input) input.value = round(s.thickness);
  });
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
