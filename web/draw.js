/**
 * Canvas rendering for the layout and spot diagrams.
 *
 * Everything here reads an analysis and writes pixels. It holds no state and knows
 * nothing about the engine, so the same drawing code will serve the desktop shell.
 */

/** Per-field ray colours, chosen to stay distinguishable on both themes. */
const FIELD_COLOURS = ["#2563eb", "#0d9488", "#d97706", "#9333ea", "#dc2626", "#65a30d"];

export function fieldColour(index) {
  return FIELD_COLOURS[index % FIELD_COLOURS.length];
}

/**
 * Approximate sRGB for a visible wavelength, so a spot diagram is coloured the way the
 * light actually is. Dan Bruton's piecewise fit, with the usual falloff at the ends.
 */
export function wavelengthColour(um) {
  const nm = um * 1000;
  let r = 0;
  let g = 0;
  let b = 0;
  if (nm < 440) {
    r = -(nm - 440) / 60;
    b = 1;
  } else if (nm < 490) {
    g = (nm - 440) / 50;
    b = 1;
  } else if (nm < 510) {
    g = 1;
    b = -(nm - 510) / 20;
  } else if (nm < 580) {
    r = (nm - 510) / 70;
    g = 1;
  } else if (nm < 645) {
    r = 1;
    g = -(nm - 645) / 65;
  } else {
    r = 1;
  }
  let falloff = 1;
  if (nm < 420) falloff = 0.3 + (0.7 * (nm - 380)) / 40;
  else if (nm > 700) falloff = 0.3 + (0.7 * (780 - nm)) / 80;
  const channel = (x) => Math.round(255 * Math.pow(Math.max(0, Math.min(1, x)) * falloff, 0.8));
  return `rgb(${channel(r)}, ${channel(g)}, ${channel(b)})`;
}

/** Match the backing store to the element's real size, and return a CSS-pixel context. */
function prepare(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  const width = Math.max(1, Math.round(rect.width));
  const height = Math.max(1, Math.round(rect.height));
  if (canvas.width !== width * dpr || canvas.height !== height * dpr) {
    canvas.width = width * dpr;
    canvas.height = height * dpr;
  }
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, width, height);
  return { ctx, width, height };
}

function themeColours() {
  const style = getComputedStyle(document.documentElement);
  return {
    ink: style.getPropertyValue("--ink").trim() || "#1f2933",
    faint: style.getPropertyValue("--faint").trim() || "#c9d2dc",
    glass: style.getPropertyValue("--glass").trim() || "rgba(37, 99, 235, 0.13)",
    glassEdge: style.getPropertyValue("--glass-edge").trim() || "#5b7fb5",
  };
}

/**
 * Draw the optical layout: glass elements in cross-section, traced rays, the stop and
 * the image plane. The aspect ratio is preserved, because a lens drawing that is not to
 * scale is worse than no drawing.
 */
export function drawLayout(canvas, analysis) {
  const { ctx, width, height } = prepare(canvas);
  const theme = themeColours();
  const [z0, y0, z1, y1] = analysis.layout.bounds;
  const pad = 14;

  const spanZ = Math.max(z1 - z0, 1e-6);
  const spanY = Math.max(y1 - y0, 1e-6);
  const scale = Math.min((width - 2 * pad) / spanZ, (height - 2 * pad) / spanY);
  const offsetX = (width - spanZ * scale) / 2 - z0 * scale;
  const offsetY = height / 2;
  const X = (z) => offsetX + z * scale;
  const Y = (y) => offsetY - y * scale;

  const path = (points) => {
    ctx.beginPath();
    points.forEach(([z, y], i) => (i ? ctx.lineTo(X(z), Y(y)) : ctx.moveTo(X(z), Y(y))));
  };

  // Optical axis.
  ctx.save();
  ctx.strokeStyle = theme.faint;
  ctx.setLineDash([5, 5]);
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(X(z0), Y(0));
  ctx.lineTo(X(z1), Y(0));
  ctx.stroke();
  ctx.restore();

  // Glass, filled first so rays draw over it.
  ctx.lineWidth = 1.2;
  for (const element of analysis.layout.elements) {
    path(element.points);
    ctx.closePath();
    ctx.fillStyle = theme.glass;
    ctx.fill();
  }

  // Rays.
  ctx.lineWidth = 0.9;
  for (const ray of analysis.layout.rays) {
    const index = ray.field_index;
    ctx.strokeStyle = fieldColour(index);
    ctx.globalAlpha = ray.complete ? 0.75 : 0.35;
    path(ray.points);
    ctx.stroke();
    if (!ray.complete && ray.points.length) {
      // Mark where a vignetted ray was stopped.
      const [z, y] = ray.points[ray.points.length - 1];
      ctx.beginPath();
      ctx.arc(X(z), Y(y), 2, 0, Math.PI * 2);
      ctx.fillStyle = fieldColour(index);
      ctx.fill();
    }
  }
  ctx.globalAlpha = 1;

  // Glass outlines over the rays.
  ctx.strokeStyle = theme.glassEdge;
  ctx.lineWidth = 1.3;
  for (const element of analysis.layout.elements) {
    path(element.points);
    ctx.closePath();
    ctx.stroke();
  }

  // The stop, drawn as blades, and the image plane.
  const specs = analysis.system.surfaces;
  ctx.strokeStyle = theme.ink;
  ctx.lineWidth = 2.5;
  analysis.layout.surfaces.forEach((surface, i) => {
    if (!specs[i] || !surface.points.length) return;
    const top = surface.points[surface.points.length - 1];
    const bottom = surface.points[0];
    if (specs[i].stop) {
      const blade = Math.max(spanY * 0.06, Math.abs(top[1]) * 0.45);
      ctx.beginPath();
      ctx.moveTo(X(top[0]), Y(top[1]));
      ctx.lineTo(X(top[0]), Y(top[1] + blade));
      ctx.moveTo(X(bottom[0]), Y(bottom[1]));
      ctx.lineTo(X(bottom[0]), Y(bottom[1] - blade));
      ctx.stroke();
    }
  });

  ctx.save();
  ctx.strokeStyle = theme.ink;
  ctx.globalAlpha = 0.55;
  ctx.lineWidth = 1.6;
  ctx.beginPath();
  ctx.moveTo(X(analysis.layout.image_z), Y(y0 * 0.85));
  ctx.lineTo(X(analysis.layout.image_z), Y(y1 * 0.85));
  ctx.stroke();
  ctx.restore();

  // Scale bar, so the drawing can be read as a measurement.
  const barMm = niceStep(spanZ / 5);
  ctx.save();
  ctx.strokeStyle = theme.ink;
  ctx.fillStyle = theme.ink;
  ctx.globalAlpha = 0.7;
  ctx.lineWidth = 1.5;
  const barY = height - 10;
  const barX = 14;
  ctx.beginPath();
  ctx.moveTo(barX, barY);
  ctx.lineTo(barX + barMm * scale, barY);
  ctx.moveTo(barX, barY - 3);
  ctx.lineTo(barX, barY + 3);
  ctx.moveTo(barX + barMm * scale, barY - 3);
  ctx.lineTo(barX + barMm * scale, barY + 3);
  ctx.stroke();
  ctx.font = "11px ui-sans-serif, system-ui, sans-serif";
  ctx.fillText(`${barMm} mm`, barX + barMm * scale + 6, barY + 4);
  ctx.restore();
}

function niceStep(raw) {
  const exponent = Math.floor(Math.log10(Math.max(raw, 1e-9)));
  const base = Math.pow(10, exponent);
  for (const m of [1, 2, 5]) {
    if (raw <= m * base) return m * base;
  }
  return 10 * base;
}

/**
 * Draw one spot diagram. `halfWidthUm` is shared across fields so the panels can be
 * compared directly, which is the whole point of showing them side by side.
 */
export function drawSpot(canvas, spot, halfWidthUm, wavelengths, showAiry = true) {
  const { ctx, width, height } = prepare(canvas);
  const theme = themeColours();
  const size = Math.min(width, height);
  const cx = width / 2;
  const cy = size / 2 + 2;
  const scale = (size / 2 - 16) / halfWidthUm;

  ctx.strokeStyle = theme.faint;
  ctx.lineWidth = 1;
  ctx.strokeRect(cx - size / 2 + 8, cy - size / 2 + 8, size - 16, size - 16);

  ctx.save();
  ctx.setLineDash([3, 3]);
  ctx.beginPath();
  ctx.moveTo(cx - size / 2 + 8, cy);
  ctx.lineTo(cx + size / 2 - 8, cy);
  ctx.moveTo(cx, cy - size / 2 + 8);
  ctx.lineTo(cx, cy + size / 2 - 8);
  ctx.stroke();
  ctx.restore();

  if (!spot.points.length) {
    ctx.fillStyle = theme.ink;
    ctx.font = "12px ui-sans-serif, system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("no rays reach the image", cx, cy);
    return;
  }

  // The Airy disc, for scale: a spot inside this circle is diffraction limited.
  if (showAiry && spot.airy * scale > 2) {
    ctx.save();
    ctx.strokeStyle = theme.ink;
    ctx.globalAlpha = 0.45;
    ctx.setLineDash([2, 3]);
    ctx.beginPath();
    ctx.arc(cx, cy, spot.airy * scale, 0, Math.PI * 2);
    ctx.stroke();
    ctx.restore();
  }

  for (const p of spot.points) {
    ctx.fillStyle = wavelengthColour(wavelengths[p.w] ?? 0.55);
    ctx.globalAlpha = 0.72;
    ctx.beginPath();
    ctx.arc(cx + p.x * scale, cy - p.y * scale, 1.15, 0, Math.PI * 2);
    ctx.fill();
  }
  ctx.globalAlpha = 1;
}
