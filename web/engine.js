/**
 * The seam between the interface and whatever is doing the physics.
 *
 * Today there is one implementation, backed by `optic-core` compiled to WebAssembly.
 * The desktop shell will add a second one backed by native Rust over Tauri's IPC. Both
 * wrap the same crate, so they cannot disagree about optics; keeping the UI on this
 * side of the seam is what lets one interface serve both.
 *
 * Strings cross the boundary as UTF-8 JSON, length-prefixed on the way back. Every
 * pointer handed out is returned in the same call that used it.
 */

const ENCODER = new TextEncoder();
const DECODER = new TextDecoder();

export class OpticEngine {
  #exports;

  constructor(instance) {
    this.#exports = instance.exports;
  }

  static async fromBytes(bytes) {
    const { instance } = await WebAssembly.instantiate(bytes, {});
    return new OpticEngine(instance);
  }

  static async fromUrl(url) {
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`could not load ${url}: ${response.status} ${response.statusText}`);
    }
    return OpticEngine.fromBytes(await response.arrayBuffer());
  }

  /** The built-in prescriptions and the glass list. */
  presets() {
    return this.#readResult(this.#exports.optic_presets());
  }

  /**
   * Trace a system and return everything needed to draw it.
   * A malformed request comes back as `{ ok: false, error }` rather than throwing.
   */
  analyze(request) {
    const bytes = ENCODER.encode(JSON.stringify(request));
    const ptr = this.#exports.optic_alloc(bytes.length);
    // The view must be taken after the allocation: growing memory detaches any buffer
    // captured before it.
    new Uint8Array(this.#exports.memory.buffer, ptr, bytes.length).set(bytes);
    const result = this.#exports.optic_analyze(ptr, bytes.length);
    this.#exports.optic_free(ptr, bytes.length);
    return this.#readResult(result);
  }

  #readResult(ptr) {
    const memory = this.#exports.memory.buffer;
    const length = new DataView(memory).getUint32(ptr, true);
    const text = DECODER.decode(new Uint8Array(memory, ptr + 4, length));
    this.#exports.optic_free_result(ptr);
    return JSON.parse(text);
  }
}
