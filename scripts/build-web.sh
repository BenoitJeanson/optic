#!/usr/bin/env bash
# Build the browser demo: compile the kernel to WebAssembly and stage it next to the page.
set -euo pipefail

cd "$(dirname "$0")/.."

# Non-login shells often miss the Rust toolchain; pick it up where rustup puts it.
if ! command -v cargo >/dev/null 2>&1 && [ -f "$HOME/.cargo/env" ]; then
  # shellcheck disable=SC1091
  . "$HOME/.cargo/env"
fi

rustup target add wasm32-unknown-unknown >/dev/null 2>&1 || true
cargo build -p optic-wasm --target wasm32-unknown-unknown --release
cp target/wasm32-unknown-unknown/release/optic_wasm.wasm web/

printf '\nBuilt web/optic_wasm.wasm (%s bytes)\n' "$(wc -c < web/optic_wasm.wasm | tr -d ' ')"
printf 'Serve it with:  python3 -m http.server -d web 8080\n'
