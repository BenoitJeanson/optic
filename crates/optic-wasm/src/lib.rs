//! WebAssembly entry points for the browser demo.
//!
//! The real work lives in [`api`], which is ordinary Rust and tested natively. What is
//! left here is a hand-written C ABI: a few functions moving UTF-8 JSON across the
//! boundary in linear memory.
//!
//! This is deliberately not `wasm-bindgen`. The generated-bindings route needs a CLI
//! whose version must match the crate exactly, which is a standing source of breakage in
//! a deploy pipeline that nobody looks at for months. For an interface this small —
//! JSON in, JSON out — sixty lines of pointer handling is a better trade. Nothing here
//! is on a hot path; the ray tracing all happens on the Rust side of one call.
//!
//! ## Protocol
//!
//! Strings crossing the boundary are UTF-8. Values *returned* to the caller are prefixed
//! with their byte length as a little-endian `u32`, so the caller can find the end
//! without a second call. Every returned pointer must be handed back to
//! [`optic_free_result`].

pub mod api;

pub use api::{analyze, analyze_json, presets_json, Analysis, Request, SystemSpec};

/// Hand a buffer's memory to the caller, who becomes responsible for returning it.
fn leak(bytes: Vec<u8>) -> *mut u8 {
    // Boxing the slice trims the capacity to the length, so the pointer and length the
    // caller holds are enough to reconstruct the exact allocation later.
    Box::into_raw(bytes.into_boxed_slice()) as *mut u8
}

/// Length-prefixed UTF-8, ready to return across the boundary.
fn leak_prefixed(text: String) -> *mut u8 {
    let bytes = text.into_bytes();
    let mut out = Vec::with_capacity(4 + bytes.len());
    out.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
    out.extend_from_slice(&bytes);
    leak(out)
}

/// Reserve `len` bytes for the caller to write a request into.
///
/// Returns null if `len` is zero.
#[no_mangle]
pub extern "C" fn optic_alloc(len: usize) -> *mut u8 {
    if len == 0 {
        return core::ptr::null_mut();
    }
    leak(vec![0u8; len])
}

/// Release a buffer obtained from [`optic_alloc`].
///
/// # Safety
/// `ptr` must come from [`optic_alloc`] with exactly this `len`, and must not be used
/// afterwards.
#[no_mangle]
pub unsafe extern "C" fn optic_free(ptr: *mut u8, len: usize) {
    if ptr.is_null() || len == 0 {
        return;
    }
    drop(Box::from_raw(core::ptr::slice_from_raw_parts_mut(ptr, len)));
}

/// Release a length-prefixed result returned by [`optic_analyze`] or [`optic_presets`].
///
/// # Safety
/// `ptr` must be a pointer this module returned and not yet freed.
#[no_mangle]
pub unsafe extern "C" fn optic_free_result(ptr: *mut u8) {
    if ptr.is_null() {
        return;
    }
    let mut header = [0u8; 4];
    header.copy_from_slice(core::slice::from_raw_parts(ptr, 4));
    let len = u32::from_le_bytes(header) as usize;
    drop(Box::from_raw(core::ptr::slice_from_raw_parts_mut(
        ptr,
        4 + len,
    )));
}

/// Analyse the system described by `len` bytes of UTF-8 JSON at `ptr`.
///
/// Returns a length-prefixed JSON response. A malformed request produces a response with
/// `ok: false` rather than a failure, so the caller always has something to display.
///
/// # Safety
/// `ptr` must point to `len` readable bytes.
#[no_mangle]
pub unsafe extern "C" fn optic_analyze(ptr: *const u8, len: usize) -> *mut u8 {
    let request = if ptr.is_null() {
        String::new()
    } else {
        String::from_utf8_lossy(core::slice::from_raw_parts(ptr, len)).into_owned()
    };
    leak_prefixed(analyze_json(&request))
}

/// The built-in prescriptions and the glass list, as length-prefixed JSON.
#[no_mangle]
pub extern "C" fn optic_presets() -> *mut u8 {
    leak_prefixed(presets_json())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Drive the C ABI exactly as the browser does: allocate, write, call, read, free.
    #[test]
    fn the_abi_round_trips_a_request() {
        let presets: serde_json::Value = serde_json::from_str(&presets_json()).unwrap();
        let request = serde_json::json!({ "system": presets["presets"][0] }).to_string();

        unsafe {
            let buf = optic_alloc(request.len());
            assert!(!buf.is_null());
            core::ptr::copy_nonoverlapping(request.as_ptr(), buf, request.len());

            let result = optic_analyze(buf, request.len());
            optic_free(buf, request.len());

            let mut header = [0u8; 4];
            header.copy_from_slice(core::slice::from_raw_parts(result, 4));
            let len = u32::from_le_bytes(header) as usize;
            let body = core::slice::from_raw_parts(result.add(4), len);
            let parsed: serde_json::Value = serde_json::from_slice(body).unwrap();

            assert_eq!(parsed["ok"], true);
            assert!((parsed["first_order"]["efl"].as_f64().unwrap() - 50.02).abs() < 0.1);

            optic_free_result(result);
        }
    }

    #[test]
    fn a_zero_length_allocation_is_null() {
        assert!(optic_alloc(0).is_null());
        // Freeing null, or a zero-length buffer, must be harmless.
        unsafe {
            optic_free(core::ptr::null_mut(), 0);
            optic_free_result(core::ptr::null_mut());
        }
    }

    #[test]
    fn a_null_request_still_produces_a_response() {
        unsafe {
            let result = optic_analyze(core::ptr::null(), 0);
            let mut header = [0u8; 4];
            header.copy_from_slice(core::slice::from_raw_parts(result, 4));
            let len = u32::from_le_bytes(header) as usize;
            let body = core::slice::from_raw_parts(result.add(4), len);
            let parsed: serde_json::Value = serde_json::from_slice(body).unwrap();
            assert_eq!(parsed["ok"], false);
            assert!(parsed["error"].is_string());
            optic_free_result(result);
        }
    }

    #[test]
    fn presets_cross_the_boundary_intact() {
        unsafe {
            let result = optic_presets();
            let mut header = [0u8; 4];
            header.copy_from_slice(core::slice::from_raw_parts(result, 4));
            let len = u32::from_le_bytes(header) as usize;
            let body = core::slice::from_raw_parts(result.add(4), len);
            let parsed: serde_json::Value = serde_json::from_slice(body).unwrap();
            assert!(parsed["presets"].as_array().unwrap().len() >= 2);
            assert!(parsed["glasses"]
                .as_array()
                .unwrap()
                .contains(&"N-BK7".into()));
            optic_free_result(result);
        }
    }
}
