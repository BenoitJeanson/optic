//! Reading and writing optical prescriptions.
//!
//! Interchange is not a side feature. A design tool that cannot read the files its users
//! already have asks them to start from nothing, and a tool whose results cannot be
//! checked in another program asks them to take its word. Both are fatal to adoption.

#![deny(rust_2018_idioms)]
#![warn(missing_debug_implementations)]

pub mod zmx;

pub use zmx::{decode, Import};
