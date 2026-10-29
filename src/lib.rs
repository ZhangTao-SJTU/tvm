#![forbid(unsafe_code)]
//! Three-dimensional vertex model for periodic cellular tissues.
pub type Result<T> = std::result::Result<T, String>;
pub mod io;
pub mod mechanics;
pub mod model;
pub mod reconnection;
pub mod rng;
pub mod simulation;
