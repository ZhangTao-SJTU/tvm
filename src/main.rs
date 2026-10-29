// SPDX-License-Identifier: GPL-3.0-or-later
use std::path::PathBuf;

fn main() {
    let arguments: Vec<_> = std::env::args_os().skip(1).collect();
    if arguments.len() == 1 && (arguments[0] == "--help" || arguments[0] == "-h") {
        println!("Usage: tvm [--dir DIRECTORY]\nReads conf and sample.topo; writes outputs into the same directory. Existing output files are never overwritten.");
        return;
    }
    let directory = match arguments.as_slice() {
        [] => PathBuf::from("."),
        [flag, path] if flag == "--dir" => PathBuf::from(path),
        _ => {
            eprintln!("TVM error: Usage: tvm [--dir DIRECTORY]");
            std::process::exit(1);
        }
    };
    if let Err(error) = vm3d::simulation::run(&directory) {
        eprintln!("TVM error: {error}");
        std::process::exit(1);
    }
}
