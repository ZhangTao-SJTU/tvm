// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2021-2026 Tao Zhang; coauthor Jennifer Schwarz.
//! Euler--Maruyama integration and checked, non-overwriting simulation outputs.
use crate::io::{load_config, load_topology, topology_text, vtk_text, Config};
use crate::mechanics;
use crate::model::{add, scale, sub, Mesh};
use crate::reconnection::Reconnection;
use crate::rng::Rng;
use crate::Result;
use std::fmt::Write as _;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

struct Output {
    path: PathBuf,
    writer: BufWriter<File>,
}

impl Output {
    fn create(path: PathBuf) -> Result<Self> {
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(|e| {
                format!(
                    "Cannot create output file {} (existing files are never overwritten): {e}",
                    path.display()
                )
            })?;
        Ok(Self {
            path,
            writer: BufWriter::new(file),
        })
    }
    fn append(&mut self, text: &str) -> Result<()> {
        self.writer
            .write_all(text.as_bytes())
            .and_then(|_| self.writer.flush())
            .map_err(|e| {
                format!(
                    "Cannot write/flush output file {}: {e}",
                    self.path.display()
                )
            })
    }
}

fn new_text(directory: &Path, name: &str, text: &str) -> Result<()> {
    Output::create(directory.join(name))?.append(text)
}

/// Overdamped forces with center-of-mass drift removed.
pub fn update_velocities(mesh: &mut Mesh) -> Result<()> {
    if mesh.vertices.is_empty() {
        return Err("Cannot evolve an empty network".into());
    }
    let mut average = [0.; 3];
    for vertex in mesh.vertices.values_mut() {
        vertex.velocity = add(vertex.volume_force, vertex.interface_force);
        average = add(average, vertex.velocity);
    }
    average = scale(average, 1. / mesh.vertices.len() as f64);
    for vertex in mesh.vertices.values_mut() {
        vertex.velocity = sub(vertex.velocity, average);
        if !vertex.velocity.iter().all(|x| x.is_finite()) {
            return Err("Nonfinite vertex velocity".into());
        }
    }
    Ok(())
}

/// Advance positions; zero temperature deliberately consumes no random samples.
pub fn advance_positions(mesh: &mut Mesh, dt: f64, temperature: f64, rng: &mut Rng) -> Result<()> {
    if !(dt > 0. && dt.is_finite() && temperature >= 0. && temperature.is_finite()) {
        return Err("Invalid integration step or temperature".into());
    }
    let noise = (2. * temperature * dt).sqrt();
    if !noise.is_finite() {
        return Err("Nonfinite thermal displacement scale".into());
    }
    let size = mesh.box_size;
    for vertex in mesh.vertices.values_mut() {
        for (k, length) in size.iter().copied().enumerate() {
            let mut x = vertex.position[k] + vertex.velocity[k] * dt;
            if noise != 0. {
                x += noise * rng.normal();
            }
            if !x.is_finite() {
                return Err("Nonfinite vertex position".into());
            }
            x %= length;
            if x < 0. {
                x += length;
            }
            vertex.position[k] = x;
        }
    }
    Ok(())
}

fn metadata(conf: &Config, mesh: &Mesh) -> String {
    format!(concat!("model periodic-bulk-target-area\nimplementation tvm {}\n",
        "repository https://github.com/ZhangTao-SJTU/tvm\n",
        "integrator Euler-Maruyama (Euler when T=0)\nrng mt19937_64 + Marsaglia polar normal\n",
        "seed {}\ntime {:.17e} {:.17e} {:.17e}\ndtr {:.17e}\nT {:.17e}\nkv {:.17e}\ns0 {:.17e}\nLth {:.17e}\n",
        "box {:.17e} {:.17e} {:.17e}\nvalidate_topology {}\ninitial_counts {} {} {} {}\n",
        "Reproducibility requires the same executable, configuration and topology.\n",
        "rng_state.txt stores the final engine and normal-distribution state; it is not a full restart checkpoint.\n"),
        env!("CARGO_PKG_VERSION"),conf.seed,conf.start,conf.end,conf.dt,conf.dtr,conf.temperature,conf.kv,conf.s0,conf.lth,
        conf.box_size[0],conf.box_size[1],conf.box_size[2],if conf.validate_topology {"on"} else {"off"},
        mesh.vertices.len(),mesh.edges.len(),mesh.polygons.len(),mesh.cells.len())
}

struct Dumps {
    topology: Output,
    centers: Output,
    shapes: Output,
    volumes: Output,
    reconnections: Output,
}

impl Dumps {
    fn create(directory: &Path) -> Result<Self> {
        Ok(Self {
            topology: Output::create(directory.join("topo.txt"))?,
            centers: Output::create(directory.join("cellCenter.txt"))?,
            shapes: Output::create(directory.join("cellShapeIndex.txt"))?,
            volumes: Output::create(directory.join("cellVolume.txt"))?,
            reconnections: Output::create(directory.join("reconnections.txt"))?,
        })
    }
    fn events(&mut self, reconnection: &mut Reconnection) -> Result<()> {
        if !reconnection.events.is_empty() {
            let mut text = String::new();
            for event in &reconnection.events {
                text.push_str(event);
                if !event.ends_with('\n') {
                    text.push('\n');
                }
            }
            self.reconnections.append(&text)?;
            reconnection.events.clear();
        }
        Ok(())
    }
    fn frame(&mut self, mesh: &Mesh, time: f64, reconnection: &mut Reconnection) -> Result<()> {
        self.topology
            .append(&format!("time {time:.17e}\n{}\n", topology_text(mesh)))?;
        let mut centers = format!("time {time:.17e}\n");
        let mut volumes = format!(
            "time {time:.17e} {:.17e}\n",
            mesh.cells.values().map(|c| c.volume).sum::<f64>()
        );
        let mut shapes = String::new();
        let mut mean_shape = 0.;
        for cell in mesh.cells.values() {
            let reference = mesh.vertices[&cell.vertices[0]].position;
            let mut center = [0.; 3];
            for v in &cell.vertices {
                center = add(
                    center,
                    mesh.distance(sub(mesh.vertices[v].position, reference)),
                );
            }
            center = mesh.wrap(add(
                scale(center, 1. / cell.vertices.len() as f64),
                reference,
            ));
            writeln!(
                centers,
                "{} {:.17e} {:.17e} {:.17e}",
                cell.id, center[0], center[1], center[2]
            )
            .unwrap();
            writeln!(volumes, "{} {:.17e}", cell.id, cell.volume).unwrap();
            let area: f64 = cell.polygons.iter().map(|p| mesh.polygons[p].area).sum();
            let shape = area * cell.volume.powf(-2. / 3.);
            mean_shape += shape;
            writeln!(shapes, "{} {:.17e}", cell.id, shape).unwrap();
        }
        mean_shape /= mesh.cells.len() as f64;
        centers.push('\n');
        volumes.push('\n');
        self.centers.append(&centers)?;
        self.volumes.append(&volumes)?;
        self.shapes
            .append(&format!("time {time:.17e} {mean_shape:.17e}\n{shapes}\n"))?;
        self.events(reconnection)
    }
}

/// Run using `conf` and `sample.topo` in a directory. Never changes the process cwd.
pub fn run(directory: &Path) -> Result<()> {
    let conf = load_config(&directory.join("conf"))?;
    let mut mesh = load_topology(&directory.join("sample.topo"), conf.box_size)?;
    println!("Initialized {} vertices, {} edges, {} faces, {} cells\nseed {}, dt {}, dtr {}, topology checks {}",
        mesh.vertices.len(),mesh.edges.len(),mesh.polygons.len(),mesh.cells.len(),conf.seed,conf.dt,conf.dtr,
        if conf.validate_topology {"on"} else {"off"});
    for name in [
        "run_metadata.txt",
        "topo.txt",
        "cellCenter.txt",
        "cellShapeIndex.txt",
        "cellVolume.txt",
        "reconnections.txt",
        "rng_state.txt",
        "run_summary.txt",
    ] {
        match fs::symlink_metadata(directory.join(name)) {
            Ok(_) => {
                return Err(
                    "This directory already contains a run; use a fresh output directory".into(),
                )
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => (),
            Err(e) => return Err(format!("Cannot inspect output path {name}: {e}")),
        }
    }
    new_text(directory, "run_metadata.txt", &metadata(&conf, &mesh))?;
    let mut dumps = Dumps::create(directory)?;
    let mut rng = Rng::new(conf.seed);
    let mut reconnect = Reconnection::new(conf.lth, conf.verbose);
    reconnect.diagnostic_directory = Some(directory.to_path_buf());
    let mut time = conf.start;
    let tolerance =
        (0.01 * conf.dt).min(32. * f64::EPSILON * 1_f64.max(conf.start.abs()).max(conf.end.abs()));
    let (mut step, mut count_dump, mut count_log, mut count_reconnect) =
        (0_u64, 0_u64, 0_u64, 0_u64);
    let timer = Instant::now();
    println!("Time Elapsed Volume I_to_H H_to_I E_volume E_interface Energy");
    let final_energy = loop {
        let energy = mechanics::forces(&mut mesh, conf.kv, conf.s0)?;
        mesh.validate_geometry()?;
        update_velocities(&mut mesh)?;
        let final_frame = time >= conf.end;
        let elapsed = time - conf.start + tolerance;
        if elapsed >= count_log as f64 * conf.log_period || final_frame {
            println!(
                "{:.17e} {:.9} {:.17e} {} {} {:.17e} {:.17e} {:.17e}",
                time,
                timer.elapsed().as_secs_f64(),
                energy.total_volume,
                reconnect.count_ih,
                reconnect.count_hi,
                energy.volume,
                energy.interface,
                energy.total()
            );
            reconnect.count_ih = 0;
            reconnect.count_hi = 0;
            count_log = (elapsed / conf.log_period).floor() as u64 + 1;
        }
        let dump_due = elapsed >= count_dump as f64 * conf.dump_period || final_frame;
        if dump_due && (time >= -tolerance || final_frame) {
            dumps.frame(&mesh, time, &mut reconnect)?;
            if conf.dump_vtk {
                new_text(
                    directory,
                    &format!("{count_dump:07}.sample.vtk"),
                    &vtk_text(&mesh),
                )?;
            }
        }
        if dump_due {
            count_dump = (elapsed / conf.dump_period).floor() as u64 + 1;
        }
        if final_frame {
            break energy;
        }
        step += 1;
        let next_time = conf.end.min(conf.start + step as f64 * conf.dt);
        let dt = next_time - time;
        if dt <= 0. {
            return Err("Time step cannot advance at this precision".into());
        }
        advance_positions(&mut mesh, dt, conf.temperature, &mut rng)?;
        time = next_time;
        if time - conf.start + tolerance >= count_reconnect as f64 * conf.dtr {
            reconnect.time = time;
            reconnect.start(&mut mesh)?;
            if conf.validate_topology {
                mesh.validate_topology()?;
            }
            count_reconnect = ((time - conf.start + tolerance) / conf.dtr).floor() as u64 + 1;
        }
    };
    dumps.events(&mut reconnect)?;
    new_text(directory, "rng_state.txt", &format!("{}\n", rng.state()))?;
    new_text(
        directory,
        "run_summary.txt",
        &format!(
            concat!(
                "completed 1\nfinal_time {:.17e}\nenergy {:.17e}\nvolume {:.17e}\n",
                "rejected_geometry {}\nrejected_topology {}\n"
            ),
            time,
            final_energy.total(),
            final_energy.total_volume,
            reconnect.rejected_geometry,
            reconnect.rejected_topology
        ),
    )?;
    Ok(())
}
