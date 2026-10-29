// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2021-2026 Tao Zhang; coauthor Jennifer Schwarz.
//! Strict configuration/topology readers and portable text serializers.
use crate::model::{Cell, Edge, Id, Mesh, Polygon, Vec3, Vertex};
use crate::Result;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct Config {
    pub start: f64,
    pub end: f64,
    pub dt: f64,
    pub dtr: f64,
    pub dump_period: f64,
    pub log_period: f64,
    pub s0: f64,
    pub lth: f64,
    pub temperature: f64,
    pub kv: f64,
    pub box_size: Vec3,
    pub seed: u64,
    pub validate_topology: bool,
    pub dump_vtk: bool,
    pub verbose: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            start: 0.,
            end: 0.,
            dt: 0.001,
            dtr: 0.01,
            dump_period: 1.,
            log_period: 1.,
            s0: 0.,
            lth: 0.,
            temperature: 0.,
            kv: 0.,
            box_size: [1.; 3],
            seed: 5489,
            validate_topology: true,
            dump_vtk: false,
            verbose: false,
        }
    }
}

fn rows(path: &Path) -> Result<Vec<(usize, Vec<String>)>> {
    let text =
        fs::read_to_string(path).map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
    Ok(text
        .lines()
        .enumerate()
        .filter_map(|(line, s)| {
            let fields: Vec<String> = s
                .split('#')
                .next()
                .unwrap_or("")
                .split_whitespace()
                .map(str::to_owned)
                .collect();
            (!fields.is_empty()).then_some((line + 1, fields))
        })
        .collect())
}

fn require(ok: bool, message: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(message.into())
    }
}

fn number(text: &str) -> Result<f64> {
    let value = text
        .parse::<f64>()
        .map_err(|_| format!("Invalid number: {text}"))?;
    require(value.is_finite(), format!("Expected finite number: {text}"))?;
    Ok(value)
}

fn identifier(text: &str) -> Result<Id> {
    let value = text
        .parse::<usize>()
        .map_err(|_| format!("Invalid or out-of-range ID: {text}"))?;
    require(
        value < usize::MAX && (value as u128) < i64::MAX as u128,
        format!("ID must be nonnegative and leave room for the next ID: {text}"),
    )?;
    Ok(value)
}

fn field_count(t: &[String], n: usize) -> Result<()> {
    require(t.len() == n, format!("Wrong field count for {}", t[0]))
}

/// Read the `conf` simulation configuration. Unknown/duplicate keys are errors.
pub fn load_config(path: &Path) -> Result<Config> {
    let mut conf = Config::default();
    let mut seen = BTreeSet::new();
    for (line, t) in rows(path)? {
        let result = (|| -> Result<()> {
            let key = &t[0];
            require(
                seen.insert(key.clone()),
                format!("Duplicate parameter: {key}"),
            )?;
            match key.as_str() {
                "time" => {
                    field_count(&t, 4)?;
                    conf.start = number(&t[1])?;
                    conf.end = number(&t[2])?;
                    conf.dt = number(&t[3])?;
                }
                "dump" => {
                    field_count(&t, 3)?;
                    require(
                        t[1] == "vtk" || t[1] == "topo",
                        "dump format must be vtk or topo",
                    )?;
                    conf.dump_vtk = t[1] == "vtk";
                    conf.dump_period = number(&t[2])?;
                }
                "box" => {
                    field_count(&t, 7)?;
                    for k in 0..3 {
                        conf.box_size[k] = number(&t[k + 1])?;
                        require(
                            t[k + 4] == "p",
                            "This branch supports periodic bulk boundaries only (p p p)",
                        )?;
                    }
                }
                "Lth" => {
                    require(
                        t.len() == 2 || (t.len() == 3 && t[2] == "verbose"),
                        "Lth requires length and optional verbose",
                    )?;
                    conf.lth = number(&t[1])?;
                    conf.verbose = t.len() == 3;
                }
                "log" | "s0" | "T" | "kv" | "dtr" => {
                    field_count(&t, 2)?;
                    let x = number(&t[1])?;
                    match key.as_str() {
                        "log" => conf.log_period = x,
                        "s0" => conf.s0 = x,
                        "T" => conf.temperature = x,
                        "kv" => conf.kv = x,
                        _ => conf.dtr = x,
                    }
                }
                "seed" => {
                    field_count(&t, 2)?;
                    require(
                        !t[1].starts_with('-'),
                        "Seed must be an unsigned 64-bit integer",
                    )?;
                    conf.seed = t[1].parse().map_err(|_| "Invalid seed".to_owned())?;
                }
                "validate_topology" => {
                    field_count(&t, 2)?;
                    require(
                        t[1] == "on" || t[1] == "off",
                        "validate_topology must be on or off",
                    )?;
                    conf.validate_topology = t[1] == "on";
                }
                _ => return Err(format!("Unknown parameter: {key}")),
            }
            Ok(())
        })();
        result.map_err(|e| format!("{}:{line}: {e}", path.display()))?;
    }
    for key in ["time", "dump", "log", "s0", "Lth", "T", "kv", "box"] {
        require(
            seen.contains(key),
            format!("Missing configuration parameter: {key}"),
        )?;
    }
    if !seen.contains("dtr") {
        conf.dtr = 10. * conf.dt;
    }
    require(
        conf.dt > 0. && conf.end >= conf.start,
        "Require dt>0 and end>=start",
    )?;
    let duration = conf.end - conf.start;
    require(
        duration.is_finite() && conf.start + conf.dt > conf.start,
        "Time range/step is outside representable precision",
    )?;
    require(
        conf.dtr.is_finite() && conf.dtr >= conf.dt,
        "Require dtr>=dt",
    )?;
    require(
        conf.dump_period >= conf.dt && conf.log_period >= conf.dt,
        "Output periods must be >= dt",
    )?;
    let limit = (i64::MAX as f64) - 2.;
    require(
        duration / conf.dump_period < limit
            && duration / conf.log_period < limit
            && duration / conf.dtr < limit,
        "Scheduled event count exceeds supported integer range",
    )?;
    require(
        duration / conf.dt < 9_007_199_254_740_991.,
        "Too many time steps for exact step indexing",
    )?;
    require(
        conf.temperature >= 0. && conf.kv >= 0. && conf.lth >= 0.,
        "T, kv and Lth must be nonnegative",
    )?;
    require(conf.s0 >= 0., "s0 must be nonnegative")?;
    require(
        conf.box_size.iter().all(|x| *x > 0.),
        "Box lengths must be positive",
    )?;
    let volume = conf.box_size.iter().product::<f64>();
    require(
        volume.is_finite() && volume > 0.,
        "Box volume must be finite and positive",
    )?;
    Ok(conf)
}

fn referenced<T>(map: &BTreeMap<Id, T>, text: &str) -> Result<Id> {
    let id = identifier(text)?;
    require(
        map.contains_key(&id),
        format!("Unknown referenced ID: {text}"),
    )?;
    Ok(id)
}

/// Load and validate a complete four-section periodic topology.
pub fn load_topology(path: &Path, box_size: Vec3) -> Result<Mesh> {
    let mut mesh = Mesh::new(box_size);
    let sections = ["vertices", "edges", "polygons", "cells"];
    let mut section = None;
    let mut remaining = 0;
    for (line, t) in rows(path)? {
        let result = (|| -> Result<()> {
            if let Some(index) = sections.iter().position(|s| *s == t[0]) {
                require(
                    remaining == 0 && index == section.map_or(0, |s| s + 1),
                    "Wrong topology section order/count",
                )?;
                require(t.len() == 2, "Section header requires its count")?;
                remaining = identifier(&t[1])?;
                require(remaining > 0, "Topology sections must not be empty")?;
                section = Some(index);
                return Ok(());
            }
            require(
                section.is_some() && remaining > 0,
                "Unexpected topology row",
            )?;
            remaining -= 1;
            let id = identifier(&t[0])?;
            match section.unwrap() {
                0 => {
                    require(
                        t.len() == 4,
                        "Vertex requires ID and exactly three coordinates",
                    )?;
                    require(!mesh.vertices.contains_key(&id), "Duplicate vertex ID")?;
                    mesh.vertices.insert(
                        id,
                        Vertex {
                            id,
                            position: [number(&t[1])?, number(&t[2])?, number(&t[3])?],
                            ..Default::default()
                        },
                    );
                    mesh.next_vertex = mesh.next_vertex.max(id + 1);
                }
                1 => {
                    require(t.len() == 3, "Edge requires exactly two endpoint IDs")?;
                    require(!mesh.edges.contains_key(&id), "Duplicate edge ID")?;
                    let vertices = [
                        referenced(&mesh.vertices, &t[1])?,
                        referenced(&mesh.vertices, &t[2])?,
                    ];
                    require(vertices[0] != vertices[1], "Self edge is invalid")?;
                    mesh.edges.insert(
                        id,
                        Edge {
                            id,
                            vertices,
                            ..Default::default()
                        },
                    );
                    mesh.next_edge = mesh.next_edge.max(id + 1);
                }
                2 => {
                    require(t.len() >= 4, "Face requires at least three edges")?;
                    require(!mesh.polygons.contains_key(&id), "Duplicate polygon ID")?;
                    let edges = t[1..]
                        .iter()
                        .map(|s| referenced(&mesh.edges, s))
                        .collect::<Result<Vec<_>>>()?;
                    mesh.polygons.insert(
                        id,
                        Polygon {
                            id,
                            edges,
                            ..Default::default()
                        },
                    );
                    mesh.next_polygon = mesh.next_polygon.max(id + 1);
                }
                _ => {
                    require(t.len() >= 5, "Cell requires at least four face IDs")?;
                    require(!mesh.cells.contains_key(&id), "Duplicate cell ID")?;
                    let polygons = t[1..]
                        .iter()
                        .map(|s| referenced(&mesh.polygons, s))
                        .collect::<Result<Vec<_>>>()?;
                    mesh.cells.insert(
                        id,
                        Cell {
                            id,
                            polygons,
                            ..Default::default()
                        },
                    );
                }
            }
            Ok(())
        })();
        result.map_err(|e| format!("{}:{line}: {e}", path.display()))?;
    }
    require(
        section == Some(3) && remaining == 0,
        "Incomplete topology sections/counts",
    )?;
    mesh.rebuild()?;
    mesh.validate_topology()?;
    mesh.validate_geometry()?;
    Ok(mesh)
}

/// Serialize an input-compatible snapshot (without a time header).
pub fn topology_text(mesh: &Mesh) -> String {
    let mut text = String::new();
    writeln!(text, "vertices {}", mesh.vertices.len()).unwrap();
    for v in mesh.vertices.values() {
        writeln!(
            text,
            "{} {:.17e} {:.17e} {:.17e}",
            v.id, v.position[0], v.position[1], v.position[2]
        )
        .unwrap();
    }
    writeln!(text, "edges {}", mesh.edges.len()).unwrap();
    for e in mesh.edges.values() {
        writeln!(text, "{} {} {}", e.id, e.vertices[0], e.vertices[1]).unwrap();
    }
    writeln!(text, "polygons {}", mesh.polygons.len()).unwrap();
    for p in mesh.polygons.values() {
        write!(text, "{}", p.id).unwrap();
        for e in &p.edges {
            write!(text, " {e}").unwrap();
        }
        text.push('\n');
    }
    writeln!(text, "cells {}", mesh.cells.len()).unwrap();
    for c in mesh.cells.values() {
        write!(text, "{}", c.id).unwrap();
        for p in &c.polygons {
            write!(text, " {p}").unwrap();
        }
        text.push('\n');
    }
    text
}

/// Legacy VTK polydata; faces crossing the periodic seam are omitted.
pub fn vtk_text(mesh: &Mesh) -> String {
    let mut text = format!(
        "# vtk DataFile Version 2.0\npolydata\nASCII\nDATASET POLYDATA\nPOINTS {} double\n",
        mesh.vertices.len()
    );
    let mut indices = BTreeMap::new();
    for (index, v) in mesh.vertices.values().enumerate() {
        indices.insert(v.id, index);
        writeln!(
            text,
            "{:.17e} {:.17e} {:.17e}",
            v.position[0], v.position[1], v.position[2]
        )
        .unwrap();
    }
    let faces: Vec<_> = mesh
        .polygons
        .values()
        .filter(|p| {
            p.edges.iter().all(|id| {
                let e = &mesh.edges[id];
                let a = mesh.vertices[&e.vertices[0]].position;
                let b = mesh.vertices[&e.vertices[1]].position;
                (0..3).all(|k| (b[k] - a[k]).abs() <= mesh.box_size[k] / 2.)
            })
        })
        .collect();
    writeln!(
        text,
        "\nPOLYGONS {} {}",
        faces.len(),
        faces.iter().map(|p| p.vertices.len() + 1).sum::<usize>()
    )
    .unwrap();
    for p in faces {
        write!(text, "{}", p.vertices.len()).unwrap();
        for v in &p.vertices {
            write!(text, " {}", indices[v]).unwrap();
        }
        text.push('\n');
    }
    text
}
