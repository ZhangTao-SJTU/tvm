// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2021-2023 Tao Zhang; coauthor Jennifer Schwarz.
//! Stable-ID mesh, periodic geometry and manifold validation.
use crate::Result;
use std::collections::{BTreeMap, BTreeSet};

pub type Id = usize;
pub type Vec3 = [f64; 3];
pub fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn scale(a: Vec3, b: f64) -> Vec3 {
    [a[0] * b, a[1] * b, a[2] * b]
}
pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}
fn pair(a: Id, b: Id) -> (Id, Id) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}
fn require(ok: bool, message: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(message.into())
    }
}

#[derive(Clone, Debug, Default)]
pub struct Vertex {
    pub id: Id,
    pub position: Vec3,
    pub volume_force: Vec3,
    pub interface_force: Vec3,
    pub velocity: Vec3,
    pub edges: Vec<Id>,
    pub cells: Vec<Id>,
}
#[derive(Clone, Debug, Default)]
pub struct Edge {
    pub id: Id,
    pub vertices: [Id; 2],
    pub vector: Vec3,
    pub center: Vec3,
    pub length: f64,
    pub candidate: bool,
    pub triangle_count: usize,
}
#[derive(Clone, Debug, Default)]
pub struct Polygon {
    pub id: Id,
    pub edges: Vec<Id>,
    pub vertices: Vec<Id>,
    pub cells: Vec<Id>,
    pub center: Vec3,
    pub area: f64,
    pub tension: f64,
}
#[derive(Clone, Debug, Default)]
pub struct Cell {
    pub id: Id,
    pub polygons: Vec<Id>,
    pub vertices: Vec<Id>,
    pub directions: BTreeMap<Id, bool>,
    pub volume: f64,
    pub pressure: f64,
    pub shape_index: f64,
}
#[derive(Clone, Debug)]
pub struct Mesh {
    pub vertices: BTreeMap<Id, Vertex>,
    pub edges: BTreeMap<Id, Edge>,
    pub polygons: BTreeMap<Id, Polygon>,
    pub cells: BTreeMap<Id, Cell>,
    pub box_size: Vec3,
    pub next_vertex: Id,
    pub next_edge: Id,
    pub next_polygon: Id,
}
impl Mesh {
    pub fn new(box_size: Vec3) -> Self {
        Self {
            vertices: BTreeMap::new(),
            edges: BTreeMap::new(),
            polygons: BTreeMap::new(),
            cells: BTreeMap::new(),
            box_size,
            next_vertex: 0,
            next_edge: 0,
            next_polygon: 0,
        }
    }
    pub fn wrap(&self, mut p: Vec3) -> Vec3 {
        for (x, l) in p.iter_mut().zip(self.box_size) {
            *x %= l;
            if *x < 0. {
                *x += l;
            }
        }
        p
    }
    /// Retain the sign of distances equal to half the box length.
    pub fn distance(&self, mut p: Vec3) -> Vec3 {
        for (x, l) in p.iter_mut().zip(self.box_size) {
            *x %= l;
            if *x > 0.5 * l {
                *x -= l;
            }
            if *x < -0.5 * l {
                *x += l;
            }
        }
        p
    }
    fn check_box(&self) -> Result<()> {
        require(
            self.box_size.iter().all(|l| l.is_finite() && *l > 0.),
            "Invalid periodic box length",
        )?;
        let volume = self.box_size.iter().product::<f64>();
        require(
            volume.is_finite() && volume > 0.,
            "Invalid periodic box volume",
        )
    }
    /// Reconstruct polygon cycles starting from the last listed edge.
    fn ordered_vertices(&self, p: &Polygon) -> Result<Vec<Id>> {
        require(
            p.edges.len() >= 3,
            format!("Polygon {} has fewer than three edges", p.id),
        )?;
        let mut incidence: BTreeMap<Id, Vec<Id>> = BTreeMap::new();
        let mut unique = BTreeSet::new();
        for eid in &p.edges {
            let e = self
                .edges
                .get(eid)
                .ok_or_else(|| format!("Missing edge {eid} in polygon {}", p.id))?;
            require(
                unique.insert(*eid) && e.vertices[0] != e.vertices[1],
                "Invalid or repeated face edge",
            )?;
            for v in e.vertices {
                require(
                    self.vertices.contains_key(&v),
                    format!("Missing endpoint {v}"),
                )?;
                incidence.entry(v).or_default().push(*eid);
            }
        }
        require(
            incidence.values().all(|v| v.len() == 2),
            format!("Polygon {} is not a simple cycle", p.id),
        )?;
        let mut next = *p.edges.last().unwrap();
        let start = self.edges[&next].vertices[0];
        let mut current = start;
        let mut ordered = Vec::new();
        let mut used = BTreeSet::new();
        loop {
            require(used.insert(next), "Polygon repeats an edge in its cycle")?;
            ordered.push(current);
            let e = &self.edges[&next];
            current = if e.vertices[0] == current {
                e.vertices[1]
            } else {
                e.vertices[0]
            };
            let neighbors = &incidence[&current];
            next = if neighbors[0] == next {
                neighbors[1]
            } else {
                neighbors[0]
            };
            if current == start {
                break;
            }
        }
        require(used.len() == p.edges.len(), "Disconnected polygon cycles")?;
        Ok(ordered)
    }
    /// Rebuild topology caches and outward winding after loading or a topology change.
    /// Do not use this to repair winding during continuous dynamics; use update_geometry.
    pub fn rebuild(&mut self) -> Result<()> {
        self.check_box()?;
        for v in self.vertices.values_mut() {
            v.edges.clear();
            v.cells.clear();
        }
        for (id, e) in &mut self.edges {
            require(
                *id == e.id && e.vertices[0] != e.vertices[1],
                "Invalid edge ID or self-edge",
            )?;
            e.triangle_count = 0;
            for v in e.vertices {
                self.vertices
                    .get_mut(&v)
                    .ok_or_else(|| format!("Missing edge endpoint {v}"))?
                    .edges
                    .push(*id);
            }
        }
        for id in self.polygons.keys().copied().collect::<Vec<_>>() {
            let vertices = self.ordered_vertices(&self.polygons[&id])?;
            let p = self.polygons.get_mut(&id).unwrap();
            p.vertices = vertices;
            p.cells.clear();
            if p.edges.len() == 3 {
                for e in &p.edges {
                    self.edges.get_mut(e).unwrap().triangle_count += 1;
                }
            }
        }
        for (id, c) in &mut self.cells {
            c.vertices.clear();
            for pid in &c.polygons {
                let p = self
                    .polygons
                    .get_mut(pid)
                    .ok_or_else(|| format!("Missing cell face {pid}"))?;
                p.cells.push(*id);
                for vid in &p.vertices {
                    if !c.vertices.contains(vid) {
                        c.vertices.push(*vid);
                    }
                }
            }
            for vid in &c.vertices {
                self.vertices.get_mut(vid).unwrap().cells.push(*id);
            }
        }
        self.update_face_geometry()?;
        for id in self.cells.keys().copied().collect::<Vec<_>>() {
            let directions = self.cell_directions(&self.cells[&id])?;
            self.cells.get_mut(&id).unwrap().directions = directions;
            let mut volume = self.cell_volume(&self.cells[&id])?;
            if volume < 0. {
                for d in self.cells.get_mut(&id).unwrap().directions.values_mut() {
                    *d = !*d;
                }
                volume = self.cell_volume(&self.cells[&id])?;
            }
            self.set_cell_geometry(id, volume);
        }
        // Keep monotonically increasing IDs after deletions and across rebuilds.
        if let Some(id) = self.vertices.keys().next_back() {
            self.next_vertex = self
                .next_vertex
                .max(id.checked_add(1).ok_or("Vertex ID exhausted")?);
        }
        if let Some(id) = self.edges.keys().next_back() {
            self.next_edge = self
                .next_edge
                .max(id.checked_add(1).ok_or("Edge ID exhausted")?);
        }
        if let Some(id) = self.polygons.keys().next_back() {
            self.next_polygon = self
                .next_polygon
                .max(id.checked_add(1).ok_or("Polygon ID exhausted")?);
        }
        Ok(())
    }
    fn cell_directions(&self, cell: &Cell) -> Result<BTreeMap<Id, bool>> {
        require(
            cell.polygons.len() >= 4,
            format!("Cell {} has fewer than four faces", cell.id),
        )?;
        let mut winding = BTreeMap::new();
        let mut direction: BTreeMap<(Id, Id), bool> = BTreeMap::new();
        let mut count: BTreeMap<(Id, Id), usize> = BTreeMap::new();
        let mut remaining = cell.polygons.clone();
        while !remaining.is_empty() {
            let mut progressed = false;
            for ri in 0..remaining.len() {
                let id = remaining[ri];
                let p = &self.polygons[&id];
                require(p.vertices.len() >= 3, "Invalid face in cell")?;
                let mut known = direction.is_empty();
                let mut forward = true;
                for i in 0..p.vertices.len() {
                    let a = p.vertices[i];
                    let b = p.vertices[(i + 1) % p.vertices.len()];
                    if let Some(found) = direction.get(&pair(a, b)) {
                        let required = (a < b) != *found;
                        require(
                            !known || forward == required,
                            format!("Inconsistent face winding in cell {}", cell.id),
                        )?;
                        known = true;
                        forward = required;
                    }
                }
                if !known {
                    continue;
                }
                require(winding.insert(id, forward).is_none(), "Duplicate cell face")?;
                for i in 0..p.vertices.len() {
                    let a = p.vertices[i];
                    let b = p.vertices[(i + 1) % p.vertices.len()];
                    let key = pair(a, b);
                    direction
                        .entry(key)
                        .or_insert(if forward { a < b } else { a > b });
                    let n = count.entry(key).or_default();
                    *n += 1;
                    require(*n <= 2, "Non-manifold cell shell")?;
                }
                remaining.remove(ri);
                progressed = true;
                break;
            }
            require(
                progressed,
                format!("Disconnected face shell in cell {}", cell.id),
            )?;
        }
        require(count.values().all(|n| *n == 2), "Open face shell")?;
        Ok(winding)
    }
    fn update_face_geometry(&mut self) -> Result<()> {
        self.check_box()?;
        for v in self.vertices.values() {
            require(
                v.position.iter().all(|x| x.is_finite()),
                "Nonfinite vertex coordinate",
            )?;
        }
        for id in self.edges.keys().copied().collect::<Vec<_>>() {
            let e = &self.edges[&id];
            let a = self
                .vertices
                .get(&e.vertices[0])
                .ok_or("Missing edge endpoint")?
                .position;
            let b = self
                .vertices
                .get(&e.vertices[1])
                .ok_or("Missing edge endpoint")?
                .position;
            let vector = self.distance(sub(b, a));
            let length = vector[0].hypot(vector[1]).hypot(vector[2]);
            require(length.is_finite(), "Nonfinite edge length")?;
            let e = self.edges.get_mut(&id).unwrap();
            e.vector = vector;
            e.length = length;
            e.center = add(a, scale(vector, 0.5));
        }
        for id in self.polygons.keys().copied().collect::<Vec<_>>() {
            let p = &self.polygons[&id];
            let first = p.edges.first().ok_or("Empty polygon")?;
            let e = self.edges.get(first).ok_or("Missing polygon edge")?;
            let origin = self.vertices[&e.vertices[0]].position;
            let mut weighted = [0.; 3];
            let mut perimeter = 0.;
            for eid in &p.edges {
                let e = self.edges.get(eid).ok_or("Missing polygon edge")?;
                weighted = add(
                    weighted,
                    scale(self.distance(sub(e.center, origin)), e.length),
                );
                perimeter += e.length;
            }
            require(
                perimeter.is_finite() && perimeter > 0.,
                format!("Zero/nonfinite perimeter at polygon {id}"),
            )?;
            let center = add(scale(weighted, 1. / perimeter), origin);
            require(
                center.iter().all(|x| x.is_finite()),
                "Nonfinite polygon center",
            )?;
            let mut area = 0.;
            for eid in &p.edges {
                let e = &self.edges[eid];
                let a = self.distance(sub(self.vertices[&e.vertices[0]].position, center));
                let b = self.distance(sub(self.vertices[&e.vertices[1]].position, center));
                area += 0.5 * norm(cross(a, b));
            }
            require(area.is_finite(), "Nonfinite polygon area")?;
            let p = self.polygons.get_mut(&id).unwrap();
            p.center = center;
            p.area = area;
        }
        Ok(())
    }
    fn cell_volume(&self, c: &Cell) -> Result<f64> {
        let first = c.polygons.first().ok_or("Empty cell geometry")?;
        let vid = self.polygons[first]
            .vertices
            .first()
            .ok_or("Empty face geometry")?;
        let origin = self.vertices[vid].position;
        let mut volume = 0.;
        for pid in &c.polygons {
            let p = &self.polygons[pid];
            let cc = self.distance(sub(p.center, origin));
            let forward = *c
                .directions
                .get(pid)
                .ok_or("Missing cell face winding; call rebuild first")?;
            for i in 0..p.vertices.len() {
                let a = self.distance(sub(self.vertices[&p.vertices[i]].position, p.center));
                let b = self.distance(sub(
                    self.vertices[&p.vertices[(i + 1) % p.vertices.len()]].position,
                    p.center,
                ));
                let tetrahedron = (1. / 6.) * dot(cc, cross(a, b));
                if forward {
                    volume += tetrahedron;
                } else {
                    volume -= tetrahedron;
                }
            }
        }
        require(
            volume.is_finite(),
            format!("Nonfinite volume in cell {}", c.id),
        )?;
        Ok(volume)
    }
    fn set_cell_geometry(&mut self, id: Id, volume: f64) {
        let area = self.cells[&id]
            .polygons
            .iter()
            .map(|p| self.polygons[p].area)
            .sum::<f64>();
        let c = self.cells.get_mut(&id).unwrap();
        c.volume = volume;
        c.shape_index = area * volume.powf(-2. / 3.);
    }
    pub fn update_geometry(&mut self) -> Result<()> {
        self.update_face_geometry()?;
        for id in self.cells.keys().copied().collect::<Vec<_>>() {
            let v = self.cell_volume(&self.cells[&id])?;
            self.set_cell_geometry(id, v);
        }
        Ok(())
    }
    pub fn validate_geometry(&self) -> Result<()> {
        self.check_box()?;
        for c in self.cells.values() {
            require(
                c.volume.is_finite() && c.volume > 0.,
                format!(
                    "Nonpositive/nonfinite cell volume at ID {}; reduce dt or inspect geometry",
                    c.id
                ),
            )?;
        }
        let total = self.cells.values().map(|c| c.volume).sum::<f64>();
        let expected = self.box_size.iter().product::<f64>();
        require(
            total.is_finite() && (total / expected - 1.).abs() <= 1e-8,
            "Periodic total volume differs from box volume; reduce dt or inspect geometry",
        )
    }
    /// Validate the closed, trivalent-face, fourfold-vertex periodic bulk complex.
    /// Derived caches are deliberately ignored: this also validates a proposed topology before rebuilding.
    pub fn validate_topology(&self) -> Result<()> {
        require(
            !self.vertices.is_empty()
                && !self.edges.is_empty()
                && !self.polygons.is_empty()
                && !self.cells.is_empty(),
            "Topology validation: empty network",
        )?;
        for (id, v) in &self.vertices {
            require(
                *id == v.id && v.position.iter().all(|x| x.is_finite()),
                "Invalid vertex ID or coordinate",
            )?;
        }
        let mut endpoint_pairs = BTreeSet::new();
        let mut vertex_edges: BTreeMap<Id, BTreeSet<Id>> = BTreeMap::new();
        let mut edge_faces: BTreeMap<Id, Vec<Id>> = BTreeMap::new();
        let mut face_vertices: BTreeMap<Id, BTreeSet<Id>> = BTreeMap::new();
        let mut face_cells: BTreeMap<Id, Vec<Id>> = BTreeMap::new();
        let mut vertex_cells: BTreeMap<Id, BTreeSet<Id>> = BTreeMap::new();
        for (id, e) in &self.edges {
            let [a, b] = e.vertices;
            require(
                *id == e.id && self.vertices.contains_key(&a) && self.vertices.contains_key(&b),
                "Invalid edge ID or endpoints",
            )?;
            require(
                a != b && endpoint_pairs.insert(pair(a, b)),
                "Self edge or duplicate endpoint pair",
            )?;
            vertex_edges.entry(a).or_default().insert(*id);
            vertex_edges.entry(b).or_default().insert(*id);
        }
        let mut reached = BTreeSet::new();
        let mut pending = vec![*self.vertices.keys().next().unwrap()];
        while let Some(v) = pending.pop() {
            if !reached.insert(v) {
                continue;
            }
            if let Some(es) = vertex_edges.get(&v) {
                for e in es {
                    pending.extend(self.edges[e].vertices);
                }
            }
        }
        require(
            reached.len() == self.vertices.len(),
            "Global network is disconnected",
        )?;
        for (id, p) in &self.polygons {
            require(*id == p.id, "Invalid polygon ID")?;
            face_vertices.insert(*id, self.ordered_vertices(p)?.into_iter().collect());
            for e in &p.edges {
                edge_faces.entry(*e).or_default().push(*id);
            }
        }
        for (id, c) in &self.cells {
            require(
                *id == c.id && c.polygons.len() >= 4,
                "Invalid cell ID or fewer than four faces",
            )?;
            let mut unique = BTreeSet::new();
            let mut cv = BTreeSet::new();
            let mut ce: BTreeMap<Id, Vec<Id>> = BTreeMap::new();
            for pid in &c.polygons {
                require(
                    self.polygons.contains_key(pid) && unique.insert(*pid),
                    "Missing or repeated cell face",
                )?;
                face_cells.entry(*pid).or_default().push(*id);
                for v in &face_vertices[pid] {
                    cv.insert(*v);
                    vertex_cells.entry(*v).or_default().insert(*id);
                }
                for e in &self.polygons[pid].edges {
                    ce.entry(*e).or_default().push(*pid);
                }
            }
            let mut adjacent: BTreeMap<Id, Vec<Id>> = BTreeMap::new();
            for faces in ce.values() {
                require(
                    faces.len() == 2,
                    "Cell shell edge is not incident to two faces",
                )?;
                adjacent.entry(faces[0]).or_default().push(faces[1]);
                adjacent.entry(faces[1]).or_default().push(faces[0]);
            }
            let mut visited = BTreeSet::new();
            let mut pending = vec![*unique.first().unwrap()];
            while let Some(p) = pending.pop() {
                if !visited.insert(p) {
                    continue;
                }
                if let Some(ns) = adjacent.get(&p) {
                    pending.extend(ns);
                }
            }
            require(visited.len() == unique.len(), "Disconnected cell shell")?;
            require(
                cv.len() as i128 - ce.len() as i128 + unique.len() as i128 == 2,
                "Cell Euler relation != 2",
            )?;
        }
        for id in self.vertices.keys() {
            require(
                vertex_edges.get(id).map_or(0, |v| v.len()) == 4,
                format!("Vertex degree != 4 at {id}"),
            )?;
            require(
                vertex_cells.get(id).map_or(0, |v| v.len()) == 4,
                format!("Vertex cell incidence != 4 at {id}"),
            )?;
        }
        let mut face_pairs = BTreeSet::new();
        for id in self.edges.keys() {
            let list = edge_faces.get(id).ok_or("Edge has no faces")?;
            require(list.len() == 3, "Edge face incidence != 3")?;
            for i in 0..list.len() {
                for j in i + 1..list.len() {
                    require(
                        face_pairs.insert(pair(list[i], list[j])),
                        "Two faces share multiple edges",
                    )?;
                }
            }
        }
        let mut cell_pairs = BTreeSet::new();
        for id in self.polygons.keys() {
            let list = face_cells.get(id).ok_or("Face has no cells")?;
            require(list.len() == 2, "Face cell incidence != 2")?;
            require(
                cell_pairs.insert(pair(list[0], list[1])),
                "Two cells share multiple faces",
            )?;
        }
        require(
            self.vertices.len() as i128 - self.edges.len() as i128 + self.polygons.len() as i128
                - self.cells.len() as i128
                == 0,
            "Periodic global Euler relation != 0",
        )
    }
}
