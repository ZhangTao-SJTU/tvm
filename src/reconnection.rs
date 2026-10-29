// Copyright 2021-2026 Tao Zhang; GPL-3.0-or-later.
//! Okuda I-H / H-I rearrangements for a periodic bulk mesh.
//!
//! Edits are constructed in a private mesh clone. A rejected candidate leaves
//! coordinates, connectivity and the monotonically increasing object IDs intact.
use crate::model::{add, cross, dot, norm, scale, sub, Edge, Id, Mesh, Polygon, Vec3, Vertex};
use crate::Result;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as FmtWrite;
use std::fs::OpenOptions;
use std::io::Write;
use std::path::PathBuf;

#[derive(Clone, Debug)]
pub struct Reconnection {
    pub lth: f64,
    pub verbose: bool,
    pub count_ih: u64,
    pub count_hi: u64,
    pub count_events: u64,
    pub rejected_geometry: u64,
    pub rejected_topology: u64,
    pub last_rejection: String,
    pub events: Vec<String>,
    pub time: f64,
    /// Directory for verbose event VTK files. None uses the working directory.
    pub diagnostic_directory: Option<PathBuf>,
}

#[derive(Debug)]
enum Rejection {
    Geometry(String),
    Topology(String),
}
type Attempt<T> = std::result::Result<T, Rejection>;
fn topology(message: &str) -> Rejection {
    Rejection::Topology(message.into())
}
fn geometry(message: &str) -> Rejection {
    Rejection::Geometry(message.into())
}
fn require(ok: bool, message: &str) -> Attempt<()> {
    if ok {
        Ok(())
    } else {
        Err(topology(message))
    }
}
fn geometric(ok: bool, message: &str) -> Attempt<()> {
    if ok {
        Ok(())
    } else {
        Err(geometry(message))
    }
}
fn distinct(ids: &[Id]) -> bool {
    ids.iter().copied().collect::<BTreeSet<_>>().len() == ids.len()
}
fn finite(v: Vec3) -> bool {
    v.iter().all(|x| x.is_finite())
}
fn epsilon(m: &Mesh) -> f64 {
    1e-12 * m.box_size.iter().copied().fold(1., f64::max)
}
fn common(a: &[Id], b: &[Id]) -> Attempt<Id> {
    let ids: Vec<_> = a.iter().copied().filter(|v| b.contains(v)).collect();
    require(ids.len() == 1, "missing or nonunique shared object")?;
    Ok(ids[0])
}
fn any_common(a: &[Id], b: &[Id]) -> bool {
    a.iter().any(|v| b.contains(v))
}
fn common_face(m: &Mesh, a: Id, b: Id) -> Attempt<Id> {
    let a = m.cells.get(&a).ok_or_else(|| topology("missing cell"))?;
    let b = m.cells.get(&b).ok_or_else(|| topology("missing cell"))?;
    common(&a.polygons, &b.polygons)
}
fn common_edge(m: &Mesh, a: Id, b: Id) -> Attempt<Id> {
    let a = m.polygons.get(&a).ok_or_else(|| topology("missing face"))?;
    let b = m.polygons.get(&b).ok_or_else(|| topology("missing face"))?;
    common(&a.edges, &b.edges)
}
fn other(m: &Mesh, e: Id, v: Id) -> Attempt<Id> {
    let e = m.edges.get(&e).ok_or_else(|| topology("missing edge"))?;
    if e.vertices[0] == v && e.vertices[1] != v {
        Ok(e.vertices[1])
    } else if e.vertices[1] == v && e.vertices[0] != v {
        Ok(e.vertices[0])
    } else {
        Err(topology("invalid edge endpoint incidence"))
    }
}
fn has_edges(m: &Mesh, p: Id, edges: &[Id]) -> Attempt<()> {
    let p = m.polygons.get(&p).ok_or_else(|| topology("missing face"))?;
    require(
        edges
            .iter()
            .all(|e| m.edges.contains_key(e) && p.edges.iter().filter(|x| *x == e).count() == 1),
        "inconsistent face incidence",
    )
}
fn direction(m: &Mesh, a: Vec3, b: Vec3) -> Attempt<Vec3> {
    geometric(finite(a) && finite(b), "nonfinite neighboring position")?;
    let d = m.distance(sub(b, a));
    let length = norm(d);
    geometric(
        finite(d) && length.is_finite() && length > epsilon(m),
        "unresolved neighboring direction",
    )?;
    Ok(scale(d, 1. / length))
}
fn add_vertex(m: &mut Mesh, position: Vec3) -> Id {
    let id = m.next_vertex;
    m.next_vertex += 1;
    m.vertices.insert(
        id,
        Vertex {
            id,
            position,
            ..Default::default()
        },
    );
    id
}
fn add_edge(m: &mut Mesh, a: Id, b: Id) -> Id {
    let id = m.next_edge;
    m.next_edge += 1;
    let vertices = if a < b { [a, b] } else { [b, a] };
    m.edges.insert(
        id,
        Edge {
            id,
            vertices,
            ..Default::default()
        },
    );
    id
}
fn edit_face(m: &mut Mesh, p: Id, remove: &[Id], append: &[Id]) {
    let p = m.polygons.get_mut(&p).expect("validated face");
    p.edges.retain(|e| !remove.contains(e));
    p.edges.extend_from_slice(append);
}
fn triangle_counts(m: &mut Mesh) {
    for e in m.edges.values_mut() {
        e.triangle_count = 0;
    }
    for p in m.polygons.values() {
        if p.edges.len() == 3 {
            for e in &p.edges {
                if let Some(e) = m.edges.get_mut(e) {
                    e.triangle_count += 1;
                }
            }
        }
    }
}
fn h_candidate(m: &Mesh, id: Id) -> bool {
    m.polygons.get(&id).is_some_and(|p| {
        p.edges.len() == 3
            && p.edges.iter().all(|e| {
                m.edges
                    .get(e)
                    .is_some_and(|e| e.candidate && e.triangle_count == 1)
            })
    })
}

/// Nine adjacent faces in the local neighborhood:
/// p14, p25, p36, p12, p23, p13, p45, p56, p46.
fn faces(m: &Mesh, c: [Id; 5]) -> Attempt<[Id; 9]> {
    let [top, bottom, a, b, d] = c;
    let p = [
        common_face(m, a, d)?,
        common_face(m, b, a)?,
        common_face(m, d, b)?,
        common_face(m, a, top)?,
        common_face(m, b, top)?,
        common_face(m, d, top)?,
        common_face(m, a, bottom)?,
        common_face(m, b, bottom)?,
        common_face(m, d, bottom)?,
    ];
    require(distinct(&p), "repeated local face")?;
    Ok(p)
}

struct Proposal {
    mesh: Mesh,
    cells: [Id; 5],
    before_edges: Vec<Id>,
    after_edges: Vec<Id>,
    new_vertex_cells: Vec<(Id, Vec<Id>)>,
}

impl Reconnection {
    pub fn new(lth: f64, verbose: bool) -> Self {
        Self {
            lth,
            verbose,
            count_ih: 0,
            count_hi: 0,
            count_events: 0,
            rejected_geometry: 0,
            rejected_topology: 0,
            last_rejection: String::new(),
            events: Vec::new(),
            time: 0.,
            diagnostic_directory: None,
        }
    }

    pub fn start(&mut self, m: &mut Mesh) -> Result<()> {
        if !self.lth.is_finite() || self.lth < 0. {
            return Err("Invalid reconnection threshold".into());
        }
        if self.lth == 0. {
            return Ok(());
        }
        m.update_geometry()?;
        m.validate_geometry()?;
        for e in m.edges.values_mut() {
            e.candidate = e.length < self.lth;
        }
        triangle_counts(m);
        // Snapshot the initial candidates: new edges/faces wait until next sweep.
        let triangles: Vec<_> = m
            .polygons
            .keys()
            .copied()
            .filter(|p| h_candidate(m, *p))
            .collect();
        let edges: Vec<_> = m.edges.keys().copied().collect();
        for e in edges {
            if m.edges
                .get(&e)
                .is_some_and(|e| e.candidate && e.triangle_count == 0)
            {
                self.i_h(m, e)?;
            }
        }
        for p in triangles {
            if h_candidate(m, p) {
                self.h_i(m, p)?;
            }
        }
        m.rebuild()?;
        Ok(())
    }

    /// Attempt an I-H event directly. Threshold eligibility is enforced by start;
    /// this entry point also supports controlled transition and regression tests.
    pub fn i_h(&mut self, m: &mut Mesh, e: Id) -> Result<bool> {
        let proposal = self.propose_ih(m, e);
        self.finish(m, proposal, true)
    }
    pub fn h_i(&mut self, m: &mut Mesh, p: Id) -> Result<bool> {
        let proposal = self.propose_hi(m, p);
        self.finish(m, proposal, false)
    }

    fn finish(&mut self, m: &mut Mesh, proposal: Attempt<Proposal>, ih: bool) -> Result<bool> {
        let mut proposal = match proposal {
            Ok(p) => p,
            Err(reason) => {
                let (kind, text) = match reason {
                    Rejection::Geometry(t) => {
                        self.rejected_geometry += 1;
                        ("geometry", t)
                    }
                    Rejection::Topology(t) => {
                        self.rejected_topology += 1;
                        ("topology", t)
                    }
                };
                self.last_rejection = text.clone();
                if self.verbose {
                    self.events.push(format!("rejected {kind}: {text}\n"));
                }
                return Ok(false);
            }
        };
        // Additional transactional guarantee: never publish an invalid mesh.
        if let Err(reason) = proposal.mesh.validate_topology() {
            return self.finish(
                m,
                Err(topology(&format!("proposed topology: {reason}"))),
                ih,
            );
        }
        if let Err(reason) = proposal
            .mesh
            .rebuild()
            .and_then(|_| proposal.mesh.validate_geometry())
        {
            return self.finish(
                m,
                Err(geometry(&format!("proposed geometry: {reason}"))),
                ih,
            );
        }
        triangle_counts(&mut proposal.mesh);
        // Retain incidence-list insertion order until the end of a sweep.
        // This order names side cells, hence it determines new object IDs during
        // subsequent interacting events. Preserve it across transactional rebuilds.
        for (id, vertex) in &m.vertices {
            if let Some(survivor) = proposal.mesh.vertices.get_mut(id) {
                survivor.cells = vertex.cells.clone();
            }
        }
        for (id, cells) in &proposal.new_vertex_cells {
            proposal.mesh.vertices.get_mut(id).unwrap().cells = cells.clone();
        }
        if self.verbose {
            self.write_diagnostics(m, &proposal, ih)?;
        }
        let label = if ih { "I->H" } else { "H->I" };
        let mut record = format!(
            "{label} time {:.17e} event {}\n",
            self.time, self.count_events
        );
        if self.verbose {
            for cell in proposal.cells {
                write!(record, "{cell} ").unwrap();
            }
            record.push('\n');
            for cell in proposal.cells {
                let cell = &m.cells[&cell];
                write!(record, "{}", cell.polygons.len()).unwrap();
                for face in &cell.polygons {
                    write!(record, " {}", m.polygons[face].edges.len()).unwrap();
                }
                record.push('\n');
            }
            // Include both removed and created endpoint coordinates for auditing.
            let mut positions = BTreeMap::new();
            for (mesh, edges) in [
                (m as &Mesh, &proposal.before_edges),
                (&proposal.mesh, &proposal.after_edges),
            ] {
                for edge in edges {
                    for vertex in mesh.edges[edge].vertices {
                        positions.insert(vertex, mesh.vertices[&vertex].position);
                    }
                }
            }
            for id in positions.keys() {
                write!(record, "{id} ").unwrap();
            }
            record.push('\n');
            for p in positions.values() {
                writeln!(record, "{:.17e} {:.17e} {:.17e}", p[0], p[1], p[2]).unwrap();
            }
            for cell in proposal.cells {
                let cell = &proposal.mesh.cells[&cell];
                write!(record, "{}", cell.polygons.len()).unwrap();
                for face in &cell.polygons {
                    write!(record, " {}", proposal.mesh.polygons[face].edges.len()).unwrap();
                }
                record.push('\n');
            }
            record.push('\n');
        }
        self.events.push(record);
        *m = proposal.mesh;
        if ih {
            self.count_ih += 1;
        } else {
            self.count_hi += 1;
        }
        self.count_events += 1;
        self.last_rejection.clear();
        Ok(true)
    }

    fn propose_ih(&self, m: &Mesh, eid: Id) -> Attempt<Proposal> {
        let e = m
            .edges
            .get(&eid)
            .ok_or_else(|| topology("I->H missing edge"))?;
        let [v10, v11] = e.vertices;
        require(v10 != v11, "I->H invalid edge endpoints")?;
        let v10ref = m
            .vertices
            .get(&v10)
            .ok_or_else(|| topology("I->H missing endpoint"))?;
        let v11ref = m
            .vertices
            .get(&v11)
            .ok_or_else(|| topology("I->H missing endpoint"))?;
        let eps = epsilon(m);
        geometric(
            self.lth.is_finite()
                && self.lth > eps
                && e.length.is_finite()
                && e.length > eps
                && finite(e.center)
                && finite(e.vector)
                && finite(v10ref.position)
                && finite(v11ref.position),
            "I->H unresolved edge length or threshold",
        )?;
        let a = &v10ref.cells;
        let b = &v11ref.cells;
        require(
            a.len() == 4 && b.len() == 4 && distinct(a) && distinct(b),
            "I->H invalid local incidence",
        )?;
        let side: Vec<_> = a.iter().copied().filter(|c| b.contains(c)).collect();
        require(side.len() == 3, "I->H requires three side cells")?;
        let top = *a
            .iter()
            .find(|c| !b.contains(c))
            .ok_or_else(|| topology("I->H missing top cell"))?;
        let bottom = *b
            .iter()
            .find(|c| !a.contains(c))
            .ok_or_else(|| topology("I->H missing bottom cell"))?;
        let cells = [top, bottom, side[0], side[1], side[2]];
        require(distinct(&cells), "I->H repeated neighboring cell")?;
        require(
            !any_common(&m.cells[&top].polygons, &m.cells[&bottom].polygons),
            "I->H top and bottom already share a face",
        )?;
        let p = faces(m, cells)?;
        for i in 3..6 {
            require(
                !any_common(&m.polygons[&p[i]].edges, &m.polygons[&p[i + 3]].edges),
                "I->H top and bottom faces already share an edge",
            )?;
        }
        let old = [
            common_edge(m, p[3], p[5])?,
            common_edge(m, p[4], p[3])?,
            common_edge(m, p[5], p[4])?,
            common_edge(m, p[6], p[8])?,
            common_edge(m, p[7], p[6])?,
            common_edge(m, p[8], p[7])?,
        ];
        let mut outer = [0; 6];
        for i in 0..6 {
            outer[i] = other(m, old[i], if i < 3 { v10 } else { v11 })?;
        }
        for i in 0..3 {
            require(
                m.polygons[&p[i]].edges.len() >= 4,
                "I->H proposed two-edge face",
            )?;
            has_edges(m, p[i], &[eid, old[i], old[i + 3]])?;
        }
        for (i, pair) in [[0, 1], [1, 2], [0, 2]].iter().enumerate() {
            has_edges(m, p[i + 3], &[old[pair[0]], old[pair[1]]])?;
            has_edges(m, p[i + 6], &[old[pair[0] + 3], old[pair[1] + 3]])?;
        }
        require(
            distinct(&old) && !old.contains(&eid),
            "I->H repeated local edge",
        )?;
        require(
            (0..3).all(|i| outer[i] != outer[i + 3])
                && outer.iter().all(|v| *v != v10 && *v != v11),
            "I->H repeated new endpoint pair or removed exterior vertex",
        )?;
        let u = direction(m, [0.; 3], e.vector)?;
        let mut projected = [[0.; 3]; 3];
        for i in 0..3 {
            let w = scale(
                add(
                    direction(m, e.center, m.vertices[&outer[i]].position)?,
                    direction(m, e.center, m.vertices[&outer[i + 3]].position)?,
                ),
                0.5,
            );
            projected[i] = sub(w, scale(u, dot(w, u)));
        }
        let longest = (0..3)
            .map(|i| norm(sub(projected[i], projected[(i + 1) % 3])))
            .fold(0., f64::max);
        geometric(
            longest.is_finite() && longest > 1e-12,
            "I->H collapsed projected triangle",
        )?;
        let mut next = [[0.; 3]; 3];
        for i in 0..3 {
            let r = add(e.center, scale(projected[i], self.lth / longest));
            geometric(finite(r), "I->H nonfinite proposed coordinate")?;
            next[i] = m.wrap(r);
        }
        let normal = cross(
            m.distance(sub(next[1], next[0])),
            m.distance(sub(next[2], next[0])),
        );
        geometric(
            finite(normal) && norm(normal) > eps * self.lth,
            "I->H collinear proposed triangle",
        )?;
        for i in 0..6 {
            let d = m.distance(sub(m.vertices[&outer[i]].position, next[i % 3]));
            geometric(
                finite(d) && norm(d) > eps,
                "I->H collapsed proposed outer edge",
            )?;
        }
        require(
            m.next_vertex <= usize::MAX - 3
                && m.next_edge <= usize::MAX - 9
                && m.next_polygon < usize::MAX
                && self.count_events < u64::MAX,
            "I->H object/event ID space exhausted",
        )?;
        let mut n = m.clone();
        let v = [
            add_vertex(&mut n, next[0]),
            add_vertex(&mut n, next[1]),
            add_vertex(&mut n, next[2]),
        ];
        let tri = [
            add_edge(&mut n, v[0], v[1]),
            add_edge(&mut n, v[0], v[2]),
            add_edge(&mut n, v[1], v[2]),
        ];
        let face = n.next_polygon;
        n.next_polygon += 1;
        n.polygons.insert(
            face,
            Polygon {
                id: face,
                edges: tri.to_vec(),
                ..Default::default()
            },
        );
        n.cells.get_mut(&top).unwrap().polygons.push(face);
        n.cells.get_mut(&bottom).unwrap().polygons.push(face);
        let mut new = [0; 6];
        for i in 0..6 {
            new[i] = add_edge(&mut n, v[i % 3], outer[i]);
        }
        for i in 0..3 {
            edit_face(
                &mut n,
                p[i],
                &[eid, old[i], old[i + 3]],
                &[new[i], new[i + 3]],
            );
        }
        for (i, pair) in [[0, 1], [1, 2], [0, 2]].iter().enumerate() {
            let t = tri[[0, 2, 1][i]];
            edit_face(
                &mut n,
                p[i + 3],
                &[old[pair[0]], old[pair[1]]],
                &[t, new[pair[0]], new[pair[1]]],
            );
            edit_face(
                &mut n,
                p[i + 6],
                &[old[pair[0] + 3], old[pair[1] + 3]],
                &[t, new[pair[0] + 3], new[pair[1] + 3]],
            );
        }
        n.vertices.remove(&v10);
        n.vertices.remove(&v11);
        n.edges.remove(&eid);
        for e in old {
            n.edges.remove(&e);
        }
        let mut before_edges = vec![eid];
        before_edges.extend_from_slice(&old);
        let mut after_edges = tri.to_vec();
        after_edges.extend_from_slice(&new);
        let [top, bottom, a, b, d] = cells;
        let new_vertex_cells = vec![
            (v[0], vec![a, d, top, bottom]),
            (v[1], vec![a, b, top, bottom]),
            (v[2], vec![b, d, top, bottom]),
        ];
        Ok(Proposal {
            mesh: n,
            cells,
            before_edges,
            after_edges,
            new_vertex_cells,
        })
    }

    fn propose_hi(&self, m: &Mesh, pid: Id) -> Attempt<Proposal> {
        let triangle = m
            .polygons
            .get(&pid)
            .ok_or_else(|| topology("H->I missing triangle"))?;
        require(
            triangle.edges.len() == 3
                && triangle.vertices.len() == 3
                && distinct(&triangle.vertices),
            "H->I requires a triangular face cycle",
        )?;
        let eps = epsilon(m);
        geometric(
            self.lth.is_finite() && self.lth > eps,
            "H->I unresolved threshold",
        )?;
        let v = [
            triangle.vertices[0],
            triangle.vertices[1],
            triangle.vertices[2],
        ];
        for id in v {
            let vertex = &m.vertices[&id];
            geometric(finite(vertex.position), "H->I nonfinite input coordinate")?;
            require(
                vertex.cells.len() == 4 && distinct(&vertex.cells),
                "H->I invalid local incidence",
            )?;
        }
        let tb: Vec<_> = m.vertices[&v[0]]
            .cells
            .iter()
            .copied()
            .filter(|c| m.cells[c].polygons.contains(&pid))
            .collect();
        require(tb.len() == 2, "H->I requires two adjacent cells")?;
        let mut side = [0; 3];
        for i in 0..3 {
            let a: Vec<_> = m.vertices[&v[i]]
                .cells
                .iter()
                .copied()
                .filter(|c| !tb.contains(c))
                .collect();
            side[i] = common(&a, &m.vertices[&v[(i + 1) % 3]].cells)?;
        }
        let cells = [tb[0], tb[1], side[0], side[1], side[2]];
        require(distinct(&cells), "H->I repeated neighboring cell")?;
        let p = faces(m, cells)?;
        for i in 0..3 {
            require(
                !any_common(&m.polygons[&p[i]].edges, &m.polygons[&p[(i + 1) % 3]].edges),
                "H->I side faces already share an edge",
            )?;
        }
        let tri = [
            common_edge(m, p[3], p[6])?,
            common_edge(m, p[5], p[8])?,
            common_edge(m, p[4], p[7])?,
        ];
        let old = [
            common_edge(m, p[3], p[5])?,
            common_edge(m, p[4], p[3])?,
            common_edge(m, p[5], p[4])?,
            common_edge(m, p[6], p[8])?,
            common_edge(m, p[7], p[6])?,
            common_edge(m, p[8], p[7])?,
        ];
        let mut outer = [0; 6];
        for i in 0..6 {
            outer[i] = other(m, old[i], v[i % 3])?;
        }
        has_edges(m, pid, &tri)?;
        for i in 0..3 {
            has_edges(m, p[i], &[old[i], old[i + 3]])?;
        }
        for (i, pair) in [[0, 1], [1, 2], [0, 2]].iter().enumerate() {
            let t = tri[[0, 2, 1][i]];
            require(
                m.polygons[&p[i + 3]].edges.len() >= 4 && m.polygons[&p[i + 6]].edges.len() >= 4,
                "H->I proposed two-edge face",
            )?;
            has_edges(m, p[i + 3], &[t, old[pair[0]], old[pair[1]]])?;
            has_edges(m, p[i + 6], &[t, old[pair[0] + 3], old[pair[1] + 3]])?;
        }
        let r78 = m.distance(sub(m.vertices[&v[1]].position, m.vertices[&v[0]].position));
        let r79 = m.distance(sub(m.vertices[&v[2]].position, m.vertices[&v[0]].position));
        let r0 = add(m.vertices[&v[0]].position, scale(add(r78, r79), 1. / 3.));
        let normal = cross(r78, r79);
        let length = norm(normal);
        geometric(
            finite(r0) && length.is_finite() && length > eps * norm(r78).max(norm(r79)),
            "H->I unresolved triangle normal",
        )?;
        let normal = scale(normal, 1. / length);
        let d0 = *m.cells[&tb[0]]
            .directions
            .get(&pid)
            .ok_or_else(|| topology("H->I missing face orientation"))?;
        let d1 = *m.cells[&tb[1]]
            .directions
            .get(&pid)
            .ok_or_else(|| topology("H->I missing face orientation"))?;
        require(d0 != d1, "H->I inconsistent adjacent face orientations")?;
        let d = if d0 { -1. } else { 1. };
        let proposed = [
            add(r0, scale(normal, d * 0.5 * self.lth)),
            sub(r0, scale(normal, d * 0.5 * self.lth)),
        ];
        geometric(
            proposed.iter().all(|r| finite(*r)),
            "H->I nonfinite proposed coordinate",
        )?;
        let next = [m.wrap(proposed[0]), m.wrap(proposed[1])];
        for i in 0..6 {
            let distance = m.distance(sub(m.vertices[&outer[i]].position, next[i / 3]));
            geometric(
                finite(distance) && norm(distance) > eps,
                "H->I collapsed proposed outer edge",
            )?;
            require(
                !v.contains(&outer[i]),
                "H->I exterior endpoint is removed vertex",
            )?;
        }
        let mut old_all = tri.to_vec();
        old_all.extend_from_slice(&old);
        require(
            distinct(&outer[..3])
                && distinct(&outer[3..])
                && distinct(&old_all)
                && !p.contains(&pid),
            "H->I repeated local edge, face or endpoint pair",
        )?;
        require(
            m.next_vertex <= usize::MAX - 2
                && m.next_edge <= usize::MAX - 7
                && self.count_events < u64::MAX,
            "H->I object/event ID space exhausted",
        )?;
        let mut n = m.clone();
        let newv = [add_vertex(&mut n, next[0]), add_vertex(&mut n, next[1])];
        let central = add_edge(&mut n, newv[0], newv[1]);
        let mut new = [0; 6];
        for i in 0..6 {
            new[i] = add_edge(&mut n, newv[i / 3], outer[i]);
        }
        for i in 0..3 {
            edit_face(
                &mut n,
                p[i],
                &[old[i], old[i + 3]],
                &[central, new[i], new[i + 3]],
            );
        }
        for (i, pair) in [[0, 1], [1, 2], [0, 2]].iter().enumerate() {
            let t = tri[[0, 2, 1][i]];
            edit_face(
                &mut n,
                p[i + 3],
                &[t, old[pair[0]], old[pair[1]]],
                &[new[pair[0]], new[pair[1]]],
            );
            edit_face(
                &mut n,
                p[i + 6],
                &[t, old[pair[0] + 3], old[pair[1] + 3]],
                &[new[pair[0] + 3], new[pair[1] + 3]],
            );
        }
        for id in v {
            n.vertices.remove(&id);
        }
        for id in &old_all {
            n.edges.remove(id);
        }
        for c in tb {
            n.cells
                .get_mut(&c)
                .unwrap()
                .polygons
                .retain(|id| *id != pid);
        }
        n.polygons.remove(&pid);
        let mut after_edges = vec![central];
        after_edges.extend_from_slice(&new);
        let [top, bottom, a, b, d] = cells;
        let new_vertex_cells = vec![
            (newv[0], vec![a, d, b, top]),
            (newv[1], vec![a, d, b, bottom]),
        ];
        Ok(Proposal {
            mesh: n,
            cells,
            before_edges: old_all,
            after_edges,
            new_vertex_cells,
        })
    }

    fn write_diagnostics(&self, before: &Mesh, p: &Proposal, ih: bool) -> Result<()> {
        let dir = self
            .diagnostic_directory
            .clone()
            .unwrap_or_else(|| PathBuf::from("."));
        let prefix = if ih { "IH" } else { "HI" };
        let mut files = Vec::new();
        for (stage, m, edges) in [(0, before, &p.before_edges), (1, &p.mesh, &p.after_edges)] {
            for kind in ["cells", "edges"] {
                let name = format!("{prefix}_{:07}_{stage}.{kind}.vtk", self.count_events);
                let data = diagnostic_vtk(m, &p.cells, edges, kind == "cells");
                files.push((dir.join(name), data));
            }
        }
        // A collision in an after-event file must be caught before any writes.
        for (path, _) in &files {
            match path.try_exists() {
                Ok(true) => {
                    return Err(format!(
                        "Refusing to overwrite existing output file {}",
                        path.display()
                    ))
                }
                Err(e) => return Err(format!("Cannot inspect {}: {e}", path.display())),
                _ => {}
            }
        }
        for (path, data) in files {
            let mut f = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
                .map_err(|e| format!("Cannot create {}: {e}", path.display()))?;
            f.write_all(data.as_bytes())
                .and_then(|_| f.flush())
                .map_err(|e| format!("Cannot write {}: {e}", path.display()))?;
        }
        Ok(())
    }
}

fn diagnostic_vtk(m: &Mesh, cells: &[Id], edges: &[Id], faces: bool) -> String {
    let rows: Vec<Vec<Id>> = if faces {
        let mut seen = BTreeSet::new();
        cells
            .iter()
            .flat_map(|c| &m.cells[c].polygons)
            .filter(|p| seen.insert(**p))
            .map(|p| m.polygons[p].vertices.clone())
            .collect()
    } else {
        edges.iter().map(|e| m.edges[e].vertices.to_vec()).collect()
    };
    let mut vertices = Vec::new();
    let mut indices = BTreeMap::new();
    for v in rows.iter().flatten() {
        if !indices.contains_key(v) {
            indices.insert(*v, vertices.len());
            vertices.push(*v);
        }
    }
    let mut out = String::from("# vtk DataFile Version 2.0\npolydata\nASCII\nDATASET POLYDATA\n");
    writeln!(out, "POINTS {} double", vertices.len()).unwrap();
    for v in vertices {
        let p = m.vertices[&v].position;
        writeln!(out, "{:.17e} {:.17e} {:.17e}", p[0], p[1], p[2]).unwrap();
    }
    writeln!(
        out,
        "\n{} {} {}",
        if faces { "POLYGONS" } else { "LINES" },
        rows.len(),
        rows.iter().map(|r| r.len() + 1).sum::<usize>()
    )
    .unwrap();
    for row in rows {
        write!(out, "{}", row.len()).unwrap();
        for v in row {
            write!(out, " {}", indices[&v]).unwrap();
        }
        out.push('\n');
    }
    out
}
