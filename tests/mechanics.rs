// SPDX-License-Identifier: GPL-3.0-or-later
//! Independent energy-gradient and geometry tests.
use std::collections::BTreeMap;
use std::path::Path;
use vm3d::mechanics::forces;
use vm3d::model::*;

fn cube(warp: f64, shift: f64, high_ids: bool) -> Mesh {
    let mut mesh = Mesh::new([10.; 3]);
    let mut positions = [
        [0., 0., 0.],
        [1., 0., 0.],
        [1., 1., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [1., 0., 1.],
        [1., 1., 1.],
        [0., 1., 1.],
    ];
    positions[6][2] += warp;
    positions[6][0] += 0.37 * warp;
    positions[6][1] -= 0.23 * warp;
    positions[5][1] += 0.18 * warp;
    positions[5][2] -= 0.31 * warp;
    positions[7][0] -= 0.24 * warp;
    positions[7][2] += 0.11 * warp;
    let ids = if high_ids {
        [1, 20000, 42950, 87296, 90000, 90001, 90002, 90003]
    } else {
        [0, 1, 2, 3, 4, 5, 6, 7]
    };
    for (i, p) in positions.into_iter().enumerate() {
        mesh.vertices.insert(
            ids[i],
            Vertex {
                id: ids[i],
                position: mesh.wrap(add(p, [shift; 3])),
                ..Default::default()
            },
        );
    }
    let faces = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ];
    let mut pairs = BTreeMap::new();
    for (f, vs) in faces.iter().enumerate() {
        let mut p = Polygon {
            id: f,
            ..Default::default()
        };
        for i in 0..4 {
            let a = ids[vs[i]];
            let b = ids[vs[(i + 1) % 4]];
            let key = (a.min(b), a.max(b));
            let next = mesh.edges.len();
            let eid = *pairs.entry(key).or_insert_with(|| {
                mesh.edges.insert(
                    next,
                    Edge {
                        id: next,
                        vertices: [a, b],
                        ..Default::default()
                    },
                );
                next
            });
            p.edges.push(eid);
        }
        mesh.polygons.insert(f, p);
    }
    mesh.cells.insert(
        0,
        Cell {
            id: 0,
            polygons: (0..6).collect(),
            ..Default::default()
        },
    );
    mesh.rebuild().unwrap();
    mesh
}

fn component_energy(mesh: &mut Mesh, component: usize) -> f64 {
    mesh.update_geometry().unwrap();
    mesh.cells
        .values()
        .map(|c| {
            if component == 0 {
                10. * (c.volume - 1.).powi(2)
            } else {
                (c.polygons
                    .iter()
                    .map(|p| mesh.polygons[p].area)
                    .sum::<f64>()
                    - 5.6)
                    .powi(2)
            }
        })
        .sum()
}

fn gradient(mesh: &mut Mesh, ids: &[Id], h: f64, tolerance: f64) {
    forces(mesh, 10., 5.6).unwrap();
    let expected: BTreeMap<_, _> = mesh
        .vertices
        .iter()
        .map(|(id, v)| (*id, [v.volume_force, v.interface_force]))
        .collect();
    for mode in 0..2 {
        let net = mesh.vertices.values().fold([0.; 3], |sum, v| {
            add(
                sum,
                if mode == 0 {
                    v.volume_force
                } else {
                    v.interface_force
                },
            )
        });
        assert!(norm(net) < 1e-9, "Net force {net:?}");
        for id in ids {
            for k in 0..3 {
                let old = mesh.vertices[id].position[k];
                mesh.vertices.get_mut(id).unwrap().position[k] = old + h;
                let plus = component_energy(mesh, mode);
                mesh.vertices.get_mut(id).unwrap().position[k] = old - h;
                let minus = component_energy(mesh, mode);
                mesh.vertices.get_mut(id).unwrap().position[k] = old;
                let numerical = -(plus - minus) / (2. * h);
                let analytical = expected[id][mode][k];
                assert!((analytical-numerical).abs()<tolerance,"gradient mode {mode}, vertex {id}, coordinate {k}: analytic {analytical}, finite difference {numerical}");
            }
        }
    }
}

#[test]
fn planar_warped_wrapped_and_high_id_energy_gradients() {
    for (warp, shift, high) in [
        (0., 2., false),
        (0.7, 2., false),
        (0.7, 9.5, false),
        (0.7, 2., true),
    ] {
        let mut mesh = cube(warp, shift, high);
        let ids = mesh.vertices.keys().copied().collect::<Vec<_>>();
        gradient(&mut mesh, &ids, 1e-6, 2e-7);
    }
    let mut regular = cube(0.7, 2., false);
    let mut high = cube(0.7, 2., true);
    assert!(
        (forces(&mut regular, 10., 5.6).unwrap().total()
            - forces(&mut high, 10., 5.6).unwrap().total())
        .abs()
            < 1e-12
    );
}

fn bulk() -> Mesh {
    vm3d::io::load_topology(
        &Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/periodic64/sample.topo"),
        [4.; 3],
    )
    .unwrap()
}

#[test]
fn periodic_bulk_gradients_volume_and_translation() {
    let mut mesh = bulk();
    mesh.validate_topology().unwrap();
    mesh.validate_geometry().unwrap();
    let ids = mesh
        .vertices
        .keys()
        .copied()
        .step_by(31)
        .collect::<Vec<_>>();
    gradient(&mut mesh, &ids, 1e-7, 2e-6);
    let pi = std::f64::consts::PI;
    for v in mesh.vertices.values_mut() {
        let [x, y, z] = scale(v.position, pi / 2.);
        v.position = add(
            v.position,
            [
                0.04 * y.sin() * z.sin(),
                0.03 * z.sin() * x.sin(),
                0.035 * x.sin() * y.sin(),
            ],
        );
    }
    gradient(&mut mesh, &ids, 1e-7, 2e-6);
    let before = forces(&mut mesh, 10., 5.6).unwrap();
    mesh.validate_geometry().unwrap();
    assert!((before.total_volume - 64.).abs() < 1e-10);
    for id in mesh.vertices.keys().copied().collect::<Vec<_>>() {
        let p = mesh.wrap(add(mesh.vertices[&id].position, [2.4, 1.7, 3.1]));
        mesh.vertices.get_mut(&id).unwrap().position = p;
    }
    let after = forces(&mut mesh, 10., 5.6).unwrap();
    assert!((before.volume - after.volume).abs() < 1e-10);
    assert!((before.interface - after.interface).abs() < 1e-10);
    mesh.validate_geometry().unwrap();
    gradient(&mut mesh, &ids, 1e-7, 2e-6);
}

#[test]
fn degenerate_face_fails_before_nan_normalization() {
    let mut mesh = cube(0., 2., false);
    for v in mesh.vertices.values_mut() {
        v.position = [v.id as f64, 0., 0.];
    }
    assert!(forces(&mut mesh, 10., 5.6).is_err());
    assert!(mesh
        .vertices
        .values()
        .all(|v| v.interface_force.iter().all(|x| x.is_finite())));
}

#[test]
fn topology_rejects_duplicate_edges_and_disconnected_cycles() {
    let mut mesh = bulk();
    let id = mesh.next_edge;
    let mut e = mesh.edges.values().next().unwrap().clone();
    e.id = id;
    mesh.edges.insert(id, e);
    assert!(mesh.validate_topology().is_err());
    let mut mesh = bulk();
    let p = mesh.polygons.values_mut().next().unwrap();
    p.edges.push(p.edges[0]);
    assert!(mesh.validate_topology().is_err());
    assert!(mesh.rebuild().is_err());
    // Two disjoint four-edge loops satisfy degree-two local incidence but are
    // not one face. This catches the subtler disconnected-cycle failure.
    let mut mesh = cube(0., 2., false);
    let opposite = mesh.polygons[&1].edges.clone();
    mesh.polygons.get_mut(&0).unwrap().edges.extend(opposite);
    assert!(mesh
        .rebuild()
        .unwrap_err()
        .contains("Disconnected polygon cycles"));
}

#[test]
fn geometry_refresh_keeps_winding_to_detect_inversion() {
    let mut mesh = cube(0., 2., false);
    for v in mesh.vertices.values_mut() {
        v.position[2] = 5. - v.position[2];
    }
    mesh.update_geometry().unwrap();
    assert!(
        mesh.cells[&0].volume < 0.,
        "Continuous dynamics must not silently reorient inverted cells"
    );
    assert!(mesh.validate_geometry().is_err());
}

#[test]
fn periodic_arithmetic_preserves_half_box_ties_and_large_wraps() {
    let mesh = Mesh::new([4.; 3]);
    assert_eq!(mesh.distance([2., -2., 14.]), [2., -2., 2.]);
    assert_eq!(mesh.wrap([-8.25, 12.25, 4.]), [3.75, 0.25, 0.]);
}
