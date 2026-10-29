// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2021-2023 Tao Zhang; coauthor Jennifer Schwarz.
//! Analytical gradients including the full perimeter-weighted face-center Jacobian.
use crate::model::*;
use crate::Result;

#[derive(Clone, Copy, Debug, Default)]
pub struct Energy {
    pub volume: f64,
    pub interface: f64,
    pub total_volume: f64,
}
impl Energy {
    pub fn total(&self) -> f64 {
        self.volume + self.interface
    }
}

fn center_contribution(
    mesh: &Mesh,
    p: &Polygon,
    e: &Edge,
    force: Vec3,
    perimeter: f64,
) -> Result<[Vec3; 2]> {
    if !(e.length.is_finite() && e.length > 0. && perimeter.is_finite() && perimeter > 0.) {
        return Err(format!(
            "Degenerate center-force geometry at polygon {}",
            p.id
        ));
    }
    let offset = mesh.distance(sub(e.center, p.center));
    let derivative = dot(force, offset) / perimeter;
    let midpoint = scale(force, 0.5 * e.length / perimeter);
    let weight = scale(e.vector, derivative / e.length);
    let answer = [sub(midpoint, weight), add(midpoint, weight)];
    if answer.iter().flatten().any(|x| !x.is_finite()) {
        return Err("Nonfinite center force".into());
    }
    Ok(answer)
}

/// Refresh geometry and calculate both energy components and vertex forces.
/// Winding is retained through continuous dynamics so inverted cells remain detectable.
pub fn forces(mesh: &mut Mesh, kv: f64, s0: f64) -> Result<Energy> {
    if !kv.is_finite() || kv < 0. || !s0.is_finite() || s0 < 0. {
        return Err("Invalid mechanics parameters".into());
    }
    mesh.update_geometry()?;
    for v in mesh.vertices.values_mut() {
        v.volume_force = [0.; 3];
        v.interface_force = [0.; 3];
    }
    for p in mesh.polygons.values_mut() {
        p.tension = 0.;
    }
    let mut energy = Energy::default();
    for cell in mesh.cells.values_mut() {
        cell.pressure = -2. * kv * (cell.volume - 1.);
        energy.volume += kv * (cell.volume - 1.).powi(2);
        energy.total_volume += cell.volume;
        let area = cell
            .polygons
            .iter()
            .map(|p| mesh.polygons[p].area)
            .sum::<f64>();
        energy.interface += (area - s0).powi(2);
        for pid in &cell.polygons {
            mesh.polygons.get_mut(pid).unwrap().tension += 2. * (area - s0);
        }
    }
    // Evaluate the cell-pressure contribution in cell and face iteration order.
    for cell in mesh.cells.values() {
        for pid in &cell.polygons {
            let p = &mesh.polygons[pid];
            let mut center_force = [0.; 3];
            let direction = if cell.directions[pid] { 1. } else { -1. };
            for i in 0..p.vertices.len() {
                let a = p.vertices[i];
                let b = p.vertices[(i + 1) % p.vertices.len()];
                let av = mesh.distance(sub(mesh.vertices[&a].position, p.center));
                let bv = mesh.distance(sub(mesh.vertices[&b].position, p.center));
                let interface = scale(cross(av, bv), 0.5 * direction);
                let force = scale(interface, (1. / 3.) * cell.pressure);
                for vid in [a, b] {
                    let v = mesh.vertices.get_mut(&vid).unwrap();
                    v.volume_force = add(v.volume_force, force);
                }
                center_force = add(center_force, force);
            }
            let perimeter = p.edges.iter().map(|e| mesh.edges[e].length).sum();
            for eid in &p.edges {
                let e = &mesh.edges[eid];
                let contribution = center_contribution(mesh, p, e, center_force, perimeter)?;
                for (vid, f) in e.vertices.into_iter().zip(contribution) {
                    let v = mesh.vertices.get_mut(&vid).unwrap();
                    v.volume_force = add(v.volume_force, f);
                }
            }
        }
    }
    for p in mesh.polygons.values() {
        let mut center_force = [0.; 3];
        for eid in &p.edges {
            let e = &mesh.edges[eid];
            let a = mesh.distance(sub(mesh.vertices[&e.vertices[0]].position, p.center));
            let b = mesh.distance(sub(mesh.vertices[&e.vertices[1]].position, p.center));
            let vv = sub(b, a);
            let normal = cross(a, b);
            let magnitude = norm(normal);
            if !magnitude.is_finite() || magnitude <= 64. * f64::EPSILON * norm(a) * norm(b) {
                return Err(format!("Degenerate area triangle at polygon {}", p.id));
            }
            let normal = scale(normal, 1. / magnitude);
            let fvv = scale(cross(normal, vv), p.tension);
            let fa = scale(cross(normal, a), p.tension);
            let fb = scale(cross(b, normal), p.tension);
            for (vid, f) in e
                .vertices
                .into_iter()
                .zip([scale(add(fa, fvv), 0.5), scale(add(fb, fvv), 0.5)])
            {
                let v = mesh.vertices.get_mut(&vid).unwrap();
                v.interface_force = add(v.interface_force, f);
            }
            center_force = add(center_force, scale(add(fa, fb), 0.5));
        }
        let perimeter = p.edges.iter().map(|e| mesh.edges[e].length).sum();
        for eid in &p.edges {
            let e = &mesh.edges[eid];
            let contribution = center_contribution(mesh, p, e, center_force, perimeter)?;
            for (vid, f) in e.vertices.into_iter().zip(contribution) {
                let v = mesh.vertices.get_mut(&vid).unwrap();
                v.interface_force = add(v.interface_force, f);
            }
        }
    }
    if !energy.total().is_finite()
        || mesh.vertices.values().any(|v| {
            v.volume_force
                .iter()
                .chain(v.interface_force.iter())
                .any(|x| !x.is_finite())
        })
    {
        return Err("Nonfinite energy or force".into());
    }
    Ok(energy)
}
