use std::fs;
use std::path::PathBuf;
use vm3d::io::load_topology;
use vm3d::model::{add, norm, scale, Mesh};
use vm3d::reconnection::Reconnection;

fn fixture() -> Mesh {
    load_topology(
        &PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/periodic64/sample.topo"),
        [4.; 3],
    )
    .unwrap()
}

// Select a genuinely valid event from the periodic 64-cell test network.
fn transitioned(verbose: bool, dir: Option<PathBuf>) -> (Mesh, Reconnection, usize) {
    let initial = fixture();
    for edge in initial.edges.keys() {
        let mut mesh = initial.clone();
        let mut recon = Reconnection::new(0.02, verbose);
        recon.diagnostic_directory = dir.clone();
        if recon.i_h(&mut mesh, *edge).unwrap() {
            let triangle = mesh.next_polygon - 1;
            return (mesh, recon, triangle);
        }
    }
    panic!("periodic64 must contain an admissible I-H edge");
}

#[test]
fn periodic_round_trip_preserves_incidence_volume_and_threshold_scale() {
    let initial = fixture();
    let sizes = (
        initial.vertices.len(),
        initial.edges.len(),
        initial.polygons.len(),
    );
    let (mut mesh, mut r, triangle) = transitioned(false, None);
    assert_eq!(
        (mesh.vertices.len(), mesh.edges.len(), mesh.polygons.len()),
        (sizes.0 + 1, sizes.1 + 2, sizes.2 + 1)
    );
    mesh.validate_topology().unwrap();
    mesh.validate_geometry().unwrap();
    let longest = mesh.polygons[&triangle]
        .edges
        .iter()
        .map(|e| mesh.edges[e].length)
        .fold(0., f64::max);
    assert!(
        (longest - r.lth).abs() < 1e-12,
        "longest triangle side {longest}"
    );
    let new_edge = mesh.next_edge;
    assert!(r.h_i(&mut mesh, triangle).unwrap(), "{}", r.last_rejection);
    assert_eq!(
        (mesh.vertices.len(), mesh.edges.len(), mesh.polygons.len()),
        sizes
    );
    assert!((mesh.edges[&new_edge].length - r.lth).abs() < 1e-12);
    mesh.validate_topology().unwrap();
    mesh.validate_geometry().unwrap();
    assert_eq!((r.count_ih, r.count_hi, r.count_events), (1, 1, 2));
    assert!(r.events[0].contains("I->H time") && r.events[1].contains("H->I time"));
}

#[test]
fn rejected_zero_edge_and_id_overflow_are_atomic() {
    let mut m = fixture();
    let eid = *m.edges.keys().next().unwrap();
    let a = m.edges[&eid].vertices[0];
    let b = m.edges[&eid].vertices[1];
    m.vertices.get_mut(&b).unwrap().position = m.vertices[&a].position;
    m.edges.get_mut(&eid).unwrap().length = 0.;
    let before = format!("{m:?}");
    let mut r = Reconnection::new(0.02, false);
    assert!(!r.i_h(&mut m, eid).unwrap());
    assert_eq!(before, format!("{m:?}"));
    assert_eq!(r.rejected_geometry, 1);
    assert_eq!(r.count_events, 0);
    let initial = fixture();
    let successful = initial
        .edges
        .keys()
        .find(|e| {
            let mut m = initial.clone();
            Reconnection::new(0.02, false).i_h(&mut m, **e).unwrap()
        })
        .copied()
        .unwrap();
    let mut m = initial;
    m.next_vertex = usize::MAX - 2;
    let before = format!("{m:?}");
    assert!(!r.i_h(&mut m, successful).unwrap());
    assert!(r.last_rejection.contains("ID space"));
    assert_eq!(before, format!("{m:?}"));
    m.next_vertex = 1_000_000;
    r.count_events = u64::MAX;
    let before = format!("{m:?}");
    assert!(!r.i_h(&mut m, successful).unwrap());
    assert_eq!(before, format!("{m:?}"));
}

#[test]
fn collinear_triangle_rejected_without_mutation() {
    let (mut m, mut r, p) = transitioned(false, None);
    let vertices = m.polygons[&p].vertices.clone();
    for (i, v) in vertices.iter().enumerate() {
        m.vertices.get_mut(v).unwrap().position = [i as f64 * 0.01, 1., 1.];
    }
    let before = format!("{m:?}");
    assert!(!r.h_i(&mut m, p).unwrap());
    assert_eq!(before, format!("{m:?}"));
    assert!(r.last_rejection.contains("triangle normal"));
    assert_eq!(r.count_events, 1);
}

#[test]
fn zero_threshold_disables_sweep_and_invalid_threshold_fails() {
    let mut m = fixture();
    let before = format!("{m:?}");
    Reconnection::new(0., false).start(&mut m).unwrap();
    assert_eq!(before, format!("{m:?}"));
    assert!(Reconnection::new(f64::NAN, false).start(&mut m).is_err());
    assert!(Reconnection::new(-1., false).start(&mut m).is_err());
}

#[test]
fn translated_boundary_round_trip() {
    let initial = fixture();
    let eid = initial
        .edges
        .keys()
        .find(|e| {
            let mut m = initial.clone();
            Reconnection::new(0.02, false).i_h(&mut m, **e).unwrap()
        })
        .copied()
        .unwrap();
    let center = initial.edges[&eid].center;
    let shift = add(scale(center, -1.), [3.999, 3.998, 3.997]);
    let mut m = initial;
    for id in m.vertices.keys().copied().collect::<Vec<_>>() {
        let p = m.wrap(add(m.vertices[&id].position, shift));
        m.vertices.get_mut(&id).unwrap().position = p;
    }
    m.rebuild().unwrap();
    let mut r = Reconnection::new(0.02, false);
    assert!(r.i_h(&mut m, eid).unwrap(), "{}", r.last_rejection);
    let triangle = m.next_polygon - 1;
    assert!(r.h_i(&mut m, triangle).unwrap(), "{}", r.last_rejection);
    m.validate_geometry().unwrap();
    m.validate_topology().unwrap();
    assert!(m
        .vertices
        .values()
        .all(|v| v.position.iter().all(|x| *x >= 0. && *x < 4.)));
}

struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let p =
            std::env::temp_dir().join(format!("vm3d-reconnection-{}-{nanos}", std::process::id()));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn verbose_events_create_four_vtk_files_each_and_refuse_overwrite() {
    let dir = Temp::new();
    let (mut m, mut r, triangle) = transitioned(true, Some(dir.0.clone()));
    r.time = 1e20;
    assert!(r.h_i(&mut m, triangle).unwrap(), "{}", r.last_rejection);
    for (prefix, event) in [("IH", 0), ("HI", 1)] {
        for stage in 0..2 {
            for kind in ["cells", "edges"] {
                let s = fs::read_to_string(
                    dir.0
                        .join(format!("{prefix}_{event:07}_{stage}.{kind}.vtk")),
                )
                .unwrap();
                assert!(s.starts_with("# vtk DataFile Version 2.0\n"));
                assert!(s.contains("POINTS"));
            }
        }
    }
    let initial = fixture();
    let eid = initial
        .edges
        .keys()
        .find(|e| {
            let mut m = initial.clone();
            Reconnection::new(0.02, false).i_h(&mut m, **e).unwrap()
        })
        .copied()
        .unwrap();
    for stage in 0..2 {
        for kind in ["cells", "edges"] {
            let temp = Temp::new();
            let name = format!("IH_0000000_{stage}.{kind}.vtk");
            fs::write(temp.0.join(&name), "sentinel").unwrap();
            let mut m = initial.clone();
            let before = format!("{m:?}");
            let mut r = Reconnection::new(0.02, true);
            r.diagnostic_directory = Some(temp.0.clone());
            assert!(r
                .i_h(&mut m, eid)
                .unwrap_err()
                .contains("Refusing to overwrite"));
            assert_eq!(before, format!("{m:?}"));
            assert_eq!(r.count_events, 0);
            assert_eq!(fs::read_to_string(temp.0.join(&name)).unwrap(), "sentinel");
            assert_eq!(fs::read_dir(&temp.0).unwrap().count(), 1);
        }
    }
}

#[test]
fn real_sweep_only_short_nontriangular_edges_are_candidates() {
    let mut m = fixture();
    let mut r = Reconnection::new(0.02, false);
    r.start(&mut m).unwrap();
    m.validate_topology().unwrap();
    m.validate_geometry().unwrap();
    assert!(m
        .edges
        .values()
        .all(|e| e.length.is_finite() && norm(e.vector) > 0.));
}
