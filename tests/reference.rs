//! Numerical reference values for the periodic 64-cell test configuration.
use std::path::Path;
use vm3d::{io, mechanics, rng::Rng};

#[test]
fn matches_reference_forces_volumes_energy_and_noise() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/periodic64");
    let mut mesh = io::load_topology(&root.join("sample.topo"), [4.; 3]).unwrap();
    let e = mechanics::forces(&mut mesh, 10., 5.6).unwrap();
    let mut rng = Rng::new(20261028);
    let text = std::fs::read_to_string(root.join("reference.tsv")).unwrap();
    let mut counts = [0; 4];
    let close = |actual: f64, expected: f64| {
        assert!(
            (actual - expected).abs() < 2e-10 * expected.abs().max(1.),
            "actual {actual:.17e}, expected {expected:.17e}"
        );
    };
    for line in text.lines() {
        let t: Vec<_> = line.split_whitespace().collect();
        match t[0] {
            "energy" => {
                counts[0] += 1;
                close(e.volume, t[1].parse().unwrap());
                close(e.interface, t[2].parse().unwrap());
                close(e.total_volume, t[3].parse().unwrap());
            }
            "vertex" => {
                counts[1] += 1;
                let v = &mesh.vertices[&t[1].parse::<usize>().unwrap()];
                for k in 0..3 {
                    close(v.volume_force[k], t[2 + k].parse().unwrap());
                    close(v.interface_force[k], t[5 + k].parse().unwrap());
                }
            }
            "cell" => {
                counts[2] += 1;
                close(
                    mesh.cells[&t[1].parse::<usize>().unwrap()].volume,
                    t[2].parse().unwrap(),
                );
            }
            "normal" => {
                counts[3] += 1;
                close(rng.normal(), t[1].parse().unwrap());
            }
            _ => panic!("unexpected golden row"),
        }
    }
    assert_eq!(counts, [1, 445, 64, 12]);
}
