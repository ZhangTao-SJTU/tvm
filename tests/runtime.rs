// SPDX-License-Identifier: GPL-3.0-or-later
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use vm3d::io::{load_config, load_topology, topology_text, vtk_text};
use vm3d::model::{Mesh, Vertex};
use vm3d::rng::Rng;
use vm3d::simulation::{advance_positions, run};

static COUNTER: AtomicU64 = AtomicU64::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        loop {
            let path = std::env::temp_dir().join(format!(
                "vm3d-runtime-{}-{}",
                std::process::id(),
                COUNTER.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&path) {
                Ok(()) => return Self(path),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => panic!("Cannot create test directory: {e}"),
            }
        }
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
    fn write(&self, name: &str, text: &str) {
        fs::write(self.path(name), text).unwrap();
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        // Never remove an arbitrary path: resolve both sides and restrict deletion
        // to this generated test-directory name directly inside the temp root.
        if let (Ok(path), Ok(temp)) = (self.0.canonicalize(), std::env::temp_dir().canonicalize()) {
            if path.parent() == Some(temp.as_path())
                && path
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with("vm3d-runtime-")
            {
                let _ = fs::remove_dir_all(path);
            }
        }
    }
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/periodic64/sample.topo")
}
fn config(start: f64, end: f64, dt: f64) -> String {
    format!("time {start} {end} {dt}\ndtr 0.01\ndump vtk {dt}\nlog {dt}\ns0 5.6\nLth 0\nT 0\nkv 10\nbox 4 4 4 p p p\nseed 42\nvalidate_topology on\n")
}
fn prepare(conf: &str) -> Scratch {
    let directory = Scratch::new();
    directory.write("conf", conf);
    fs::copy(fixture(), directory.path("sample.topo")).unwrap();
    directory
}
fn times(path: &Path) -> Vec<f64> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|s| s.starts_with("time "))
        .map(|s| s.split_whitespace().nth(1).unwrap().parse().unwrap())
        .collect()
}

#[test]
fn config_rejects_malformed_nonfinite_and_unsafe_ranges() {
    let directory = Scratch::new();
    let base = config(0., 0.01, 0.001);
    for (old, new, expected) in [
        ("time 0 0.01 0.001", "time 0 0.01 0", "Require dt>0"),
        (
            "time 0 0.01 0.001",
            "time 1e20 1e20 0.001",
            "representable precision",
        ),
        ("T 0\n", "T -1\n", "nonnegative"),
        ("box 4 4 4", "box 0 4 4", "positive"),
        ("p p p", "p p f", "periodic bulk"),
        ("seed 42", "seed 12junk", "Invalid seed"),
        ("seed 42", "seed -1", "unsigned 64-bit"),
        ("s0 5.6", "s0 NaN", "finite number"),
        ("dtr 0.01", "dtr 0", "dtr>=dt"),
        ("kv 10", "unknown 10", "Unknown parameter"),
        ("kv 10", "kv 10\nkv 10", "Duplicate parameter"),
        ("dump vtk 0.001", "dump vtk 1e-300", "periods must be >= dt"),
        ("log 0.001", "log 1e-300", "periods must be >= dt"),
        (
            "box 4 4 4",
            "box 1e300 1e300 1e300",
            "volume must be finite",
        ),
        ("Lth 0", "Lth 0 stray", "optional verbose"),
    ] {
        assert!(base.contains(old), "{old}");
        directory.write("conf", &base.replace(old, new));
        let error = load_config(&directory.path("conf")).unwrap_err();
        assert!(error.contains(expected), "{new}: {error}");
    }
    directory.write("conf", &base.replace("kv 10\n", ""));
    assert!(load_config(&directory.path("conf"))
        .unwrap_err()
        .contains("Missing configuration parameter: kv"));
    directory.write(
        "conf",
        &format!(
            "# comment\n{}",
            base.replace("dtr 0.01\n", "")
                .replace("kv 10", "kv 10 # comment")
        ),
    );
    let parsed = load_config(&directory.path("conf")).unwrap();
    assert_eq!(parsed.dtr, 10. * parsed.dt);
}

#[test]
fn topology_rejects_corrupt_references_and_section_counts() {
    let directory = Scratch::new();
    let original = fs::read_to_string(fixture()).unwrap();
    let lines: Vec<_> = original.lines().map(str::to_owned).collect();
    let edges = lines.iter().position(|s| s.starts_with("edges ")).unwrap();
    for kind in 0..5 {
        let mut bad = lines.clone();
        match kind {
            0 => {
                let mut fields: Vec<_> = bad[1].split_whitespace().map(str::to_owned).collect();
                fields[1] = "nan".into();
                bad[1] = fields.join(" ");
            }
            1 => {
                let id = bad[1].split_whitespace().next().unwrap().to_owned();
                let mut fields: Vec<_> = bad[2].split_whitespace().map(str::to_owned).collect();
                fields[0] = id;
                bad[2] = fields.join(" ");
            }
            2 => {
                let mut fields: Vec<_> = bad[edges + 1]
                    .split_whitespace()
                    .map(str::to_owned)
                    .collect();
                fields[2] = "99999999".into();
                bad[edges + 1] = fields.join(" ");
            }
            3 => {
                let mut fields: Vec<_> = bad[edges + 1]
                    .split_whitespace()
                    .map(str::to_owned)
                    .collect();
                fields[2] = fields[1].clone();
                bad[edges + 1] = fields.join(" ");
            }
            _ => {
                while bad.last().is_some_and(|x| x.trim().is_empty()) {
                    bad.pop();
                }
                bad.pop();
            }
        }
        directory.write("sample.topo", &bad.join("\n"));
        let error = load_topology(&directory.path("sample.topo"), [4.; 3]).unwrap_err();
        let expected = [
            "finite number",
            "Duplicate vertex ID",
            "Unknown referenced ID",
            "Self edge",
            "Incomplete topology",
        ][kind];
        assert!(error.contains(expected), "{kind}: {error}");
    }
}

#[test]
fn input_compatible_snapshot_roundtrip_and_vtk_indices() {
    let mesh = load_topology(&fixture(), [4.; 3]).unwrap();
    let directory = Scratch::new();
    directory.write("snapshot", &topology_text(&mesh));
    let restored = load_topology(&directory.path("snapshot"), [4.; 3]).unwrap();
    assert_eq!(topology_text(&mesh), topology_text(&restored));
    let vtk = vtk_text(&mesh);
    let words: Vec<_> = vtk.split_whitespace().collect();
    let start = words.iter().position(|x| *x == "POINTS").unwrap();
    let points: usize = words[start + 1].parse().unwrap();
    assert_eq!(points, mesh.vertices.len());
    let mut cursor = start + 3 + 3 * points;
    assert_eq!(words[cursor], "POLYGONS");
    let count: usize = words[cursor + 1].parse().unwrap();
    let total: usize = words[cursor + 2].parse().unwrap();
    cursor += 3;
    let begin = cursor;
    assert!(count > 0);
    for _ in 0..count {
        let n: usize = words[cursor].parse().unwrap();
        cursor += 1;
        assert!(n >= 3);
        for _ in 0..n {
            assert!(words[cursor].parse::<usize>().unwrap() < points);
            cursor += 1;
        }
    }
    assert_eq!(cursor - begin, total);
    assert_eq!(cursor, words.len());
}

#[test]
fn partial_final_step_negative_times_and_zero_duration() {
    for (start, end, expected) in [
        (0., 0.013, vec![0., 0.005, 0.01, 0.013]),
        (-0.02, -0.007, vec![-0.007]),
        (-0.01, 0.007, vec![0., 0.005, 0.007]),
        (2., 2., vec![2.]),
    ] {
        let directory = prepare(&config(start, end, 0.005));
        run(&directory.0).unwrap();
        let actual = times(&directory.path("topo.txt"));
        assert_eq!(actual.len(), expected.len(), "{actual:?}");
        for (a, b) in actual.iter().zip(expected) {
            assert!((a - b).abs() < 1e-15);
        }
        assert_eq!(actual.last().copied(), Some(end));
        for name in ["cellVolume.txt", "cellCenter.txt", "cellShapeIndex.txt"] {
            assert_eq!(times(&directory.path(name)), actual);
        }
        let summary = fs::read_to_string(directory.path("run_summary.txt")).unwrap();
        assert!(summary.starts_with("completed 1\n"));
        let metadata = fs::read_to_string(directory.path("run_metadata.txt")).unwrap();
        assert!(metadata.contains(&format!("implementation tvm {}", env!("CARGO_PKG_VERSION"))));
        assert!(metadata.contains("repository https://github.com/ZhangTao-SJTU/tvm"));
    }
}

#[test]
fn same_seed_repeats_noise_and_zero_temperature_consumes_nothing() {
    let mut mesh = Mesh::new([100.; 3]);
    mesh.vertices.insert(
        0,
        Vertex {
            id: 0,
            position: [50.; 3],
            velocity: [0.2, 0., 0.],
            ..Default::default()
        },
    );
    let mut repeat = mesh.clone();
    let mut other = mesh.clone();
    let (mut a, mut b, mut c) = (Rng::new(42), Rng::new(42), Rng::new(43));
    for _ in 0..7 {
        advance_positions(&mut mesh, 0.01, 0.1, &mut a).unwrap();
        advance_positions(&mut repeat, 0.01, 0.1, &mut b).unwrap();
        advance_positions(&mut other, 0.01, 0.1, &mut c).unwrap();
    }
    assert_eq!(mesh.vertices[&0].position, repeat.vertices[&0].position);
    assert_ne!(mesh.vertices[&0].position, other.vertices[&0].position);
    let saved = a.state();
    advance_positions(&mut mesh, 0.01, 0., &mut a).unwrap();
    assert_eq!(a.state(), saved);
}

#[test]
fn full_run_reproducible_and_existing_outputs_preserved() {
    let conf = config(0., 0.002, 0.001).replace("T 0\n", "T 0.0001\n");
    let a = prepare(&conf);
    let b = prepare(&conf);
    run(&a.0).unwrap();
    run(&b.0).unwrap();
    for name in [
        "topo.txt",
        "cellVolume.txt",
        "run_summary.txt",
        "rng_state.txt",
        "reconnections.txt",
    ] {
        assert_eq!(
            fs::read(a.path(name)).unwrap(),
            fs::read(b.path(name)).unwrap()
        );
    }
    let before = fs::read(a.path("topo.txt")).unwrap();
    assert!(run(&a.0).unwrap_err().contains("already contains a run"));
    assert_eq!(before, fs::read(a.path("topo.txt")).unwrap());
}

#[test]
fn output_collisions_and_invalid_input_return_errors() {
    let directory = prepare(&config(0., 0., 0.001));
    directory.write("0000000.sample.vtk", "sentinel");
    assert!(run(&directory.0)
        .unwrap_err()
        .contains("existing files are never overwritten"));
    assert_eq!(
        fs::read_to_string(directory.path("0000000.sample.vtk")).unwrap(),
        "sentinel"
    );
    assert!(!directory.path("run_summary.txt").exists());
    let bad = prepare(&config(0., 0., 0.001).replace("kv 10", "kv nan"));
    assert!(run(&bad.0).is_err());
    assert!(!bad.path("run_metadata.txt").exists());
}

#[test]
fn binary_accepts_directory_and_rejects_unknown_arguments() {
    let directory = prepare(&config(0., 0., 0.001));
    let binary = env!("CARGO_BIN_EXE_tvm");
    let output = std::process::Command::new(binary)
        .arg("--dir")
        .arg(&directory.0)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let bad = std::process::Command::new(binary)
        .arg("--unknown")
        .output()
        .unwrap();
    assert_eq!(bad.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&bad.stderr).contains("Usage:"));
}

#[test]
fn extreme_absolute_time_and_tiny_duration_reach_exact_end() {
    for (start, end, dt) in [(1e12, 1e12 + 0.001, 0.0002), (0., 1e-16, 1.)] {
        let directory = prepare(&config(start, end, dt).replace("dtr 0.01", "dtr 10"));
        run(&directory.0).unwrap();
        let frames = times(&directory.path("topo.txt"));
        assert!(frames.len() >= 2, "{frames:?}");
        assert_eq!(frames[0], start);
        assert_eq!(frames.last().copied(), Some(end));
        assert!(frames.windows(2).all(|w| w[1] > w[0]));
    }
}

#[test]
fn unstable_geometry_stops_without_completion_or_nan_output() {
    let directory = prepare(&config(0., 0.5, 0.5).replace("dtr 0.01", "dtr 1"));
    let error = run(&directory.0).unwrap_err().to_lowercase();
    assert!(
        error.contains("geometry") || error.contains("volume"),
        "{error}"
    );
    assert!(!directory.path("run_summary.txt").exists());
    let shapes = fs::read_to_string(directory.path("cellShapeIndex.txt")).unwrap();
    assert!(!shapes.to_lowercase().contains("nan"));
}
