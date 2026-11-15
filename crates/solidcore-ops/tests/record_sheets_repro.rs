//! Reproductions for two limits met when records of freeform sheets are decoded: a record of
//! separate sheets, and the cost of healing a large shell. Kernel only: records made with
//! `solidcore-record` types, read with `solidcore-ops`; no scan data, no application code.
//!
//! Run: `cargo test -p solidcore-ops --test record_sheets_repro -- --ignored --nocapture`
#![allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::cast_possible_truncation
)]

use solidcore_ops::{Body, stitch::stitch};
use solidcore_record::brep::{
    Brep, Curve, Curve2, Edge, EdgeUseRef, Face, Loop, OrientedEdge, OrientedFace, Region, Shell,
    Surface, VERSION, Vertex,
};
use std::collections::BTreeMap;
use std::time::Instant;

// ---------------------------------------------------------------------------------------------
// A closed network of bicubic Bezier patches (a cube sphere), edges shared.
// ---------------------------------------------------------------------------------------------

type Net = Vec<Vec<[f64; 3]>>;

/// `6 n^2` bicubic patches, control points on a sphere of `radius`, and their corner ids.
fn cube_sphere(n: usize, radius: f64) -> (Vec<Net>, Vec<[u32; 4]>) {
    let steps = 3 * n;
    let point = |lat: [usize; 3]| -> [f64; 3] {
        let c = lat.map(|a| 2.0 * a as f64 / steps as f64 - 1.0);
        let len = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        c.map(|x| x / len * radius)
    };
    let mut ids: BTreeMap<[usize; 3], u32> = BTreeMap::new();
    let mut patches = Vec::new();
    let mut corners = Vec::new();
    for axis in 0..3 {
        for side in 0..2 {
            let (p, q) = ((axis + 1) % 3, (axis + 2) % 3);
            let (au, av) = if side == 1 { (p, q) } else { (q, p) };
            for i in 0..n {
                for j in 0..n {
                    let lattice = |m: usize, k: usize| -> [usize; 3] {
                        let mut l = [0; 3];
                        l[axis] = side * steps;
                        l[au] = 3 * i + m;
                        l[av] = 3 * j + k;
                        l
                    };
                    let net: Net = (0..4)
                        .map(|m| (0..4).map(|k| point(lattice(m, k))).collect())
                        .collect();
                    let mut id = |l: [usize; 3]| {
                        let next = ids.len() as u32;
                        *ids.entry(l).or_insert(next)
                    };
                    corners.push([
                        id(lattice(0, 0)),
                        id(lattice(3, 0)),
                        id(lattice(3, 3)),
                        id(lattice(0, 3)),
                    ]);
                    patches.push(net);
                }
            }
        }
    }
    (patches, corners)
}

/// The network as a record: one `Surface::Spline` face per patch, edges and vertices shared,
/// one closed shell, in the canonical order of first use.
fn record_of(patches: &[Net], corners: &[[u32; 4]]) -> Brep {
    let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    let uv = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let mut vertex_of: BTreeMap<u32, u32> = BTreeMap::new();
    let mut vertices: Vec<Vertex> = Vec::new();
    let mut edge_of: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    let mut edges: Vec<Edge> = Vec::new();
    let mut faces: Vec<Face> = Vec::new();
    for (f, (net, c)) in patches.iter().zip(corners).enumerate() {
        let side = |s: usize| -> Vec<[f64; 3]> {
            (0..4)
                .map(|k| match s {
                    0 => net[k][0],
                    1 => net[3][k],
                    2 => net[3 - k][3],
                    _ => net[0][3 - k],
                })
                .collect()
        };
        let mut loop_edges = Vec::new();
        for s in 0..4 {
            let (a, b) = (c[s], c[(s + 1) % 4]);
            let (lo, hi, forward) = (a.min(b), a.max(b), a < b);
            let eid = *edge_of.entry((lo, hi)).or_insert_with(|| {
                let mut pts = side(s);
                if !forward {
                    pts.reverse();
                }
                let mut ids = [0u32; 2];
                for (slot, v) in ids.iter_mut().zip([lo, hi]) {
                    let next = vertices.len() as u32;
                    *slot = *vertex_of.entry(v).or_insert_with(|| {
                        vertices.push(Vertex {
                            point: if v == lo { pts[0] } else { pts[3] },
                            tolerance: 1e-6,
                        });
                        next
                    });
                }
                edges.push(Edge {
                    curve: Curve::Spline {
                        knots: knots.clone(),
                        points: pts,
                    },
                    nurbs: None,
                    interval: [0.0, 1.0],
                    vertices: ids,
                    tolerance: 1e-6,
                    uses: Vec::new(),
                });
                (edges.len() - 1) as u32
            });
            edges[eid as usize].uses.push(EdgeUseRef {
                face: f as u32,
                r#loop: 0,
                position: s as u32,
                upper: !forward,
            });
            let (p0, p1) = (uv[s], uv[(s + 1) % 4]);
            let (start, end) = if forward { (p0, p1) } else { (p1, p0) };
            let pcurve = Curve2::Spline {
                knots: vec![0.0, 0.0, 1.0, 1.0],
                points: vec![start, end],
            };
            loop_edges.push(OrientedEdge {
                edge: eid,
                forward,
                pcurve,
            });
        }
        faces.push(Face {
            surface: Surface::Spline {
                u_knots: knots.clone(),
                v_knots: knots.clone(),
                points: net.clone(),
            },
            nurbs: None,
            domain: [[0.0, 1.0], [0.0, 1.0]],
            loops: vec![Loop::Edges { edges: loop_edges }],
        });
    }
    let oriented = (0..faces.len() as u32)
        .map(|face| OrientedFace {
            face,
            forward: true,
        })
        .collect();
    Brep {
        version: VERSION,
        tolerance: 1e-6,
        regions: vec![
            Region {
                shells: vec![0],
                void: false,
            },
            Region {
                shells: vec![0],
                void: false,
            },
        ],
        shells: vec![Shell::Faces {
            faces: oriented,
            wires: vec![],
        }],
        faces,
        edges,
        vertices,
    }
}

fn seconds_to_decode(n: usize, trusted: bool) -> (usize, f64) {
    let (patches, corners) = cube_sphere(n, 50.0);
    let record = record_of(&patches, &corners);
    record.validate_built().unwrap();
    let faces = record.faces.len();
    let start = Instant::now();
    Body::from_record(&record, trusted).unwrap();
    (faces, start.elapsed().as_secs_f64())
}

/// A valid closed record of bicubic patches is healed in time near linear in its face count.
/// Four times the faces would take four times as long; the factor is printed for each pair.
#[test]
#[ignore = "reproduction: from_record(.., false) grows faster than linearly with the face count"]
fn healer_time_is_close_to_linear_in_the_face_count() {
    let mut times = Vec::new();
    for n in [16, 32] {
        let (faces, healed) = seconds_to_decode(n, false);
        let (_, trusted) = seconds_to_decode(n, true);
        println!(
            "n={n} faces={faces} from_record(false)={healed:.2}s from_record(true)={trusted:.3}s"
        );
        times.push(healed);
    }
    let factor = times[1] / times[0];
    println!("4x the faces took x{factor:.1} as long (linear: x4)");
    assert!(factor <= 6.0, "x{factor:.1} for 4x the faces");
}

// ---------------------------------------------------------------------------------------------
// Separate sheets: six bilinear faces of a block, each with its own edges and vertices.
// ---------------------------------------------------------------------------------------------

fn unsewn_block_record(size: [f64; 3]) -> Brep {
    let p = |x: f64, y: f64, z: f64| [x * size[0], y * size[1], z * size[2]];
    // each face as its four corners (u, v) = (0,0) (1,0) (1,1) (0,1)
    let quads: [[[f64; 3]; 4]; 6] = [
        [p(0., 0., 0.), p(0., 1., 0.), p(1., 1., 0.), p(1., 0., 0.)],
        [p(0., 0., 1.), p(1., 0., 1.), p(1., 1., 1.), p(0., 1., 1.)],
        [p(0., 0., 0.), p(1., 0., 0.), p(1., 0., 1.), p(0., 0., 1.)],
        [p(0., 1., 0.), p(0., 1., 1.), p(1., 1., 1.), p(1., 1., 0.)],
        [p(0., 0., 0.), p(0., 0., 1.), p(0., 1., 1.), p(0., 1., 0.)],
        [p(1., 0., 0.), p(1., 1., 0.), p(1., 1., 1.), p(1., 0., 1.)],
    ];
    let uv = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let line = vec![0.0, 0.0, 1.0, 1.0];
    let mut faces = Vec::new();
    let mut edges = Vec::new();
    let mut vertices = Vec::new();
    for (f, q) in quads.iter().enumerate() {
        for s in 0..4 {
            vertices.push(Vertex {
                point: q[s],
                tolerance: 1e-6,
            });
            edges.push(Edge {
                curve: Curve::Spline {
                    knots: line.clone(),
                    points: vec![q[s], q[(s + 1) % 4]],
                },
                nurbs: None,
                interval: [0.0, 1.0],
                vertices: [(4 * f + s) as u32, (4 * f + (s + 1) % 4) as u32],
                tolerance: 1e-6,
                uses: vec![EdgeUseRef {
                    face: f as u32,
                    r#loop: 0,
                    position: s as u32,
                    upper: false,
                }],
            });
        }
        let loop_edges = (0..4)
            .map(|s| OrientedEdge {
                edge: (4 * f + s) as u32,
                forward: true,
                pcurve: Curve2::Spline {
                    knots: line.clone(),
                    points: vec![uv[s], uv[(s + 1) % 4]],
                },
            })
            .collect();
        faces.push(Face {
            surface: Surface::Spline {
                u_knots: line.clone(),
                v_knots: line.clone(),
                points: vec![vec![q[0], q[3]], vec![q[1], q[2]]],
            },
            nurbs: None,
            domain: [[0.0, 1.0], [0.0, 1.0]],
            loops: vec![Loop::Edges { edges: loop_edges }],
        });
    }
    let shells = (0..6u32)
        .map(|face| Shell::Faces {
            faces: vec![OrientedFace {
                face,
                forward: true,
            }],
            wires: vec![],
        })
        .collect();
    Brep {
        version: VERSION,
        tolerance: 1e-6,
        regions: vec![Region {
            shells: (0..6).collect(),
            void: false,
        }],
        shells,
        faces,
        edges,
        vertices,
    }
}

/// A record of six separate sheets that share their corners exactly (the coordinates are the
/// same numbers) is a valid record. Decoding it, and sewing the result, should give the block.
#[test]
#[ignore = "reproduction: separate sheets of one record"]
fn unsewn_sheets_of_a_block_sew_into_a_solid() {
    let size = [40.0, 30.0, 20.0];
    let record = unsewn_block_record(size);
    println!(
        "validate: {:?}  validate_built: {:?}",
        record.validate(),
        record.validate_built()
    );
    for trusted in [false, true] {
        match Body::from_record(&record, trusted) {
            Err(e) => println!("from_record(trusted={trusted}): Err({e:?})"),
            Ok((body, _)) => {
                println!("from_record(trusted={trusted}): Ok");
                for tolerance in [1e-6, 1e-2] {
                    let copy = body.copy().unwrap();
                    match stitch(copy, tolerance) {
                        Ok((sewn, closed)) => {
                            let volume = sewn.volume(1e-6);
                            println!("  stitch({tolerance:e}): closed={closed} volume={volume:?}");
                        }
                        Err(e) => println!("  stitch({tolerance:e}): Err({e:?})"),
                    }
                }
            }
        }
    }
    let (body, _) = Body::from_record(&record, false).expect("a record of separate sheets decodes");
    let (sewn, closed) = stitch(body, 1e-6).unwrap();
    assert!(closed);
    let volume = sewn.volume(1e-6).unwrap();
    assert!(
        (volume - size[0] * size[1] * size[2]).abs() < 1e-3,
        "volume {volume}"
    );
}

/// The same block with every face a bicubic Bezier patch (flat, but with curved-edge machinery:
/// degree 3 side curves), the top domed by lifting its four inner control points.
fn unsewn_bicubic_block_record(size: [f64; 3], dome: f64) -> Brep {
    let p = |x: f64, y: f64, z: f64| [x * size[0], y * size[1], z * size[2]];
    let quads: [[[f64; 3]; 4]; 6] = [
        [p(0., 0., 0.), p(0., 1., 0.), p(1., 1., 0.), p(1., 0., 0.)],
        [p(0., 0., 1.), p(1., 0., 1.), p(1., 1., 1.), p(0., 1., 1.)],
        [p(0., 0., 0.), p(1., 0., 0.), p(1., 0., 1.), p(0., 0., 1.)],
        [p(0., 1., 0.), p(0., 1., 1.), p(1., 1., 1.), p(1., 1., 0.)],
        [p(0., 0., 0.), p(0., 0., 1.), p(0., 1., 1.), p(0., 1., 0.)],
        [p(1., 0., 0.), p(1., 1., 0.), p(1., 1., 1.), p(1., 0., 1.)],
    ];
    let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    let uv = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let mut faces = Vec::new();
    let mut edges = Vec::new();
    let mut vertices = Vec::new();
    for (f, q) in quads.iter().enumerate() {
        let mut net: Net = (0..4)
            .map(|i| {
                (0..4)
                    .map(|j| {
                        let (u, v) = (i as f64 / 3.0, j as f64 / 3.0);
                        std::array::from_fn(|k| {
                            (1.0 - u) * (1.0 - v) * q[0][k]
                                + u * (1.0 - v) * q[1][k]
                                + u * v * q[2][k]
                                + (1.0 - u) * v * q[3][k]
                        })
                    })
                    .collect()
            })
            .collect();
        if f == 1 {
            for i in 1..3 {
                for j in 1..3 {
                    net[i][j][2] += dome;
                }
            }
        }
        let side = |s: usize| -> Vec<[f64; 3]> {
            (0..4)
                .map(|k| match s {
                    0 => net[k][0],
                    1 => net[3][k],
                    2 => net[3 - k][3],
                    _ => net[0][3 - k],
                })
                .collect()
        };
        for s in 0..4 {
            let pts = side(s);
            vertices.push(Vertex {
                point: pts[0],
                tolerance: 1e-6,
            });
            edges.push(Edge {
                curve: Curve::Spline {
                    knots: knots.clone(),
                    points: pts,
                },
                nurbs: None,
                interval: [0.0, 1.0],
                vertices: [(4 * f + s) as u32, (4 * f + (s + 1) % 4) as u32],
                tolerance: 1e-6,
                uses: vec![EdgeUseRef {
                    face: f as u32,
                    r#loop: 0,
                    position: s as u32,
                    upper: false,
                }],
            });
        }
        let loop_edges = (0..4)
            .map(|s| OrientedEdge {
                edge: (4 * f + s) as u32,
                forward: true,
                pcurve: Curve2::Spline {
                    knots: vec![0.0, 0.0, 1.0, 1.0],
                    points: vec![uv[s], uv[(s + 1) % 4]],
                },
            })
            .collect();
        faces.push(Face {
            surface: Surface::Spline {
                u_knots: knots.clone(),
                v_knots: knots.clone(),
                points: net,
            },
            nurbs: None,
            domain: [[0.0, 1.0], [0.0, 1.0]],
            loops: vec![Loop::Edges { edges: loop_edges }],
        });
    }
    let shells = (0..6u32)
        .map(|face| Shell::Faces {
            faces: vec![OrientedFace {
                face,
                forward: true,
            }],
            wires: vec![],
        })
        .collect();
    Brep {
        version: VERSION,
        tolerance: 1e-6,
        regions: vec![Region {
            shells: (0..6).collect(),
            void: false,
        }],
        shells,
        faces,
        edges,
        vertices,
    }
}

#[test]
#[ignore = "reproduction: separate bicubic sheets of one record, flat and domed top"]
fn unsewn_bicubic_sheets_sew_into_a_solid() {
    for dome in [0.0, 6.0] {
        let record = unsewn_bicubic_block_record([40.0, 30.0, 20.0], dome);
        println!("dome {dome}: validate {:?}", record.validate());
        for trusted in [false, true] {
            match Body::from_record(&record, trusted) {
                Err(e) => println!("  from_record(trusted={trusted}): Err({e:?})"),
                Ok((body, _)) => {
                    println!("  from_record(trusted={trusted}): Ok");
                    for tolerance in [1e-6, 1e-2] {
                        match stitch(body.copy().unwrap(), tolerance) {
                            Ok((sewn, closed)) => {
                                println!(
                                    "    stitch({tolerance:e}): closed={closed} volume={:?}",
                                    sewn.volume(1e-6)
                                );
                            }
                            Err(e) => println!("    stitch({tolerance:e}): Err({e:?})"),
                        }
                    }
                }
            }
        }
    }
}
