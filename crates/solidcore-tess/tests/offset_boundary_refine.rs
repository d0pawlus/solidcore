//! A cylinder face whose boundary edges stand off the cylinder by up to the edge's tolerance:
//! the sag refinement should settle at a chord a little coarser than the offset, and at finer
//! chords keep the mesh small and its area right. Kernel only: one face made with
//! `solidcore-record` types, meshed with `solidcore_tess::plan`; no application code or scan data.
//!
//! Run: `cargo test -p solidcore-tess --test offset_boundary_refine -- --ignored --nocapture`
#![allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::float_arithmetic,
    clippy::cast_possible_truncation,
    clippy::print_stdout
)] // test helpers

use solidcore_record::brep::{
    Brep, Curve, Curve2, Edge, EdgeUseRef, Face, Frame, Loop, OrientedEdge, OrientedFace, Region,
    Shell, Surface, VERSION, Vertex,
};
use solidcore_tess::{Mesh, MeshTolerance, plan};

const RADIUS: f64 = 4.0;
const START: f64 = 0.5;
const SPAN: f64 = 2.0;
const HEIGHT: f64 = 6.0;

/// The part of the cylinder `r = 4` between angles 0.5 and 2.5 and heights 0 and 6, as one face
/// whose four boundary edges lie on the cylinder of radius `4 + offset` (circle arcs and lines),
/// each edge carrying `tolerance = offset`.
fn cylinder_face(offset: f64) -> Brep {
    let (a0, a1) = (START, START + SPAN);
    let rr = RADIUS + offset;
    let at = |a: f64, z: f64| [rr * a.cos(), rr * a.sin(), z];
    let frame = Frame {
        origin: [0.0; 3],
        x_axis: [1.0, 0.0, 0.0],
        y_axis: [0.0, 1.0, 0.0],
    };
    let tol = offset.max(1e-6);
    let line = |s: [f64; 2], e: [f64; 2]| Curve2::Spline {
        knots: vec![0.0, 0.0, 1.0, 1.0],
        points: vec![s, e],
    };
    // the arc from START to START + SPAN as four quadratic rational Bezier arcs, exactly, on [0, 1]
    let arcs = 4;
    let phi = SPAN / (2.0 * arcs as f64);
    let arc_knots: Vec<f64> = [0.0, 0.0, 0.0]
        .into_iter()
        .chain((1..arcs).flat_map(|i| [i as f64 / arcs as f64; 2]))
        .chain([1.0, 1.0, 1.0])
        .collect();
    let ring = |radius: f64, z: f64| -> Vec<[f64; 4]> {
        (0..=2 * arcs)
            .map(|i| {
                let (r, w) = if i % 2 == 1 {
                    (radius / phi.cos(), phi.cos())
                } else {
                    (radius, 1.0)
                };
                let a = a0 + phi * i as f64;
                [r * a.cos() * w, r * a.sin() * w, z * w, w]
            })
            .collect()
    };
    let arc = |z: f64| Curve::RationalSpline {
        knots: arc_knots.clone(),
        homogeneous: ring(rr, z),
    };
    let vertices = vec![
        Vertex {
            point: at(a0, 0.0),
            tolerance: tol,
        },
        Vertex {
            point: at(a1, 0.0),
            tolerance: tol,
        },
        Vertex {
            point: at(a1, HEIGHT),
            tolerance: tol,
        },
        Vertex {
            point: at(a0, HEIGHT),
            tolerance: tol,
        },
    ];
    let edge = |curve: Curve, vertices: [u32; 2], position: u32, upper: bool| Edge {
        curve,
        nurbs: None,
        interval: [0.0, 1.0],
        vertices,
        tolerance: tol,
        uses: vec![EdgeUseRef {
            face: 0,
            r#loop: 0,
            position,
            upper,
        }],
    };
    let edges = vec![
        edge(arc(0.0), [0, 1], 0, false),
        edge(
            Curve::Line {
                start: at(a1, 0.0),
                end: at(a1, HEIGHT),
            },
            [1, 2],
            1,
            false,
        ),
        edge(arc(HEIGHT), [3, 2], 2, true),
        edge(
            Curve::Line {
                start: at(a0, 0.0),
                end: at(a0, HEIGHT),
            },
            [0, 3],
            3,
            true,
        ),
    ];
    let loop_edges = vec![
        OrientedEdge {
            edge: 0,
            forward: true,
            pcurve: line([0.0, 0.0], [1.0, 0.0]),
        },
        OrientedEdge {
            edge: 1,
            forward: true,
            pcurve: line([1.0, 0.0], [1.0, 1.0]),
        },
        OrientedEdge {
            edge: 2,
            forward: false,
            pcurve: line([0.0, 1.0], [1.0, 1.0]),
        },
        OrientedEdge {
            edge: 3,
            forward: false,
            pcurve: line([0.0, 0.0], [0.0, 1.0]),
        },
    ];
    // the face is on the cylinder of RADIUS, in the NURBS form of that cylinder
    let form = Surface::RationalSpline {
        u_knots: arc_knots,
        v_knots: vec![0.0, 0.0, 1.0, 1.0],
        homogeneous: ring(RADIUS, 0.0)
            .into_iter()
            .zip(ring(RADIUS, HEIGHT))
            .map(|(a, b)| vec![a, b])
            .collect(),
    };
    Brep {
        version: VERSION,
        tolerance: tol,
        regions: vec![Region {
            shells: vec![0],
            void: false,
        }],
        shells: vec![Shell::Faces {
            faces: vec![OrientedFace {
                face: 0,
                forward: true,
            }],
            wires: vec![],
        }],
        faces: vec![Face {
            surface: Surface::Cylinder {
                frame,
                radius: RADIUS,
            },
            nurbs: Some(form),
            domain: [[0.0, 1.0], [0.0, 1.0]],
            loops: vec![Loop::Edges { edges: loop_edges }],
        }],
        edges,
        vertices,
    }
}

fn area(m: &Mesh) -> f64 {
    m.triangles
        .iter()
        .map(|t| {
            let p = |i: u32| m.positions[i as usize];
            let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
            let (u, v) = (
                [b.x - a.x, b.y - a.y, b.z - a.z],
                [c.x - a.x, c.y - a.y, c.z - a.z],
            );
            let n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
        })
        .sum()
}

#[test]
#[ignore = "reproduction: refinement runs to its cap when the boundary stands off the surface"]
fn a_boundary_off_the_surface_does_not_run_the_refinement_to_its_cap() {
    let exact = RADIUS * SPAN * HEIGHT;
    let mut worst = (0, 0.0_f64);
    for offset in [0.0, 0.05, 0.1] {
        let record = cylinder_face(offset);
        println!("offset {offset}: validate {:?}", record.validate());
        for chord in [0.4, 0.1, 0.05, 0.03, 0.02, 0.015, 0.01, 0.008] {
            let (p, _) = plan(&record, &MeshTolerance::chord(chord), 8.0).unwrap();
            let (mesh, refused) = p.evaluate();
            let a = area(&mesh) / exact;
            println!(
                "  chord {chord}: {} triangles, area {a:.4} of the cylinder's, refused {}",
                mesh.triangles.len(),
                refused.len()
            );
            if offset > 0.0 {
                worst = (
                    worst.0.max(mesh.triangles.len()),
                    worst.1.max((a - 1.0).abs()),
                );
            }
        }
    }
    assert!(
        worst.0 <= 5_000,
        "up to {} triangles for a face of 4 mm radius",
        worst.0
    );
    assert!(worst.1 <= 0.02, "area off by {:.3}", worst.1);
}
