//! `BSplineSurface::copy_and_add_analytics` on an exact rational-spline cylinder (two quadratic
//! Bezier arcs in u, one linear span in v): the surface is a cylinder wherever it is, so it should
//! come back as the cone/cylinder class. Kernel only: no records, no application code, no scan data.
//!
//! Run: `cargo test -p solidcore --test rational_cylinder_analytics -- --ignored --nocapture`
#![allow(clippy::unwrap_used, clippy::float_arithmetic, clippy::print_stdout)] // test helpers

use solidcore::Kernel;
use solidcore::nurbs::b_spline_surface::BSplineSurface;
use solidcore::nurbs::math::CPoint;
use solidcore::nurbs::srf_basic::Surface as Nurb;
use solidcore::surface::surface::p2::SurfaceKind;

/// The cylinder of `radius` about the z axis through `origin`, from angle `start` to `start + arc`,
/// between heights 0 and `height`, exactly: two quadratic rational Bezier arcs.
fn cylinder(origin: [f64; 3], radius: f64, start: f64, arc: f64, height: f64) -> BSplineSurface {
    let phi = arc / 4.0;
    let w = phi.cos();
    let ring = |z: f64| -> Vec<CPoint> {
        (0..5)
            .map(|i| {
                let a = start + phi * i as f64;
                let (r, wi) = if i % 2 == 1 {
                    (radius / w, w)
                } else {
                    (radius, 1.0)
                };
                CPoint::new(
                    (origin[0] + r * a.cos()) * wi,
                    (origin[1] + r * a.sin()) * wi,
                    (origin[2] + z) * wi,
                    wi,
                )
            })
            .collect()
    };
    let net: Vec<Vec<CPoint>> = ring(0.0)
        .into_iter()
        .zip(ring(height))
        .map(|(a, b)| vec![a, b])
        .collect();
    let u = vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0];
    let v = vec![0.0, 0.0, 1.0, 1.0];
    BSplineSurface::from_nurb(Nurb::new(net, 2, 1, u, v).unwrap(), false)
}

fn is_cylinder(s: &BSplineSurface) -> bool {
    let k = Kernel::default();
    let analytic = s.copy_and_add_analytics(&k).unwrap();
    matches!(analytic.kind(), SurfaceKind::Cone(_))
}

#[test]
#[ignore = "reproduction: an exact rational cylinder is not recognised at some placements"]
fn an_exact_rational_cylinder_is_recognised_wherever_it_is() {
    let cases = [
        (
            "r 2.4, arc 0.2, h 3.1, at the origin",
            [0.0, 0.0, 0.0],
            2.4,
            0.0,
            0.2,
            3.1,
        ),
        (
            "same, moved 1e-6 in x",
            [1e-6, 0.0, 0.0],
            2.4,
            0.0,
            0.2,
            3.1,
        ),
        (
            "same, moved (7.3, 1.1, 0)",
            [7.3, 1.1, 0.0],
            2.4,
            0.0,
            0.2,
            3.1,
        ),
        (
            "same, starting at angle 0.3",
            [0.0, 0.0, 0.0],
            2.4,
            0.3,
            0.2,
            3.1,
        ),
        ("arc 0.15", [0.0; 3], 2.4, 0.0, 0.15, 3.1),
        ("arc 0.25", [0.0; 3], 2.4, 0.0, 0.25, 3.1),
        ("arc 0.3", [0.0; 3], 2.4, 0.0, 0.3, 3.1),
        ("r 5", [0.0; 3], 5.0, 0.0, 0.2, 3.1),
        ("h 8", [0.0; 3], 2.4, 0.0, 0.2, 8.0),
        (
            "quarter turn, r 2.4, h 3.1",
            [0.0; 3],
            2.4,
            0.0,
            std::f64::consts::FRAC_PI_2,
            3.1,
        ),
    ];
    let mut missed = Vec::new();
    for (name, origin, r, start, arc, h) in cases {
        let ok = is_cylinder(&cylinder(origin, r, start, arc, h));
        println!("{name}: {}", if ok { "cylinder" } else { "spline" });
        if !ok {
            missed.push(name);
        }
    }
    assert!(
        missed.is_empty(),
        "not recognised as a cylinder: {missed:?}"
    );
}
