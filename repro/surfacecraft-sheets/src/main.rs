//! Reproductions for upstream issues (SurfaceCraft scan-to-CAD sheets).
mod nurbs;
mod err {
    #[derive(Debug)]
    pub enum KernelError { Argument(String), Op { op: &'static str, message: String } }
    pub type Result<T> = std::result::Result<T, KernelError>;
    pub fn op_err(op: &'static str, e: impl std::fmt::Display) -> KernelError { KernelError::Op { op, message: e.to_string() } }
}
use solidcore::kernel::{Kernel, KernelOptions};
use solidcore::topology::Brep;

fn wave() -> Vec<Vec<[f64; 3]>> {
    (0..30).map(|i| (0..20).map(|j| {
        let (x, y) = (60.0 * f64::from(i) / 29.0, 40.0 * f64::from(j) / 19.0);
        [x, y, 3.0 * (x / 12.0).sin() * (y / 9.0).cos()]
    }).collect()).collect()
}

fn main() {
    let s = nurbs::fit_grid(&wave(), 7, 6, 3).expect("fit");
    // 1. sheet from surface through the kernel
    let solidcore_record::brep::Surface::Spline { .. } = &s else { panic!() };
    let mut k = Kernel::new(&KernelOptions::default());
    let ks = {
        let n = nurbs::to_nurb_pub(&s).expect("nurb");
        solidcore::surface::surface::Surface::new(solidcore::surface::surface::SurfaceKind::BSpline(
            solidcore::nurbs::b_spline_surface::BSplineSurface::from_nurb(n, false)))
    };
    let mut sheet = Brep::new();
    let r = ks.natural_uv_domain_with_nurb(&k).and_then(|d| solidcore::topology::brep::sheet::face_on(&mut k, &mut sheet, ks, &d));
    println!("1. kernel sheet::face_on on a fitted cubic B-spline: {:?}", r.map(|_| ()));
    // 2. thicken: one-face sheet vs the same spline as its Bezier patch network
    let one = nurbs::sheet_from_spline(&s).expect("one-face sheet");
    let b1 = solidcore_ops::Body::from_record(&one, false).expect("decode").0;
    println!("2a. thicken one-face sheet 2 mm: {:?} volume {:?}", solidcore_ops::offset::thicken(&b1, 2.0, false), b1.volume(1e-6));
    let net = nurbs::sheet_from_patches(&nurbs::decompose(&s).expect("decompose")).expect("network");
    let b2 = solidcore_ops::Body::from_record(&net, false).expect("decode").0;
    println!("2b. same geometry as {} faces sharing edges, counts {:?}: thicken 2 mm: {:?}", net.faces.len(), b2.counts(), solidcore_ops::offset::thicken(&b2, 2.0, false));
}
mod solid {
    use crate::err::{Result, op_err};
    pub fn decode(r: &solidcore_record::brep::Brep) -> Result<solidcore_ops::Body> {
        solidcore_ops::Body::from_record(r, false).map(|b| b.0).map_err(|e| op_err("from_record", e))
    }
    pub fn encode(b: &solidcore_ops::Body) -> Result<solidcore_record::brep::Brep> {
        b.to_record().map(|r| r.0).map_err(|e| op_err("to_record", e))
    }
}
