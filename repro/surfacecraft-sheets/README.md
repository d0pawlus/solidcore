# Reproduction: sheets from fitted surfaces (SurfaceCraft)

`cargo run --release` from this folder (standalone workspace, path dependencies on this checkout).

1. `topology::brep::sheet::face_on` on a least-squares cubic B-spline (`srf_fit::lst_sq_approx`, 7 × 6 net) returns `Upstream { operation: "topology::brep::p8", code: 1001 }`.
2. `ops::offset::thicken(2 mm)` succeeds on the spline as a one-face sheet (record built by `sheet_from_spline`, volume 4893.99) but returns `NotSolid` on the same geometry as its 12 Bézier patches sharing edges (`decompose_to_bez` + `sheet_from_patches`, counts 12/31/20).
