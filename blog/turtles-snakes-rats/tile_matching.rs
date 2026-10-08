#!/usr/bin/env rust-script
//! ```cargo
//! [package]
//! edition = "2021"
//!
//! [dependencies]
//! tilezz = "=0.1.2"
//! ```
//!
//! Standalone generator for the "Finding matches between tiles" figure:
//! the threefold junction-angle case distinction when gluing two tiles
//! along a boundary match. Run with just `rust-script`:
//!
//!   rust-script tile_matching.rs [OUTDIR]
//!
//! It writes one SVG (`tile_matching.svg`) into OUTDIR (default `.`) with
//! three panels, all using two spectre tiles glued along a chosen match
//! interval.
//!
//! ## Setup
//!
//! Two boundaries are cyclic angle sequences. A match glues A's edges
//! `a1..al` to B's edges `bl..b1` (anti-parallel reverse-complement: a1
//! pairs with bl, al with b1). The cheap pre-filter checks the two
//! junction sums where the match meets the surviving boundary -- exactly
//! `tilezz::geom::glue::junction_sums`:
//!
//!   cw  = a[ns]              + b[ne]
//!   ccw = a[(ns + mlen) % n] + b[(ne + m - mlen) % m]
//!
//! with the threefold meaning (junction angle = sum - hturn):
//!
//!   sum > 0 : valid, non-degenerate corner; the match is MAXIMAL here.
//!   sum = 0 : the revcomp identity still holds one step out -> the match
//!             is NOT maximal, it could be EXTENDED past this end.
//!   sum < 0 : the junction angle drops below -hturn -> the surviving
//!             edges fold back and INTERSECT; no simple glue exists.
//!
//! ## Tiles
//!
//! The spectre over ZZ12 is the 14-edge turn sequence below (units of
//! 1/12 turn, turn = 12, hturn = 6). The figure renders it through
//! tilezz's own `vis` pipeline; tile B is placed by the rigid transform
//! that maps its matched edge onto A's (self-verified to coincide).
//!
//! ## Self-checks (panic before writing anything if violated)
//!
//! - the placed B's matched run coincides with A's vertex-for-vertex;
//! - `junction_sums` reproduces the documented sign for each panel;
//! - `>0` panel: A and B share only the matched edges (no other crossing);
//! - `<0` panel: a surviving A-edge and B-edge genuinely cross (the fold);
//! - `=0` panel: the next edge out coincides (the match could extend).

use std::path::PathBuf;

use tilezz::cyclotomic::ZZ12;
use tilezz::geom::snake::Turtle;
use tilezz::geom::tiles::spectre;
use tilezz::vis::plotutils::{P64, R64};
use tilezz::vis::scene::{Color, Fill, Item, Scene, Stroke, TextStyle, Viewport};

/// Spectre boundary turn sequence over ZZ12 (see `tiles::spectre`).
const SPECTRE: [i32; 14] = [3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
const N: usize = 14;

// ---- small vector helpers ----
fn sub(a: P64, b: P64) -> P64 {
    (a.0 - b.0, a.1 - b.1)
}
fn dist(a: P64, b: P64) -> f64 {
    (sub(a, b).0.powi(2) + sub(a, b).1.powi(2)).sqrt()
}
fn lerp(a: P64, b: P64, t: f64) -> P64 {
    (a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t)
}
/// Unit left normal of the directed edge a->b (points into a CCW polygon).
fn left_normal(a: P64, b: P64) -> P64 {
    let d = sub(b, a);
    let l = (d.0 * d.0 + d.1 * d.1).sqrt();
    (-d.1 / l, d.0 / l)
}

/// Ray-cast point-in-polygon test (f64). Boundary cases are not relied on.
fn point_in_poly(p: P64, poly: &[P64]) -> bool {
    let mut inside = false;
    let n = poly.len();
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = poly[i];
        let (xj, yj) = poly[j];
        if ((yi > p.1) != (yj > p.1)) && (p.0 < (xj - xi) * (p.1 - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn centroid(pts: &[P64]) -> P64 {
    let s = pts.iter().fold((0.0, 0.0), |a, &q| (a.0 + q.0, a.1 + q.1));
    (s.0 / pts.len() as f64, s.1 / pts.len() as f64)
}

/// The 14 spectre vertices (no trailing duplicate), in math coords.
fn spectre_verts() -> Vec<P64> {
    let mut v = spectre::<ZZ12>().to_polyline_f64(Turtle::default());
    if v.len() > 1 && v.first() == v.last() {
        v.pop();
    }
    assert_eq!(v.len(), N, "spectre should have {N} vertices");
    v
}

// ---- junction sums (the pre-filter being illustrated) ----
fn junction_sums(ns: usize, mlen: usize, ne: usize) -> (i32, i32) {
    let cw = SPECTRE[ns % N] + SPECTRE[ne % N];
    let ccw = SPECTRE[(ns + mlen) % N] + SPECTRE[(ne + N - mlen % N) % N];
    (cw, ccw)
}

/// Place tile B (a second spectre, same vertices `va`) so its matched run
/// `[ne-mlen, ne)` lands anti-parallel onto A's matched run `[ns, ns+mlen)`.
/// Returns B's 14 placed vertices. Panics if the rigid placement does not
/// reproduce the match (which would mean a bogus interval / reflection).
fn place_b(va: &[P64], ns: usize, mlen: usize, ne: usize) -> Vec<P64> {
    let pa_s = va[ns % N];
    let pa_e = va[(ns + mlen) % N];
    let pb_s = va[ne % N];
    let pb_e = va[(ne + N - mlen % N) % N];
    let ua = sub(pa_e, pa_s);
    let ub = sub(pb_e, pb_s);
    assert!(
        ((ua.0 * ua.0 + ua.1 * ua.1).sqrt() - (ub.0 * ub.0 + ub.1 * ub.1).sqrt()).abs() < 1e-9,
        "matched runs differ in length"
    );
    let ang = ua.1.atan2(ua.0) - ub.1.atan2(ub.0);
    let (c, s) = (ang.cos(), ang.sin());
    let place = |p: P64| {
        let q = sub(p, pb_s);
        (q.0 * c - q.1 * s + pa_s.0, q.0 * s + q.1 * c + pa_s.1)
    };
    let placed: Vec<P64> = va.iter().map(|&p| place(p)).collect();
    for k in 0..=mlen {
        let d = dist(va[(ns + k) % N], placed[(ne + N - k % N) % N]);
        assert!(d < 1e-6, "placement mismatch at k={k}: d={d}");
    }
    placed
}

// ---- proper segment crossing (for the fold / clean-glue self-checks) ----
fn orient(a: P64, b: P64, c: P64) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}
/// Intersection point if `ab` and `cd` cross in their interiors, else None.
fn proper_cross(a: P64, b: P64, c: P64, d: P64) -> Option<P64> {
    let (o1, o2, o3, o4) = (
        orient(a, b, c),
        orient(a, b, d),
        orient(c, d, a),
        orient(c, d, b),
    );
    if o1 * o2 < -1e-9 && o3 * o4 < -1e-9 {
        let t = o1 / (o1 - o2);
        Some(lerp(c, d, t))
    } else {
        None
    }
}

/// Every proper crossing between an A-edge and a placed-B-edge, excluding
/// the shared matched edges on each side.
fn cross_points(va: &[P64], pb: &[P64], ns: usize, mlen: usize, ne: usize) -> Vec<P64> {
    let a_matched: Vec<usize> = (0..mlen).map(|k| (ns + k) % N).collect();
    let b_matched: Vec<usize> = (0..mlen).map(|k| (ne + N - 1 - k) % N).collect();
    let mut hits = Vec::new();
    for i in 0..N {
        if a_matched.contains(&i) {
            continue;
        }
        for j in 0..N {
            if b_matched.contains(&j) {
                continue;
            }
            if let Some(p) = proper_cross(va[i], va[(i + 1) % N], pb[j], pb[(j + 1) % N]) {
                hits.push(p);
            }
        }
    }
    hits
}

// ============================================================
// rendering
// ============================================================

const A_FILL: Color = Color::rgba(255, 210, 90, 80);
const B_FILL: Color = Color::rgba(150, 152, 168, 70);
const OUTLINE: Color = Color::rgb(60, 60, 60);
const MATCH_COL: Color = Color::rgb(20, 95, 205); // the match itself: blue
const A_LBL: Color = Color::rgb(150, 110, 0);
const B_LBL: Color = Color::rgb(70, 80, 100);
const RED: Color = Color::rgb(210, 35, 45); // the two edges of a sum<0 junction
const GREEN: Color = Color::rgb(20, 150, 70); // the could-extend edge
const SUMPOS: Color = Color::rgb(90, 90, 90); // a junction sum > 0

/// Which special junction behaviour to draw on top of a panel.
enum Case {
    Valid,
    Extend,
    Invalid,
}

struct Panel {
    ns: usize,
    mlen: usize,
    ne: usize,
    case: Case,
    title: String,
    note: String,
}

/// Push every primitive of one panel into `scene`, shifted by `off`.
/// `cx` is the panel's centre x (its slot centre); `y_top` / `y_bot` are
/// the shared title / caption baselines, so all three panels line up.
fn add_panel(scene: &mut Scene, off: P64, cx: f64, va: &[P64], p: &Panel, y_top: f64, y_bot: f64) {
    let (ns, mlen, ne) = (p.ns, p.mlen, p.ne);
    let pb = place_b(va, ns, mlen, ne);
    let (cw, ccw) = junction_sums(ns, mlen, ne);

    let sh = |q: P64| (q.0 + off.0, q.1 + off.1);
    let shv = |pts: &[P64]| -> Vec<P64> { pts.iter().map(|&q| sh(q)).collect() };

    // Combined centroid: side labels / sum annotations are pushed away
    // from it so they sit outside the tile bodies.
    let cen = {
        let s = va.iter().chain(&pb).fold((0.0, 0.0), |a, &q| (a.0 + q.0, a.1 + q.1));
        let nn = (va.len() + pb.len()) as f64;
        (s.0 / nn, s.1 / nn)
    };
    let outward = |q: P64, d: f64| {
        let v = (q.0 - cen.0, q.1 - cen.1);
        let l = (v.0 * v.0 + v.1 * v.1).sqrt().max(1e-9);
        sh((q.0 + v.0 / l * d, q.1 + v.1 / l * d))
    };

    // Tiles A (pale yellow) and B (neutral grey), translucent so the
    // overlap in the invalid case shows through.
    scene.push(Item::Polygon {
        points: shv(va),
        fill: Some(Fill::solid(A_FILL)),
        stroke: Some(Stroke::solid(OUTLINE, 0.03)),
        arrow: None,
    });
    scene.push(Item::Polygon {
        points: shv(&pb),
        fill: Some(Fill::solid(B_FILL)),
        stroke: Some(Stroke::solid(OUTLINE, 0.03)),
        arrow: None,
    });

    // Tile letters, placed on each tile's *own* part: if a tile's centroid
    // falls inside the other (the invalid overlap case), use the centroid
    // of just its non-overlapping vertices instead.
    let own_spot = |mine: &[P64], other: &[P64]| -> P64 {
        let c = centroid(mine);
        if point_in_poly(c, other) {
            let out: Vec<P64> = mine.iter().copied().filter(|&v| !point_in_poly(v, other)).collect();
            if out.is_empty() { c } else { centroid(&out) }
        } else {
            c
        }
    };
    // In the invalid (overlap) case A's non-overlap centroid still hugs the
    // intersection; shift it horizontally onto A's own (left) spike, keeping
    // its height.
    let a_pos = own_spot(va, &pb);
    let a_pos = if point_in_poly(centroid(va), &pb) {
        (a_pos.0 - 1.0, a_pos.1)
    } else {
        a_pos
    };
    for (name, pos) in [("A", a_pos), ("B", own_spot(&pb, va))] {
        scene.push(Item::Text {
            pos: sh(pos),
            text: name.into(),
            style: TextStyle::new(0.5, OUTLINE).bold(),
        });
    }

    // The matched run, drawn blue; per-edge labels straddle it -- a_i on
    // A's side, b_{l+1-i} on B's side, so a1 pairs with bl.
    let run: Vec<P64> = (0..=mlen).map(|k| sh(va[(ns + k) % N])).collect();
    scene.push(Item::Polyline {
        points: run,
        stroke: Stroke::solid(MATCH_COL, 0.085),
        arrow: None,
    });
    for k in 0..mlen {
        let a = va[(ns + k) % N];
        let b = va[(ns + k + 1) % N];
        let mid = lerp(a, b, 0.5);
        let nl = left_normal(a, b);
        // Keep the offset SMALL so each label tracks its own edge midpoint
        // (well separated in x along the run); a large offset pulls
        // neighbours together toward the corner bisector and they collide.
        let d = 0.22;
        scene.push(Item::Text {
            pos: sh((mid.0 + nl.0 * d, mid.1 + nl.1 * d)),
            text: format!("a{}", k + 1),
            style: TextStyle::new(0.20, A_LBL).bold(),
        });
        scene.push(Item::Text {
            pos: sh((mid.0 - nl.0 * d, mid.1 - nl.1 * d)),
            text: format!("b{}", mlen - k),
            style: TextStyle::new(0.20, B_LBL).bold(),
        });
    }

    // Junction dots, each annotated with its signed sum (coloured by
    // sign): the threefold check, shown where it applies. The label sits
    // just outside the corner along the exterior bisector of the two
    // surviving edges -- using the direct bisector for a reflex corner
    // (sum < hturn) and its negation for a convex one, and falling back
    // to "away from the match" when those edges are antiparallel (sum 0).
    let norm = |v: P64| {
        let l = (v.0 * v.0 + v.1 * v.1).sqrt().max(1e-9);
        (v.0 / l, v.1 / l)
    };
    let junctions = [
        (va[ns % N], va[(ns + N - 1) % N], pb[(ne + 1) % N], va[(ns + 1) % N], cw),
        (
            va[(ns + mlen) % N],
            va[(ns + mlen + 1) % N],
            pb[(ne + N - mlen - 1) % N],
            va[(ns + mlen - 1) % N],
            ccw,
        ),
    ];
    for (j, p1, p2, madj, val) in junctions {
        scene.push(Item::Marker {
            center: sh(j),
            shape: tilezz::vis::scene::MarkerShape::Circle,
            size: 0.14,
            fill: Some(Fill::solid(OUTLINE)),
            stroke: None,
        });
        let u1 = norm((p1.0 - j.0, p1.1 - j.1));
        let u2 = norm((p2.0 - j.0, p2.1 - j.1));
        let s = (u1.0 + u2.0, u1.1 + u2.1);
        let ls = (s.0 * s.0 + s.1 * s.1).sqrt();
        let dir0 = if ls < 0.4 {
            norm((j.0 - madj.0, j.1 - madj.1)) // antiparallel: away from the match
        } else if val < 6 {
            (s.0 / ls, s.1 / ls) // reflex corner
        } else {
            (-s.0 / ls, -s.1 / ls) // convex corner
        };
        // Nudge upward so the number clears the (often low) edge; the
        // antiparallel (sum 0) case sits right on the edge, so lift it more.
        let up = if val == 0 { 0.9 } else { 0.3 };
        let dir = norm((dir0.0, dir0.1 + up));
        let col = if val < 0 {
            RED
        } else if val == 0 {
            GREEN
        } else {
            SUMPOS
        };
        let txt = if val > 0 {
            format!("+{val}")
        } else {
            format!("{val}")
        };
        scene.push(Item::Text {
            pos: sh((j.0 + dir.0 * 0.5, j.1 + dir.1 * 0.5)),
            text: txt,
            style: TextStyle::new(0.34, col).bold(),
        });
    }

    // Case-specific overlay.
    match p.case {
        Case::Valid => {}
        Case::Extend => {
            // Not maximal: the next edge out coincides with B's next edge.
            // Draw it green/dashed; put its label off to the side.
            let nxt = va[(ns + N - 1) % N];
            let here = va[ns % N];
            assert!(
                dist(pb[(ne + 1) % N], nxt) < 1e-6,
                "extend panel: next edge does not coincide (cw should be 0)"
            );
            scene.push(Item::Segment {
                a: sh(nxt),
                b: sh(here),
                stroke: Stroke::dashed(GREEN, 0.085, vec![0.13, 0.10]),
                arrow: None,
            });
            scene.push(Item::Text {
                pos: outward(lerp(nxt, here, 0.5), 1.75),
                text: "could extend".into(),
                style: TextStyle::new(0.30, GREEN).bold(),
            });
        }
        Case::Invalid => {
            // The two surviving edges meeting at the sum<0 junction -- the
            // pair whose angle folds below -hturn and would intersect.
            assert!(
                !cross_points(va, &pb, ns, mlen, ne).is_empty(),
                "invalid panel: expected a fold crossing"
            );
            let (ea, eb) = if cw < 0 {
                (
                    (va[(ns + N - 1) % N], va[ns % N]),
                    (pb[ne % N], pb[(ne + 1) % N]),
                )
            } else {
                (
                    (va[(ns + mlen) % N], va[(ns + mlen + 1) % N]),
                    (pb[(ne + N - mlen - 1) % N], pb[(ne + N - mlen) % N]),
                )
            };
            for (q1, q2) in [ea, eb] {
                scene.push(Item::Segment {
                    a: sh(q1),
                    b: sh(q2),
                    stroke: Stroke::solid(RED, 0.09),
                    arrow: None,
                });
            }
        }
    }

    // Aligned title (top) and caption (bottom).
    scene.push(Item::Text {
        pos: (cx, y_top),
        text: p.title.clone(),
        style: TextStyle::new(0.44, OUTLINE).bold(),
    });
    scene.push(Item::Text {
        pos: (cx, y_bot),
        text: p.note.clone(),
        style: TextStyle::new(0.34, OUTLINE),
    });
}

fn bbox(a: &[P64], b: &[P64]) -> R64 {
    let mut lo = (f64::INFINITY, f64::INFINITY);
    let mut hi = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &(x, y) in a.iter().chain(b) {
        lo = (lo.0.min(x), lo.1.min(y));
        hi = (hi.0.max(x), hi.1.max(y));
    }
    (lo, hi)
}

fn main() {
    let dir: PathBuf = std::env::args().nth(1).unwrap_or_else(|| ".".to_string()).into();
    let va = spectre_verts();

    let panels = [
        Panel {
            ns: 9, mlen: 3, ne: 5, case: Case::Valid,
            title: "sum > 0: valid".into(),
            note: "both junctions are proper corners".into(),
        },
        Panel {
            ns: 11, mlen: 3, ne: 5, case: Case::Extend,
            title: "sum = 0: not maximal".into(),
            note: "next edges coincide: extend the match".into(),
        },
        Panel {
            ns: 13, mlen: 1, ne: 8, case: Case::Invalid,
            title: "sum < 0: intersection".into(),
            note: "adjacent edges cross: reject".into(),
        },
    ];

    // Measure each panel's glued-pair bbox to size the slots and the
    // shared title / caption baselines.
    let mut widths = Vec::new();
    let mut heights = Vec::new();
    for p in &panels {
        let pb = place_b(&va, p.ns, p.mlen, p.ne);
        let bb = bbox(&va, &pb);
        widths.push(bb.1 .0 - bb.0 .0);
        heights.push(bb.1 .1 - bb.0 .1);
    }
    let slot = widths.iter().cloned().fold(0.0_f64, f64::max) + 3.0;
    let maxh = heights.iter().cloned().fold(0.0_f64, f64::max);
    let y_top = maxh / 2.0 + 1.0;
    let y_bot = -maxh / 2.0 - 1.0;

    // Each panel is centred (in both axes) within its slot at x = j*slot.
    let mut scene = Scene::new().with_background(Color::WHITE);
    for (j, p) in panels.iter().enumerate() {
        let pb = place_b(&va, p.ns, p.mlen, p.ne);
        let bb = bbox(&va, &pb);
        let cxl = (bb.0 .0 + bb.1 .0) * 0.5;
        let cyl = (bb.0 .1 + bb.1 .1) * 0.5;
        let cx = j as f64 * slot;
        add_panel(&mut scene, (cx - cxl, -cyl), cx, &va, p, y_top, y_bot);
    }

    // Explicit bounds so the titles / captions (which auto_bounds ignores)
    // are always framed.
    let bounds = (
        (-slot / 2.0, y_bot - 0.6),
        (2.0 * slot + slot / 2.0, y_top + 0.6),
    );
    let w = bounds.1 .0 - bounds.0 .0;
    let h = bounds.1 .1 - bounds.0 .1;
    let px_h = 900u32;
    let px_w = (px_h as f64 * w / h).round() as u32;
    let vp = Viewport::rect_for(px_w, px_h, bounds, 8);

    let path = dir.join("tile_matching.svg");
    std::fs::write(&path, scene.to_svg(&vp)).expect("write svg");
    println!("wrote {}", path.display());
    for p in &panels {
        let (cw, ccw) = junction_sums(p.ns, p.mlen, p.ne);
        println!("  {} (ns={} mlen={} ne={}): cw={cw} ccw={ccw}", p.title, p.ns, p.mlen, p.ne);
    }
}
