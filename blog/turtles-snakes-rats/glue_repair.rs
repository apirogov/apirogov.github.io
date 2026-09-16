#!/usr/bin/env rust-script
//! ```cargo
//! [package]
//! edition = "2021"
//!
//! [dependencies]
//! tilezz = "=0.1.2"
//! ```
//!
//! Standalone generator for the "Matching Boundaries" glue/angle-repair
//! figure. Run with just `rust-script`:
//!
//!   rust-script glue_repair.rs [OUTDIR]
//!
//! It writes one SVG (`glue_repair.svg`) into OUTDIR (default `.`) with two
//! panels (BEFORE -> AFTER) of the same two spectres glued along the same
//! maximal match used in the matching figure (ns=9, mlen=3, ne=5):
//!
//! - BEFORE: tiles A and B glued along the match. The matched run a1..al is
//!   drawn solid (it glues to bl..b1, a1<->bl); the four surviving turns
//!   adjacent to the match -- a0, a_{l+1} on A and b0, b_{l+1} on B -- are
//!   labelled, since those are the ones the glue merges.
//! - AFTER: the merged tile. The matched edges are gone (drawn dashed: they
//!   are now interior). The two junction vertices carry the repaired angles
//!     x = a_{l+1} + b0 - H      (where surviving a_{l+1} meets surviving b0)
//!     y = b_{l+1} + a0 - H      (where surviving b_{l+1} meets surviving a0)
//!   with H the half-turn. The resulting boundary sequence is
//!   x a_{l+2}..a_{-1} y b_{l+2}..b_{-1}.
//!
//! The placement of B is the rigid transform that maps its matched edge onto
//! A's; it is self-verified to coincide before anything is drawn.

use std::path::PathBuf;

use tilezz::cyclotomic::ZZ12;
use tilezz::geom::snake::Turtle;
use tilezz::geom::tiles::spectre;
use tilezz::vis::plotutils::{P64, R64};
use tilezz::vis::scene::{Color, Fill, Item, MarkerShape, Scene, Stroke, TextStyle, Viewport};

const N: usize = 14;
const NS: usize = 9;
const MLEN: usize = 3;
const NE: usize = 5;

// ---- vector helpers ----
fn sub(a: P64, b: P64) -> P64 {
    (a.0 - b.0, a.1 - b.1)
}
fn dist(a: P64, b: P64) -> f64 {
    (sub(a, b).0.powi(2) + sub(a, b).1.powi(2)).sqrt()
}
fn lerp(a: P64, b: P64, t: f64) -> P64 {
    (a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t)
}
fn left_normal(a: P64, b: P64) -> P64 {
    let d = sub(b, a);
    let l = (d.0 * d.0 + d.1 * d.1).sqrt();
    (-d.1 / l, d.0 / l)
}

fn spectre_verts() -> Vec<P64> {
    let mut v = spectre::<ZZ12>().to_polyline_f64(Turtle::default());
    if v.len() > 1 && v.first() == v.last() {
        v.pop();
    }
    assert_eq!(v.len(), N, "spectre should have {N} vertices");
    v
}

/// Place tile B onto A's matched run (anti-parallel); panics if the rigid
/// placement does not reproduce the match.
fn place_b(va: &[P64]) -> Vec<P64> {
    let pa_s = va[NS % N];
    let pa_e = va[(NS + MLEN) % N];
    let pb_s = va[NE % N];
    let pb_e = va[(NE + N - MLEN) % N];
    let ua = sub(pa_e, pa_s);
    let ub = sub(pb_e, pb_s);
    let ang = ua.1.atan2(ua.0) - ub.1.atan2(ub.0);
    let (c, s) = (ang.cos(), ang.sin());
    let place = |p: P64| {
        let q = sub(p, pb_s);
        (q.0 * c - q.1 * s + pa_s.0, q.0 * s + q.1 * c + pa_s.1)
    };
    let placed: Vec<P64> = va.iter().map(|&p| place(p)).collect();
    for k in 0..=MLEN {
        assert!(
            dist(va[(NS + k) % N], placed[(NE + N - k) % N]) < 1e-6,
            "placement mismatch at k={k}"
        );
    }
    placed
}

// ---- colours ----
const A_FILL: Color = Color::rgba(255, 210, 90, 80);
const B_FILL: Color = Color::rgba(150, 152, 168, 70);
const MERGED_FILL: Color = Color::rgba(120, 190, 130, 80);
const OUTLINE: Color = Color::rgb(60, 60, 60);
const MATCH_COL: Color = Color::rgb(20, 95, 205);
const SEAM_COL: Color = Color::rgb(150, 150, 150);
const A_LBL: Color = Color::rgb(150, 110, 0);
const B_LBL: Color = Color::rgb(70, 80, 100);
const XY_COL: Color = Color::rgb(170, 0, 120);

fn text(scene: &mut Scene, pos: P64, s: &str, size: f64, color: Color, bold: bool) {
    let mut st = TextStyle::new(size, color);
    if bold {
        st = st.bold();
    }
    scene.push(Item::Text { pos, text: s.into(), style: st });
}

/// Label edge `i` of polyline `pts`, offset `off` into the tile interior
/// (left of the CCW edge), already shifted by `sh`.
fn elabel(scene: &mut Scene, sh: &dyn Fn(P64) -> P64, pts: &[P64], i: usize, s: &str, color: Color) {
    let a = pts[i % N];
    let b = pts[(i + 1) % N];
    let mid = lerp(a, b, 0.5);
    let nl = left_normal(a, b);
    text(scene, sh((mid.0 + nl.0 * 0.24, mid.1 + nl.1 * 0.24)), s, 0.2, color, true);
}

fn centroid(pts: &[P64]) -> P64 {
    let s = pts.iter().fold((0.0, 0.0), |a, &q| (a.0 + q.0, a.1 + q.1));
    (s.0 / pts.len() as f64, s.1 / pts.len() as f64)
}

/// Draw the polyline of the matched run (A vertices NS..NS+MLEN), shifted.
fn matched_run(va: &[P64], sh: &dyn Fn(P64) -> P64) -> Vec<P64> {
    (0..=MLEN).map(|k| sh(va[(NS + k) % N])).collect()
}

fn poly(scene: &mut Scene, pts: Vec<P64>, fill: Color) {
    scene.push(Item::Polygon {
        points: pts,
        fill: Some(Fill::solid(fill)),
        stroke: Some(Stroke::solid(OUTLINE, 0.03)),
        arrow: None,
    });
}

fn dot(scene: &mut Scene, p: P64) {
    scene.push(Item::Marker {
        center: p,
        shape: MarkerShape::Circle,
        size: 0.14,
        fill: Some(Fill::solid(OUTLINE)),
        stroke: None,
    });
}

fn add_before(scene: &mut Scene, off: P64, va: &[P64], pb: &[P64], y_top: f64, y_bot: f64) {
    let sh = move |q: P64| (q.0 + off.0, q.1 + off.1);
    poly(scene, va.iter().map(|&q| sh(q)).collect(), A_FILL);
    poly(scene, pb.iter().map(|&q| sh(q)).collect(), B_FILL);

    text(scene, sh(centroid(va)), "A", 0.5, OUTLINE, true);
    text(scene, sh(centroid(pb)), "B", 0.5, OUTLINE, true);

    // Matched run, solid blue.
    scene.push(Item::Polyline {
        points: matched_run(va, &sh),
        stroke: Stroke::solid(MATCH_COL, 0.085),
        arrow: None,
    });

    // Matched-edge labels: a1..al on A, bl..b1 on B (a1 pairs with bl).
    for k in 0..MLEN {
        elabel(scene, &sh, va, NS + k, &format!("a{}", k + 1), A_LBL);
        elabel(scene, &sh, pb, NE - 1 - k, &format!("b{}", MLEN - k), B_LBL);
    }
    // Surviving turns adjacent to the match -- the ones the glue merges.
    elabel(scene, &sh, va, NS + N - 1, "a0", A_LBL); // a_0  (before a1)
    elabel(scene, &sh, va, NS + MLEN, "a4", A_LBL); // a_{l+1} (after al)
    elabel(scene, &sh, pb, NE, "b4", B_LBL); // b_{l+1}
    elabel(scene, &sh, pb, NE + N - MLEN - 1, "b0", B_LBL); // b_0

    dot(scene, sh(va[NS % N]));
    dot(scene, sh(va[(NS + MLEN) % N]));

    let cx = (sh(centroid(va)).0 + sh(centroid(pb)).0) * 0.5;
    text(scene, (cx, y_top), "before: maximal match", 0.42, OUTLINE, true);
    text(scene, (cx, y_bot), "glue removes a1..al / bl..b1", 0.34, OUTLINE, false);
}

fn add_after(scene: &mut Scene, off: P64, va: &[P64], pb: &[P64], y_top: f64, y_bot: f64) {
    let sh = move |q: P64| (q.0 + off.0, q.1 + off.1);
    // One merged tile: both halves in the same fill.
    poly(scene, va.iter().map(|&q| sh(q)).collect(), MERGED_FILL);
    poly(scene, pb.iter().map(|&q| sh(q)).collect(), MERGED_FILL);

    text(scene, sh(centroid(va)), "A+B", 0.46, OUTLINE, true);

    // Matched run, dashed: it is now interior (removed from the boundary).
    scene.push(Item::Polyline {
        points: matched_run(va, &sh),
        stroke: Stroke::dashed(SEAM_COL, 0.05, vec![0.14, 0.10]),
        arrow: None,
    });

    // The two junction vertices carry the repaired angles.
    //   ccw junction (NS+MLEN): x, where a_{l+1} meets b0
    //   cw  junction (NS):      y, where b_{l+1} meets a0
    let jx = va[(NS + MLEN) % N];
    let jy = va[NS % N];
    dot(scene, sh(jx));
    dot(scene, sh(jy));
    // Push the x / y label outward from the merged centroid.
    let cen = centroid(va);
    let out = |p: P64, d: f64| {
        let v = sub(p, cen);
        let l = (v.0 * v.0 + v.1 * v.1).sqrt().max(1e-9);
        sh((p.0 + v.0 / l * d, p.1 + v.1 / l * d))
    };
    text(scene, out(jx, 0.55), "x", 0.5, XY_COL, true);
    text(scene, out(jy, 0.55), "y", 0.5, XY_COL, true);

    let cx = (sh(centroid(va)).0 + sh(centroid(pb)).0) * 0.5;
    text(scene, (cx, y_top), "after: match is interior", 0.42, OUTLINE, true);
    text(scene, (cx, y_bot + 0.45), "x = a4 + b0 - H", 0.34, XY_COL, true);
    text(scene, (cx, y_bot - 0.1), "y = b4 + a0 - H", 0.34, XY_COL, true);
}

fn bbox(pts: &[P64]) -> R64 {
    let mut lo = (f64::INFINITY, f64::INFINITY);
    let mut hi = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &(x, y) in pts {
        lo = (lo.0.min(x), lo.1.min(y));
        hi = (hi.0.max(x), hi.1.max(y));
    }
    (lo, hi)
}

fn main() {
    let dir: PathBuf = std::env::args().nth(1).unwrap_or_else(|| ".".to_string()).into();
    let va = spectre_verts();
    let pb = place_b(&va);

    let all: Vec<P64> = va.iter().chain(&pb).copied().collect();
    let bb = bbox(&all);
    let (w, h) = (bb.1 .0 - bb.0 .0, bb.1 .1 - bb.0 .1);
    let (cxl, cyl) = ((bb.0 .0 + bb.1 .0) * 0.5, (bb.0 .1 + bb.1 .1) * 0.5);
    let slot = w + 3.2;
    let y_top = h / 2.0 + 1.0;
    let y_bot = -h / 2.0 - 1.0;

    let mut scene = Scene::new().with_background(Color::WHITE);
    // BEFORE at x=0, AFTER at x=slot, each centred in its slot.
    add_before(&mut scene, (-cxl, -cyl), &va, &pb, y_top, y_bot);
    add_after(&mut scene, (slot - cxl, -cyl), &va, &pb, y_top, y_bot);

    // "glue" arrow between the panels.
    let gap_x = slot * 0.5;
    scene.draw_arrow((gap_x - 0.7, 0.0), (gap_x + 0.7, 0.0), Stroke::solid(OUTLINE, 0.05), 0.35);
    text(&mut scene, (gap_x, 0.5), "glue", 0.4, OUTLINE, true);

    let bounds = ((-slot / 2.0, y_bot - 0.7), (slot + slot / 2.0, y_top + 0.7));
    let (bw, bh) = (bounds.1 .0 - bounds.0 .0, bounds.1 .1 - bounds.0 .1);
    let px_h = 900u32;
    let px_w = (px_h as f64 * bw / bh).round() as u32;
    let vp = Viewport::rect_for(px_w, px_h, bounds, 8);

    let path = dir.join("glue_repair.svg");
    std::fs::write(&path, scene.to_svg(&vp)).expect("write svg");
    println!("wrote {}", path.display());
}
