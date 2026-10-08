#!/usr/bin/env rust-script
//! ```cargo
//! [package]
//! edition = "2021"
//!
//! [dependencies]
//! # Pinned: the figures are rendered through tilezz's own `vis` pipeline,
//! # so they stay in the project's visual language with zero duplicated
//! # drawing logic. `vis` is available with default features (no `raster`
//! # needed for SVG output).
//! tilezz = "=0.1.2"
//! ```
//!
//! Standalone generator for the "Maintaining Valid Snakes" figures. Run it
//! with no repo checkout, just `rust-script`:
//!
//!   rust-script grid_figures.rs [OUTDIR]
//!
//! It writes two SVGs into OUTDIR (default `.`):
//!
//! `unit_grid_locality.svg` -- a 3x3 block of unit cells. One endpoint of
//! every drawn segment lives in the *center* cell (the tinted one); each
//! segment has unit length and reaches into one of the nine cells of the
//! block: eight fan out into the eight neighbour cells, and one stays
//! entirely inside the center cell (the rare "both endpoints share a cell"
//! case). A unit step from the center cell lands no further than an adjacent
//! cell, so the far endpoints never escape the surrounding 3x3 neighbourhood
//! -- the locality the self-intersection check exploits.
//!
//! `unit_grid_neighborhood.svg` -- a single unit segment (red) crossing
//! diagonally between two cells that differ by (1,1), missing their shared
//! corner. Shaded on top is exactly the cell set the snake's
//! `UnitSquareGrid::edge_neighborhood_of` walks for that edge: the union of
//! the two 5-cell `+` crosses around the endpoint cells (8 cells), with the
//! two endpoint cells tinted darker. Gray "existing" edges show what the
//! grid does and does not collect: edges with a vertex in a blue cell are
//! candidates (amber vertices) -- some cross the red query, some do not
//! (collected is necessary, not sufficient) -- while edges wholly in white
//! cells are never even tested. A self-check asserts that no edge crossing
//! the red query can lie entirely outside the blue cells.
//!
//! Both figures share the same visual language (white ground, grey grid,
//! blue cell tints, dark vertex dots).

use std::f64::consts::PI;
use std::path::{Path, PathBuf};

use tilezz::vis::draw::{MarkerStyle, rainbow};
use tilezz::vis::plotutils::{P64, R64};
use tilezz::vis::scene::{Color, Fill, Item, MarkerShape, Scene, Stroke, Viewport};

/// Side of a grid cell, in scene units.
const CELL: f64 = 1.0;
/// Light tint for a cell that belongs to the search region.
const TINT_LIGHT: Color = Color::rgba(70, 130, 210, 26);
/// Stronger tint for a cell that holds a segment endpoint.
const TINT_STRONG: Color = Color::rgba(70, 130, 210, 80);
/// Grid line stroke shared by both figures.
fn grid_stroke() -> Stroke {
    Stroke::solid(Color::rgb(176, 176, 176), 0.012)
}
/// Vertex-dot marker shared by both figures.
fn vertex_dot() -> MarkerStyle {
    MarkerStyle::filled_circle(0.075, Color::rgb(40, 40, 40))
}

fn add(a: P64, b: P64) -> P64 {
    (a.0 + b.0, a.1 + b.1)
}
fn scale(a: P64, s: f64) -> P64 {
    (a.0 * s, a.1 * s)
}
fn dir(theta: f64) -> P64 {
    (theta.cos(), theta.sin())
}

// ---- shared scene helpers ----

/// Fill the unit cell with lower-left integer corner `(cx, cy)`.
fn fill_cell(scene: &mut Scene, cx: i64, cy: i64, color: Color) {
    let (x, y) = (cx as f64, cy as f64);
    scene.push(Item::Polygon {
        points: vec![(x, y), (x + 1.0, y), (x + 1.0, y + 1.0), (x, y + 1.0)],
        fill: Some(Fill::solid(color)),
        stroke: None,
        arrow: None,
    });
}

/// Draw the integer grid lines spanning `[lo, hi]^2`.
fn draw_grid(scene: &mut Scene, lo: f64, hi: f64) {
    let s = grid_stroke();
    let mut x = lo;
    while x <= hi + 1e-9 {
        scene.draw_segment((x, lo), (x, hi), s.clone());
        x += CELL;
    }
    let mut y = lo;
    while y <= hi + 1e-9 {
        scene.draw_segment((lo, y), (hi, y), s.clone());
        y += CELL;
    }
}

// ============================================================
// Figure 1: locality fan
// ============================================================

/// Radius of the hub the eight fan segments start from. Small enough
/// that all eight near endpoints stay well inside the center cell
/// (`0.5 + HUB_R < 1`) and large enough that the hub reads as a ring
/// of distinct anchor dots rather than a single blob.
const HUB_R: f64 = 0.3;
/// Center of the center cell `[0,1]^2`.
const FAN_C: P64 = (0.5, 0.5);

/// One drawn segment: its center-cell endpoint `tail` and its head
/// `head`, exactly one unit apart.
struct Seg {
    tail: P64,
    head: P64,
}

/// Build the nine unit segments of figure 1. Index 0..8 are the eight
/// fan spokes (compass order, starting due east, going CCW); index 8
/// is the in-cell segment.
fn fan_segments() -> Vec<Seg> {
    let mut segs = Vec::with_capacity(9);

    // Eight spokes: tail on the hub circle, head one unit further out
    // along the same radial direction, so the head lands deep inside
    // the neighbour cell in that direction.
    for k in 0..8 {
        let theta = (k as f64) * PI / 4.0;
        let u = dir(theta);
        let tail = add(FAN_C, scale(u, HUB_R));
        let head = add(tail, scale(u, CELL));
        segs.push(Seg { tail, head });
    }

    // In-cell segment: a unit chord centered on the cell center, tilted
    // to 22.5 degrees so it bisects the gap between the 0- and
    // 45-degree spokes and is colinear with none of them. Both
    // endpoints stay inside `[0,1]^2`.
    let u = dir(PI / 8.0);
    segs.push(Seg {
        tail: add(FAN_C, scale(u, -0.5 * CELL)),
        head: add(FAN_C, scale(u, 0.5 * CELL)),
    });

    segs
}

fn fig_locality() -> (Scene, R64) {
    let segs = fan_segments();
    assert_non_crossing(&segs);

    let mut scene = Scene::new().with_background(Color::WHITE);

    // Tint the center cell so "the cell our anchors live in" reads at a
    // glance.
    fill_cell(&mut scene, 0, 0, TINT_LIGHT);
    draw_grid(&mut scene, -1.0, 2.0);

    // The nine unit segments: a color-wheel of rainbow strokes.
    let palette = rainbow(segs.len(), 0.68, 0.46);
    for (i, s) in segs.iter().enumerate() {
        scene.draw_segment(s.tail, s.head, Stroke::solid(palette[i], 0.034));
    }

    // Mark both endpoints of every segment: these are the
    // representative points the snake buckets into grid cells.
    let mut verts: Vec<P64> = Vec::with_capacity(segs.len() * 2);
    for s in &segs {
        verts.push(s.tail);
        verts.push(s.head);
    }
    scene.draw_points(&verts, &vertex_dot());

    (scene, ((-1.0, -1.0), (2.0, 2.0)))
}

// ============================================================
// Figure 2: edge neighborhood
// ============================================================

/// Tail endpoint of the diagonal segment, inside cell `(0,0)`.
const EDGE_A: P64 = (0.8, 0.4);
/// Direction of the diagonal segment. 54 degrees is clearly off the
/// 45-degree cell diagonal, so the unit step lands in cell `(1,1)`
/// while passing well clear of the shared corner `(1,1)` rather than
/// through it.
const EDGE_THETA: f64 = 54.0 * PI / 180.0;

/// The eight cells `UnitSquareGrid::edge_neighborhood_of` walks for a
/// `(dx,dy) = (1,1)` unit edge: the union of the 5-cell `+` crosses
/// around cell `(0,0)` and cell `(1,1)`. (Kept in sync by
/// `assert_neighborhood_cells` below.)
const EDGE_REGION: [(i64, i64); 8] = [
    (-1, 0),
    (0, -1),
    (0, 0),
    (0, 1),
    (1, 0),
    (1, 1),
    (1, 2),
    (2, 1),
];
/// The two cells holding the segment endpoints.
const EDGE_ENDPOINT_CELLS: [(i64, i64); 2] = [(0, 0), (1, 1)];

/// Existing snake edges drawn as context, each as `(tail, theta_deg)`
/// with the head one unit away at that heading. They split into three
/// groups, demonstrating exactly what the grid lookup does and does
/// not collect:
///
/// - **crossers** -- unit edges that actually intersect the red query
///   edge. By the locality guarantee each one has a vertex in a blue
///   cell, so the grid finds it (asserted below). One has both ends in
///   blue; the other is the key edge case -- only *one* vertex in a
///   blue cell (the light cell (1,0), which red merely passes through,
///   not an endpoint cell), the other vertex out in white (2,0). It is
///   still caught, because that single shared cell is enough.
/// - **near misses** -- at least one vertex in a blue cell, so the grid
///   collects them, yet they do *not* cross the query: landing a vertex
///   in a searched cell is necessary, not sufficient. Includes cases
///   with *both* vertices in blue (even one in an endpoint cell) to
///   make the point that proximity alone is not a hit.
/// - **irrelevant** -- both vertices in white cells, never filed in any
///   searched cell, never tested. Mixed orientations, and a couple
///   that span two white cells, so they aren't all near-diagonal
///   look-alikes.
const GRAY_CROSSERS: [(P64, f64); 2] = [
    ((0.5605, 0.9199), -36.0), // steep cross low on the query, both ends blue
    ((1.07, 0.86), 0.0),       // edge case: tail in light cell (1,0), head in white (2,0)
];
const GRAY_NEAR_MISS: [(P64, f64); 4] = [
    ((-0.60, 0.55), 35.0),  // both blue (-1,0)->(0,1), up-left of the query
    ((0.50, -0.70), 100.0), // both blue (0,-1)->(0,0), below the query
    ((2.40, 1.30), 150.0),  // both blue (2,1)->(1,1): in an endpoint cell, still no hit
    ((1.50, 2.30), -10.0),  // one blue (1,2), head out into white (2,2)
];
const GRAY_IRRELEVANT: [(P64, f64); 5] = [
    ((1.15, -0.50), 0.0),   // horizontal, spans two white cells (1,-1)->(2,-1)
    ((-0.50, 1.20), 90.0),  // vertical, spans two white cells (-1,1)->(-1,2)
    ((-0.90, -0.90), 45.0), // diagonal inside white (-1,-1)
    ((2.85, 2.15), 135.0),  // anti-diagonal inside white (2,2)
    ((0.15, 2.15), 50.0),   // oblique inside white (0,2)
];

/// Is `p` filed in one of the blue (searched) cells?
fn in_blue(p: P64) -> bool {
    EDGE_REGION.contains(&cell_of(p))
}

/// Materialize a `(tail, theta_deg)` spec into a unit-length segment.
fn unit_seg((tail, theta_deg): (P64, f64)) -> (P64, P64) {
    (tail, add(tail, dir(theta_deg.to_radians())))
}

fn fig_neighborhood() -> (Scene, R64) {
    let red = unit_seg((EDGE_A, EDGE_THETA.to_degrees()));
    let (a, b) = red;
    assert_neighborhood_geometry(a, b);

    let crossers: Vec<(P64, P64)> = GRAY_CROSSERS.iter().map(|&s| unit_seg(s)).collect();
    let near_miss: Vec<(P64, P64)> = GRAY_NEAR_MISS.iter().map(|&s| unit_seg(s)).collect();
    let irrelevant: Vec<(P64, P64)> = GRAY_IRRELEVANT.iter().map(|&s| unit_seg(s)).collect();
    assert_gray_classification(&crossers, &near_miss, &irrelevant, red);

    let mut scene = Scene::new().with_background(Color::WHITE);

    // Shade the search region: light tint for every neighbourhood cell,
    // a stronger tint for the two cells that actually hold an endpoint.
    for &(cx, cy) in &EDGE_REGION {
        let tint = if EDGE_ENDPOINT_CELLS.contains(&(cx, cy)) {
            TINT_STRONG
        } else {
            TINT_LIGHT
        };
        fill_cell(&mut scene, cx, cy, tint);
    }

    // The region's cells span x,y in {-1,0,1,2}, so the 4x4 block
    // [-1,3]^2 frames it exactly -- no empty padding cells.
    draw_grid(&mut scene, -1.0, 3.0);

    // Existing edges in gray; the red query edge on top so it stays the
    // focus.
    let gray = Stroke::solid(Color::rgb(125, 125, 125), 0.022);
    for &(p, q) in crossers.iter().chain(&near_miss).chain(&irrelevant) {
        scene.draw_segment(p, q, gray.clone());
    }
    scene.draw_segment(a, b, Stroke::solid(Color::rgb(205, 35, 45), 0.045));

    // Vertex markers. A vertex in a blue cell is "collected" -> amber
    // highlight; a vertex in a white cell is never looked at -> small
    // plain grey dot. The red query's own endpoints stay black.
    let mut highlighted: Vec<P64> = Vec::new();
    let mut plain: Vec<P64> = Vec::new();
    for &(p, q) in crossers.iter().chain(&near_miss).chain(&irrelevant) {
        for v in [p, q] {
            if in_blue(v) {
                highlighted.push(v);
            } else {
                plain.push(v);
            }
        }
    }
    scene.draw_points(&plain, &MarkerStyle::filled_circle(0.05, Color::rgb(155, 155, 155)));
    scene.draw_points(&highlighted, &highlight_dot());
    scene.draw_points(&[a, b], &vertex_dot());

    (scene, ((-1.0, -1.0), (3.0, 3.0)))
}

/// Marker for a "collected" vertex: an amber dot ringed in dark so it
/// reads clearly on the blue cell tint.
fn highlight_dot() -> MarkerStyle {
    MarkerStyle {
        shape: MarkerShape::Circle,
        size: 0.09,
        fill: Some(Fill::solid(Color::rgb(240, 150, 25))),
        stroke: Some(Stroke::solid(Color::rgb(60, 40, 0), 0.014)),
    }
}

// ---- self-verification ----

fn orient(a: P64, b: P64, c: P64) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}

fn on_seg(a: P64, b: P64, p: P64) -> bool {
    p.0 <= a.0.max(b.0) + 1e-12
        && p.0 >= a.0.min(b.0) - 1e-12
        && p.1 <= a.1.max(b.1) + 1e-12
        && p.1 >= a.1.min(b.1) - 1e-12
}

/// Closed-segment intersection on f64 with a small epsilon. Returns
/// true if `ab` and `cd` share any point at all (including touching).
fn segments_meet(a: P64, b: P64, c: P64, d: P64) -> bool {
    let eps = 1e-9;
    let (o1, o2, o3, o4) = (
        orient(a, b, c),
        orient(a, b, d),
        orient(c, d, a),
        orient(c, d, b),
    );
    if o1 * o2 < -eps && o3 * o4 < -eps {
        return true; // proper crossing
    }
    (o1.abs() <= eps && on_seg(a, b, c))
        || (o2.abs() <= eps && on_seg(a, b, d))
        || (o3.abs() <= eps && on_seg(c, d, a))
        || (o4.abs() <= eps && on_seg(c, d, b))
}

/// Panic if any two segments of figure 1 touch -- the figure's whole
/// claim is that these unit steps are mutually non-crossing.
fn assert_non_crossing(segs: &[Seg]) {
    for i in 0..segs.len() {
        for j in (i + 1)..segs.len() {
            assert!(
                !segments_meet(segs[i].tail, segs[i].head, segs[j].tail, segs[j].head),
                "fig1: segments {i} and {j} cross -- illustration is not clean"
            );
        }
    }
    println!(
        "fig1 verified: {} unit segments are pairwise non-crossing",
        segs.len()
    );
}

/// Distance from point `p` to segment `ab`.
fn point_seg_dist(a: P64, b: P64, p: P64) -> f64 {
    let ab = (b.0 - a.0, b.1 - a.1);
    let len2 = ab.0 * ab.0 + ab.1 * ab.1;
    let t = if len2 <= 1e-18 {
        0.0
    } else {
        (((p.0 - a.0) * ab.0 + (p.1 - a.1) * ab.1) / len2).clamp(0.0, 1.0)
    };
    let proj = (a.0 + t * ab.0, a.1 + t * ab.1);
    ((p.0 - proj.0).powi(2) + (p.1 - proj.1).powi(2)).sqrt()
}

fn cell_of(p: P64) -> (i64, i64) {
    (p.0.floor() as i64, p.1.floor() as i64)
}

/// Verify figure 2's segment really is a unit edge from cell `(0,0)` to
/// cell `(1,1)` that clears the shared corner, and that `EDGE_REGION`
/// is exactly the union of the two 5-cell `+` crosses around the
/// endpoint cells.
fn assert_neighborhood_geometry(a: P64, b: P64) {
    assert_eq!(cell_of(a), (0, 0), "fig2: A not in cell (0,0)");
    assert_eq!(cell_of(b), (1, 1), "fig2: B not in cell (1,1)");
    let len = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
    assert!((len - 1.0).abs() < 1e-9, "fig2: segment not unit length: {len}");
    let clearance = point_seg_dist(a, b, (1.0, 1.0));
    assert!(
        clearance > 0.1,
        "fig2: segment passes too close to shared corner (1,1): {clearance}"
    );
    assert_neighborhood_cells();
    println!("fig2 verified: unit edge (0,0)->(1,1), corner clearance {clearance:.3}, 8-cell region");
}

/// Recompute the union of the two `+` crosses around the endpoint cells
/// and confirm it matches `EDGE_REGION`, so the shaded cells can never
/// silently drift from `edge_neighborhood_of`'s actual behaviour.
fn assert_neighborhood_cells() {
    let cross = |(cx, cy): (i64, i64)| {
        [
            (cx - 1, cy),
            (cx, cy - 1),
            (cx, cy),
            (cx, cy + 1),
            (cx + 1, cy),
        ]
    };
    let mut want: Vec<(i64, i64)> = EDGE_ENDPOINT_CELLS
        .iter()
        .flat_map(|&c| cross(c))
        .collect();
    want.sort_unstable();
    want.dedup();
    let mut have = EDGE_REGION.to_vec();
    have.sort_unstable();
    assert_eq!(want, have, "fig2: EDGE_REGION != union of the two + crosses");
}

/// Verify the gray context edges actually demonstrate what they claim,
/// so the figure can never tell a lie about the locality guarantee:
///
/// - the *core* invariant: **every** gray edge that intersects the red
///   query edge has a vertex in a blue cell -- i.e. searching the blue
///   cells alone misses no real intersection;
/// - crossers cross and are collected; near misses are collected but
///   do not cross; irrelevant edges neither cross nor land any vertex
///   in a blue cell.
fn assert_gray_classification(
    crossers: &[(P64, P64)],
    near_miss: &[(P64, P64)],
    irrelevant: &[(P64, P64)],
    red: (P64, P64),
) {
    let (a, b) = red;
    let crosses = |&(p, q): &(P64, P64)| segments_meet(p, q, a, b);
    let collected = |&(p, q): &(P64, P64)| in_blue(p) || in_blue(q);

    // Core invariant across *all* gray edges.
    for s in crossers.iter().chain(near_miss).chain(irrelevant) {
        if crosses(s) {
            assert!(
                collected(s),
                "fig2: gray edge {s:?} crosses the query but lands no vertex \
                 in a blue cell -- the blue-box search would miss it"
            );
        }
    }
    for s in crossers {
        assert!(crosses(s), "fig2: a 'crosser' {s:?} does not cross the query");
        assert!(collected(s), "fig2: a 'crosser' {s:?} is not collected");
    }
    for s in near_miss {
        assert!(!crosses(s), "fig2: a 'near miss' {s:?} actually crosses");
        assert!(collected(s), "fig2: a 'near miss' {s:?} is not collected");
    }
    for s in irrelevant {
        assert!(!crosses(s), "fig2: an 'irrelevant' edge {s:?} crosses the query");
        assert!(
            !collected(s),
            "fig2: an 'irrelevant' edge {s:?} has a vertex in a blue cell"
        );
    }
    println!(
        "fig2 verified: {} crossers (all collected), {} near misses (collected, no hit), \
         {} irrelevant (never collected); no crossing edge escapes the blue cells",
        crossers.len(),
        near_miss.len(),
        irrelevant.len()
    );
}

// ---- output ----

fn write_figure(dir: &Path, name: &str, scene: &Scene, bounds: R64) {
    let vp = Viewport::square_for(800, bounds, 12);
    let svg_path = dir.join(format!("{name}.svg"));
    std::fs::write(&svg_path, scene.to_svg(&vp)).expect("write svg");
    println!("wrote {}", svg_path.display());
}

fn main() {
    let dir: PathBuf = std::env::args().nth(1).unwrap_or_else(|| ".".to_string()).into();

    let (loc, loc_bounds) = fig_locality();
    write_figure(&dir, "unit_grid_locality", &loc, loc_bounds);

    let (nbr, nbr_bounds) = fig_neighborhood();
    write_figure(&dir, "unit_grid_neighborhood", &nbr, nbr_bounds);
}
