#!/usr/bin/env rust-script
//! ```cargo
//! [dependencies]
//! tilezz = { path = "..", features = ["cli"] }
//! ```
//!
//! Render one periodic-tiling patch from a freshly minted certificate.
//!
//!   fig_tilings.rs <out.svg|-> <word> <conway|bn|iso|aniso|torus> [radius] [cap] [q|RRGGBB]
//!
//! The detector is chosen explicitly rather than letting the classification
//! cascade pick: the figure set wants one panel per criterion, and most tiles
//! satisfy several. `radius` and `cap` bound how far the orbit is grown, and the
//! last argument overrides colouring (see the block near the bottom).
//!
//! Also runnable as a cargo example: drop it in examples/ and
//!   cargo run --release --example fig_tilings --features cli -- <args>
//!
//! Requires three items to be `pub` in the library (they are peers of the
//! already-public cert constructors): mint::torus_cert, mint::translation_cert
//! and tiling::DETECT_ORBIT_CAP.
use tilezz::classify::aniso::tiles_anisohedral_restart;
use tilezz::classify::cascade::{AcceptBounds, PeriodicVia};
use tilezz::classify::cert::PeriodicCert;
use tilezz::classify::conway::{bn_criterion, build_tiling};
use tilezz::classify::grow::replay_placements;
use tilezz::classify::isohedral::isohedral_tiling;
use tilezz::classify::mint::{cert_from_cluster, cert_from_tiling, torus_cert, translation_cert};
use tilezz::classify::render::{orientation_colors, scene_from_placements};
use tilezz::vis::draw::rainbow;
use tilezz::classify::tiling::{DETECT_ORBIT_CAP, ORBIT_RADIUS_FACTOR};
use tilezz::cyclotomic::ZZ12;
use tilezz::geom::iso::Iso;
use tilezz::geom::patch::boundary_vertices;
use tilezz::geom::rat::Rat;
use tilezz::vis::scene::{Color, Viewport};

type Z = ZZ12;

fn mint(seq: &[i8], which: &str) -> Option<PeriodicCert> {
    let base = Rat::<Z>::from_slice_trusted(seq);
    let b = AcceptBounds::default();
    let radius = ORBIT_RADIUS_FACTOR * seq.len() as f64;
    match which {
        "conway" => cert_from_tiling(&base, &build_tiling::<Z>(seq, radius, DETECT_ORBIT_CAP)?, PeriodicVia::Conway),
        "bn" => bn_criterion::<Z>(seq).into_iter().find_map(|(v1, v2)| translation_cert(&base, v1, v2)),
        "iso" => cert_from_tiling(&base, &isohedral_tiling::<Z>(seq, radius, DETECT_ORBIT_CAP, b.iso_builds)?, PeriodicVia::Isohedral),
        "aniso" => {
            let w = tiles_anisohedral_restart::<Z>(seq, b.aniso_kmax, b.aniso_cap, b.aniso_budget, b.aniso_restarts)?;
            cert_from_cluster::<Z>(&base, &w.build, &w.tiling, PeriodicVia::Anisohedral(w.k))
        }
        "torus" => (2..=7).rev().find_map(|c| torus_cert::<Z>(&base, 32, c)),
        _ => panic!("unknown detector {which}"),
    }
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let seq: Vec<i8> = a[1].split(',').map(|s| s.trim().parse().unwrap()).collect();
    let radius: f64 = a.get(3).map_or(5.5, |s| s.parse().unwrap());
    let cap: usize = a.get(4).map_or(96, |s| s.parse().unwrap());
    let base = Rat::<Z>::from_slice_trusted(&seq);
    let Some(cert) = mint(&seq, &a[2]) else { eprintln!("no cert for {}", a[2]); return };
    assert!(cert.verify(&base), "minted cert must verify");

    // grow() orbits the META-tile; the meta-tile itself is `k` base copies,
    // replayed from cert.build. The tiling is every (meta o inner) composite.
    let inner: Vec<Iso<Z>> = if cert.build.is_empty() {
        vec![Iso::id()]
    } else {
        replay_placements(&base, &cert.build).expect("replay failed")
    };
    let metas = cert.grow(&base, radius, cap).expect("grow failed");
    let placements: Vec<Iso<Z>> = metas
        .iter()
        .flat_map(|m| inner.iter().map(move |b| m.after(b)))
        .collect();

    let verts = boundary_vertices::<Z>(base.seq());
    // arg 6: "q" = qualitative remap (for panels where the rainbow puts two
    // near-identical hues side by side), RRGGBB = one flat fill, else the rainbow.
    // Same 12-hue rainbow as the other panels, walked so that consecutive picks
    // land on opposite sides of the wheel: red, green, blue, yellow, cyan, magenta.
    const WHEEL: [usize; 12] = [0, 4, 8, 2, 6, 10, 1, 5, 9, 3, 7, 11];
    let pal = rainbow(12, 0.60, 0.55);
    let raw = orientation_colors::<Z>(&placements);
    let colors: Vec<Color> = match a.get(5).map(String::as_str) {
        Some("q") => {
            let mut seen: Vec<Color> = Vec::new();
            for c in &raw {
                if !seen.contains(c) {
                    seen.push(*c);
                }
            }
            raw.iter()
                .map(|c| pal[WHEEL[seen.iter().position(|s| s == c).unwrap() % 12]])
                .collect()
        }
        Some(hex) => {
            let v = u32::from_str_radix(hex, 16).expect("fill must be RRGGBB hex");
            vec![Color::rgb((v >> 16) as u8, (v >> 8) as u8, v as u8); placements.len()]
        }
        None => raw,
    };
    let scene = scene_from_placements::<Z>(&verts, &placements, &colors, None, 0.05);
    let bb = scene.auto_bounds().expect("empty scene");
    if a[0] != "-" {
        std::fs::write(&a[0], scene.to_svg(&Viewport::square_for(720, bb, 14))).unwrap();
    }
    eprintln!("via {:?}: {} meta x {} inner = {} tiles -> {}", cert.via, metas.len(), inner.len(), placements.len(), a[0]);
}
