#!/usr/bin/env python3
"""Figures for "Geometry as a finite automaton".

Two panels. The first draws a walk over the rhombic cell grid: every vertex is
labelled with the integer coordinates of the cell it falls in, and coloured by
its *state* -- where inside that cell it sits. The second zooms one cell to show
the three distinct states as three dots.

Drawn over Z[zeta_10] rather than Z[zeta_12], because a rhombus with only three
states per cell makes the point in one picture; the square cell of Z[zeta_12]
would need a bigger grid to show the same thing.

Config below: N/K choose the ring and the second basis direction, WORD is the
walk (absolute step directions), COLS is one colour per state.
"""
import cmath, math

N, K = 10, 3
u, v = complex(1, 0), cmath.exp(2j*cmath.pi*K/N)
WORD = [1, 9, 7, 0, 1, 9, 7]
COLS = ["#d14747", "#47d147", "#4747d1"]          # same rainbow hues as the other figures
UNIT, STROKE = 74.0, 2.6

det = u.real*v.imag - u.imag*v.real
def coords(p): return ((p.real*v.imag - p.imag*v.real)/det, (u.real*p.imag - u.imag*p.real)/det)
def fold(p):
    a, b = coords(p); ca, cb = round(a), round(b)
    return (ca, cb), (round(a-ca, 6), round(b-cb, 6))

pts, p = [0j], 0j
for d in WORD:
    p += cmath.exp(2j*cmath.pi*d/N); pts.append(p)
cells  = [fold(q)[0] for q in pts]
states = [fold(q)[1] for q in pts]
order  = sorted(set(states))
cidx   = [order.index(s) for s in states]

xs = [q.real for q in pts]; ys = [q.imag for q in pts]
x0, x1 = min(xs)-1.15, max(xs)+1.15
y0, y1 = min(ys)-1.15, max(ys)+1.15
W, H = (x1-x0)*UNIT, (y1-y0)*UNIT
to = lambda q: ((q.real-x0)*UNIT, (y1-q.imag)*UNIT)

out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.0f} {H:.0f}" width="{W:.0f}" '
       f'height="{H:.0f}" fill="none" role="img" aria-label="A unit-step walk on the rhombic cell '
       f'grid of Z[zeta_10]; each vertex is coloured by its position within its own cell.">',
       '<title>Cells and states</title>']
# cell walls: the two line families at half-integer coordinates
rng = range(-12, 13)
for m in rng:
    for (fixed, along) in ((True, False), (False, True)):
        a0, b0 = (m+0.5, -12) if fixed else (-12, m+0.5)
        a1, b1 = (m+0.5, 12) if fixed else (12, m+0.5)
        p0, p1 = to(a0*u + b0*v), to(a1*u + b1*v)
        out.append(f'<line x1="{p0[0]:.1f}" y1="{p0[1]:.1f}" x2="{p1[0]:.1f}" y2="{p1[1]:.1f}" '
                   f'stroke="currentColor" stroke-width="0.7" opacity="0.28"/>')
out.append(f'<rect x="0" y="0" width="{W:.0f}" height="{H:.0f}" fill="none"/>')
# the walk
for i in range(len(pts)-1):
    a, b = to(pts[i]), to(pts[i+1])
    dx, dy = b[0]-a[0], b[1]-a[1]; L = math.hypot(dx, dy); ux, uy = dx/L, dy/L
    bx, by = b[0]-ux*13, b[1]-uy*13; px, py = -uy, ux
    out.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{bx:.1f}" y2="{by:.1f}" '
               f'stroke="currentColor" stroke-width="{STROKE}" stroke-linecap="round" opacity="0.75"/>')
    out.append(f'<polygon points="{b[0]:.1f},{b[1]:.1f} {bx+px*5:.1f},{by+py*5:.1f} '
               f'{bx-px*5:.1f},{by-py*5:.1f}" fill="currentColor" opacity="0.75"/>')
# vertices coloured by state, cell coordinates in small grey text
for q, c, ci in zip(pts, cells, cidx):
    x, y = to(q)
    out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="{COLS[ci]}" stroke="currentColor" stroke-width="1.4"/>')
    cc = to(c[0]*u + c[1]*v)
    # the label belongs to the cell, but must clear the vertex dot, which sits
    # off-centre by the state -- for a state below the centre they collide.
    ly = max(cc[1], y) + 24
    out.append(f'<text x="{cc[0]:.1f}" y="{ly:.1f}" fill="currentColor" opacity="0.55" '
               f'font-family="system-ui, sans-serif" font-size="13" text-anchor="middle">'
               f'({c[0]},{c[1]})</text>')
out.append('</svg>')
open("fig-automaton-walk.svg", "w").write("\n".join(out)+"\n")

# panel B: the base cell alone, with the distinct states as dots
S = 260.0; pad = 34.0
cx, cy = S/2, S/2
sc = (S - 2*pad) / max(abs(u+v), abs(u-v))
tob = lambda q: (cx + q.real*sc, cy - q.imag*sc)
corners = [(-0.5*u - 0.5*v), (0.5*u - 0.5*v), (0.5*u + 0.5*v), (-0.5*u + 0.5*v)]
poly = " ".join(f"{tob(c)[0]:.1f},{tob(c)[1]:.1f}" for c in corners)
b = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {S:.0f} {S:.0f}" width="{S:.0f}" '
     f'height="{S:.0f}" fill="none" role="img" aria-label="One rhombic cell containing the three '
     f'distinct states of the walk as coloured dots.">', '<title>The states inside one cell</title>',
     f'<polygon points="{poly}" fill="none" stroke="currentColor" stroke-width="1.6" opacity="0.55"/>']
for i, s in enumerate(order):
    q = s[0]*u + s[1]*v; x, y = tob(q)
    b.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="9" fill="{COLS[i]}" stroke="currentColor" stroke-width="1.4"/>')
b.append('</svg>')
open("fig-automaton-cell.svg", "w").write("\n".join(b)+"\n")

print(f"walk {WORD}: {len(pts)} vertices in {len(set(cells))} cells, {len(order)} distinct states")
for i, s in enumerate(order):
    who = [j for j, c in enumerate(cidx) if c == i]
    print(f"   {COLS[i]}  state {s}  at vertices {who}  cells {[cells[j] for j in who]}")
