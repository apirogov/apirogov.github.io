"""Shared SVG drawing for the figures in this post.

Everything is a walk over a cyclotomic ring: a list of turning angles (or of
absolute directions), turned into points in the complex plane, turned into an
SVG polyline with per-edge rainbow colouring and arrowheads.

Two emitters, because the figures want two different framings:

  emit()  -- crops tightly to the walk plus MARGIN, so a pair of panels can be
             laid side by side and each fills its own box. Optionally draws the
             origin-to-head dashed line with a length label (the seesaw pair).
  boxed() -- fixed square canvas with the walk centred, for a figure that
             stands alone or is deliberately magnified.

Colours are HLS-derived so the edge hue sweeps the wheel once per walk; a shared
hue means "same edge" when two panels show the same word under different
embeddings.
"""
import cmath, colorsys, math

# --- defaults, overridable per call -----------------------------------------
N = 12          # ring: Z[zeta_N]
UNIT = 46.0     # px per unit edge; shared by panels meant to look the same scale
MARGIN = 0.30   # blank border, in units, around a tightly-cropped panel
STROKE = 3.0    # edge width in px
ARROW = (12.0, 5.0)   # arrowhead (length, half-width) in px
FONT = 15       # label size in px


def hexcol(hue):
    """Edge colour at `hue` degrees around the wheel."""
    r, g, b = colorsys.hls_to_rgb((hue % 360) / 360.0, 0.48, 0.72)
    return "#%02x%02x%02x" % (round(r * 255), round(g * 255), round(b * 255))


def walk(turns, g=1, n=None):
    """Points of a turning-angle word, embedded via zeta -> zeta^g."""
    n = n or N
    d, p, pts = 0, 0j, [0j]
    for a in turns:
        d = (d + a) % n
        p += cmath.exp(2j * cmath.pi * ((g * d) % n) / n)
        pts.append(p)
    return pts


def walk_dirs(dirs, g=1, n=None):
    """Points of a word given as absolute step directions, not turns."""
    n = n or N
    p, pts = 0j, [0j]
    for k in dirs:
        p += cmath.exp(2j * cmath.pi * ((g * k) % n) / n)
        pts.append(p)
    return pts


def _edges(out, pts, to, stroke, arrow):
    n = len(pts) - 1
    for i in range(n):
        p0, p1 = to(pts[i]), to(pts[i + 1])
        dx, dy = p1[0] - p0[0], p1[1] - p0[1]
        L = math.hypot(dx, dy)
        ux, uy = dx / L, dy / L
        bx, by = p1[0] - ux * arrow[0], p1[1] - uy * arrow[0]   # arrowhead base
        px, py = -uy, ux                                        # perpendicular
        col = hexcol(i * 360.0 / n)
        out.append(f'<line x1="{p0[0]:.2f}" y1="{p0[1]:.2f}" '
                   f'x2="{bx + ux * 1.5:.2f}" y2="{by + uy * 1.5:.2f}" '
                   f'stroke="{col}" stroke-width="{stroke}" stroke-linecap="round"/>')
        out.append(f'<polygon points="{p1[0]:.2f},{p1[1]:.2f} '
                   f'{bx + px * arrow[1]:.2f},{by + py * arrow[1]:.2f} '
                   f'{bx - px * arrow[1]:.2f},{by - py * arrow[1]:.2f}" fill="{col}"/>')


def _open(W, H, title, fmt=".1f"):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:{fmt}} {H:{fmt}}" '
            f'width="{W:.0f}" height="{H:.0f}" fill="none" role="img" '
            f'aria-label="{title}">', f'<title>{title}</title>']


def emit(pts, path, title, label=None, unit=None, stroke=None, arrow=None, font=None):
    """Tightly-cropped panel. With `label`, draw origin->head dashed and annotate."""
    unit = unit or UNIT
    stroke = STROKE if stroke is None else stroke
    arrow = arrow or ARROW
    font = font or FONT
    xs = [z.real for z in pts]
    ys = [z.imag for z in pts]
    x0, x1 = min(xs) - MARGIN, max(xs) + MARGIN
    y0, y1 = min(ys) - MARGIN, max(ys) + MARGIN
    if label:                      # reserve room for the in-figure text
        x1 += 1.05
        y1 += 0.12
    W, H = (x1 - x0) * unit, (y1 - y0) * unit
    to = lambda z: ((z.real - x0) * unit, (y1 - z.imag) * unit)

    out = _open(W, H, title)
    if label:
        o, head = to(0j), to(pts[-1])
        out.append(f'<line x1="{o[0]:.2f}" y1="{o[1]:.2f}" x2="{head[0]:.2f}" y2="{head[1]:.2f}" '
                   f'stroke="currentColor" stroke-width="1.6" stroke-dasharray="5 4" opacity="0.65"/>')
    _edges(out, pts, to, stroke, arrow)
    if label:
        o, head = to(0j), to(pts[-1])
        out.append(f'<circle cx="{o[0]:.2f}" cy="{o[1]:.2f}" r="4.5" fill="none" '
                   f'stroke="currentColor" stroke-width="1.8"/>')
        out.append(f'<circle cx="{head[0]:.2f}" cy="{head[1]:.2f}" r="4.5" fill="currentColor"/>')
        mx, my = (o[0] + head[0]) / 2, (o[1] + head[1]) / 2
        out.append(f'<text x="{mx + 9:.2f}" y="{my - 9:.2f}" fill="currentColor" '
                   f'font-family="system-ui, sans-serif" font-size="{font}">'
                   f'{label} = {abs(pts[-1]):.3f}</text>')
    else:
        _dots(out, pts, to)
    out.append('</svg>')
    _write(path, out)
    return W


def boxed(pts, path, title, box, unit=None, stroke=None, arrow=None):
    """Fixed square canvas, walk centred. `unit` sets the magnification."""
    unit = unit or UNIT
    stroke = STROKE if stroke is None else stroke
    arrow = arrow or ARROW
    cx = (max(z.real for z in pts) + min(z.real for z in pts)) / 2
    cy = (max(z.imag for z in pts) + min(z.imag for z in pts)) / 2
    W = H = box
    to = lambda z: (W / 2 + (z.real - cx) * unit, H / 2 - (z.imag - cy) * unit)

    out = _open(W, H, title, fmt=".0f")
    _edges(out, pts, to, stroke, arrow)
    _dots(out, pts, to)
    out.append('</svg>')
    _write(path, out)
    return W


def _dots(out, pts, to):
    """Vertex dots, plus a ring on the starting vertex."""
    for z in pts[1:]:
        x, y = to(z)
        out.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="1.9" fill="currentColor" opacity="0.55"/>')
    x, y = to(pts[0])
    out.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.2" fill="none" '
               f'stroke="currentColor" stroke-width="1.8"/>')


def _write(path, out):
    open(path, "w").write("\n".join(out) + "\n")
    print(f"  wrote {path}")
