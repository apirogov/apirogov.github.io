from dataclasses import dataclass

@dataclass
class Point:
    x: float
    y: float

def wedge(a: Point, b: Point) -> float:
    """Return wedge product of two vectors."""
    return (a.x * b.y) - (a.y * b.x)

def dot(a: Point, b: Point) -> float:
    """Return dot product of two vectors."""
    return (a.x * b.x) + (a.y * b.y)

def line_through(p: Point, q: Point) -> tuple[float, float, float]:
    """Return (a, b, c) for line ax + by + c = 0 though points P and Q."""
    return (p.y - q.y, q.x - p.x, p.x * q.y - q.x * p.y)

def is_on_line(p: Point, l: tuple[float, float, float]) -> bool:
    """Return whether point P is on the line given by the coefficients."""
    (a, b, c) = l
    return a * p.x + b * p.y + c == 0

def is_between(a: Point, b: Point, c: Point) -> bool:
    """Return whether point B is strictly between A and C."""
    ba = Point(a.x-b.x, a.y-b.y)
    cb = Point(b.x-c.x, b.y-c.y)
    return wedge(ba, cb) == 0 and dot(ba, cb) > 0

def is_ccw(a: Point, b: Point, c: Point) -> bool:
    """Return whether line segment AB is counter-clockwise from AC."""
    ab = Point(b.x-a.x, b.y-a.y)
    ac = Point(c.x-a.x, c.y-a.y)
    return wedge(ab, ac) > 0

def intersect(a: Point, b: Point, c: Point, d: Point) -> bool:
    """Return true if line segments AB and CD intersect."""
    if a == c or a == d or b == c or b == d:
        return True
    l_ab = line_through(a, b)
    if is_on_line(c, l_ab) and is_on_line(d, l_ab):
        return is_between(a, c, b) or is_between(a, d, b)
    return is_ccw(a,c,d) != is_ccw(b,c,d) and is_ccw(a,b,c) != is_ccw(a,b,d)

# Some test cases using this point layout: E F
#                                          A B C D
A, B, C, D, E, F = Point(0, 0), Point(1, 0), Point(2, 0), Point(3, 0), Point(0, 1), Point(1, 1)
print(f"intersect ({A}, {B}) and ({E}, {F}): {intersect(A, B, E, F)} (exp: False)")  # case 1
print(f"intersect ({A}, {B}) and ({C}, {D}): {intersect(A, B, C, D)} (exp: False)")  # case 2.1
print(f"intersect ({A}, {B}) and ({B}, {C}): {intersect(A, B, B, C)} (exp: True)")  # case 2.2
print(f"intersect ({A}, {C}) and ({B}, {D}): {intersect(A, C, B, D)} (exp: True)")  # case 2.3a
print(f"intersect ({A}, {D}) and ({B}, {C}): {intersect(A, D, B, C)} (exp: True)")  # case 2.3b
print(f"intersect ({A}, {B}) and ({F}, {D}): {intersect(A, B, F, D)} (exp: False)")  # case 3.1
print(f"intersect ({A}, {F}) and ({B}, {E}): {intersect(A, F, B, E)} (exp: True)")  # case 3.2
print(f"intersect ({A}, {C}) and ({B}, {F}): {intersect(A, C, B, F)} (exp: True)")  # case 3.3
print(f"intersect ({A}, {E}) and ({A}, {B}): {intersect(A, E, A, B)} (exp: True)")  # case 3.4
