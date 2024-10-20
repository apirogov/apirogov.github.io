"""Naive implementation of complex integer rings in Python."""
# Copyright (C) 2024 Anton Pirogov, licensed under the MIT License
from __future__ import annotations

from typing import TypeVar, Generic, Type
from fractions import Fraction

from typing import Any
Self = Any

T = TypeVar("T")

class SymNum(Generic[T]):
    """Number from a complex integer ring with discrete rotation steps.

    This base class is just a wrapper around an integer, providing machinery for non-trivial subclasses.

    A number is represented as a vector of rational numbers,
    each component denoting the coefficient multiplied by a symbolic root term,
    such as i, sqrt(2), etc. (technically, we use it only for square roots
    of positive and negative integers).

    Note that there is no meaningful notion of a classical angle or distance between two
    points, because the result usually is not itself part of the ring or not uniquely
    defined inside the ring. Those values can be computed after conversion to the
    conventional complex numbers.

    See https://pirogov.de/blog/perfect-precision-2d-geometry-complex-integers/
    """
    NUM_TYPE: Type[T] = Fraction

    D_SQ: list[T] = [1]
    """Squares of units of the symbolic root dimensions that are represented by the number."""

    R: int = 1
    """Rotation unit in the ring is 1/R of a full turn."""

    L: int = 1
    """Dimension unit scaling factor is 1/L (needed to make all coefficients integers)."""

    _CCW: list[int] = [1]
    """Coefficients of normal form representation of the rotation unit."""

    # ----

    vec: list[T]
    """Vector with coefficients of the number."""

    @staticmethod
    def _mul(x: list[T], y: list[T]) -> list[T]:
        """Symbolic dimension multiplication table."""
        (a,), (b,) = x, y
        return [a * b]

    @classmethod
    def _ccw(cls) -> list[Fraction]:
        """Smallest unit rotation vector (generates all rotations)."""
        return list(map(lambda x: x * Fraction(1, cls.L), cls._CCW))

    def __init__(self, val: list[T] | int | Self[T] | None = None):
        if val is None:
            self.vec = [self.NUM_TYPE() for _ in range(len(self.D_SQ))]
        elif isinstance(val, SymNum):
            self.vec = list(val.vec)
        elif isinstance(val, int):
            one = self.NUM_TYPE(1) if not issubclass(self.NUM_TYPE, Fraction) else self.NUM_TYPE(1,1)
            self.vec = ([val*one] + [self.NUM_TYPE() for _ in range(len(self.D_SQ)-1)])
        else:
            assert len(val) == len(self.D_SQ)
            self.vec = list(val)

    @classmethod
    def zero(cls) -> Self:
        """Zero point."""
        return cls()

    @classmethod
    def one(cls) -> Self:
        """Unit pointing toward positive side of the real axis."""
        return cls(1)

    @classmethod
    def ccw(cls) -> Self:
        """Vector of smallest supported rotation."""
        return cls(cls._ccw())

    @classmethod
    def unit(cls, k: int = 0) -> Self:
        """Get a unit vector oriented into the desired direction."""
        return cls.one() * (cls.ccw()**k)

    # ----

    def __neg__(self) -> Self:
        return type(self)(list(map(lambda x: -x, self.vec)))

    def __add__(self, other: Self) -> Self:
        return type(self)(list(map(lambda a, b: a+b, self.vec, other.vec)))

    def __radd__(self, other: Self):
        return self + other 

    def __sub__(self, other: Self) -> Self:
        return self + (-other)

    def __mul__(self, other: int | Self) -> Self:
        if isinstance(other, int):
            # scalar multiplication
            return type(self)(list(map(lambda x: other*x, self.vec)))
        # special vector multiplication
        return type(self)(list(map(lambda x: x, self._mul(self.vec, other.vec))))

    def __rmul__(self, other):
        return self * other

    def __pow__(self, k: int) -> Self:
        if k < 0:
            k %= self.R
        result = type(self).one()
        for _ in range(k):
            result *= self
        return result


    def __iadd__(self, other: Self) -> Self:
        self.vec = (self+other).vec
        return self

    def __isub__(self, other: Self) -> Self:
        self.vec = (self-other).vec
        return self

    def __imul__(self, other: Self) -> Self:
        self.vec = (self*other).vec
        return self

    def __ipow__(self, k: int):
        self.vec = (self**k).vec
        return self

    # ----

    def __complex__(self) -> complex:
        """Convert to complex number, i.e. multiply and sum the encoded symbolic terms."""
        vec = self.vec
        if isinstance(vec[0], SymNum):
            vec = list(map(complex, self.vec))
        return sum(map(lambda c,d: c * complex(d)**(1/2), vec, self.D_SQ))

    def xy(self) -> tuple[float, float]:
        """Returns (Re(val), Im(val))."""
        c = complex(self)
        return (c.real, c.imag)

    def x(self) -> float:
        """Return Re(val)."""
        return self.xy()[0]

    def y(self) -> float:
        """Return Im(val)."""
        return self.xy()[1]

    # ----

    def __getitem__(self, ix):
        return self.vec[ix]

    def __iter__(self):
        return iter(self.vec)

    def __eq__(self, other) -> bool:
        for x, y in zip(self, other):
            if x != y:
                return False
        return True


    def __hash__(self) -> int:
        return hash(tuple(self.vec))


    def __bool__(self) -> bool:
        return any(map(bool, self.vec))


    def __repr__(self) -> str:
        return repr(self.vec).replace("array", type(self).__name__)

    def __str__(self) -> str:
        """Pretty-print number in human-readable algebraic form."""

        def fmt_d(d: int) -> str:
            if d == 1:
                return ""
            elif d == -1:
                return "1j"
            else:
                d_str = f"sqrt({abs(d)})"
                i_str ='1j*' if  d < 0 else ''
                return f"{i_str}{d_str}"

        def fmt_v(v) -> str:
            scaled = v * self.L
            if isinstance(v, (int, Fraction)):
                return str(int(scaled))
            else:
                tmp = str(scaled)
                return f"({tmp})" if "*" in tmp or "-" in tmp else tmp

        def fmt_term(v, d):
            v_str, d_str = fmt_v(v), fmt_d(d)
            if not d_str:
                return v_str
            if v_str == "1":
                return d_str
            return f"{v_str}*{d_str}"

        terms =[ fmt_term(v, d) for v, d in zip(self.vec, self.D_SQ) if v and d ]
        ret = "+".join(terms)

        if not ret:
            return "0"

        if "+" in ret:
            ret = f"({ret})"

        if self.L == 1:
            return ret
        return f"{ret}/{self.L}"


class GaussInt(SymNum):
    """Numbers with 90 degree rotations (Gaussian integers).

    This is the only instance where the resulting dimensions are orthogonal
    and exactly match the two axes of the complex plane.
    """
    R = 4
    D_SQ = [1, -1]
    L = 1
    _CCW = [0, 1]

    @staticmethod
    def _mul(x, y):
        (a, b), (c, d) = x, y
        return [
            a*c - b*d,
            a*d + b*c
        ]

class EisensteinInt(SymNum):
    """Numbers with 60 degree rotations (Eisenstein).

    Note that this is equivalent to using 120 degree rotations as basis, but
    the 120 degree rotation basis is improper, because it is not minimal
    (one can produce 60 deg. rotations from linear combinations of 120 degree rotations).

    Can also be represented in symbolic root normal form (see EisenClockInt).
    """
    R = 6
    D_SQ = [1, -3]
    L = 2
    _CCW = [1, 1]

    @staticmethod
    def _mul(x, y):
        (a, b), (c, d) = x, y
        return [
            a*c - 3*b*d,
            a*d +   b*c
        ]


class CSymNum(SymNum[GaussInt]):
    """A complex symbolic root number with Gaussian Integers as coefficients.

    This base class is just a wrapper around a GaussInt, providing machinery
    to represent different complex integer rings in symbolic root normal form, i.e

    1/l(c_i*sqrt(d_i)) with c_i being Gaussian integers, l in Z, d_i in R+.
    """
    NUM_TYPE = GaussInt
    _CCW = [(0, 1)]

    @classmethod
    def _ccw(cls) -> list[Fraction]:
        def to_sym(xy: tuple[int, int]):
            x, y = xy
            return cls.NUM_TYPE([Fraction(x, cls.L), Fraction(y, cls.L)])

        return list(map(to_sym, cls._CCW))



class CompassInt(CSymNum):
    """Numbers with 45 degree rotations (compass integers)."""
    R = 8
    D_SQ = [1, 2]
    L = 2
    _CCW = [(0,0), (1,1)]

    @staticmethod
    def _mul(x, y):
        # Dimension multiplication matrix (for Z[i]-valued vectors):
        #    c d
        # a [1 s]
        # b [s 2]
        # where s = sqrt(2)
        (a, b), (c, d) = x, y
        return [
            a*c + 2*b*d,
            a*d +   b*c,
        ]


class PenroseInt(CSymNum):
    """Numbers with 36 degree rotations (Penrose integers)."""
    PENTA = 2*(5-5**(1/2))  # common term for needed symbolic root dimensions

    R = 10
    D_SQ = [1, 5, PENTA, 5*PENTA]
    L = 8
    _CCW = [(2,0), (2,0), (0,2), (0,0)]

    # condensed version (if vectors were Z[i]-valued):
    # x=sq5, y=2(5-x)
    #    e     f      g        h
    # a [1   ,  x   ,   y    ,   xy   ]
    # b [ x  , 5    ,  xy    ,  5 y   ]
    # c [  y ,  xy  , 10-2x  , 10(x-1)]
    # d [ xy , 5 y  , 10(x-1), 10(5-x)]
    # where x = sqrt(5), y = sqrt(2*(5-sqrt(5)))
    @staticmethod
    def _mul(x, y):
        (a, b, c, d), (e, f, g, h) = x, y
        return [
            a*e + 5*b*f + 10*(c*g + 5*d*h - c*h - d*g),
            a*f + b*e - 2*c*g + 10*(c*h + d*g - d*h),
            a*g + 5*(b*h + d*f) + c*e,
            a*h + b*g + c*f + d*e,
        ]


class ClockInt(CSymNum):
    """Numbers with 30 degree rotations (clock integers)"""
    R = 12
    D_SQ = [1, 3]
    L = 2
    _CCW = [(0, 1),(1, 0)]

    @staticmethod
    def _mul(x, y):
        # Dimension multiplication matrix (for Z[i]-valued vectors):
        #    c d
        # a [1 s]
        # b [s 3]
        # where s = sqrt(3)
        (a, b), (c, d) = x, y
        return [
            a*c + 3*b*d,
            a*d +   b*c,
        ]


class EisenClockInt(ClockInt):
    """Normal-form representation of Eisenstein numbers (just for illustration purposes)."""
    R = 6

    @classmethod
    def _ccw(cls):
        return (ClockInt.ccw()**2).vec


class DigiClockInt(CSymNum):
    """Numbers with 15 degree rotations (digiclock numbers)."""
    R = 24
    D_SQ = [1, 2, 3, 6]
    L = 4
    _CCW = [(0,0), (1, -1), (0,0), (1,1)]

    @staticmethod
    def _mul(x, y):
        # Dimension multiplication matrix (for Z[i]-valued vectors):
        #     e    f     g    h
        # a [1   ,  x  ,   y ,  xy ]
        # b [ x  , 2   ,  xy , 2 y ]
        # c [  y ,  xy , 3   , 3x  ]
        # d [ xy , 2 y , 3x  , 6   ]
        (a, b, c, d), (e, f, g, h) = (x, y)
        return [
            a*e + 2*b*f + 3*c*g + 6*d*h,
            a*f + b*e + 3*(c*h + d*g),
            a*g + c*e + 2*(b*h + d*f),
            a*h + b*g + c*f + d*e,
        ]

ZZ = {
    4: GaussInt,
    6: EisensteinInt,
    8: CompassInt,
    10: PenroseInt,
    12: ClockInt,
    24: DigiClockInt,
}
"""Convenience access to different complex integer classes."""


def test_rings():
    from math import e, pi

    for cl in ZZ.values():
        print("Testing units of ring", cl.__name__, f"= ZZ_{cl.R}")

        # half cycle = -1
        assert cl.unit(cl.R//2) == -cl.one()
        assert complex(-cl.one()) == -1

        # full cycle = 1
        assert cl.unit(cl.R) == cl.one()
        assert complex(cl.one()) == 1

        # check deviation from reference of complex rotation vs matrix
        # and correctness of direct multiplication table vs rotation matrix
        for k in range(-cl.R, cl.R + 1):
            ref = e**(2j*pi*k/cl.R)
            unit = cl.unit(k)
            assert abs(complex(unit)-ref) < 1e-15

        # check unit multiplication (should rotate correctly)
        for m in range(-cl.R, cl.R + 1):
            for n in range(-cl.R, cl.R + 1):
                mul_result = cl.unit(m) * cl.unit(n)
                exp_add_result = cl.unit(m + n)
                assert mul_result == exp_add_result
