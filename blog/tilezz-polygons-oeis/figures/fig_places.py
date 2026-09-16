#!/usr/bin/env python3
"""Figures for "Reachability, at every place of the field".

Each figure is one call with its parameters spelled out -- change a word, a ring
relabelling `g`, or a scale here and re-run. Nothing else reads these values.
"""
import common as C

# --- the words being drawn ---------------------------------------------------
DODECAGON = [1] * 12        # twelve unit steps, each +30 degrees: the regular 12-gon
HEXAGON = [2] * 6           # a rat: unit-edge hexagon
SEESAW = (0, 5, 5, 10)      # four ABSOLUTE step directions; head lands ~0.268 from origin
SHADOW_G = 5                # the relabelling zeta -> zeta^5, the other archimedean place

# --- the 12-gon and its shadow, side by side (38% / 50% display width) -------
# Both at the base scale, so the star's canvas comes out much smaller in px and
# therefore renders magnified when the two are shown at similar display widths.
print("12-gon pair:")
C.emit(C.walk(DODECAGON, 1), "fig-dodecagon.svg",
       "Regular 12-gon over Z[zeta_12]")
C.emit(C.walk(DODECAGON, SHADOW_G), "fig-dodecagon-shadow.svg",
       "The same walk under zeta -> zeta^5: the {12/5} star")

# --- the star again, alone and deliberately blown up ------------------------
print("star, magnified:")
C.boxed(C.walk(DODECAGON, SHADOW_G), "fig-dodecagon-shadow-zoom.svg",
        "The {12/5} star: the regular 12-gon under the relabelling zeta -> zeta^5",
        box=200, unit=170.0, stroke=5.0, arrow=(26.0, 11.0))

# --- a plain rat, for the "what is a rat" slot -------------------------------
print("hexagon:")
C.boxed(C.walk(HEXAGON, 1), "fig-hexagon.svg",
        "A rat: unit-edge hexagon over Z[zeta_12]",
        box=130, stroke=2.6, arrow=(11.0, 4.6))

# --- the seesaw: same four steps, two places, lengths multiplying to 1 -------
print("seesaw pair:")
C.emit(C.walk_dirs(SEESAW, 1), "fig-seesaw-physical.svg",
       "A four-step walk whose head lands 0.268 from the origin", label="|h|")
C.emit(C.walk_dirs(SEESAW, SHADOW_G), "fig-seesaw-shadow.svg",
       "The same four steps relabelled: the head lands 3.732 from the origin",
       label="|σ₅(h)|")
