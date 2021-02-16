import numpy as np
import phy4910

# build a nonrelativistic white dwarf, just for fun
# these numbers match worksheet 3
phy4910.build_white_dwarf_polytrope(1.5, 3.166e12, 4.045e6)

# now do a relativistic one for Assignment 1

phy4910.build_white_dwarf_polytrope(3, 4.936e14, 53.31e6)
