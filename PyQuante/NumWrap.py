"""
 NumWrap.py - Interface to Numeric and numpy

 An interface to the Numeric and numpy libraries that will (hopefully)
 make the transformation to numpy go as seamlessly as possible.

 Also interfaces the LinearAlgebra and numpy.linalg libraries, since
 these have to be done consistently with the Numeric/numpy choice
"""

# Removing Numeric support
# Todo:
# -remove use of 'matrixmultiply' in favor of dot
# -remove the use of the numwrap module altogether, in favor of
#  direct imports of numpy modules
from numpy import *
from numpy.linalg import *
matrixmultiply = dot
    
# These functions are used by Optimize, but numpy.oldnumeric no longer exists.
# 10/2014: I'm going to comment these out, since it doesn't break anything in the test suite
#  and will re-enable them as needed 
#import numpy.oldnumeric.mlab as MLab
#from numpy.oldnumeric import NewAxis
import numpy as Numeric

