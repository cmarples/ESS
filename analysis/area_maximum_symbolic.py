# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 18:12:48 2021

@author: Callum Marples

Symbolic calculation used in Appendix A to help derive the maximum surface area
element on a triaixal ellipsoid.
"""

from sympy import simplify
from sympy.abc import a, b, c

f = a**2 / (2*(a**2-c**2)) * (b**2*c**2*a**2/((2*(a**2-c**2))) - a**2*b**2*a**2/(2*(a**2-c**2)) + a**2*b**2)

g = b**2 / (2*(b**2-c**2)) * (a**2*c**2*b**2/((2*(b**2-c**2))) - a**2*b**2*b**2/(2*(b**2-c**2)) + a**2*b**2)

d = f - g
D = simplify(d)

j = g - (a*c)**2
J = simplify(j)