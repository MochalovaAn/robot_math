import Parts.Constans as co
import Parts.Visual as vi
import Parts.MainFunction as mf
import Parts.differ as di

import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.integrate
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from sklearn import preprocessing


def getDiscr(ksi10, ksi20):
    Sist = [co.psi0, co.w0, ksi10, ksi20]
    def eevent(t, Sist): return mf.getY(t, Sist) / (mf.getX(t, Sist[1]) + 0.0000001) + 10000000000
    step = 0.001
    res = di.resolveSist(mf.SoprSist, 10, Sist, eevent, step)
    l = len(res.t) - 1
    T = res.t[l]
    p = res.y[0][l]
    pt = res.y[1][l]

    return di.discrepancy(T, p ,pt)



a = np.linspace(-1,1)
b = np.linspace(-1,1)

X, Y = np.meshgrid(a,b, indexing='ij')
Z = []

for i in range (1, len(a)):
    for j in range (1, len(a)):
        Z = getDiscr(X[i,j], Y[i,j]);


ybnd = np.linspace(-3, 3);
xbnd = -0.25 * np.power(ybnd,2)


plt.contour(X, Y, Z, 'ShowText', 'on')
plt.plot(xbnd, ybnd, 'r', 'LineWidth', 2.0)

