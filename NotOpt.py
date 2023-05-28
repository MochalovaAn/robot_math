import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from prettytable import PrettyTable

#x0 = 0
#y0 = 0
#xt0 = 0
#yt0 = 0
psi0 = 0
psit0 = 0
tetta0 = 0
m = 1
R = 5
M = 10
g = 9.80665
H = 10
I = 0.5
#C1 = -m * R * np.sin(psi0)
#C2 = m * R * np.cos(psi0) * psit0

#C3 = -m * R * np.cos(psi0)
#C4 = m * R * np.sin(psi0) * psit0

Sist = [psi0, psit0]

#S = []
#Y = []

ts = [0.0001, 0, 0]
L = (m+M)/(R**2 * m**2) *(I + m*M*R**2/(m+M))

def getX(t, p): return (C1 * t + C2 + m * R * np.sin(p)) / (m + M)


def getY(t, Sist): return ((getL(t, Sist) - (m + M) * g) / 2 * t ** 2 + C3 * t + C4 + m * R * np.cos(Sist[0])) / (m + M)


def u(t): return 6


def getL(t, Sist): return getW(t,Sist)[1]*np.sin(Sist[0]) + Sist[1]**2*np.cos(Sist[0]) + 1


def getPsi(Sist): return Sist[1]


def getW(t, Sist): return [Sist[1],
                           (u(t) - Sist[1] ** 2 * np.sin(Sist[0]) * np.cos(Sist[0]) - np.sin(Sist[0])) / (
                                   L + np.sin(Sist[0]) ** 2 + 1)]


getW(0, [0, 0])
ttiks = np.linspace(0, 1, 1001)


def eevent(t, Sist): return getL(t, Sist)


tmp = []

eevent.terminal = True
eevent.direction = -1

step = 0.01
solv0 = solve_ivp(getW, [ts[0], 10], y0=[Sist[0], Sist[1]], events=eevent,
                  max_step=step)
tmp.clear()
La = []
La.append(0)
for j in range(1, len(solv0.t) - 1):
    #S.append(getX(solv0.t[j], solv0.y[0][j]))
    #Y.append(getY(solv0.t[j], [solv0.y[0][j], solv0.y[1][j]]))
    La.append(getL(solv0.t[j], [solv0.y[0][j], solv0.y[1][j]]))
La.append(0)
ax = plt.subplot()
#plt.plot(solv0.t, S, 'o', ms=1.2, label='Without F, not opt')
plt.plot(solv0.t, La, 'o', ms=1.2, label='Without F, not opt')
plt.xlabel('X')
plt.ylabel('Y')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.show()
