import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# Константы
fi0 = sig0 = x0 = ksi0 = psi0 = 0
om0 = 1.1
J = 0.1
m = 0.2
M = 0.3
L = 4
R = 6
mu = 0.8
T = 10
k = 4

kappa = (m * L) / ((m + M) * R)
nu = (m * L * R) / (J + m * L ** 2)
ro = J / (J + m * L ** 2)


# Система
def u(): return -1


Sist = [fi0, sig0, x0, ksi0, psi0, om0]


# def s1(Sist): return Sist[1]
#
#
# def s2(Sist): return u() / ro
#
#
# def s3(Sist): return Sist[3]
#
# def s4(Sist): return s6(Sist) * (kappa * tmp1(Sist) * np.sin(Sist[4]) - kappa * np.cos(Sist[4])) \
#                      + tmp1(Sist) * kappa * Sist[5] ** 2 * np.cos(Sist[4]) + tmp1(Sist) \
#                      + kappa * Sist[5] ** 2 * np.sin(Sist[4])
#
#
# def s5(Sist): return Sist[5]
#
#
# def s6(Sist): return (u() - nu * np.sin(Sist[4]) - nu * np.cos(Sist[4]) * (
#             tmp1(Sist) * kappa * Sist[5] ** 2 * np.cos(Sist[4]) + tmp1(Sist) + kappa * Sist[5] ** 2 * np.sin(Sist[4]))) \
#                      / (1 + nu * np.cos(Sist[4]) * (kappa * tmp1(Sist)) * np.sin(Sist[4]) - kappa * np.cos(Sist[4]))
#
#
# def getSist(t, Sist): return [s1(Sist), s2(Sist), s3(Sist), s4(Sist), s5(Sist), s6(Sist)]



def tmp1(Sist): return 2 * mu * np.arctan(k * (Sist[3] + Sist[5] - Sist[1])) / np.pi


def A(Sist): return [[1, kappa * (np.cos(Sist[4]) - tmp1(Sist) * np.sin(Sist[4]))],
                     [nu * np.cos(Sist[4]), 1]]


def b(Sist): return [tmp1(Sist) * (1 + kappa * Sist[5] ** 2 * np.cos(Sist[4])) + kappa * Sist[5] ** 2 * np.sin(Sist[4]),
                     u() - nu * np.sin(Sist[4])]


def tmpSist(Sist): return np.linalg.inv(A(Sist)).dot(b(Sist))


def N(Sist): return kappa*(tmpSist(Sist)[1]*np.sin(Sist[4]) + Sist[5]*np.cos(Sist[4])) + 1

# Sist = [fi0, sig0, x0, ksi0, psi0, om0]
def getSist2(t, Sist): return [Sist[1], u() / ro, Sist[3], tmpSist(Sist)[0], Sist[5], tmpSist(Sist)[1]]

def eevent(t, Sist): return N(Sist)
eevent.terminal = True
eevent.direction = -1

Tmp = getSist2(1, Sist)
# Основные вычисления
step = 0.001
res = solve_ivp(getSist2, [0, 20], y0=Sist, max_step=step, events = eevent)

# Визуализация результатов
X = []
X.append(x0)

Fi = []
Fi.append(fi0)

Nm = []
Nm.append(N([res.y[0][0], res.y[1][0], res.y[2][0], res.y[3][0], res.y[4][0], res.y[5][0]]))

for i in range(0, len(res.t) - 1):
    X.append(res.y[2][i])
    Fi.append(res.y[4][i])
    Nm.append(N([res.y[0][i], res.y[1][i], res.y[2][i], res.y[3][i], res.y[4][i], res.y[5][i]]))

plt.figure(1)
plt.plot(res.t, X, 'o', ms=1.2, label='X/T')
plt.xlabel('t')
plt.ylabel('X')
plt.figure(2)
plt.plot(res.t, Fi, 'o', ms=1.2, label='Fi/T')
plt.xlabel('t')
plt.ylabel('Fi')
plt.figure(3)
plt.plot(res.t, Nm, 'o', ms=1.2, label='N/T')
plt.xlabel('t')
plt.ylabel('N')
plt.grid()
plt.show()
