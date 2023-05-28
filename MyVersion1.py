import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from prettytable import PrettyTable

ksi0 = 0
eta0 = 0
alf0 = 1
beta0 = 1
w0 = 0.5
mu = 0.5
k = 0.5
delt = 0.01
kappa = 10
Eps = 12

# Xsn - X simple, not optimal, ind = 0
Xsn = [ksi0, alf0, eta0, beta0, w0]
# Xso - X simple, optimal, ind = 1
Xso = [ksi0, alf0, eta0, beta0, w0]
# Xn - X with forse,not optimal, ind = 2
Xn = [ksi0, alf0, eta0, beta0, w0]
# Xo - X with forse, more optimal, ind = 2
Xo = [ksi0, alf0, eta0, beta0, w0]


S = [[], [], [], []]
Y = [[], [], [], []]

eps = [mu * (1 + k) * beta0 / (1 - k), mu * (1 + k) * beta0 / (1 - k)]
Q = [eps[0] - 200 * delt, eps[1] + 200 * delt, 0]
ts = [0, 0, 0]

psi3 = 0
u = []
ttiks = np.linspace(0, 1, 1001)

def getXYsn(t, X): return [X[1], 0, X[3], -1, Q[0]]  # without forse, not optimal


def getXYso(t, X): return [X[1], 0, X[3], -1, Q[1]]  # without forse, optimal


def mysqrt(X): return (X[1] ** 2 + X[3] ** 2) ** 0.5


def getAlf(X): return -(Eps * X[1] * mysqrt(X) + X[3] * X[4])


def getBeta(X): return -(Eps * X[3] * mysqrt(X) + 1 - X[1] * X[4])


def getXYn(t, X): return [X[1], getAlf(X), X[3], getBeta(X),
                          Q[2] - kappa * np.abs(X[4]) * X[4]]  # with forse, not optimal


def getXYo(t, X): return [X[1] * T, T * (-Eps * X[1] * mysqrt(X) - X[4] * X[3]), X[3] * T,
                         T * (-Eps * X[3] * mysqrt(X) + X[4] * X[1] - 1),
                         (Ufun(1 - t) - kappa * np.sign(X[4]) * X[4] ** 2) * T]


# X = [X0=psi2, X1=psi4, X2=psi5]
# Сопряженная система
def getPsi2(t, X): return 1 + X[0] * Eps * (
        np.sqrt(VXfun(1 - t) ** 2 + VYfun(1 - t) ** 2) + (VXfun(1 - t) ** 2) / np.sqrt(
    VXfun(1 - t) ** 2 + VYfun(1 - t) ** 2)) + X[1] * (
                                  Eps * VXfun(1 - t) * VYfun(1 - t) / np.sqrt(
                              VXfun(1 - t) ** 2 + VYfun(1 - t) ** 2) - Wfun(1 - t))


def getPsi4(t, X): return X[0] * (
        Eps * VXfun(1 - t) * VYfun(1 - t) / np.sqrt(VXfun(1 - t) ** 2 + VYfun(1 - t) ** 2) + Wfun(1 - t)) - psi3 + X[
                              1] * Eps * (
                                  np.sqrt(VXfun(1 - t) ** 2 + VYfun(1 - t) ** 2) + VYfun(1 - t) ** 2 / np.sqrt(
                              VXfun(1 - t) ** 2 + VYfun(1 - t) ** 2))


def getPsi5(t, X): return X[0] * VYfun(1 - t) - X[1] * VXfun(1 - t) + X[2] * 2 * kappa * np.abs(Wfun(1 - t))


def oneStep(t, X): return [-getPsi2(t, X) * T, -getPsi4(t, X) * T, -getPsi5(t, X) * T]


# Функция задания новых начальных условий
def getNewX(solv, ind, sind,  q):
    eps = solv.y[1][sind] - solv.y[4][sind]
    return [solv.y[0][ind],
            solv.y[1][sind] - mu * (1 + k) * solv.y[3][sind]*np.sign(eps),
            solv.y[2][ind],
            -k * solv.y[3][sind],
            solv.y[4][sind] - mu * (1 + k) * solv.y[3][sind]*np.sign(eps) + q]


def eevent(t, X): return X[2]


eevent.terminal = True
eevent.direction = -1

for i in range(0, 10):
    step = 0.0001
    solv0 = solve_ivp(getXYsn, [ts[0], float("inf")], y0=[Xsn[0], Xsn[1], Xsn[2], Xsn[3], Xsn[4]], events=eevent,
                      max_step=step)
    solv1 = solve_ivp(getXYso, [ts[1], float("inf")], y0=[Xso[0], Xso[1], Xso[2], Xso[3], Xso[4]], events=eevent,
                      max_step=step)
    solv2 = solve_ivp(getXYn, [ts[2], float("inf")], y0=[Xn[0], Xn[1], Xn[2], Xn[3], Xn[4]], events=eevent,
                      max_step=step)

    u.clear()

    for i in range(0, 1001):
        u.append(1)
    Ufun = CubicSpline(ttiks, u)
    psi3 = 0
    for ind1 in range(0, 4):

        Ta = 0.1
        Tb = 99.9
        T = 0.5 * (Ta + Tb)
        for ind2 in range(20):
            if abs(Tb - Ta) < 10 ** (-3):
                break
            solv31 = solve_ivp(getXYo, [0, 1], y0=[Xo[0], Xo[1], Xo[2], Xo[3], Xo[4]], events=eevent, max_step=0.001, t_eval=ttiks)
            if len(solv31.t_events[0]) == 0:
                Ta = T
            else:
                Tb = T
            T = 0.5 * (Ta + Tb)
        psi3 = solv31.y[1][len(solv31.t) - 2] / solv31.y[3][len(solv31.t) - 2]
        Xfun = CubicSpline(solv31.t, solv31.y[0])
        VXfun = CubicSpline(solv31.t, solv31.y[1])
        Yfun = CubicSpline(solv31.t, solv31.y[2])
        VYfun = CubicSpline(solv31.t, solv31.y[3])
        Wfun = CubicSpline(solv31.t, solv31.y[4])
        solv32 = solve_ivp(oneStep, [0, 1], y0=[0, 0, 0], max_step=0.001)
        u.clear()
        for f in range(0, len(solv32.t)):
            Psi5 = solv32.y[2][f]
            u.append(np.sign(Psi5))
        Ufun = CubicSpline(solv32.t, u)

    for j in range(1, len(solv0.t) - 1):
        S[0].append(solv0.y[0][j])
        Y[0].append(solv0.y[2][j])
    for j in range(1, len(solv1.t) - 1):
        S[1].append(solv1.y[0][j])
        Y[1].append(solv1.y[2][j])
    for j in range(1, len(solv2.t) - 1):
        S[2].append(solv2.y[0][j])
        Y[2].append(solv2.y[2][j])
    for j in range(1, len(solv31.t) - 1):
        S[3].append(solv31.y[0][j])
        Y[3].append(solv31.y[2][j])

    last_ind = [len(solv0.t) - 2, len(solv1.t) - 2, len(solv2.t) - 2, len(solv31.t) - 2]
    ts = [solv0.t[last_ind[0]], solv1.t[last_ind[1]], solv2.t[last_ind[2]], solv31.t[last_ind[3]]]
    Q = [solv0.y[1][last_ind[0]] - solv0.y[4][last_ind[0]] - 200 * delt,
         solv1.y[1][last_ind[1]] - solv1.y[4][last_ind[1]] + 200 * delt, 0, Ufun(solv32.t[len(solv32.t) -1 ])]
    Xsn = getNewX(solv0, last_ind[0], 0, Q[0])
    Xso = getNewX(solv1, last_ind[1], 0, Q[1])
    Xn = getNewX(solv2, last_ind[2], 0, Q[2])
    Xo = getNewX(solv31, last_ind[3], last_ind[3], Q[3])

ax = plt.subplot()
plt.plot(S[0], Y[0], 'o', ms =1.2, label='Without F, not opt')
plt.plot(S[1], Y[1], label='Without F, opt')
plt.plot(S[2], Y[2], label='With F, not opt')
plt.plot(S[3], Y[3], label='With F, more opt')
plt.xlabel('X')
plt.ylabel('Y')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.show()
