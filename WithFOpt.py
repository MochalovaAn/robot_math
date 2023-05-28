import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, pi
from numpy import arctan as atan
from scipy.integrate import solve_ivp

# Константы
fi0 = sig0 = x0 = ksi0 = psi0  =l20=l30=l40=l50=l60=0
l10 = 0.2
om0 = 1.1
J = 0.1
m = 0.1
M = 1
L = 4.5
R = 4
mu = 0.8
T = 10
k = 4
u0 = 1

kappa = (m * L) / ((m + M) * R)
nu = (m * L * R) / (J + m * L ** 2)
ro = J / (J + m * L ** 2)





Sist0 = [fi0, sig0, x0, ksi0, psi0, om0]
Sist1 = [fi0, sig0, x0, ksi0, psi0, om0, l10, l20, l30, l40, l50, l60]

# Система
def u(b, d, e, f, l2, l4, l6): return u0* np.sign(l2/ro - l4*(kappa*np.cos(e) - kappa * np.sin(e) * 2 * mu * np.arctan(k * (-b + d + f)) / np.pi)/
                            (1 - nu *  kappa * cos(e)*(cos(e) - (2*mu*atan(k*(d + f -b))/pi) * sin(e))) +
                            l6/(1 - nu *  kappa * cos(e)*(cos(e) - (2*mu*atan(k*(d + f -b))/pi) * sin(e))))
def tmp1(Sist): return 2 * mu * np.arctan(k * (Sist[3] + Sist[5] - Sist[1])) / np.pi


def A(Sist): return [[1, kappa * (np.cos(Sist[4]) - tmp1(Sist) * np.sin(Sist[4]))],
                     [nu * np.cos(Sist[4]), 1]]


def b(Sist): return [tmp1(Sist) * (1 + kappa * Sist[5] ** 2 * np.cos(Sist[4])) + kappa * Sist[5] ** 2 * np.sin(Sist[4]),
                     u(Sist[1], Sist[3], Sist[4], Sist[5], Sist[7], Sist[9], Sist[11]) - nu * np.sin(Sist[4])]

def bb(Sist): return [tmp1(Sist) * (1 + kappa * Sist[5] ** 2 * np.cos(Sist[4])) + kappa * Sist[5] ** 2 * np.sin(Sist[4]),
                     u0 - nu * np.sin(Sist[4])]

def tmpSist(Sist): return np.linalg.inv(A(Sist)).dot(b(Sist))

def tmpSist2(Sist): return np.linalg.inv(A(Sist)).dot(bb(Sist))

def Nnoopt(Sist): return kappa*(tmpSist2(Sist)[1]*np.sin(Sist[4]) + Sist[5]*np.cos(Sist[4])) + 1
def Nopt(Sist): return kappa*(tmpSist(Sist)[1]*np.sin(Sist[4]) + Sist[5]*np.cos(Sist[4])) + 1


def getL1(): return 0

def getL2(a, b, c, d, e, f, l1, l2, l3, l4, l5, l6): return -1 * (l1 + l4*(-2*k*kappa*mu*nu*(-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*sin(e)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) - 2*k*kappa*mu*(-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*sin(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)) - 2*k*mu*(f**2*kappa*cos(e) + 1)*(nu*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1) + 1)/(pi*(k**2*(-b + d + f)**2 + 1)) + (f**2*kappa*sin(e) + 2*mu*(f**2*kappa*cos(e) + 1)*atan(k*(-b + d + f))/pi)*(2*k*kappa*mu*nu**2*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*sin(e)*cos(e)**2/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) + 2*k*kappa*mu*nu*sin(e)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)))) + l6*(-2*k*kappa*mu*nu**2*(f**2*kappa*sin(e) + 2*mu*(f**2*kappa*cos(e) + 1)*atan(k*(-b + d + f))/pi)*sin(e)*cos(e)**2/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) + 2*k*kappa*mu*nu*(-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*sin(e)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) + 2*k*mu*nu*(f**2*kappa*cos(e) + 1)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)))
)

def getL3(): return 0

def getL4(a, b, c, d, e, f, l1, l2, l3, l4, l5, l6): return -1 * (l3 + l4*(2*k*kappa*mu*nu*(-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*sin(e)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) + 2*k*kappa*mu*(-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*sin(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)) + 2*k*mu*(f**2*kappa*cos(e) + 1)*(nu*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1) + 1)/(pi*(k**2*(-b + d + f)**2 + 1)) + (f**2*kappa*sin(e) + 2*mu*(f**2*kappa*cos(e) + 1)*atan(k*(-b + d + f))/pi)*(-2*k*kappa*mu*nu**2*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*sin(e)*cos(e)**2/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) - 2*k*kappa*mu*nu*sin(e)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)))) + l6*(2*k*kappa*mu*nu**2*(f**2*kappa*sin(e) + 2*mu*(f**2*kappa*cos(e) + 1)*atan(k*(-b + d + f))/pi)*sin(e)*cos(e)**2/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) - 2*k*kappa*mu*nu*(-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*sin(e)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) - 2*k*mu*nu*(f**2*kappa*cos(e) + 1)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)))
)

def getL5(a, b, c, d, e, f, l1, l2, l3, l4, l5, l6): return -1 *(l4*(nu*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1) - (-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*sin(e) + kappa*nu*(-2*mu*cos(e)*atan(k*(-b + d + f))/pi - sin(e))*cos(e))*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2 - (-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*(-2*kappa*mu*cos(e)*atan(k*(-b + d + f))/pi - kappa*sin(e))/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1) + (f**2*kappa*sin(e) + 2*mu*(f**2*kappa*cos(e) + 1)*atan(k*(-b + d + f))/pi)*(nu*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*sin(e) + kappa*nu*(-2*mu*cos(e)*atan(k*(-b + d + f))/pi - sin(e))*cos(e))*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2 - nu*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*sin(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1) + nu*(-2*kappa*mu*cos(e)*atan(k*(-b + d + f))/pi - kappa*sin(e))*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)) + (nu*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1) + 1)*(-2*f**2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + f**2*kappa*cos(e))) + l6*(-nu*(f**2*kappa*sin(e) + 2*mu*(f**2*kappa*cos(e) + 1)*atan(k*(-b + d + f))/pi)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*sin(e) + kappa*nu*(-2*mu*cos(e)*atan(k*(-b + d + f))/pi - sin(e))*cos(e))*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2 + nu*(f**2*kappa*sin(e) + 2*mu*(f**2*kappa*cos(e) + 1)*atan(k*(-b + d + f))/pi)*sin(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1) - nu*(-2*f**2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + f**2*kappa*cos(e))*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1) - nu*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1) + (-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*sin(e) + kappa*nu*(-2*mu*cos(e)*atan(k*(-b + d + f))/pi - sin(e))*cos(e))/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2)
)

def getL6(a, b, c, d, e, f, l1, l2, l3, l4, l5, l6): return -1 * (l4*(2*k*kappa*mu*nu*(-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*sin(e)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) + 2*k*kappa*mu*(-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*sin(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)) + (f**2*kappa*sin(e) + 2*mu*(f**2*kappa*cos(e) + 1)*atan(k*(-b + d + f))/pi)*(-2*k*kappa*mu*nu**2*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*sin(e)*cos(e)**2/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) - 2*k*kappa*mu*nu*sin(e)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1))) + (nu*(-2*kappa*mu*sin(e)*atan(k*(-b + d + f))/pi + kappa*cos(e))*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1) + 1)*(4*f*kappa*mu*cos(e)*atan(k*(-b + d + f))/pi + 2*f*kappa*sin(e) + 2*k*mu*(f**2*kappa*cos(e) + 1)/(pi*(k**2*(-b + d + f)**2 + 1)))) + l5 + l6*(2*k*kappa*mu*nu**2*(f**2*kappa*sin(e) + 2*mu*(f**2*kappa*cos(e) + 1)*atan(k*(-b + d + f))/pi)*sin(e)*cos(e)**2/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) - 2*k*kappa*mu*nu*(-nu*sin(e) + u(b, d, e, f, l2, l4, l6))*sin(e)*cos(e)/(pi*(k**2*(-b + d + f)**2 + 1)*(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1)**2) - nu*(4*f*kappa*mu*cos(e)*atan(k*(-b + d + f))/pi + 2*f*kappa*sin(e) + 2*k*mu*(f**2*kappa*cos(e) + 1)/(pi*(k**2*(-b + d + f)**2 + 1)))*cos(e)/(-kappa*nu*(-2*mu*sin(e)*atan(k*(-b + d + f))/pi + cos(e))*cos(e) + 1))
)



#Sist0 = [fi0, sig0, x0, ksi0, psi0, om0]
#Sist1 = [fi0 a, sig0 b, x0 c, ksi0 d, psi0 e, om0 f, l10, l20, l30, l40, l50, l60]
def getSist2(t, Sist): return [Sist[1], u(Sist[1], Sist[3], Sist[4], Sist[5], Sist[7], Sist[9], Sist[11]) / ro, Sist[3], tmpSist(Sist)[0], Sist[5], tmpSist(Sist)[1],
                               getL1(),
                               getL2(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11]),
                               getL3(),
                               getL4(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11]),
                               getL5(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11]),
                               getL6(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11])]
def notOpt(t, Sist): return [Sist[1], u0 / ro, Sist[3], tmpSist2(Sist)[0], Sist[5], tmpSist2(Sist)[1]]


def soprSist(t, Sist): return [-getL1(),
                               -getL2(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11]),
                               -getL3(),
                               -getL4(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11]),
                               -getL5(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11]),
                               -getL6(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11])]


def eevent1(t, Sist): return Nopt(Sist)
eevent1.terminal = True
eevent1.direction = -1


def eevent2(t, Sist): return Nnoopt(Sist)
eevent2.terminal = True
eevent2.direction = -1

# Основные вычисления
step = 0.001
res = solve_ivp(getSist2, [0, 10    ], y0=Sist1, max_step=step, events=eevent1)
res2 = solve_ivp(notOpt, [0, 20], y0=Sist0, max_step=step, events=eevent2)

#Sist0 = [fi0, sig0, x0, ksi0, psi0, om0]
Fi = []
Sig = []
X = []
Ksi = []
Psi = []
Om = []

#Метод Крылова-Черноусько
for ind in range (0, 10):
    Ta = 0.1
    Tb = 99.9
    T = 0.5 * (Ta + Tb)
    for ind2 in range(20):
        if abs(Tb - Ta) < 10 ** (-3):
            break
        solv31 = solve_ivp(getSist2, [0, 1], y0=[Sist1[0], Sist1[1], Sist1[2], Sist1[3], Sist1[4]], events=[eevent1, eevent2], max_step=0.001,
                          )
        if len(solv31.t_events[0]) == 0:
            Ta = T
        else:
            Tb = T
        T = 0.5 * (Ta + Tb)
    #Считаем граничные условия из системы
    #Применяем CubicSpline(solv31.t, solv31.y[]) для основной системы?
    Fi = CubicSpline(solv31.t, solv31.y[0])
    Sig = CubicSpline(solv31.t, solv31.y[1])
    X = CubicSpline(solv31.t, solv31.y[2])
    Ksi = CubicSpline(solv31.t, solv31.y[3])
    Psi = CubicSpline(solv31.t, solv31.y[4])
    Om = CubicSpline(solv31.t, solv31.y[5])
    #Решаем сопряженную систему в обратном времени
    solv32 = solve_ivp(soprSist, [0, 1], y0=[0, 0, 0, 0, 0, 0], max_step=0.001)
    u.clear()
    for f in range(0, len(solv32.t)):
        u.append(u(solv31.y[0][ind3], ))
    Ufun = CubicSpline(solv32.t, u)



# Визуализация результатов
X = []
#X.append(x0)

X2 = []
#X2.append(x0)

Fi = []
Fiopt = []
#Fi.append(fi0)


Nm1 = []
#Nm1.append(0)

Nm2 = []
#Nm2.append(0)

U = []

for i in range(0, len(res.t)):
   X.append(res.y[2][i])
   Fiopt.append(res.y[4][i])
   Nm2.append(Nopt(
       [res.y[0][i], res.y[1][i], res.y[2][i], res.y[3][i], res.y[4][i], res.y[5][i], res.y[6][i], res.y[7][i],
        res.y[8][i], res.y[9][i], res.y[10][i], res.y[11][i]]))
   U.append(u(res.y[1][i], res.y[3][i], res.y[4][i], res.y[5][i], res.y[7][i], res.y[9][i], res.y[11][i]))

#Sist1 = [fi0 0, sig0 1, x0 2, ksi0 3, psi0 4, om0 5, l10, l20, l30, l40, l50, l60]
#Vx = xt + kappa*psit*cos(psi) = ksi + kappa * omega * cos(psi)
Vx = res.y[3][len(res.t)-1] + kappa * res.y[5][len(res.t)-1] * np.cos(res.y[4][len(res.t)-1])
#Vy = kappa * omega * sin(psi)
Vy = kappa * res.y[5][len(res.t)-1] * np.sin(res.y[4][len(res.t)-1])
print("Vx = " + str(Vx))
print("Vy = " + str(Vy))
print("X = " + str(res.y[2][len(res.t)-1]))
print("Xt = " + str(res.y[3][len(res.t)-1]))

i=0
for i in range(0, len(res2.t)):
    X2.append(res2.y[2][i])
    Fi.append(res2.y[4][i])
    Nm1.append(Nnoopt([res2.y[0][i], res2.y[1][i], res2.y[2][i], res2.y[3][i], res2.y[4][i], res2.y[5][i]]))

plt.figure(1)
plt.plot(res.t, X, 'o', ms=1.2, label='X/T опт', color="g")
plt.plot(res2.t, X2, 'o', ms=1.2, label='X/T не опт', color = "b")
plt.xlabel('t')
plt.ylabel('X')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)

plt.figure(2)
plt.plot(res2.t, Fi, 'o', ms=1.2, label='Fi no opt/T', color = "b")
plt.plot(res.t, Fiopt, 'o', ms=1.2, label='Fi opt/T', color = "g")
plt.xlabel('t')
plt.ylabel('Fi')

plt.figure(3)
plt.plot(res2.t, Nm1, 'o', ms=1.2, label='N no opt/T', color="b")
plt.xlabel('t')
plt.ylabel('N')
plt.plot(res.t, Nm2, 'o', ms=1.2, label='N opt/T', color="g")
plt.xlabel('t')
plt.ylabel('N')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid()

plt.figure(4)
plt.plot(res.t, U, 'o', ms=1.2, label='U/T')
plt.xlabel('t')
plt.ylabel('U')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid()
plt.show()
