import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, pi
from numpy import arctan as atan
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

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


Fi = []
Sig = []
X = []
Ksi = []
Psi = []
Om = []
U = []
ttiks = np.linspace(0, 1, 1001)
ResX = []
ResFi = []
Vx = []
Vy = []
N = []
L1 = []
L2 = []
L3 = []
L4 = []
L5 = []
L6 = []

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

def tmpSist(Sist): return np.linalg.inv(A(Sist)).dot(b(Sist))


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
#Fi(t), Sig(t), X(t), Ksi(t), Psi(t), Om(t), Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5]
def getSist2(t, Sist): return [T*Sist[1],
                               T*u(Sist[1], Sist[3], Sist[4], Sist[5], Sist[7], Sist[9], Sist[11]) / ro,
                               T*Sist[3],
                               T*tmpSist(Sist)[0],
                               T*Sist[5],
                               T*tmpSist(Sist)[1],
                               T*getL1(),
                               T*getL2(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11]),
                               T*getL3(),
                               T*getL4(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11]),
                               T*getL5(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11]),
                               T*getL6(Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5], Sist[6], Sist[7], Sist[8], Sist[9], Sist[10], Sist[11])]


def soprSist(t, Sist): return [-T*getL1(),
                               -T*getL2(Fi(t), Sig(t), X(t), Ksi(t), Psi(t), Om(t), Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5]),
                               -T*getL3(),
                               -T*getL4(Fi(t), Sig(t), X(t), Ksi(t), Psi(t), Om(t), Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5]),
                               -T*getL5(Fi(t), Sig(t), X(t), Ksi(t), Psi(t), Om(t), Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5]),
                               -T*getL6(Fi(t), Sig(t), X(t), Ksi(t), Psi(t), Om(t), Sist[0], Sist[1], Sist[2], Sist[3], Sist[4], Sist[5])]


def eevent1(t, Sist): return Nopt(Sist)
eevent1.terminal = False
eevent1.direction = -1



# Основные вычисления
step = 0.001

#Sist0 = [fi0, sig0, x0, ksi0, psi0, om0]


#Метод Крылова-Черноусько
for ind in range (0, 2):
    Ta = 0.1
    Tb = 99.9
    T = 0.5 * (Ta + Tb)
    for ind2 in range(20):
        if abs(Tb - Ta) < 10 ** (-3):
            break
        solv31 = solve_ivp(getSist2, [0, 1], y0=Sist1, events=eevent1, max_step=0.001,
                           t_eval=ttiks)
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
    solv32 = solve_ivp(soprSist, [0, 1], y0=[0, 0, 0, 0, 0, 0], max_step=0.001, t_eval=ttiks)
    U.clear()
    L1.clear()
    L2.clear()
    L3.clear()
    L4.clear()
    L5.clear()
    L6.clear()
    for ind5 in range (0, len(solv32.t)):
       L1.append(solv32.y[0][len(solv32.t) -1 -ind5])
       L2.append(solv32.y[1][len(solv32.t) -1 -ind5])
       L3.append(solv32.y[2][len(solv32.t) -1 -ind5])
       L4.append(solv32.y[3][len(solv32.t) -1 -ind5])
       L5.append(solv32.y[4][len(solv32.t) -1 -ind5])
       L6.append(solv32.y[5][len(solv32.t) -1 -ind5])
    
    #def u(b, d, e, f, l2, l4, l6)
    # Sist1 = [fi0 a 0, sig0 b 1, x0 c 2, ksi0 d 3, psi0 e 4, om0 f 5, l10, l20, l30, l40, l50, l60]
    for ind3 in range(0, len(solv31.t)):
        U.append(u(solv31.y[1][ind3],solv31.y[3][ind3], solv31.y[4][ind3], solv31.y[5][ind3],L2[ind3], L4[ind3], L6[ind3]))
       # Vx = solv31.y[3][ind3] + kappa * solv31.y[5][ind3] * np.cos(solv31.y[4][ind3])
       # Vy = kappa * solv31.y[5][ind3] * sin(solv31.y[4][ind3])
       # N = Nopt([solv31.y[0][ind3], solv31.y[1][ind3], solv31.y[2][ind3], solv31.y[3][ind3], solv31.y[4][ind3], solv31.y[5][ind3],L1[ind3], L2[ind3], L3[ind3], L4[ind3], L5[ind3], L6[ind3]])
    Ufun = CubicSpline(solv31.t, U)
    ResX = solv31.y[2]
    ResFi = solv31.y[0]
    print("Step " + str(ind) + " T = " + str(T))



for ind3 in range(0, len(solv31.t)):
    Vx .append(solv31.y[3][ind3] + kappa * solv31.y[5][ind3] * np.cos(solv31.y[4][ind3]))
    Vy.append(kappa * solv31.y[5][ind3] * sin(solv31.y[4][ind3]))
    N.append(Nopt([solv31.y[0][ind3], solv31.y[1][ind3], solv31.y[2][ind3], solv31.y[3][ind3], solv31.y[4][ind3], solv31.y[5][ind3],L1[ind3], L2[ind3], L3[ind3], L4[ind3], L5[ind3], L6[ind3]]))





plt.figure(1)
plt.plot(ttiks, ResX, 'o', ms=1.2, label='X/T опт', color="g")
plt.xlabel('t')
plt.ylabel('X')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)

plt.figure(2)
plt.plot(ttiks, ResFi, 'o', ms=1.2, label='Fi no opt/T', color = "b")
plt.xlabel('t')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.ylabel('Fi')

plt.figure(3)
plt.plot(ttiks, N, 'o', ms=1.2, label='N opt/T', color="b")
plt.xlabel('t')
plt.ylabel('N')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid()

plt.figure(4)
plt.plot(ttiks, U, 'o', ms=1.2, label='U/T')
plt.xlabel('t')
plt.ylabel('U')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid()

plt.figure(5)
plt.plot(ttiks, Vx, 'o', ms=1.2, label='Vx/T', color="b")
plt.xlabel('t')
plt.ylabel('Vx')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid()


plt.plot(ttiks, Vy, 'o', ms=1.2, label='Vy/T', color="g")
plt.xlabel('t')
plt.ylabel('Vy')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid()


plt.show()
