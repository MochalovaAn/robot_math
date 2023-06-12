
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# Константы
fi0 = x0 = ksi0 = psi0 = 0
gr = 57.2958
sig0 = -87.890625/gr
om0 = 0

# #Real data
J = 0.0027
m = 0.270
M = 0.8
L = 0.03
R = 0.0625

#Эксп. данные
# J = 0.1
# m = 0.1
# M = 1
# L = 4.5
# R = 5



mu = 0.8
T = 10
k = 100
g = 9.80665

kappa = (m * L) / ((m + M) * R)
nu = (m * L * R) / (J + m * L ** 2)
ro = J / (J + m * L ** 2)


# Система
def u(t): return 0


Sist = [fi0, sig0, x0, ksi0, psi0, om0]
NUReal = [0, (2*np.pi*500)/2048, 0, 0, 0,(2*np.pi*500)/2048]
NUExp = [0, 0, 0, 0, 0, 0.1]


def tmp1(Sist): return 2 * mu * np.arctan(k * (Sist[3] + Sist[5] - Sist[1])) / np.pi
# def tmp1(Sist): return 1

def A(Sist): return [[1, kappa * (np.cos(Sist[4]) - tmp1(Sist) * np.sin(Sist[4]))],
                     [nu * np.cos(Sist[4]), 1]]


def b(Sist): return [tmp1(Sist) * (1 + kappa * Sist[5] ** 2 * np.cos(Sist[4])) + kappa * Sist[5] ** 2 * np.sin(Sist[4]),
                     u(Sist[2]) - nu * np.sin(Sist[4])]


def tmpSist(Sist): return np.linalg.inv(A(Sist)).dot(b(Sist))


def N(Sist): return kappa*(tmpSist(Sist)[1]*np.sin(Sist[4]) + Sist[5]*np.cos(Sist[4])) + 1

# Sist = [fi0, sig0, x0, ksi0, psi0, om0]
def getSist2(t, Sist): return [Sist[1], u(Sist[2]) / ro, Sist[3], tmpSist(Sist)[0], Sist[5], tmpSist(Sist)[1]]
#def getSist2(t, Sist): return [Sist[1], tmpSist(Sist)[0] + tmpSist(Sist)[1], Sist[3], tmpSist(Sist)[0], Sist[5], tmpSist(Sist)[1]]

def one_meter(t, Sist): return abs(Sist[2]) - 15.8730158730 #1 metr
one_meter.terminal = True
one_meter.direction = 0

def otriv(t, Sist) : return N(Sist)
otriv.terminal = True
otriv.direction = 0


Tmp = getSist2(1, Sist)
# Основные вычисления
step = 0.001
res = solve_ivp(getSist2, [0,2000], y0=  NUReal, max_step=step, events = one_meter)
# Визуализация результатов
X = []
X.append(x0)

Fi = []
Fi.append(fi0)

Psi = []
Psi.append(psi0)

Psi_1 = []
Psi_1.append(np.sign(res.y[4][0]) *((res.y[4][0]*gr + 1)%360)/360)

Nm = []
Nm.append(N([res.y[0][0], res.y[1][0], res.y[2][0], res.y[3][0], res.y[4][0], res.y[5][0]]))

T_ = np.sqrt(R/g)

Time = []
Time.append(0)

Vfi = []
Vfi.append(res.y[0][0])

Vx = []
Vx.append(res.y[3][0] + kappa * res.y[5][0] * np.cos(res.y[4][0]))

Scol = []
Scol.append(res.y[3][0] + res.y[5][0] - res.y[1][0])

Vy = []
Vy.append(kappa * res.y[5][0] * np.sin(res.y[4][0]))

# Sist = [fi0, sig0, x0, ksi0, psi0, om0]
#f = open('C:/Users/anmo1021/Desktop/testprog.txt', 'w')
for i in range(0, len(res.t) - 1):
    Time.append(res.t[i]*T_*14)
    X.append(res.y[2][i]*R)
    Fi.append(res.y[0][i]*gr)
    Scol.append(res.y[3][i] + res.y[5][i] - res.y[1][i])
    Psi_1.append((res.y[4][i] - res.y[0][i])*gr)
    Psi.append(np.sign(res.y[4][i]) *((360 - abs(res.y[4][i]*gr))%360)/360)
    Vfi.append(res.y[1][i]*gr/(360/2048)) #перевод скорости угловой в шаги в секунду
    Vx.append(res.y[3][i] + kappa * res.y[5][i] * np.cos(res.y[4][i]))
    Vy.append(kappa * res.y[5][i] * np.sin(res.y[4][i]))
    Nm.append(N([res.y[0][i], res.y[1][i], res.y[2][i], res.y[3][i], res.y[4][i], res.y[5][i]])*g*(m+M))
    #f.write("speed ")
    #f.write("acceleration 1000")
    #f.write("rotate 2048")
    #f.write("timer 2")

#for k in range (0, len(Time)):


#f.close()
# Ожидаемые результаты: для того, чтобы проехать 1 метр трубе диаметром 12.6 см надо сделать примерно 2.5 оборота
# 2.5 оборота - примерно 900 градусов или 15.707 радиан
# В реальности робот сначала слега покачивается назад, затем уверенно едет вперед

fig = 1

plt.figure(fig)
plt.plot(Time, X, 'o', ms=1.2, label='X/T')
plt.xlabel('t')
plt.ylabel('X')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid()

fig = fig+1
plt.figure(fig)
plt.plot(Time, Fi, 'o', ms=1.2, label='Fi/T')
plt.xlabel('t')
plt.ylabel('Fi')
plt.grid()
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)


fig = fig+1
plt.figure(fig)
plt.plot(Time, Vfi, 'o', ms=1.2, label='Vfi/T')
plt.xlabel('t')
plt.ylabel('Vfi')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid()

fig = fig+1
plt.figure(fig)
plt.plot(Time, Psi_1, 'o', ms=1.2, label='Угол поворота робота')
plt.xlabel('t')
plt.ylabel('Tetta')
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid()

fig = fig+1
plt.figure(fig)
plt.plot(res.t, Vx, 'o', ms=1.2, label='Vx/T')
plt.xlabel('t')
plt.ylabel('Vx')

fig = fig+1
plt.figure(fig)
plt.plot(res.t, Vy, 'o', ms=1.2, label='Vy/T')
plt.xlabel('t')
plt.ylabel('Vy')
plt.grid()
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)

fig = fig+1
plt.figure(fig)
plt.plot(Time, Scol, 'o', ms=1.2, label='Проскальзывание')
plt.xlabel('t')
plt.ylabel('Scol')
plt.grid()
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)

fig = fig+1
plt.figure(fig)
plt.plot(Time, Nm, 'o', ms=1.2, label='N')
plt.xlabel('t')
plt.ylabel('N')
plt.grid()
leg = plt.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.65, 0.9), mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)

plt.show()
