import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Población de República Dominicana
N = 10266000
# Infectados
I0 = 75660
# Recuperados
R0 = 40122 + 1222
# Suceptibles
S0 = N - I0 - R0
# Tasa de transmisión a nivel del total de evaluados
beta = 0.2717
# Tasa de recuperación en 1/días
gamma = 1 / 14
# Tiempo entre 5-Agosto-2020 & 1-Julio-2021
# Magnitud de tiempo en dias
t = np.linspace(0, 330, 330)
# Tasa de retrasmisión probable
f = 0.14
# Tasa de mortalidad
u = 0.0162


# ecuaciones diferencial del modelo SIRS
def deriv(y, t, N, beta, gamma, f, u):
    S, I, R = y
    dSdt = -(beta * S * I / N)
    dIdt = (beta * S * I / N) - (gamma * I) 
    dRdt = (gamma * I)
    return dSdt, dIdt, dRdt


# vector de condiciones iniciales
y0 = S0, I0, R0

# resolviendo las ecuaciones del metodo deriv
ret = odeint(deriv, y0, t, args=(N, beta, gamma, f, u))
S, I, R = ret.T

fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#aaaaaa', axisbelow=True)
ax.plot(t, S/1000000, 'b', alpha=0.5, lw=2, label='Suceptibles')
ax.plot(t, I/1000000, 'r', alpha=0.5, lw=2, label='Infectados')
ax.plot(t, ((R) * (1 - u))/1000000, 'g', alpha=0.5, lw=2, label='Recuperados')
ax.plot(t, ((R) * u)/1000000, '#000000', alpha=0.5, lw=2, label='Muertos')
ax.set_xlabel('Dias')
ax.set_ylabel('Población dominicana (x 1,000,000)')
ax.set_ylim(0, 11)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()


#0.02863
