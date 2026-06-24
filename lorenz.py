import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def model(t, z):
    x, y = z
    a = 1.0 + 0.5*np.sin(0.2*t)   # 外部周期扰动
    b, c, d = 0.5, 1.0, 0.5
    K, h = 40, 0.2

    dx = a*x*(1 - x/K) - (b*x*y)/(1 + h*x)
    dy = (d*x*y)/(1 + h*x) - c*y

    return [dx, dy]

t = np.linspace(0, 200, 20000)
sol = solve_ivp(model, [0, 200], [15, 8], t_eval=t)

plt.plot(sol.y[0], sol.y[1], lw=0.5)
plt.title("Forced Predator-Prey (Chaotic-like Behavior)")
plt.xlabel("Prey")
plt.ylabel("Predator")
plt.grid()
plt.show()