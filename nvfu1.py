import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def model(t, z):
    x, y = z

    r = 1.0
    K = 25

    a = 0.5
    h = 0.2

    b = 0.4
    d = 1.0   # 让系统稳定但保留振荡

    dx = r*x*(1 - x/K) - (a*x*y)/(1 + h*x)
    dy = (b*x*y)/(1 + h*x) - d*y

    return [dx, dy]

t = np.linspace(0, 100, 8000)
sol = solve_ivp(model, [0, 100], [20, 8], t_eval=t)

plt.figure()
plt.plot(sol.y[0], sol.y[1])
plt.title("Spiral Sink: Oscillation → Stable Point (λ < 0)")
plt.xlabel("Prey")
plt.ylabel("Predator")
plt.grid()
plt.show()