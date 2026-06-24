import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def model(t, z):
    x, y = z

    # 外部周期扰动（关键）
    r = 1.0 + 0.6*np.sin(0.3*t)

    K = 50
    a, h = 0.8, 0.2
    b, d = 0.4, 0.5

    dx = r*x*(1 - x/K) - (a*x*y)/(1 + h*x)
    dy = (b*x*y)/(1 + h*x) - d*y
    return [dx, dy]

t = np.linspace(0, 200, 20000)
sol = solve_ivp(model, [0, 200], [8, 3], t_eval=t)

plt.plot(sol.y[0], sol.y[1], lw=0.5)
plt.title("Chaotic-like Dynamics (λ > 0)")
plt.xlabel("Prey")
plt.ylabel("Predator")
plt.grid()
plt.show()