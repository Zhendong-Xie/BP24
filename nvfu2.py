import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def model(t, z):
    x, y = z
    r, K = 1.0, 50
    a, h = 0.5, 0.05
    b, d = 0.3, 0.4

    dx = r*x*(1 - x/K) - (a*x*y)/(1 + h*x)
    dy = (b*x*y)/(1 + h*x) - d*y
    return [dx, dy]

t = np.linspace(0, 100, 40000)
sol = solve_ivp(model, [0, 100], [10, 5], t_eval=t)

plt.plot(sol.y[0], sol.y[1])
plt.title("Limit Cycle (λ ≈ 0) - Periodic Oscillation")
plt.xlabel("Prey")
plt.ylabel("Predator")
plt.grid()
plt.show()