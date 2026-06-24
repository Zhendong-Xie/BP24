import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 参数
alpha = 1.0   # prey growth rate
beta  = 0.1   # predation rate
gamma = 1.5   # predator death rate
delta = 0.075 # predator reproduction rate

# Lotka-Volterra方程
def lv(t, z):
    x, y = z
    dxdt = alpha * x - beta * x * y
    dydt = delta * x * y - gamma * y
    return [dxdt, dydt]

# 初始条件
x0 = 40   # prey
y0 = 9    # predator
z0 = [x0, y0]

# 时间范围
t_span = (0, 200)
t_eval = np.linspace(*t_span, 10000)

# 求解
sol = solve_ivp(lv, t_span, z0, t_eval=t_eval)

x = sol.y[0]
y = sol.y[1]
t = sol.t

# -------------------------
# 1. 时间序列
# -------------------------
plt.figure(figsize=(10,4))
plt.plot(t, x, label="Prey (x)")
plt.plot(t, y, label="Predator (y)")
plt.xlabel("Time")
plt.ylabel("Population")
plt.title("Prey-Predator Dynamics (Time Series)")
plt.legend()
plt.grid()
plt.show()

# -------------------------
# 2. 相空间图
# -------------------------
plt.figure(figsize=(5,5))
plt.plot(x, y)
plt.xlabel("Prey (x)")
plt.ylabel("Predator (y)")
plt.title("Phase Space (x vs y)")
plt.grid()
plt.show()