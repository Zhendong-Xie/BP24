import numpy as np
import matplotlib.pyplot as plt

# x 轴
x = np.linspace(0, 10, 200)

# 生成随机波动（随机噪声 + 平滑趋势）
noise = np.random.randn(200) * 0.5
trend = np.sin(x)  # 基础趋势
y = trend + noise

# 绘图
plt.figure(figsize=(8, 4))
plt.plot(x, y, lw=1.5)

plt.title("Random XY Line Plot")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)

plt.show()