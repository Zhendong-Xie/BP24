import numpy as np
import matplotlib.pyplot as plt

# 生成测试数据
x = np.linspace(-5,5,100)
y = np.linspace(-5,5,100)
X,Y = np.meshgrid(x,y)
Z = np.exp(-(X**2+Y**2))

# 画等高线
plt.contourf(X, Y, Z, cmap="viridis", levels=30)

# 画色条 + 【反转色条】
cbar = plt.colorbar()
cbar.ax.invert_yaxis()  # 👈 核心代码

plt.show()