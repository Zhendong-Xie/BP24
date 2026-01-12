import numpy as np
import matplotlib.pyplot as plt 


def entd(x):
    entD = -x*np.log2(x)
    return entD

x_arr = np.linspace(0.01, 0.99, 99 )
x_arr_1 = 1-x_arr 
y = entd(x_arr)
y1 = entd(x_arr_1)
plt.figure()
plt.plot(x_arr, y)
plt.plot(x_arr, y+y1) 
plt.plot(x_arr, np.exp(x_arr))
plt.show()