from casadi import *
import matplotlib.pyplot as plt
import numpy as np

opti = Opti()

x = opti.variable()
y = opti.variable()

opti.minimize((1-x)**2 + (y-x**2)**2)

opti.solver('ipopt')
sol = opti.solve()

X = np.linspace(-0.5, 1.5, 100)
Y = np.linspace(-0.5, 1.5, 100)
XX, YY = np.meshgrid(X, Y)
ZZ = (1 - XX)**2 + (YY - XX**2)**2

print(XX.shape)
print(YY.shape)
print(ZZ.shape)

plt.contour(XX, YY, ZZ, levels=100, cmap='viridis')
plt.plot(sol.value(x), sol.value(y), 'ro', label='Optimal Solution')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-0.5, 1.5)
plt.ylim(-0.5, 1.5)
plt.legend()
plt.title('Rosenbrock Function Contour and Optimum')
plt.show()