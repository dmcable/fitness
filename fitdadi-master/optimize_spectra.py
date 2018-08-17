import cvxpy as cp
import numpy as np
import matplotlib.pyplot as plt

n = 5
p1 = np.linspace(0, 1, 30)
p2 = np.linspace(0, 1, 30) ** 2
p_arr = np.random.rand(n, 30)
z = np.linspace(0, 1, 30) * 1

for i in range(n):
    p_arr[i, :] = np.linspace(0, 1, 30) ** i
    p_arr[i, :] = p_arr[i, :] / sum(p_arr[i, :])

x = cp.Variable(n)
objective = cp.Maximize(cp.log(x * p_arr) * z)
constraints = [0 <= x, x <= 1, cp.sum(x) == 1]
prob = cp.Problem(objective, constraints)

# The optimal objective value is returned by `prob.solve()`.
result = prob.solve()
# The optimal value for x is stored in `x.value`.
print(x.value)

plt.plot(np.dot(np.expand_dims(x.value,1).T , p_arr)[0,:])
plt.hold(True)
plt.plot(z)
plt.show()
