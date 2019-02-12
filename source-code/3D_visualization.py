import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

f = open("list_optimized_acc.txt", "r")
str_optimizing_acc = f.read()
list_acc = str_optimizing_acc.split(",")
Acc_array = np.array(list_acc)


C_array = np.zeros((1, 70))
gamma_array = np.zeros((1,499))
for i in range(len(C_array[0, ])):
    C_array[0,i] = (i+1) * 0.5
for i in range(len(gamma_array[0,])):
    gamma_array[0,i] = (i+1) * 0.001

acc_array = np.zeros((len(gamma_array[0,]),len(C_array[0,])))
for i in range(len(gamma_array[0,])):
    for j in range(len(C_array[0,])):
        acc_array[i][j] = Acc_array[i+j]
        #print("acc_array[%d][%d] = %s" %(i,j, acc_array[i][j]))
#print(acc_array)

fig = plt.figure()
ax = plt.axes(projection='3d')
X = C_array[0,:]
Y = gamma_array[0,:]
plt.xlabel("Parameter C")
plt.ylabel("Parameter gamma")
ax.set_zlabel('Prediction accuracy')
X,Y = np.meshgrid(X,Y)
Z = acc_array
surf = ax.plot_surface(X, Y, Z, cmap = "coolwarm")
ax.contour(X,Y,Z,zdir='z', offset=-3,cmap="coolwarm")
# ax.scatter(X, Y, Z, c = 'blue')
#ax.plot_wireframe(X, Y, Z)
ax.view_init(elev=30, azim=200)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

# ValueError: shape mismatch: objects cannot be broadcast to a single shape
# 注意z的数据的z.shape是(len(y), len(x))

