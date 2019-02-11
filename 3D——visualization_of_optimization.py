import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

f = open("list_optimized_acc.txt", "r")
str_optimizing_acc = f.read()
list_acc = str_optimizing_acc.split(",")
acc_array = np.array(list_acc)

if(1):
    C_array = np.zeros((1, 70))
    gamma_array = np.zeros((1,499))
    for i in range(len(C_array[0, ])):
        C_array[0,i] = (i+1) * 0.5
    for i in range(len(gamma_array[0,])):
        gamma_array[0,i] = (i+1) * 0.001

    acc_array = acc_array.reshape((70,499))

if(0):
    C_array = np.zeros((1, len(acc_array)))
    gamma_array = np.zeros((1, len(acc_array)))
    temp = 0
    for i in range(len(C_array[0,])):
        if i % 499 ==0:
            temp += 0.5
        C_array[0, i] = temp
    for i in range(len(gamma_array[0,])):
        if i % 499 == 0:
            gamma_array[0, i] = 0.001
        else:
            k = i % 499
            gamma_array[0, i] = 0.001 + k * 0.001


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X = C_array
Y = gamma_array
X,Y = np.meshgrid(X,Y)
Z = acc_array
ax.plot_surface(X, Y, Z,cmap='rainbow')
# ax.scatter(X, Y, Z, c = 'blue')
plt.show()

