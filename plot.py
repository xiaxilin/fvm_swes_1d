#! /usr/bin/python
# Filename : draw.py

# A python program for drawing 1d result

import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('output.txt')

x = data[:,0]
h = data[:,1]
q = data[:,2]

plt.figure(figsize=(8,4))
plt.plot(x,h,color='blue',linewidth=2)
plt.xlabel('x(m)')
plt.ylabel('h(m)')
plt.title('water depth')
plt.ylim(-0.2,1.2)
plt.show()


plt.figure(figsize=(8,4))
plt.plot(x,q,color='red',linewidth=2)
plt.xlabel('x(m)')
plt.ylabel('q(m/s)')
plt.title('discharge')
#plt.ylim(-0.2,1.2)
plt.show()
