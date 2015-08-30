#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt('summary.txt', skiprows=1)
R_vals = data[::,0]
rhf = data[::,1]
fthf = data[::,2]
labels = ['RHF', 'FTHF']



fig = plt.figure(1)
ax1 = fig.add_subplot(111)
ax1.plot(R_vals, rhf, c = 'r',label = labels[0])
ax1.plot(R_vals, fthf,color = 'r',linestyle = '--', label = labels[1])
plt.xlabel('Radius for H_2 bond bohr')
plt.ylabel('Single point energy (au))')
plt.title('FTHF versus HF at 100 K')
plt.legend(loc = 'upper right')
plt.savefig('FTHFversusHF')
#
#plt.figure(2)
#data2=np.loadtxt('pair.dat',usecols =[0,-1,-4])
#plt.plot(data2[::,0],data2[::,1])
#plt.show()
#plt.figure(3)
#plt.plot(data2[::,0],data2[::,2])
#plt.show()

