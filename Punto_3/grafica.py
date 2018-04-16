import matplotlib.pylab as plt
import numpy as np
x=np.linspace(0,10,50);
y=np.cos(x);
data=np.array([x,y]).T
np.savetxt('First.txt',data)
plt.plot(x,y,Data[:,0],Data[:,1])
plt.show()
