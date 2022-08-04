import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure

n=1
e=1

def f(B):
    return B/(n*e);

y=np.linspace(0,1,1000);


plt.plot(y,f(y),'r-')
plt.axhline(0.2, color='blue')


plt.title("Plot of resistivity $\\rho_{xy}$")
plt.xlabel('$B$')
plt.ylabel('$\\rho$')

ax = plt.gca()
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

legend1=plt.legend(['$\\rho_{xy}$', '$\\rho_{xx}$'], loc =1)
ax = plt.gca().add_artist(legend1)


plt.show()

