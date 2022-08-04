import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure

c0=1
a=1
v1=0
v2=0.8
v3=-0.8
v4=0.001
v5=-0.001

def dx1(k):
    return(c0*(1-v1)+c0*(1+v1)*np.cos(k*a))
def dy1(k):
    return(c0*(1+v1)*np.sin(k*a))

def dx2(k):
    return(c0*(1-v2)+c0*(1+v2)*np.cos(k*a))
def dy2(k):
    return(c0*(1+v2)*np.sin(k*a))

def dx3(k):
    return(c0*(1-v3)+c0*(1+v3)*np.cos(k*a))
def dy3(k):
    return(c0*(1+v3)*np.sin(k*a))

def dx4(k):
    return(c0*(1-v4)+c0*(1+v4)*np.cos(k*a))
def dy4(k):
    return(c0*(1+v4)*np.sin(k*a))

def dx5(k):
    return(c0*(1-v5)+c0*(1+v5)*np.cos(k*a))
def dy5(k):
    return(c0*(1+v5)*np.sin(k*a))

k=np.linspace(-np.pi,np.pi,1000)



plt.plot(dx1(k),dy1(k),'r-')
plt.plot(dx2(k),dy2(k),'g-')
plt.plot(dx3(k),dy3(k),'b-')
plt.plot(dx4(k),dy3(k),'c:')
plt.plot(dx5(k),dy3(k),'m:')

plt.axvline(0, color='black', linestyle=':')
plt.axhline(0, color='black', linestyle=':')

plt.title("Plot of $d_x$ versus $d_y$")
plt.xlabel(r'$d_x$')
plt.ylabel(r'$d_y$')

legend1=plt.legend([r'$\varepsilon=0$',r'$\varepsilon=0.8$',r'$\varepsilon=-0.8$',r'$\varepsilon=0.001$',r'$\varepsilon=-0.001$'], loc =1)
ax = plt.gca().add_artist(legend1)

legend2=plt.legend(['$c_0=1$','$a=1$',r'$k:\frac{-\pi}{a}<k<\frac{\pi}{a}$'], loc =4, handlelength=0)
ax = plt.gca().add_artist(legend2)

plt.axis([-0.008, 0.008, -0.1, 0.1])
plt.show()