import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
from scipy.integrate import quad
from scipy.integrate import dblquad

a=1
v=1
c0=1
def dx(k):
    return(c0*(1-v)+c0*(1+v)*np.cos(k*a))
def dy(k):
    return(c0*(1+v)*np.sin(k*a))
def phi(k):
    return(np.arctan(dy(k)/dx(k)))

def f(k):
    return np.sqrt(0.5)*np.exp(1j*(phi(k)))

def df(k):
    return derivative(f, k, dx=1e-6)

#def df(k):
    #return -(1j*a*(v+1)*np.exp(1j*a*k)*((v-1)*np.exp(2*1j*a*k)+(-2*v-2)*np.exp(1j*a*k)+v-1))/(2*((v-1)*np.exp(1j*a*k)-v-1)**2*np.sqrt(((v+1)*np.exp(1j*a*k)-v+1)/((v+1)*np.exp(-1j*a*k)-v+1)))

def conf(k):
    return np.conj(f(k))

#def conf(k):
    #return np.sqrt(((1+v)*np.exp(-1j*k*a)+(1-v))/((1+v)*np.exp(1j*k*a)+(1-v)))

#berry connection
def confdf(k):
    return 1j*conf(k)*df(k)

#integrate
res, err = quad(confdf, np.pi, -np.pi,limit=1000)
print("The numerical result is {:f} (+-{:g})"
    .format(res, err))

def dconfdf(k):
    return derivative(confdf, k, dx=1e-6)
#print(dconfdf(y))


def c1dconfdf(k):
    return quad(dconfdf, -np.pi, np.pi,limit=1000)[0]

print("The numerical result of 1st int is {:f}".format(c1dconfdf(y)))

def c2dconfdf(k):
    return (1/(2*np.pi))*quad(c1dconfdf, -np.pi, np.pi,limit=1000)[0]

print("The numerical result of 2st int is {:f}".format(c2dconfdf(y)))



y=np.linspace(-2*np.pi/a,2*np.pi/a,1000)

plt.axvline(np.pi, color='r', linestyle=':')
plt.axvline(-np.pi, color='r', linestyle=':')
plt.text(np.pi,-0.5,'$\pi$')
plt.text(-np.pi,-0.5,'$-\pi$')

#plt.plot(y, f(y), color='purple', label=r'$f(k)$')
#plt.plot(y, df(y), color='green', label=r'$\frac{\partial f(k)}{\partial k}$')
#plt.plot(y, confdf(y), color='black', label=r'$i f(k)^*\frac{\partial f(k)}{\partial k} $')

#plt.fill_between(y, confdf(y), where=[(y > -np.pi) and (y < np.pi) for y in y],label=r'$\gamma =i\oint f(k)^*\frac{\partial f(k)}{\partial k} dk$')

#plt.plot(y, dconfdf(y), color='red', label=r'$\Omega{(k)}$')

plt.plot(y, c1dconfdf(y), color='blue', label=r'$\1{(k)}$')

#plt.fill_between(y, dconfdf(y), where=[(y > -np.pi) and (y < np.pi) for y in y],label=r'$ Berry curvature$')

legend1=plt.legend(['$a=${}'.format(a),r'$\varepsilon=${}'.format(v),r'$k:\frac{-\pi}{a}<k<\frac{\pi}{a}$'], loc =1, handlelength=0)
ax = plt.gca().add_artist(legend1)

#legend2=plt.legend(['The numerical result of integration is {}'.format(res)], loc =4, handlelength=0)
#ax = plt.gca().add_artist(legend2)

plt.title("Plot of Berry Curvature")
plt.xlabel('1st Brillouin zone $ka$')
plt.ylabel(r'$\Omega{(k)}$')

plt.legend(loc='upper left')
plt.show()