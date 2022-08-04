import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
from scipy.integrate import quad

a=1
v=-1
def f(k):
    return np.sqrt(((1+v)*np.exp(1j*k*a)+(1-v))/((1+v)*np.exp(-1j*k*a)+(1-v)))

def df(k):
    return derivative(f, k, dx=1e-6)

#def df(k):
    #return -(1j*a*(v+1)*np.exp(1j*a*k)*((v-1)*np.exp(2*1j*a*k)+(-2*v-2)*np.exp(1j*a*k)+v-1))/(2*((v-1)*np.exp(1j*a*k)-v-1)**2*np.sqrt(((v+1)*np.exp(1j*a*k)-v+1)/((v+1)*np.exp(-1j*a*k)-v+1)))

def conf(k):
    return np.conj(f(k))

#def conf(k):
    #return np.sqrt(((1+v)*np.exp(-1j*k*a)+(1-v))/((1+v)*np.exp(1j*k*a)+(1-v)))

def confdf(k):
    return 0.5*1j*conf(k)*df(k)

#integrate
res, err = quad(confdf, np.pi, -np.pi,limit=1000)
print("The numerical result is {:f} (+-{:g})"
    .format(res, err))

y=np.linspace(-2*np.pi/a,2*np.pi/a,1000)

plt.axvline(np.pi, color='r', linestyle=':')
plt.axvline(-np.pi, color='r', linestyle=':')
plt.text(np.pi,-0.5,'$\pi$')
plt.text(-np.pi,-0.5,'$-\pi$')

plt.plot(y, f(y), color='purple', label=r'$f(k)$')
plt.plot(y, df(y), color='green', label=r'$\frac{\partial f(k)}{\partial k}$')
plt.plot(y, confdf(y), color='black', label=r'$i f(k)^*\frac{\partial f(k)}{\partial k} $')
plt.fill_between(y, confdf(y), where=[(y > -np.pi) and (y < np.pi) for y in y],label=r'$\gamma =i\oint f(k)^*\frac{\partial f(k)}{\partial k} dk$')

#plt.plot(y, np.imag(f(y)), color='purple', label='Function', linestyle=':')
#plt.plot(y, np.imag(df(y)), color='green', label='Derivative', linestyle=':')
#plt.plot(y, np.imag(confdf(y)), color='green', label='phi', linestyle=':')

legend1=plt.legend(['$a=${}'.format(a),r'$\varepsilon=${}'.format(v),r'$k:\frac{-\pi}{a}<k<\frac{\pi}{a}$'], loc =1, handlelength=0)
ax = plt.gca().add_artist(legend1)

legend2=plt.legend(['The numerical result of integration is {}'.format(res)], loc =4, handlelength=0)
ax = plt.gca().add_artist(legend2)

plt.title("Plot of Zak phase")
plt.xlabel('1st Brillouin zone $ka$')
plt.ylabel(r'$arb. unit$')

plt.legend(loc='upper left')
plt.show()