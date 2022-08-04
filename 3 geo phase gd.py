import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
from scipy.integrate import quad

a=1
v=-1
c0=1
#t=(1-v)/(1+v)

def dx(k):
    return(c0*(1-v)+c0*(1+v)*np.cos(k*a))
def dy(k):
    return(c0*(1+v)*np.sin(k*a))
#def dxdy(k):
    #return((t/np.sin(k*a))+1/np.tan(k*a))

def f(k):
    return np.exp(1j*(np.arctan(dy(k)/dx(k))))
#def phi(k):
    #return(1/(np.arctan(dxdy(k))))

#def f(k):
    #return np.exp(-1j*phi(k))

def df(k):
    return derivative(f, k,dx=1e-6)

def conf(k):
    return np.conj(f(k))

def confdf(k):
    return 0.5*1j*conf(k)*df(k)

#integrate
res, err = quad(confdf, np.pi, -np.pi,limit=1000)
print("The numerical result is {:f} (+-{:g})"
    .format(res, err))

y=np.linspace(-2*np.pi/a,2*np.pi/a,1000)

plt.axvline(np.pi, color='g', linestyle=':')
plt.axvline(-np.pi, color='g', linestyle=':')
plt.text(np.pi,-0.5,'$\pi$')
plt.text(-np.pi,-0.5,'$-\pi$')

plt.plot(y, f(y), color='purple', label=r'$f(k)$')
plt.plot(y, df(y), color='green', label=r'$\frac{\partial f(k)}{\partial k}$')
plt.plot(y, confdf(y), color='black', label=r'$i f(k)^*\frac{\partial f(k)}{\partial k} $')
plt.fill_between(y, confdf(y), where=[(y > -np.pi) and (y < np.pi) for y in y],label=r'$\gamma =i\oint f(k)^*\frac{\partial f(k)}{\partial k} dk$')

legend1=plt.legend(['$a=${}'.format(a),r'$\varepsilon=${}'.format(v),r'$k:\frac{-\pi}{a}<k<\frac{\pi}{a}$'], loc =1, handlelength=0)
ax = plt.gca().add_artist(legend1)

legend2=plt.legend(['The numerical result of integration is {}'.format(res)], loc =4, handlelength=0)
ax = plt.gca().add_artist(legend2)

#plt.plot(y, np.imag(f(y)), color='purple', label='Function', linestyle=':')
#plt.plot(y, np.imag(df(y)), color='green', label='Derivative', linestyle=':')
plt.title("Plot of Zak phase using $v=e^{i\Phi}$")
plt.xlabel('1st Brillouin zone $ka$')
plt.ylabel(r'$arb. unit$')

plt.legend(loc='upper left')
plt.show()