import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
from scipy.integrate import quad
from scipy.integrate import dblquad

##constant

a=1
v=1
c0=1

##phi
def dx(k):
    return(c0*(1-v)+c0*(1+v)*np.cos(k*a))
def dy(k):
    return(c0*(1+v)*np.sin(k*a))
def phi(k):
    return(np.arctan(dy(k)/dx(k)))
##eignvector
def ephi(k):
    return(1/np.sqrt(2))*(np.exp(1j* phi(k)))
def ephi_dk(k):
    return derivative(ephi, k, dx=1e-6)
def conephi(k):
    return(np.conj(ephi(k)))
## Berry connection
def Berrycon(k):
    return 1j*ephi_dk(k)*conephi(k)
def berryphi(k):
    return quad(Berrycon, np.pi, -np.pi,limit=1000)[0]
print("The numerical result of berryphi  is {}".format(berryphi(y)))

## Berry curvature(1)
#def dphi_dy(k):
#    return derivative(phi, dy(k), dx=1e-6)
#def dy_dk(k):
#    return derivative(dy, k, dx=1e-6)
#def curA(k):
#    return -0.5*dphi_dy(k)*dy_dk(k)
#def dcurA(k):
#    return derivative(curA, dx(k), dx=1e-6)*dx_dk(k)

#def dphi_dx(k):
#    return derivative(phi, dx(k), dx=1e-6)
#def dx_dk(k):
#    return derivative(dx, k, dx=1e-6)
#def curB(k):
#    return -0.5*dphi_dx(k)*dx_dk(k)
#def dcurB(k):
#    return derivative(curB, dy(k), dx=1e-6)*dy_dk(k)

#def dcurAB(k):
#    return (dcurA(k)-dcurB(k))
#print(dcurAB(y))
## Berry curvature(2)
#def dA_dy(k):
    #return derivative(Berrycon, dy(k), dx=1e-6)
#def dA_dy_dx(k):
    #return derivative(dA_dy, dx(k), dx=1e-6)

#def dA_dx(k):
    #return derivative(Berrycon, dx(k), dx=1e-6)
#def dA_dx_dy(k):
    #return derivative(dA_dx, dx(k), dx=1e-6)

#def cur(k):
    #return dA_dy_dx(k)-dA_dx_dy(k)
## Berry curvature(3)
def dphi_dy(k):
    return derivative(phi, dy(k), dx=1e-6)
def dy_dk(k):
    return derivative(dy, k, dx=1e-6)
def curA(k):
    return -0.5*dphi_dy(k)*dy_dk(k)
def dcurA(k):
    return derivative(curA, k, dx=1e-6)

def dphi_dx(k):
    return derivative(phi, dx(k), dx=1e-6)
def dx_dk(k):
    return derivative(dx, k, dx=1e-6)
def curB(k):
    return -0.5*dphi_dx(k)*dx_dk(k)
def dcurB(k):
    return derivative(curB, k, dx=1e-6)
def dcurAB(k):
    return (dcurA(k)-dcurB(k))
#print(dcurAB(y))


## chern number
def int1(k):
    return (1/(np.pi))*quad(Berrycon, np.pi, -np.pi,limit=1000)[0]
print("The numerical result of 1st int1 is {:f}".format(int1(y)))


y=np.linspace(-2*np.pi/a,2*np.pi/a,1000)

plt.plot(y, Berrycon(y), color='blue', label=r'$Berrycon(k)$')
plt.plot(y, dcurAB(y), color='red', label=r'$dcurAB(k)$')

legend1=plt.legend(['$a=${}'.format(a),r'$\varepsilon=${}'.format(v),r'$k:\frac{-\pi}{a}<k<\frac{\pi}{a}$'], loc =1, handlelength=0)
ax = plt.gca().add_artist(legend1)

legend2=plt.legend(['berryphi is {:f}'.format(berryphi(y))], loc =4, handlelength=0)
ax = plt.gca().add_artist(legend2)

legend3=plt.legend(['Ch numberis {:f}'.format(int1(y))], loc =3, handlelength=0)
ax = plt.gca().add_artist(legend3)

plt.title("Plot of Berry Curvature")
plt.xlabel('1st Brillouin zone $ka$')
plt.ylabel(r'$\Omega{(k)}$')

plt.legend(loc='upper left')

plt.show()

plt.show()