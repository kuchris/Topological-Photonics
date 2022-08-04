import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.hermite as Herm

#Choose simple units
k=1.
h=1.
e=1.
B=1.
l=np.sqrt(h/e*B)

#Discretized space
dx = 0.05
x_lim = 12
x = np.arange(-x_lim,x_lim,dx)

def hermite(x, n):
    xi = (x-k*l**2)/l
    herm_coeffs = np.zeros(n+1)
    herm_coeffs[n] = 1
    return Herm.hermval(xi, herm_coeffs)

plt.figure()
plt.plot(x, hermite(x,0), linewidth=2,label=r"$H_0$")
plt.plot(x, hermite(x,1), linewidth=2,label=r"$H_1$")
plt.plot(x, hermite(x,2), linewidth=2,label=r"$H_2$")
plt.plot(x, hermite(x,3), linewidth=2,label=r"$H_3$")
plt.plot(x, hermite(x,4), linewidth=2,label=r"$H_4$")

#Set limits for axes
plt.xlim([-2.5,2.5])
plt.ylim([-20,20])

#Set axes labels
plt.xlabel("x")
plt.ylabel(r"$H_n(\xi)$")
plt.title(r"Hermite Polynomials, $H_n(\xi)$")
plt.legend()
plt.show()