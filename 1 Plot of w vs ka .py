import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure

m1=1
m2=1
c1=0.2
c2=1.8
a=1

def m12():
    return ((m1+m2)/(m1*m2));

def w1(k):
    return((((c1+c2)/2)*m12()+((((c1+c2)**2)/4)*(m12()**2)-(4*(c1*c2)/(m1*m2))*(np.sin(k*a/2)**2))**(1/2))**(1/2));

def w2(k):
    return((((c1+c2)/2)*m12()-((((c1+c2)**2)/4)*(m12()**2)-(4*(c1*c2)/(m1*m2))*(np.sin(k*a/2)**2))**(1/2))**(1/2));

k=np.linspace(-2*np.pi,2*np.pi,1000)

plt.figure(dpi=120)
plt.plot(k,w1(k),'r-')
plt.plot(k,w2(k),'b-')
plt.axvline(np.pi, color='g', linestyle=':')
plt.axvline(-np.pi, color='g', linestyle=':')
plt.text(np.pi,0,'$\pi$')
plt.text(-np.pi,0,'$-\pi$')

plt.title("Plot of $\omega$ versus $ka$")
plt.xlabel('1st Brillouin zone $ka$')
plt.ylabel(r'$\omega$')

#legend1=plt.legend(['Optical'r'$\sqrt{(\frac{c_1+c_2}{2})\frac{(M_1+M_2)}{M_1M_2}+\sqrt{(\frac{c_1+c_2}{2})^2(\frac{(M_1+M_2)}{M_1M_2})^2-\frac{4c_1c_2}{M_1M_2}\sin^2(\frac{ka}{2})}}$','Acoustic'r'$\sqrt{(\frac{c_1+c_2}{2})\frac{(M_1+M_2)}{M_1M_2}-\sqrt{(\frac{c_1+c_2}{2})^2(\frac{(M_1+M_2)}{M_1M_2})^2-\frac{4c_1c_2}{M_1M_2}\sin^2(\frac{ka}{2})}}$'],loc =1)

legend1=plt.legend(['Optical','Acoustic'], loc =1)
ax = plt.gca().add_artist(legend1)
legend2=plt.legend(['$M_1={}$'.format(m1),'$M_2={}$'.format(m2),'$c_1={}$'.format(c1),'$c_2={}$'.format(c2),'$a={}$'.format(a)], loc =5, handlelength=0)
ax = plt.gca().add_artist(legend2)


plt.show()

