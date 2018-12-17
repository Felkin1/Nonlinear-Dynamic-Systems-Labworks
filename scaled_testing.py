import numpy as np
from sympy import dsolve, Eq, symbols, Function
import matplotlib.pyplot as plt

def force_duffing(x):
    return x**3

def newmark(f,t,m,h,k,x0):
    '''Newmark integrator
    x0 = [3 x 1] array where the index represents the order of integration
    example - x0[0] = x, x0[1] = x', x0[2] = x''
    '''
    x = np.zeros(3) # prepare array for new values
    x[2] = (f - (h*t/2 + k*t*t/4) * x0[2] - (h+k*t) * x0[1] - k*x0[0]) / (m + h*t/2 + k*t*t/4)
    x[1] = x0[1] + (x0[2]+x[2]) * t/2
    x[0] = x0[0] + x0[1]*t + (x0[2] + x[2]) * t*t/4

    return x


a = 0.9
x = a
dx = 0
ddx = -x + x**3
t = 0
dt = 0.01
iterations = 1000
xs = []
m,h,k = 1,0,1
t_arr = []
t_arr.append(t)
xs.append([x,dx,ddx])
for j in range(iterations):
    f = force_duffing(xs[j][0])
    x_new = newmark(f,dt,m,h,k,xs[j])
    t+=dt
    xs.append(x_new)
    t_arr.append(t)
xs = np.array(xs)
plt.plot(t_arr, xs[:, 0],'g')


# x0(tau) = A*cos(tau)
# tau = w*t
# w = 1 + k1
# k1 = 3*A^2/8 (see textbook)

k1 = 3*(a**2)/8;
w = 1 + k1;

t_arr = np.array([t_arr])
tau = w*t_arr
X0 = a*np.cos(t_arr)
X0_bad = a*np.cos(tau)
print(X0[:10])
plt.plot(tau,X0,'r-',t_arr,X0_bad,'b--')


#X1 = -a^3*(cos(3*T)-cos(T))/32 ;%+ a*cos(T);
#X1_bad = -a^3*(cos(3*tau)-cos(tau))/32 ;%+ a*cos(tau);

#plot(tau,X1,'r',T,X1_bad,'r:')
plt.show()