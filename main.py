import numpy as np
import matplotlib.pyplot as plt

import solver


def f_a(t, x):
    return 1


def f_beg_smooth(x):
    return np.sin(2*3.1415*x)


def f_beg_discont(x):
    return 1 if x > 0.5 else -1


X, h, T, tau = 1, 0.01, 2, 0.01
f_beg = f_beg_smooth
print('T, tau, X, h =', T, tau, X, h)

s1 = solver.Solver(T, tau, X, h, solver.Wl, solver.Wr, f_a, f_beg)
s2 = solver.Solver(T, tau/2, X, h/2, solver.Wl, solver.Wr, f_a, f_beg)
s3 = solver.Solver(T, tau/4, X, h/4, solver.Wl, solver.Wr, f_a, f_beg)
print('courant:', s1.get_courant(1, 0.5))
s1.calculate()
print('s1: calculated')
s2.calculate()
print('s2: calculated')
s3.calculate()
print('s3: calculated')

a1 = solver.get_acc_field(T, tau, X, h, f_a, f_beg)
a2 = solver.get_acc_field(T, tau/2, X, h/2, f_a, f_beg)
a3 = solver.get_acc_field(T, tau/4, X, h/4, f_a, f_beg)

diff = a1-s1.get_field()
diff2 = a2-s2.get_field()
diff3 = a3-s3.get_field()


f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
ax1.contourf(s1.get_field(), levels=np.linspace(-1, 1, 101))
ax1.set_title('s1')
ax2.contourf(a1, levels=np.linspace(-1, 1, 101))
ax2.set_title('accurate')
ax3.contourf(diff)
ax3.set_title('difference')


print('1. max:', np.max(np.abs(diff)))
print('1. L1: ', np.sum(np.abs(diff)))
print('1. L1/cells:', np.sum(np.abs(diff)) / (1/h/tau))

print('2. max:', np.max(np.abs(diff2)))
print('2. L1: ', np.sum(np.abs(diff2)))
print('2. L1/cells:', np.sum(np.abs(diff2)) / (4/h/tau))

print('3. max:', np.max(np.abs(diff3)))
print('3. L1: ', np.sum(np.abs(diff3)))
print('3. L1/cells:', np.sum(np.abs(diff3)) / (16/h/tau))
plt.show()


x1 = np.linspace(0, 1, int(1/h))
x2 = np.linspace(0, 1, int(2/h))
x3 = np.linspace(0, 1, int(4/h))


f2, (bx1, bx2, bx3) = plt.subplots(1, 3, sharey=True)

#bx1.plot(x1, a1[-1], label='a1', color='orange')
bx1.plot(x1, (a1-s1.get_field())[-1], label='d21')
bx1.grid()
bx1.legend()

#bx2.plot(x2, a2[-1], label='a2', color='orange')
bx2.plot(x2, (a2-s2.get_field())[-1], label='d22')
bx2.grid()
bx2.legend()

#bx3.plot(x3, a3[-1], label='a3', color='orange')
bx3.plot(x3, (a3-s3.get_field())[-1], label='d23')
bx3.grid()
bx3.legend()
'''
bx3.plot(x3, (a3-s3.get_field())[-1], label='s3')
bx3.grid()
bx3.legend()
'''

plt.show()



f2, (bx1, bx2, bx3) = plt.subplots(1, 3, sharey=True)
bx1.plot(x1, s1.get_field()[-1], label='d21')
bx1.plot(x1, a1[-1], label='b21')
bx1.grid()
bx1.legend()

bx2.plot(x2, s2.get_field()[-1], label='d22')
bx2.plot(x2, a2[-1], label='b22')
bx2.grid()
bx2.legend()

bx3.plot(x3, s3.get_field()[-1], label='d23')
bx3.plot(x3, a3[-1], label='b23')
bx3.grid()
bx3.legend()
'''
bx3.plot(x3, (a3-s3.get_field())[-1], label='s3')
bx3.grid()
bx3.legend()
'''

plt.show()