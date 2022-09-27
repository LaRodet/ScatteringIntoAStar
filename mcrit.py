import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar, root, minimize_scalar

### Units:
# ap0 = 1. #initial semi-major axis of the planet
# mp = 1. #mass of the planet

### Parameters
rH = (1e-3/3)**(1/3) #dimensionless Hill radius
beta = 1 #extent of the chaotic zone accessible after scattering in unit of rH

verbose = False

a0 = 1 + 2*rH
e0 = 0.
i0 = 0.
ep0range = np.linspace(0., 0.5, 100)

ti = 1/a0 + 2*np.sqrt(a0*(1-e0*e0))*np.cos(i0) #Tisserant parameter
print(f"T = {ti}")

mcrit_smallm = np.zeros_like(ep0range)
mcrit = np.zeros_like(ep0range)

for j,ep0 in enumerate(ep0range):

    #small mass
    mcrit_smallm[j] = ep0*ep0/(ti*np.sqrt(1-ep0*ep0) - 2/(1-beta*rH))

    def f(x):
        h = np.sqrt(1-ep0*ep0) + 0.5*x*(ti-1/a0)
        e = 1 + x/a0
        return(h*h*e - 2*x/(1-beta*rH) - 1)

    sol = root_scalar(f, bracket=(0,1.), method="brentq")#x0=mcrit_smallm[j]
    if sol.converged:
        mcrit[j] = sol.root
    # if ep0>0:
    #     x = np.linspace(0, 0.5, 100)
    #     plt.plot(x, f(x))
    #     plt.axhline(y=0)
    #     plt.show()
    #     quit()

plt.plot(ep0range, mcrit, label=r"$\beta~r_{\rm H} = %.2f$" % (beta*rH))
plt.plot(ep0range, mcrit_smallm, "C0--", alpha=0.5)#, label="Small mass approximation"

### Parameters
rH = (1e-3/3)**(1/3) #dimensionless Hill radius
beta = 2 #extent of the chaotic zone accessible after scattering in unit of rH

verbose = False

a0 = 1 + 2*rH
e0 = 0.
i0 = 0.
ep0range = np.linspace(0., 0.5, 100)

ti = 1/a0 + 2*np.sqrt(a0*(1-e0*e0))*np.cos(i0) #Tisserant parameter
print(f"T = {ti}")

mcrit_smallm = np.zeros_like(ep0range)
mcrit = np.zeros_like(ep0range)

for j,ep0 in enumerate(ep0range):

    #small mass
    mcrit_smallm[j] = ep0*ep0/(ti*np.sqrt(1-ep0*ep0) - 2/(1-beta*rH))

    def f(x):
        h = np.sqrt(1-ep0*ep0) + 0.5*x*(ti-1/a0)
        e = 1 + x/a0
        return(h*h*e - 2*x/(1-beta*rH) - 1)

    sol = root_scalar(f, bracket=(0,1.), method="brentq")#x0=mcrit_smallm[j]
    if sol.converged:
        mcrit[j] = sol.root
    # if ep0>0:
    #     x = np.linspace(0, 0.5, 100)
    #     plt.plot(x, f(x))
    #     plt.axhline(y=0)
    #     plt.show()
    #     quit()

plt.plot(ep0range, mcrit, label=r"$\beta~r_{\rm H} = %.2f$" % (beta*rH))
plt.plot(ep0range, mcrit_smallm, "C1--", alpha=0.5)#, label="Small mass approximation"

### Parameters
rH = (1e-3/3)**(1/3) #dimensionless Hill radius
beta = 3 #extent of the chaotic zone accessible after scattering in unit of rH

verbose = False

a0 = 1 + 2*rH
e0 = 0.
i0 = 0.
ep0range = np.linspace(0., 0.5, 100)

ti = 1/a0 + 2*np.sqrt(a0*(1-e0*e0))*np.cos(i0) #Tisserant parameter
print(f"T = {ti}")

mcrit_smallm = np.zeros_like(ep0range)
mcrit = np.zeros_like(ep0range)

for j,ep0 in enumerate(ep0range):

    #small mass
    mcrit_smallm[j] = ep0*ep0/(ti*np.sqrt(1-ep0*ep0) - 2/(1-beta*rH))

    def f(x):
        h = np.sqrt(1-ep0*ep0) + 0.5*x*(ti-1/a0)
        e = 1 + x/a0
        return(h*h*e - 2*x/(1-beta*rH) - 1)

    sol = root_scalar(f, bracket=(0,1.), method="brentq")#x0=mcrit_smallm[j]
    if sol.converged:
        mcrit[j] = sol.root
    # if ep0>0:
    #     x = np.linspace(0, 0.5, 100)
    #     plt.plot(x, f(x))
    #     plt.axhline(y=0)
    #     plt.show()
    #     quit()

plt.plot(ep0range, mcrit, label=r"$\beta~r_{\rm H} = %.2f$" % (beta*rH))
plt.plot(ep0range, mcrit_smallm, "C2--", alpha=0.5)#, label="Small mass approximation"

plt.xlabel(r"$e_{\rm p,0}$")
plt.ylabel(r"$m_{\rm crit}/M_\star$")
plt.legend()
plt.ylim((0.,0.5))
plt.margins(x=0)

plt.tight_layout()
plt.show()
