###
# qmin as a function of m/mp for ep0 = 0 (Figure 2 of Rodet & Lai 2022)
###
from numpy import sqrt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar, root, minimize_scalar
from extrafunctions import *

### For the problem to have a solution:
# apmin < ap < ap1 or ap2 < ap

eps = 1e-10 #small quantity

### Units:
# ap0 = 1. #initial semi-major axis of the planet
# mp = 1. #mass of the planet

### Parameters
rH = (1e-3/3)**(1/3) #dimensionless Hill radius
beta = 3 #extent of the chaotic zone accessible after scattering in unit of rH

verbose = False

mrange = np.logspace(-4, -0.1, 100)
a0 = 1 + beta*rH
e0 = 0.
i0 = 0.
ep0range = [0., 0.02, 0.1, 0.2, 0.5]

### Test mass limit (Eq. 9)

ti = 1/a0 + 2*np.sqrt(a0*(1-e0*e0))*np.cos(i0) #Tisserant parameter
print(f"T = {ti}")

Qchaos = 1-beta*rH
apomin = 2/ti

if Qchaos > apomin:

    apo = Qchaos

    aa = ti*ti - 8*apo
    bb = -8*apo*apo - 4*ti + 2*apo*ti*ti
    delta = 64*apo*(2+(apo**3)-apo*ti)

    qmin = (-bb-np.sqrt(delta))/(2*aa)

else:

    qmin = 0

plt.axhline(y=qmin, ls='--', color="C0")

### Solution of B8 and B9
for ep0 in ep0range:

    qmin = np.zeros_like(mrange)
    for j,m in enumerate(mrange):

        h = sqrt(1-ep0*ep0) + m*sqrt(a0*(1-e0*e0))
        e = 1 + m/a0

        qmin0 = check_qmin0(m, beta*rH, e, h)

        if qmin0:
            qmin[j] = 0
        else:
            # print(ep0, m)

            apmin = 1/e
            apmin += eps
            apcrit = h*h #solution of ap=qp and q=0 (important for eq. B6)

            # Constraints on ap from eq. 11: mp sqrt(ap) + m sqrt(a) >= h
            def forbidden_ap(ap): # should be positive
                a = m/(e-1/ap) #eq. 10
                return(sqrt(ap)+m*sqrt(a)-h)

            # Here we use the fact that forbidden_ap is positive for ap = apmin and apcrit to find ap1 and ap2 such that forbidden_ap(ap)<0 for ap1 < ap < ap2
            sol = minimize_scalar(forbidden_ap, bracket=(apmin, apcrit), bounds=(apmin, apcrit), method='bounded')
            ap0 = sol.x
            if forbidden_ap(ap0) < 0:
                sol = root_scalar(forbidden_ap, bracket=(apmin, ap0), method="brentq")
                ap1 = sol.root
                sol = root_scalar(forbidden_ap, bracket=(ap0, apcrit), method="brentq")
                ap2 = sol.root
            else:
                ap1 = 1.
                ap2 = 1.

            apcross, qcross = find_rightapcross(m, beta*rH, e, h)
            if not np.isnan(apcross):
                if verbose:
                    print(f"qcross = {qcross}, apcross = {apcross}")

                qmin[j] = find_qmin(m, beta*rH, e, h, qcross, apcross, ap1, ap2, verbose)[1]
                # if np.isnan(qmin[j]):
                    # print(m)
                    # quit()
            else:
                if verbose:
                    print("Unable to find apcross")
                qmin[j] = np.nan

            # if verbose:
                # print(f"qmin = {qmin[j]}")

    plt.plot(mrange, qmin, label=r"$e_{\rm p,0} = %s$" % ep0)

plt.xscale('log')
plt.title(r"$\beta~r_{\rm H} = %.2f$" % (beta*rH))
plt.legend()
plt.margins(x=0)
plt.ylim(ymin=0)
plt.xlabel(r"$m/m_{\rm p}$")
plt.ylabel(r"$q_{\rm min}~[a_{\rm p,0}]$")
plt.tight_layout()

plt.show()
