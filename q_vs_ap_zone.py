###
# Figure B1 of Rodet & Lai 2022
# Scattering between a planetesimal and a planet
# Possible post-scattering values of the planetesimal periastron q and the planet semi-major axis ap
###
from numpy import sqrt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar, root, minimize_scalar
from extrafunctions import *

eps = 1e-10 #small quantity

### Units:
# ap0 = 1. #initial semi-major axis of the planet
# mp = 1. #mass of the planet

### Parameters
rH = (1e-3/3)**(1/3) #dimensionless Hill radius
beta = 2 #extent of the chaotic zone accessible after scattering in unit of rH
m = 0.1 #mass of the planetesimal

### Initial conditions
a0 = 1 + 2*rH #initial semi-major axis of the planetesimal
e0 = 0. #initial eccentricity of the planetesimal
ep0 = 0. #initial eccentricity of the planet

plt.plot(1., a0*(1-e0), "ko") #plot initial condition

e = 1 + m/a0 #energy (eq. 10)
h = sqrt(1-ep0*ep0) + m*sqrt(a0*(1-e0*e0)) #angular momentum (eq. 11)
print(f"E = {e}")
print(f"h = {h}")

### Simple constraints

apmin = 1/e #minimum possible ap (from eq. 10, corresponding to a=+inf)
print(f"ap,min = {apmin}")
plt.axvline(x=apmin, color='k')
plt.text(apmin, 1., r"$E < m_{\rm p}/a_{\rm p}$", color='#3D1C02', ha='right', va='top', rotation=90)
apmin += eps #to avoid computational divergence

apcrit = h*h #solution of ap=qp and q=0 (important for eq. B6)
print(f"ap,crit = {apcrit}")
plt.axvline(x=apcrit, color="k", ls='--')
plt.text(apcrit, 1, r"$a_{\rm p} = h^2/m_{\rm p}^2$", va="center", rotation=90)

# plotting limits
xmin = 1 - 1.5*(1-apmin)
xmax = 1. + 1.5*(apcrit-1.)
plt.xlim((xmin, xmax))
plt.axvspan(xmin, apmin, alpha=0.2, color='C8')

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
plt.axvline(x=ap1, color='k')
plt.axvline(x=ap2, color='k')
plt.text(0.5*(ap1+ap2), 1., r"$h > m_{\rm p}\sqrt{a_{\rm p}} + m\sqrt{a}$", color='#580F41', ha='center', va='top', rotation=90)
plt.axvspan(ap1, ap2, alpha=0.2, color='C4')
ap1 -= eps
ap2 += eps

#Conclusion
print(f"{apmin} < ap < {ap1}")
print(f"or ap > {ap2}")

### Further constraints

# Branch qp = ap for ap2<ap<apcrit
ap = np.linspace(ap2, apcrit, 100)
a = m/(e-1/ap)
aux = ((h-np.sqrt(ap))/m)**2 / a
q = a[aux<1]*(1-np.sqrt(1-aux[aux<1]))
plt.plot(ap[aux<1], q, "C0")
plt.fill_between(ap[aux<1], 0, q, color="C0", alpha=0.2)

# Branch qp = ap for apmin<ap<ap1
ap = np.linspace(apmin, ap1, 100)
a = m/(e-1/ap)
aux = ((h-np.sqrt(ap))/m)**2 / a
q = a[aux<1]*(1-np.sqrt(1-aux[aux<1]))
plt.plot(ap[aux<1], q, "C0")
plt.fill_between(ap[aux<1], 0, q, color="C0", alpha=0.2)

ymax = 1.8*q.max()
plt.ylim((0, ymax))

plt.text(1.01*ap2, 0, r"$q_{\rm p} > a_{\rm p}$", color="C0", va='bottom')

# Branch q = a
ap = np.linspace(apmin, xmax, 1000)
a = m/(e-1/ap)
plt.plot(ap[ap<ap1], a[ap<ap1], "k")
plt.fill_between(ap[ap<ap1], a[ap<ap1], ymax, color="k", alpha=0.2)
plt.plot(ap[ap>ap2], a[ap>ap2], "k")
plt.fill_between(ap[ap>ap2], a[ap>ap2], ymax, color="k", alpha=0.2)
plt.text(0.5*(ap2+apcrit), 1.25, r"$q > a$", color='k', ha="right", va="bottom")

# Branches qp = Q/(1-beta*rH)
apcross, qcross = find_rightapcross(m, beta*rH, e, h) #find ap such that qp = ap = Q/(1-beta*rH)
plt.plot(apcross, qcross, "C0o")

ap, q, ap2b, q2 = compute_rightbranch(m, beta*rH, e, h, apcross, qcross) #two solutions for ap>ap2

a = m/(e-1/ap)
plt.plot(ap, q, "C1", label=r"$q_{\rm p} = Q/(1-\beta r_{\rm H})$")
plt.fill_between(ap, q, a, color="C1", alpha=0.2, linewidth=0)
plt.text(ap[-1], q[-1], r"$Q < q_{\rm p}(1-\beta r_{\rm H})$", color='C1', ha="left")

if len(ap2b) > 1:
    a2 = m/(e-1/ap2b)
    plt.plot(ap2b, q2, "C1")

    aux = ((h-np.sqrt(ap2b))/m)**2 / a2
    qbis = a2[aux<1]*(1-np.sqrt(1-aux[aux<1]))

    plt.fill_between(ap2b[aux<1], qbis, q2[aux<1], color="C1", alpha=0.2, linewidth=0)

# Color the region where there is no solution to qp = Q/(1-beta*rH)
ap = np.linspace(ap[-1], xmax, 100)
a = m/(e-1/ap)

aux = ((h-np.sqrt(ap))/m)**2 / a
aux[ap>apcrit] = 0
qbis = a[aux<1]*(1-np.sqrt(1-aux[aux<1])) # q such that qp=ap

plt.fill_between(ap[aux<1], qbis, a[aux<1], color="C1", alpha=0.2, linewidth=0)

#Find qmin numerically (eqs. B8 & B9)

qmin0 = check_qmin0(m, beta*rH, e, h)
if qmin0:
    qmin = 0
    print("qmin = 0")
else:
    apqmin, qmin, success = find_qmin(m, beta*rH, e, h, qcross, apcross, ap1, ap2, True)
    if success:
        plt.plot(apqmin, qmin, "ro")

# Branches Qp = q/(1+beta*rH)
apcross, qcross = find_leftapcross(m, beta*rH, e, h, ap1) #find ap such that qp = ap and q = Qp(1+beta*rH)

if apcross>0:

    ap, q, ap2b, q2 = compute_leftbranch(m, beta*rH, e, h, apcross, qcross, ap1) #two solutions for ap<ap1

    a = m/(e-1/ap)

    plt.plot(ap, q, "C2", label=r"$q = Q_{\rm p}(1+\beta r_{\rm H})$")
    plt.fill_between(ap, q, a, color="C2", alpha=0.2, linewidth=0)
    plt.text(ap[-1], ymax, r"$q > Q_{\rm p}(1+\beta r_{\rm H})$", color='g', rotation=90, va='top')

    a = m/(e-1/ap2b)

    plt.plot(ap2b, q2, "C2")

    aux = ((h-np.sqrt(ap2b))/m)**2 / a
    qbis = a[aux<1]*(1-np.sqrt(1-aux[aux<1]))
    plt.fill_between(ap2b[aux<1], qbis, q2[aux<1], color="C2", alpha=0.2, linewidth=0)

    #Color zone when no solution
    ap = np.linspace(apmin, ap[-1]-eps, 100)
    a = m/(e-1/ap)

    aux = ((h-np.sqrt(ap))/m)**2 / a
    qbis = a[aux<1]*(1-np.sqrt(1-aux[aux<1]))

    plt.fill_between(ap[aux<1], qbis, a[aux<1], color="C2", alpha=0.2, linewidth=0)

plt.xlabel(r"$a_{\rm p}~[a_{\rm p,0}]$")
plt.ylabel(r"$q~[a_{\rm p,0}]$ ")
plt.tight_layout()

plt.show()
