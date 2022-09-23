###
# Figure B1 of Rodet & Lai 2022
# Scattering between a planetesimal and a planet
# Possible post-scattering values of the planetesimal periastron q and the planet semi-major axis ap
###
from numpy import sqrt
import numpy as np
from scipy.optimize import root_scalar, root, minimize_scalar

eps = 1e-10 #small quantity

def check_qmin0(m, brH, e, h):
    # qmin0 = True if qmin = 0
    qmin0 = (h*h*e <  (2*m/(1-brH) + 1.)) # eq B17
    if not qmin0:
        x1 = 2*m/((1-brH)*e*h*h)
        x2 = 1/(e*h*h)
        delta = ((x1+2*x2)**2 - 4*x2)*x1*x1 # eq B22
        qmin0 = (delta > 0)
        if qmin0:
            xroot = (x1*x1 + 2*x1*x2 - 2*x2 + np.sqrt(delta))/(2*(2*x1-1.))
            qmin0 = (xroot > 1.) # eq B23

    return(qmin0)

def find_qmin(m, brH, e, h, qcross, apcross, ap1, ap2, verbose):
    ''' find the minimum possible q (solution of eqs. B8 and B9)'''

    apmin = 1/e

    def f(q, ap, a, qp): # the minimum q is the root of f for fixed ap and a (eq B8)
        return(h-(sqrt(qp*(2-qp/ap))+m*sqrt(q*(2-q/a))))

    def df(q, ap, a, qp): # eq B9
        dx = (qp/ap)**2 - 4*(a*qp)/(ap*ap)/m*(1-qp/ap)/(2-q/a)
        dy = -(q*q)/(m*ap*ap)
        return(m*dy*sqrt(qp*(2-qp/ap)) + dx*sqrt(q*(2-q/a)))

    def fvec(x): #2D function to find the root of
        [q, ap] = x
        a = m/(e-1/ap)
        qp = (2*a-q)/(1-brH)
        if (q>a) or (qp>ap) or (q<0) or (qp<0) or (ap<0) or (a<0) or (ap<apmin) or ((ap>ap1) and (ap<ap2)):
            aux = [1.,1.]
        else:
            aux = [f(q, ap, a, qp), df(q, ap, a, qp)]
        return(aux)

    sol = root(fvec, [qcross, apcross])
    success = sol.success
    if success:
        if verbose:
            print(f"1, qmin = {sol.x[0]}")
        apqmin = sol.x[1]
        qmin = sol.x[0]
    else:
        apcrossmax = (2*m/(1-brH)+1)/e
        sol = root(fvec, [qcross, apcross + 0.5*(apcrossmax-apcross)])
        success = sol.success
        if success:
            if verbose:
                print(f"2, qmin = {sol.x[0]}")
            apqmin = sol.x[1]
            qmin = sol.x[0]
        else:
            apqmin = np.nan
            qmin = np.nan
            if verbose:
                print("Root failed to find qmin ")
    return(apqmin, qmin, success)

def find_rightapcross(m, brH, e, h):
    '''find ap such that qp = ap = Q/(1-brH)'''

    def cross(ap):
        qp = ap
        a = m/(e-1/ap)
        q = 2*a - qp*(1-brH)
        if (qp*(2-qp/ap) < 0) or (q*(2-q/a) < 0):
            print("Error in function cross")
            quit()
        else:
            return(h-(sqrt(qp*(2-qp/ap))+m*sqrt(q*(2-q/a))))

    apcrossmin = (m/(1-brH)+1)/e # min ap such that qp = Q/(1-brH) have a solution
    apcrossmax = (2*m/(1-brH)+1)/e # apcross < apcrossmax
    apcrossmin += eps
    apcrossmax -= eps


    if cross(apcrossmin)*cross(apcrossmax) < 0:
        sol = root_scalar(cross, bracket=(apcrossmin,apcrossmax), method="brentq")
        apcross = sol.root
        qpcross = apcross
        across = m/(e-1/apcross)
        qcross = 2*across - qpcross*(1-brH)
        # print(f"apcross = {apcross}")
    else:
        print("Unable to find apcross")
        quit()
    return(apcross, qcross)

def compute_rightbranch(m, brH, e, h, apcross, qcross):
    ''' find solution of qp = Q/(1-brH) for ap > ap2'''

    def f(q, ap, a): # the minimum q is the root of f for fixed ap and a (eq. 11 and B6)
        qp = (2*a-q)/(1-brH)
        if (qp*(2-qp/ap) < 0) or (q*(2-q/a) < 0):
            print("Error in function f")
            quit()
        else:
            return(h-(sqrt(qp*(2-qp/ap))+m*sqrt(q*(2-q/a))))

    #points of reference
    apcrossmin = (m/(1-brH)+1)/e
    apcrossmax = (2*m/(1-brH)+1)/e
    apmin = 1/e
    apcrit = h*h
    apcrossmin += eps

    #we begin the search for a solution at apcrossmin, and increment ap by dap at everystep. When a solution is found, found=True, and the algorithm will stop when no solution is found.
    stop = False
    ap = [apcrossmin]
    ap2b = [apcross]
    q = []
    q2 = [qcross]
    dap = (1-apmin)/1000.
    found = False
    while not stop:
        a = m/(e-1/ap[-1])
        qmin = max(2*a-ap[-1]*(1-brH),0)#+eps
        qmax = a#-eps
        if f(qmin, ap[-1], a)*f(qmax, ap[-1], a) < 0: ### One root
            sol = root_scalar(f, bracket=(qmin,qmax), method="brentq", args=(ap[-1], a))
            success = sol.converged
            if success:
                q += [sol.root]
        else: ### Two roots
            sol = minimize_scalar(f, bracket=(qmin,qmax), bounds=(qmin,qmax), method='bounded', args=(ap[-1], a))
            min = sol.x
            if f(qmin, ap[-1], a)*f(min, ap[-1], a) < 0:
                sol = root_scalar(f, bracket=(min, qmax), method="brentq", args=(ap[-1], a))
                success = sol.converged
                if success:
                    q += [sol.root]
                sol = root_scalar(f, bracket=(qmin, min), method="brentq", args=(ap[-1], a))
                success = (success and sol.converged)
                if success:
                    q2 += [sol.root]
                    ap2b += [ap[-1]]
            else:
                success = False

        if not success:
            if found:
                print(f"Branch qp = Q/(1-beta rH) ends at ap = {ap[-1]}")
                ap.pop()
                stop = True
            else:
                apnew = ap[-1]+dap
                ap.pop()
                ap += [apnew]
        else:
            if not found:
                found = True
            if ap[-1] > 10*apcrit:
                stop = True
                print(f"Single branche does not stop")
            else:
                apnew = ap[-1]+dap
                ap += [apnew]

    return(np.array(ap), np.array(q), np.array(ap2b), np.array(q2))


def find_leftapcross(m, brH, e, h, ap1):
    '''find ap such that qp = ap and q = Qp/(1+brH)'''

    def cross(ap):
        qp = ap
        a = m/(e-1/ap)
        q = (1+brH)*(2*ap-qp)
        if (qp*(2-qp/ap) < 0) or (q*(2-q/a) < 0):
            print("Error in function cross left")
            quit()
        else:
            return(h-(sqrt(qp*(2-qp/ap))+m*sqrt(q*(2-q/a))))

    apmin = 1/e
    apcrossmin = apmin
    apcrossmax = np.minimum((m/(1+brH)+1)/e, ap1) # max ap such that q = Qp/(1+brH) have a solution
    apcrossmin += eps
    apcrossmax -= eps

    if cross(apcrossmin)*cross(apcrossmax) < 0:
        sol = root_scalar(cross, bracket=(apcrossmin,apcrossmax), method="brentq")
        apcross = sol.root
        qpcross = apcross
        qcross = (1+brH)*(2*apcross-qpcross)
        print(f"apcross = {apcross}")
    else:
        print("Unable to find apcross left")
        apcross = -1
        qcross = -1
    return(apcross, qcross)

def compute_leftbranch(m, brH, e, h, apcross, qcross, ap1):
    ''' find solution of q = Qp/(1+brH) for ap < ap1'''

    eps = 1e-10

    def f(q, ap, a):
        qp = 2*ap-q/(1+brH)
        if (qp*(2-qp/ap) < 0) or (q*(2-q/a) < 0):
            print("Error in function f left")
            quit()
        else:
            return(h-(sqrt(qp*(2-qp/ap))+m*sqrt(q*(2-q/a))))

    #points of reference
    apmin = 1/e
    apcrossmin = apmin
    apcrossmax = np.minimum((m/(1+brH)+1)/e, ap1) # max ap such that q = Qp/(1+brH) have a solution
    apcrossmin += eps
    apcrossmax -= eps
    apmin += eps

    #we begin the search for a solution at apcrossmin, and increment ap by dap at everystep. When a solution is found, found=True, and the algorithm will stop when no solution is found.
    stop = False
    ap = [apcrossmax]
    ap2b = [apcross]
    q = []
    q2 = [qcross]
    dap = (1-apmin)/1000.
    found = False
    while not stop:
        if ap[-1] > apcrossmin:
            a = m/(e-1/ap[-1])
            qmin = ap[-1]*(1+brH)
            qmax = np.minimum(2*ap[-1]*(1+brH), a)-eps
            if f(qmin, ap[-1], a)*f(qmax, ap[-1], a) < 0: ### One root
                sol = root_scalar(f, bracket=(qmin,qmax), method="brentq", args=(ap[-1], a))
                success = sol.converged
                if success:
                    q += [sol.root]
            else: ### Two roots
                sol = minimize_scalar(f, bracket=(qmin,qmax), bounds=(qmin,qmax), method='bounded', args=(ap[-1], a))
                min = sol.x
                if f(qmin, ap[-1], a)*f(min, ap[-1], a) < 0:
                    sol = root_scalar(f, bracket=(min, qmax), method="brentq", args=(ap[-1], a))
                    success = sol.converged
                    if success:
                        q += [sol.root]
                    sol = root_scalar(f, bracket=(qmin, min), method="brentq", args=(ap[-1], a))
                    success = (success and sol.converged)
                    if success:
                        q2 += [sol.root]
                        ap2b += [ap[-1]]

                else:
                    success = False

            if not success:
                if found:
                    print(f"ap = {ap[-1]} is the min ap for which q = Qp(1+beta rH)")
                    ap.pop()
                    stop = True
                else:
                    apnew = ap[-1]-dap
                    ap.pop()
                    ap += [apnew]
            else:
                if not found:
                    found = True
                apnew = ap[-1]-dap
                ap += [apnew]
        else:
            ap.pop()
            stop = True

    return(np.array(ap), np.array(q), np.array(ap2b), np.array(q2))
