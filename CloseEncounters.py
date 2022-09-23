import numpy as np
import matplotlib.pyplot as plt
import rebound
from rebound import hash as h

dr = np.pi/180
rJ = 7e4/1.5e8 # Jupiter radius (au)

### Simulation

ms = 1.
mp = 0.001

ap = 1.
ep = 0.3

tp = np.sqrt(ap**3/(ms+mp))

radiusp = 10*rJ

mu = mp/(ms+mp)

rH = (mu/3)**(1/3)

ntp = 1000

aa = 1.15
bb = 2.4*((mp/(ms+mp))**(1/3))
cc = ep-1
dain = 1-((-bb+np.sqrt(bb*bb-4*aa*cc))/(2*aa))**2

aa = 1/(1+ep)
bb = -2.4*((mp/(ms+mp))**(1/3))
cc = -1.15
daout = ((-bb+np.sqrt(bb*bb-4*aa*cc))/(2*aa))**2-1

a = ap*(1-dain + (daout+dain)*np.random.random(ntp))
e = np.random.random(ntp)*0.05
f = np.random.random(ntp)*2*np.pi
om = np.random.random(ntp)*2*np.pi
bigom = np.random.random(ntp)*2*np.pi
ci = np.random.random(ntp)*(1-np.cos(2*dr))+np.cos(2*dr)
inc = np.arccos(ci)

dmin_est = -np.ones(ntp, dtype='float')
tmin = np.zeros(ntp)
teject = np.zeros(ntp)
tcol = np.zeros(ntp)
ce = [] # list of particles number, dmin and t where d < 0.1 ap

threshold = 0.1*ap

dmax = 100*ap

n = 2*np.pi

tinteg = 1e6*tp
nt = int(tinteg)*10
t = np.linspace(0, tinteg, nt)
# dt = np.zeros(nt)
ngr = 10
np.savetxt("ci.dat", [a, e, inc, om, bigom, f])
c = 0

hs = h('star').value
hp = h('planet').value

for k in (range(int(ntp/ngr))):

    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')

    sim.add(m=ms, hash='star')
    sim.add(m=mp, a=ap, e=ep, hash='planet', r=radiusp)
    for l in range(ngr):
        n = ngr*k+l
        sim.add(primary=sim.particles["star"], a=a[n], e=e[n], f=f[n], omega=om[n], Omega=bigom[n], inc=inc[n], hash=l)
    sim.integrator = "ias15"
    sim.collision = "line"

    def collision_resolve(sim_pointer, collision):
        global tcol, tscat
        sim = sim_pointer.contents
        if (sim.particles[collision.p1].hash.value == hp):
            n = ngr*k+sim.particles[collision.p2].hash.value
            tcol[n] = sim.t
            return 2
        else:
            n = ngr*k+sim.particles[collision.p1].hash.value
            tcol[n] = sim.t
            return 1

    sim.collision_resolve = collision_resolve #"merge"

    sim.N_active = 2

    sim.move_to_com()

    ejectionchecked = np.zeros(ngr, dtype="bool")
    dminchecked = np.zeros(ngr, dtype="bool")

    for j,tt in enumerate(t):
        if np.any(teject==0):
            sim.integrate(tt)
            # dt[j] = sim.dt
            s = sim.particles["star"]
            p = sim.particles["planet"]
            for l in range(ngr):
                n = ngr*k+l
                if (teject[n] == 0) and (tcol[n] == 0):
                    tp = sim.particles[h(l)]
                    rp = np.sqrt((p.x-tp.x)**2+(p.y-tp.y)**2+(p.z-tp.z)**2)
                    rs = np.sqrt((s.x-tp.x)**2+(s.y-tp.y)**2+(s.z-tp.z)**2)
                    if (dmin_est[n] == -1.) or (rs < dmin_est[n]):
                        dmin_est[n] = rs
                    if (not ejectionchecked[l]) and (rs > dmax):
                        orbit = tp.calculate_orbit()
                        etp = orbit.e
                        # print(etp,tt)
                        if etp > 1:
                            teject[n] = tt
                            sim.remove(hash=h(l))
                            # print(f"tp {l} ejected, time {tt} yr, min dist {dmin_est[n]}")
                        else:
                            ejectionchecked[l] = True
                    elif ejectionchecked[l] and (rs < dmax):
                        ejectionchecked[l] = False
                    if (rp > 3*rH*ap) and (rs<ap*(1-ep)):
                        orbit = tp.calculate_orbit(primary=sim.particles["star"])
                        ftp = orbit.f
                        if ftp > np.pi:
                            etp = orbit.e
                            atp = orbit.a
                            peritp = atp*(1-etp)
                            if (peritp < dmin_est[n]):
                                dmin_est[n] = peritp
                                tmin[n] = tt
                            if (not dminchecked[l]) and (peritp < threshold):
                                ce += [[n, peritp, tt]]
                                dminchecked[l] = True
                            elif dminchecked[l] and (peritp > threshold):
                                dminchecked[l] = False
                        elif dminchecked[l]:
                            dminchecked[l] = False
                    elif dminchecked[l]:
                        dminchecked[l] = False
                    # if tt == t[-1]:
                    #     print(tp.x, tp.y, tp.z)
                    #     print(tp.a, tp.e, tp.inc)


    c += ngr
    countfile = open("count.txt", "w")
    countfile.write(f"{c}/{ntp}")
    countfile.close()

np.savetxt("dmin.dat", dmin_est)
np.savetxt("tmin.dat", tmin)
np.savetxt("teject.dat", teject)
np.savetxt("tcol.dat", tcol)
np.savetxt("ce.dat", np.array(ce))
