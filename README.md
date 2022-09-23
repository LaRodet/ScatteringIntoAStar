# ScatteringIntoAStar
Material for Rodet & Lai 2022

Please contact me for a link to get all the simulations folders (lbr63@cornell.edu). The folder e=00 is the case where ep = 0, Rp = 1 RJ and mp = 1e-3 Mstar.

Each folder corresponds to the results of a simulation of >1000 test particles, the details of which are given in Section 3. The code to produce one of these folders (e=03_r=10) is given as an example in CloseEncounters.py. Particles are identified by their number, which also corresponds to their positions in the different result files. The unit of distance and time are the planet semi-major axis and orbital period (except in the BetaPictoris case, where it is in au and yr).

ci.dat has the initial orbital elements of every test particle:
  - the first line is the semi-major axis;
  - the second line is the eccentricity;
  - the thrid line is the inclination in radian;
  - the fourth line is the argument of periastron in radian;
  - the fifth line is the longitude of the node in radian;
  - the sixth line is the true anomaly in radian.

dmin.dat records the minimum stellar approach of every test particles in unit of the planet semi-major axis.

tmin.dat records the time when the test particles reached their minimum approach in unit of the planet orbital period.

teject.dat records the time at which the test particles was ejected from the system (distance > 100 x planet semi-major axis, eccentricity > 1). If a particle is not ejected, then the value is 0.

tcol.dat records the time at which the test particles collided with the planet (distance < planet radius). If a particle did not collide, then the value is 0.

ce.dat records all the close encounters that happened during the simulation. Contrary to the other files, each line is a different encounter, not a different test particles. Each close encounter is recorded as [particle #, closest approach, time of closest approach].  

The name of the folder gives the eccentricity of the planet (e=03 -> ep = 0.3) and if the mass and radius are different from the fiducial values (mp = 1e-3 Mstar and Rp = 1 RJ are the fiducial values, if different then m=01 -> mp = 1e-2 Mstar, r=10 -> Rp = 10 RJ).

Finally, qmin_vs_m.py is the code to create Fig. 2 and qmin_vs_ap_zone.py is the code to create Fig. B1. They require the numpy, matplotlib and scipy libraries, along with functions in extrafunctions.py
