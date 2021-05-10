import xraylib
energy_keV = 7
R = 200e-6
#
xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", energy_keV, 1.85)).real
F = R / (2 * xrl_delta)
p=65 #-36.0
q = 1/(1/F-1/p)
print("F (source): %g, p1: %g, q1: %g" % (F,p,q))
p=65 -36.0
q = 1/(1/F-1/p)
print("F (slit): %g, p1: %g, q1: %g" % (F,p,q))
print("R_Be [mm]= ", 1e3*R)