import xraylib

F = 10.0
xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7, 1.85)).real
R = F * (2 * xrl_delta)
print("F: %g  R_Be [m]= %g" % (F, R))


R = 200e-6
xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7, 1.85)).real
F = R / (2 * xrl_delta)
print("F: %g  R_Be [m]= %g" % (F, R))
