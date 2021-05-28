
"""

first run:
case7keV_2h_spectral_density.dat:#UFWHM  16.6572
case7keV_2v_spectral_density.dat:#UFWHM  24.0535
case7keV_4h_spectral_density.dat:#UFWHM  28.7397
case7keV_4v_spectral_density.dat:#UFWHM  31.0935
case7keV_6h_spectral_density.dat:#UFWHM  8.01583
case7keV_6v_spectral_density.dat:#UFWHM  7.92004
case7keV_8h_spectral_density.dat:#UFWHM  7.7421
case7keV_8v_spectral_density.dat:#UFWHM  10.7947
case7keV_10h_spectral_density.dat:#UFWHM 1.95508
case7keV_10v_spectral_density.dat:#UFWHM 3.05068


second run:
case7keV_2h_spectral_density.dat:#UFWHM 16.3054
case7keV_2v_spectral_density.dat:#UFWHM 23.4668
case7keV_4h_spectral_density.dat:#UFWHM 28.5441
case7keV_4v_spectral_density.dat:#UFWHM 31.0935
case7keV_6h_spectral_density.dat:#UFWHM 7.69525
case7keV_6v_spectral_density.dat:#UFWHM 8.21338
case7keV_8h_spectral_density.dat:#UFWHM 7.7421
case7keV_8v_spectral_density.dat:#UFWHM 10.7947
case7keV_10h_spectral_density.dat:#UFWHM 1.84558
case7keV_10v_spectral_density.dat:#UFWHM 2.69868
"""

marco_wofry = [
 [  20.0    ,     16.3054 ],
 [  20.0    ,     23.4668 ],
 [  30.0    ,     28.5441 ],
 [  30.1    ,     31.0935 ],
 [  11.5    ,     7.69525 ],
 [  7.5     ,     8.21338 ],
 [  7.0     ,     7.7421],
 [  12.0    ,     10.7947 ],
 [  3.2     ,      1.84558 ],
 [  3.0     ,      2.69868 ],
 ]

for element in marco_wofry:
    print("%5.1f   %5.1f   %5.1f " % (element[0], element[1], 100*(element[0]-element[1])/element[1]))


