import Shadow
import numpy
import xraylib
from orangecontrib.shadow.util.shadow_util import ShadowPhysics
from scipy.interpolate import interp1d
from scipy import optimize
import csv
import matplotlib.pyplot as plt
import pandas as pd

def lens_rad(focal_dist, photon_energy=7000, material='Be', n_lens=1):
    """ Short function to get the lens radius for a given focal_distance [m]"""
    #delta = 1 - xraylib.Refractive_Index_Re(material, photon_energy/1000, ShadowPhysics.getMaterialDensity(material))
    #radius = 2 * n_lens * focal_dist * delta
    #to use exactly what Manolo used:
    xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7, 1.85)).real
    radius = focal_dist * n_lens * (2 * xrl_delta)    

    return numpy.round(radius, 7)

def get_min(x,y):
    """Interpolate x and y points and then get the minimum and its position 
    if an error rises just get the min value of array"""
    
    f = interp1d(x, y, kind='cubic')
    #get minimum value, initial guess at minimum FWHM value #
    try:
    
        min_x = numpy.round((optimize.minimize(f, x0=x[y.index(min(y))])).get('x')[0], 3) #x0=min(y)
        value_y = numpy.round((optimize.minimize(f, x0=x[y.index(min(y))])).get('fun'), 3)
        
    except ValueError:
        print("x0 is out of limits so taking just the minimum")
        value_y = min(y)
        min_x = x[y.index(value_y)]

    return min_x, value_y

def run_shadow(beam_file, f1_h, f1_v, f2_h, f2_v, mode=170):
    """ Ray tracing that starts from a beam file, please notice that
    in this case the file is after T1-entrance (Transfocator 1 entrance)
    the input are the fosucing lengths and mode (position of f2)  """
    # Read the beamfile previously created with OASYS#
    beam = Shadow.Beam()
    beam.load(beam_file)

    #print(beam.rays)
    iwrite = 0
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()
    oe5 = Shadow.OE()
    oe6 = Shadow.OE()
    oe7 = Shadow.OE()
    oe8 = Shadow.OE()
    oe9 = Shadow.OE()
    oe10 = Shadow.OE()

    oe1.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2 * lens_rad(f1_h), 0.0])
    oe1.CIL_ANG = 90.0
    oe1.DUMMY = 100.0
    oe1.FCYL = 1
    oe1.FHIT_C = 1
    oe1.FILE_R_IND_IMA = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe1.FMIRR = 4
    oe1.FSHAPE = 2
    oe1.FWRITE = 3
    oe1.F_EXT = 1
    oe1.F_REFRAC = 1
    oe1.F_R_IND = 2
    oe1.PARAM = lens_rad(f1_h)
    oe1.RLEN2 = 0.0005
    oe1.RWIDX2 = 0.0005
    oe1.T_IMAGE = 5e-06
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 0.0

    oe2.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 2 * lens_rad(f1_h), 0.0])
    oe2.CIL_ANG = 90.0
    oe2.DUMMY = 100.0
    oe2.FCYL = 1
    oe2.FHIT_C = 1
    oe2.FILE_R_IND_OBJ = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe2.FMIRR = 4
    oe2.FSHAPE = 2
    oe2.FWRITE = 3
    oe2.F_CONVEX = 1
    oe2.F_EXT = 1
    oe2.F_REFRAC = 1
    oe2.F_R_IND = 1
    oe2.PARAM = lens_rad(f1_h)
    oe2.RLEN2 = 0.0005
    oe2.RWIDX2 = 0.0005
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 0.0
    oe2.T_REFLECTION = 180.0
    oe2.T_SOURCE = 5e-06

    oe3.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2 * lens_rad(f1_v), 0.0])
    oe3.DUMMY = 100.0
    oe3.FCYL = 1
    oe3.FHIT_C = 1
    oe3.FILE_R_IND_IMA = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe3.FMIRR = 4
    oe3.FSHAPE = 2
    oe3.FWRITE = 3
    oe3.F_EXT = 1
    oe3.F_REFRAC = 1
    oe3.F_R_IND = 2
    oe3.PARAM = lens_rad(f1_v)
    oe3.RLEN2 = 0.0005
    oe3.RWIDX2 = 0.0005
    oe3.T_IMAGE = 5e-06
    oe3.T_INCIDENCE = 0.0
    oe3.T_REFLECTION = 180.0
    oe3.T_SOURCE = 0.0

    oe4.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 2 * lens_rad(f1_v), 0.0])
    oe4.DUMMY = 100.0
    oe4.FCYL = 1
    oe4.FHIT_C = 1
    oe4.FILE_R_IND_OBJ = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe4.FMIRR = 4
    oe4.FSHAPE = 2
    oe4.FWRITE = 3
    oe4.F_CONVEX = 1
    oe4.F_EXT = 1
    oe4.F_REFRAC = 1
    oe4.F_R_IND = 1
    oe4.PARAM = lens_rad(f1_v)
    oe4.RLEN2 = 0.0005
    oe4.RWIDX2 = 0.0005
    oe4.T_IMAGE = 0.0
    oe4.T_INCIDENCE = 0.0
    oe4.T_REFLECTION = 180.0
    oe4.T_SOURCE = 5e-06
    # This is T2 aperture (Tranfocator 2 entrance) #
    oe5.DUMMY = 100.0
    oe5.FWRITE = 3
    oe5.F_REFRAC = 2
    oe5.F_SCREEN = 1
    oe5.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe5.K_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe5.N_SCREEN = 1
    oe5.RX_SLIT = numpy.array([0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe5.RZ_SLIT = numpy.array([0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe5.T_IMAGE = 0.0
    oe5.T_INCIDENCE = 0.0
    oe5.T_REFLECTION = 180.0
    # To consider the two position cases #
    if mode == 170:
        oe5.T_SOURCE = 105.0
    elif mode == 192:
        oe5.T_SOURCE = 127.0
    else:
        raise RuntimeError(f'ERROR: mode not identified {mode}')
    # f2 horizontal lens #
    oe6.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2 * lens_rad(f2_h), 0.0])
    oe6.CIL_ANG = 90.0
    oe6.DUMMY = 100.0
    oe6.FCYL = 1
    oe6.FHIT_C = 1
    oe6.FILE_R_IND_IMA = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe6.FMIRR = 4
    oe6.FSHAPE = 2
    oe6.FWRITE = 3
    oe6.F_EXT = 1
    oe6.F_REFRAC = 1
    oe6.F_R_IND = 2
    oe6.PARAM = lens_rad(f2_h)
    oe6.RLEN2 = 0.0005
    oe6.RWIDX2 = 0.0005
    oe6.T_IMAGE = 5e-06
    oe6.T_INCIDENCE = 0.0
    oe6.T_REFLECTION = 180.0
    oe6.T_SOURCE = 0.0

    oe7.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 2 * lens_rad(f2_h), 0.0])
    oe7.CIL_ANG = 90.0
    oe7.DUMMY = 100.0
    oe7.FCYL = 1
    oe7.FHIT_C = 1
    oe7.FILE_R_IND_OBJ = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe7.FMIRR = 4
    oe7.FSHAPE = 2
    oe7.FWRITE = 3
    oe7.F_CONVEX = 1
    oe7.F_EXT = 1
    oe7.F_REFRAC = 1
    oe7.F_R_IND = 1
    oe7.PARAM = lens_rad(f2_h)
    oe7.RLEN2 = 0.0005
    oe7.RWIDX2 = 0.0005
    oe7.T_IMAGE = 0.0
    oe7.T_INCIDENCE = 0.0
    oe7.T_REFLECTION = 180.0
    oe7.T_SOURCE = 5e-06

    oe8.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2 * lens_rad(f2_v), 0.0])
    oe8.DUMMY = 100.0
    oe8.FCYL = 1
    oe8.FHIT_C = 1
    oe8.FILE_R_IND_IMA = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe8.FMIRR = 4
    oe8.FSHAPE = 2
    oe8.FWRITE = 3
    oe8.F_EXT = 1
    oe8.F_REFRAC = 1
    oe8.F_R_IND = 2
    oe8.PARAM = lens_rad(f2_v)
    oe8.RLEN2 = 0.0005
    oe8.RWIDX2 = 0.0005
    oe8.T_IMAGE = 5e-06
    oe8.T_INCIDENCE = 0.0
    oe8.T_REFLECTION = 180.0
    oe8.T_SOURCE = 0.0

    oe9.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 2 * lens_rad(f2_v), 0.0])
    oe9.DUMMY = 100.0
    oe9.FCYL = 1
    oe9.FHIT_C = 1
    oe9.FILE_R_IND_OBJ = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe9.FMIRR = 4
    oe9.FSHAPE = 2
    oe9.FWRITE = 3
    oe9.F_CONVEX = 1
    oe9.F_EXT = 1
    oe9.F_REFRAC = 1
    oe9.F_R_IND = 1
    oe9.PARAM = lens_rad(f2_v)
    oe9.RLEN2 = 0.0005
    oe9.RWIDX2 = 0.0005
    oe9.T_IMAGE = 0.0
    oe9.T_INCIDENCE = 0.0
    oe9.T_REFLECTION = 180.0
    oe9.T_SOURCE = 5e-06

    oe10.DUMMY = 100.0
    oe10.FWRITE = 3
    oe10.F_REFRAC = 2
    oe10.F_SCREEN = 1
    oe10.N_SCREEN = 1
    oe10.T_IMAGE = 0.0
    oe10.T_INCIDENCE = 0.0
    oe10.T_REFLECTION = 180.0
    
    # To consider the different modes #

    if mode == 170:
        oe10.T_SOURCE = 30.0
    elif mode == 192:
        oe10.T_SOURCE = 8.0
    else:
        raise RuntimeError(f'ERROR: mode not identified {mode}')
    #
    #run optical element 1
    #
    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1,1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")
    #
    #run optical element 2
    #
    print("    Running optical element: %d"%(2))
    if iwrite:
        oe2.write("start.02")

    beam.traceOE(oe2,2)

    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")
    #
    #run optical element 3
    #
    print("    Running optical element: %d"%(3))
    if iwrite:
        oe3.write("start.03")

    beam.traceOE(oe3,3)

    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")
    #
    #run optical element 4
    #
    print("    Running optical element: %d"%(4))
    if iwrite:
        oe4.write("start.04")

    beam.traceOE(oe4,4)

    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")
    #
    #run optical element 5
    #
    print("    Running optical element: %d"%(5))
    if iwrite:
        oe5.write("start.05")

    beam.traceOE(oe5,5)

    if iwrite:
        oe5.write("end.05")
        beam.write("star.05")
    #
    #run optical element 6
    #
    print("    Running optical element: %d"%(6))
    if iwrite:
        oe6.write("start.06")

    beam.traceOE(oe6,6)

    if iwrite:
        oe6.write("end.06")
        beam.write("star.06")
    #
    #run optical element 7
    #
    print("    Running optical element: %d"%(7))
    if iwrite:
        oe7.write("start.07")

    beam.traceOE(oe7,7)

    if iwrite:
        oe7.write("end.07")
        beam.write("star.07")
    #
    #run optical element 8
    #
    print("    Running optical element: %d"%(8))
    if iwrite:
        oe8.write("start.08")

    beam.traceOE(oe8,8)

    if iwrite:
        oe8.write("end.08")
        beam.write("star.08")
    #
    #run optical element 9
    #
    print("    Running optical element: %d"%(9))
    if iwrite:
        oe9.write("start.09")

    beam.traceOE(oe9,9)

    if iwrite:
        oe9.write("end.09")
        beam.write("star.09")
    #
    #run optical element 10
    #
    print("    Running optical element: %d"%(10))
    if iwrite:
        oe10.write("start.10")

    beam.traceOE(oe10,10)

    if iwrite:
        oe10.write("end.10")
        beam.write("star.10")
    # FWHM in horizontal an vertical in microns #
    fwhm_h = 1e6 * beam.histo1(1, nbins=301, nolost=1, ref=23)["fwhm"]
    fwhm_v = 1e6 * beam.histo1(3, nbins=301, nolost=1, ref=23)["fwhm"]

    std_h = 1e6 * beam.get_standard_deviation(1, nolost=1, ref=23)
    std_v = 1e6 * beam.get_standard_deviation(3, nolost=1, ref=23)

    #Shadow.ShadowTools.plotxy(beam,1,3,nbins=301,nolost=1,title="Real space")

    return fwhm_h, fwhm_v, std_h, std_v

def scan_f1_f2(f1_i=5, f1_f=100, f1_steps = 60 , f2_i=10, f2_f=50, f2_steps=30, mode=170):

    """Short function to perform a ray tracing scanning, for a given
    f1 it scans f2, getting the minimum FWHM [um] at sample position and
    the corresponding f2 [m], it writes and output file with this info
    including the lenses curvatures"""

    f1_scan = numpy.linspace(f1_i, f1_f, num=f1_steps, endpoint=True).tolist()
    f2_scan = numpy.linspace(f2_i, f2_f, num=f2_steps, endpoint=True).tolist()

    size_h = []
    size_v = []
    f2_h = []
    f2_v = []

    #Scanning routine, scan f2 for each f1

    for f1 in f1_scan:

        tmp_h = []
        tmp_v = []       
        for f2 in f2_scan:
            print(f'Running for f1={f1} m, step {f1_scan.index(f1)+1} of {len(f1_scan)}')
            print(f'Running for f2={f2} m, step {f2_scan.index(f2)+1} of {len(f2_scan)}')        
            # run SHADOW
            fwhm_h, fwhm_v, std_h, std_v = run_shadow("T1_entrance_beam_1M.dat", f1, f1, f2, f2, mode=mode)
            #tmp_h.append(std_h); tmp_v.append(std_v) #using std
            tmp_h.append(fwhm_h); tmp_v.append(fwhm_v)  #using FWHM
        # get minimum and f2 position #
        min_x_h, value_y_h = get_min(f2_scan, tmp_h)
        min_x_v, value_y_v = get_min(f2_scan, tmp_v)
        
        f2_h.append(min_x_h)
        size_h.append(value_y_h)
        f2_v.append(min_x_v)
        size_v.append(value_y_v)

        print(f'For f1 = {f1}, h_min is {min_x_h}, at f2_h of {value_y_h}')
        print(f'For f1 = {f1}, v_min is {min_x_v}, at f2_v of {value_y_v}')

    return f1_scan, f2_h, size_h, f2_v, size_v

def get_fwhm(input_beam_file, f1_scan, f2_h, f2_v, outfile_name, mode=170):
    
    """This function gets the FWHM for a given list of f1, f2_h, f2_v,
    it returns the FWHM in both planes"""

    s_h = []
    s_v = []

    for f1, f2h, f2v in zip(f1_scan, f2_h, f2_v):
        #print(f'Running for f1={f1} m, f2_horizontal{f2h} m, f2_vertical{f2v} m, step {f1_scan.index(f1)+1} of {len(f1_scan)}')
        fwhm_h, fwhm_v, std_h, std_v = run_shadow(input_beam_file, f1, f1, f2h, f2v, mode=mode)
        s_h.append(fwhm_h)
        s_v.append(fwhm_v)

        # write output file #

    with open(outfile_name, mode="w") as csvfile:
        cwriter = csv.writer(csvfile, delimiter=",", quotechar=" ", quoting=csv.QUOTE_MINIMAL)
        cwriter.writerow(['f1[m]', 'r1[um]','f2_h[m]','r2_h[um]', 'H_FWHM[um]', 'f2_v[m]', 'r2_v[um]', 'V_FWHM[um]'])
        for i, f1 in enumerate(f1_scan):
            cwriter.writerow([f1,lens_rad(f1)*1e6, f2_h[i],lens_rad(f2_h[i])*1e6, s_h[i], f2_v[i],lens_rad(f2_v[i])*1e6, s_v[i]])

    print(f"File written to disk:{outfile_name}")

    return fwhm_h, fwhm_v

def read_f1_f2(input_file):

    df = pd.read_csv(input_file, delim_whitespace=True, header=None, skiprows=None)

    return df.iloc[:,0], df.iloc[:,1]


if __name__ == '__main__':

    #scan_f1_f2(outfile_name='f1_f2_scan_v2_1M.csv',f1_i=13, f1_f=100 f1_steps = 88, f2_i=10, f2_f=50, f2_steps=21)
    #f1_scan, f2_h, size_h, f2_v, size_v = scan_f1_f2(f1_i=13, f1_f=100, f1_steps=88, f2_i=10, f2_f=50, f2_steps=21,
    #                                                  mode=170)
    #f1_scan, f2_h, size_h, f2_v, size_v = scan_f1_f2(f1_i=20, f1_f=100, f1_steps=88, f2_i=5, f2_f=12, f2_steps=21,
    #                                                 mode=192)
    #f1_scan, f2_h, size_h, f2_v, size_v = scan_f1_f2(f1_i=20, f1_f=100, f1_steps=10, f2_i=5, f2_f=15, f2_steps=11,
    #                                                 mode=192)
    #print("f1 f2 scanning done proceeding to calculate the FWHM for each")
    
    # Manolos files f1 and f2 #
    #f1_h, f2_h = read_f1_f2('f1_vs_f2_slit40.3_h.dat')
    f1_h, f2_h = read_f1_f2('f1_vs_f2_slit85.1_h.dat')

    #f1_v, f2_v = read_f1_f2('f1_vs_f2_slit227_v.dat')
    f1_v, f2_v = read_f1_f2('f1_vs_f2_slit506.7_v.dat')

    #fwhm_h, fwhm_v = get_fwhm("T1_entrance_beam_1M.dat", f1_scan, f2_h, f2_v, "f1_f2_scan_v3_5M_test.csv")
    #fwhm_h, fwhm_v = get_fwhm("T1_entrance_beam_5M.dat", f1_scan, f2_h, f2_v, "f1_f2_scan_v3_5M_fwhm_192_std.csv",
    #                          mode=170)
    ## First scan for paper figures: this uses Manolos f1 and f2 values ##
    #fwhm_h, fwhm_v = get_fwhm("T1_entrance_beam_40.3h_227v_5M.dat", f1_h, f2_h, f2_v, "f1_f2_scan_paper_5M_fwhm_1h_1v.csv",
    #                        mode=170)
    ## Second scan for 3h_3v ##
    fwhm_h, fwhm_v = get_fwhm("T1_entrance_beam_81.5h_506.7v_5M.dat", f1_h, f2_h, f2_v, "f1_f2_scan_paper_5M_fwhm_3h_3v.csv",
                              mode=170)
