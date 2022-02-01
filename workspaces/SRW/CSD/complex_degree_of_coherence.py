"""The following lines of code are adaptation from pieces of coding by other authors. Before each segment, there is a
short description of the functionality, whom the code was copied from, who adapted it and when the inclusion to srwlibAux
was done."""
#############################################################################
# Script used for obtaining the mutual coherence plot from the SRW
# Adapted from O. Chubar IgorPro script - private correspondence
# Authors: Luca Rebuffi, Manuel Sanchez del Rio (Python version)
# Adapted by: Rafael Celestre
# 09.01.2018
#############################################################################

def _file_load(_fname, _read_labels=1):
    nLinesHead = 11
    hlp = []

    with open(_fname, 'r') as f:
        for i in range(nLinesHead):
            hlp.append(f.readline())

    ne, nx, ny = [int(hlp[i].replace('#', '').split()[0]) for i in [3, 6, 9]]
    ns = 1
    testStr = hlp[nLinesHead - 1]
    if testStr[0] == '#':
        ns = int(testStr.replace('#', '').split()[0])

    e0, e1, x0, x1, y0, y1 = [float(hlp[i].replace('#', '').split()[0]) for i in [1, 2, 4, 5, 7, 8]]

    data = numpy.squeeze(numpy.loadtxt(_fname, dtype=numpy.float64))  # get data from file (C-aligned flat)

    allrange = e0, e1, ne, x0, x1, nx, y0, y1, ny

    arLabels = ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity']
    arUnits = ['eV', 'm', 'm', 'ph/s/.1%bw/mm^2']

    if _read_labels:

        arTokens = hlp[0].split(' [')
        arLabels[3] = arTokens[0].replace('#', '')
        arUnits[3] = '';
        if len(arTokens) > 1:
            arUnits[3] = arTokens[1].split('] ')[0]

        for i in range(3):
            arTokens = hlp[i * 3 + 1].split()
            nTokens = len(arTokens)
            nTokensLabel = nTokens - 3
            nTokensLabel_mi_1 = nTokensLabel - 1
            strLabel = ''
            for j in range(nTokensLabel):
                strLabel += arTokens[j + 2]
                if j < nTokensLabel_mi_1: strLabel += ' '
            arLabels[i] = strLabel
            arUnits[i] = arTokens[nTokens - 1].replace('[', '').replace(']', '')

    return data, None, allrange, arLabels, arUnits

def _loadNumpyFormatCoh(_filename):
    data, dump, allrange, arLabels, arUnits = _file_load(_filename)

    dim_x = allrange[5]
    dim_y = allrange[8]

    dim = 1
    if dim_x > 1:
        dim = dim_x
    elif dim_y > 1:
        dim = dim_y

    np_array = data.reshape((dim, dim))
    np_array = np_array.transpose()

    if dim_x > 1:
        coordinates = numpy.linspace(allrange[3], allrange[4], dim_x)
        conj_coordinates = numpy.linspace(allrange[3], allrange[4], dim_x)
    elif dim_y > 1:
        coordinates = numpy.linspace(allrange[6], allrange[7], dim_y)
        conj_coordinates = numpy.linspace(allrange[6], allrange[7], dim_y)
    else:
        coordinates = None
        conj_coordinates = None

    return coordinates, conj_coordinates, np_array, allrange

def DegreeOfTransverseCoherence(_filename,_set_extrapolated_to_zero=True):
    """
    Converts the output files originated by srwl_wfr_emit_prop_multi_e from Mutual Intensity (Cross Spectral Density -
    CSD) to a Degree of Transverse Coherence (DoTC) file saving it in hdf5 generic file.
    :param _filename: Multual intensity array data from an ASCII file (format is defined in srwl_uti_save_intens_ascii
    :param _set_extrapolated_to_zero: zero padding of the matrix when applying rotation to it when set to TRUE
    """
    from scipy.interpolate import RectBivariateSpline
    FileName = _filename.split("/")
    print(">>>> Calculating the degree of transverse coherence: %s"%FileName[-1])

    coor, coor_conj, mutual_intensity, wfr_mesh = _loadNumpyFormatCoh(_filename)

    file_h5 = _filename.replace(".dat",".h5")

    f = h5py.File(file_h5, 'w')
    f1 = f.create_group("Cross_Spectral_Density")
    f1["CSD_method"] = "Oleg's IgorPro"
    f1["CSD_photon_energy"] = wfr_mesh[0]
    f1["CSD_photon_mesh"] = numpy.array([wfr_mesh[3::]])
    f1["CSD"] = mutual_intensity
    if wfr_mesh[5] == 1:
        f1["CSD_direction"] = "vertical"
    else:
        f1["CSD_direction"] = "horizontal"

    interpolator0 = RectBivariateSpline(coor, coor_conj, mutual_intensity, bbox=[None, None, None, None], kx=3, ky=3,s=0)

    X = numpy.outer(coor, numpy.ones_like(coor_conj))
    Y = numpy.outer(numpy.ones_like(coor), coor_conj)

    nmResDegCoh_z = numpy.abs(interpolator0(X + Y, X - Y, grid=False)) / \
                    numpy.sqrt(numpy.abs(interpolator0(X + Y, X + Y, grid=False))) / \
                    numpy.sqrt(numpy.abs(interpolator0(X - Y, X - Y, grid=False)))

    if _set_extrapolated_to_zero:
        nx, ny = nmResDegCoh_z.shape

        idx = numpy.outer(numpy.arange(nx), numpy.ones((ny)))
        idy = numpy.outer(numpy.ones((nx)), numpy.arange(ny))

        mask = numpy.ones_like(idx)

        bad = numpy.where(idy < 1. * (idx - nx / 2) * ny / nx)
        mask[bad] = 0

        bad = numpy.where(idy > ny - 1. * (idx - nx / 2) * ny / nx)
        mask[bad] = 0

        bad = numpy.where(idy < 0.5 * ny - 1. * idx * ny / nx)
        mask[bad] = 0

        bad = numpy.where(idy > 0.5 * ny + 1. * idx * ny / nx)
        mask[bad] = 0

        nmResDegCoh_z *= mask

    if _filename is not None:

        f1 = f.create_group("Degree_of_Transverse_Coherence")
        f1["DoTC_method"] = "Oleg's IgorPro"
        f1["DoTC_photon_energy"] = wfr_mesh[0]
        f1["DoTC_photon_mesh"] = numpy.array([wfr_mesh[3::]])
        f1["DoTC"] = nmResDegCoh_z
        if wfr_mesh[5] == 1:
            f1["DoTC_direction"] = "vertical"
        else:
            f1["DoTC_direction"] = "horizontal"
        f.close()

        FileName = file_h5.split("/")
        print(">>>> File %s written to disk."%FileName[-1])
