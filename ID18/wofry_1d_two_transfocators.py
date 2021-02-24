import numpy


if __name__ == "__main__":
    from srxraylib.plot.gol import plot

    directions = ['h','v']
    focii = ['large', 'small']
    energies = [7, 15, 20, 35]
    distances = [170, 192]

    key1 = []
    for energy in energies:
        for distance in distances:
            for focus in focii:
                key1.append("e%02dkeV_f2_at_%dm_%s" % (energy, distance, focus))
                # for direction in directions:
                #     print("e%02dkeV_f2_at_%dm_%s  %s" % (energy, distance, focus, direction))



    # attached the npz file, best is to read it with "datastorage": data = datastorage.read(fname)
    #
    #
    # and here the δ,β
    #
    # In [1]: from sr import materials
    #
    # In [2]: for energy in (7,10,15,20,35):
    #    ...: print(materials.get_delta_beta("Be",energy=energy))
    # ...:
    # (6.959562521724472e-06, 3.920823829597519e-09)
    # (3.4079264780162433e-06, 1.1184654031792269e-09)
    # (1.514042063277543e-06, 3.6034468029877136e-10)
    # (8.515233095307551e-07, 2.0080453378588767e-10)
    # (2.780114856104632e-07, 8.816014655748295e-11) """

    import datastorage

    data = datastorage.read("summary_of_GSM_results_for_Manuel.npz")

    # print(data["e07keV_f2_at_170m_large"]['v'])

    label = "e07keV_f2_at_170m_large"
    for direction in directions:
        tmp = data[label][direction]
        for key in tmp.keys():
            # print(key, tmp[key])
            slit = tmp["p035m_hard_aperture"]
            f1 = tmp["p066m_f1_used"]
            f2 = tmp["pf2pos_f2_used"]
            size200 = tmp["p200m_fwhm_beam"]
            sizewaist = tmp["waist_fwhm_size"]
            poswaist = tmp["waist_pos"]
        print(">>>>%s %s :  slit: %g m, f1: %g, f2: %g, size: %g um (%g at %g m)" %
              (label, direction, slit, f1, f2, 1e6*size200, 1e6*sizewaist, poswaist))

    # print(key1[0], data[key1[0]]['h'])


