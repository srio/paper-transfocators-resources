import numpy
import xraylib

# https://www.rxoptics.de/design-parameters/
# view-source:https://www.rxoptics.de/design-parameters/
# https://www.rxoptics.de/wp-content/themes/rxoptics/js/crlcalc.js

def erf(x): #  {//A&S formula 7.1.26
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911
    t = 1.0 / (1.0 + p * x)
    return 1.0 - (((((a5*t+a4)*t)+a3)*t+a2)*t+a1) * t * numpy.exp(-x*x)

def rxoptics_compute(e=10.0,        # * photon energy in keV
            presf=0,       # * 0=focal length, 1=image distance,
            desf=1000.0,   # * desired focal lengt - focal length (presf=0) in mm or image distance (presf=1) in mm
            dist=40.0,     # * distance source-lens [m]
            R=0.05,        # * radius of curvature R at the apex of the individual lenses in the CRL [mm]
            w=150.0,       # * width (FWHM) of the x-ray source [um]
            h=150.0,       # * height (FWHM) of the x-ray source [um]
            F=2.0,         # * thickness of the lens frame [mm]
            d=30,          # * web thickness [um]
            D2=0,          # * lens type: 0=1D, 1=2D
            th=1.0,        # * thickness D of the individual lenses [mm]
            MATERIAL="Be", # *
            force_n_lenses = None,
            use_xraylib=0, # internal(0) or xraylib(1) optics library
            plot_cross_sections=False):



    lambda1 = 12.398 / e
    if MATERIAL == "Be":
        rho=1.85
        tdelta=5.4e-6 * lambda1**2 * rho * 4/9.012
        muoverrhos= [[1.0,604.1],[1.5,179.7],[2.0,74.69],[3.0,21.27],[4.0,8.685],[5.0,4.369],[6.0,2.527],[8.0,1.124],[10.0,0.6466],[15.0,0.307],[20.0,0.2251],[30.0,0.1792],[40.0,0.164],[50.0,0.1554],[60.0,0.1493],[80.0,0.1401],[100.0,0.1328],[150.0,0.119],[200.0,0.1089],[300.0,0.09463],[400.0,0.08471],[500.0,0.07739],[600.0,0.07155],[800.0,0.06286],[1000.0,0.05652],[1250.0,0.05054],[1500.0,0.04597],[2000.0,0.03938],[3000.0,0.03138],[4000.0,0.02664],[5000.0,0.02347],[6000.0,0.02121],[8000.0,0.01819],[10000.0,0.01627],[15000.0,0.01361],[20000.0,0.01227]]
    elif MATERIAL == "Al":
        rho=2.7
        tdelta=5.4*1e-6 * lambda1**2 * rho * 13/26.982
        muoverrhos=[[1.0,1185.0],[1.5,402.2],[1.5596,362.1],[1.5596,3957.0],[2.0,2263.0],[3.0,788.0],[4.0,360.5],[5.0,193.4],[6.0,115.3],[8.0,50.33],[10.0,26.23],[15.0,7.955],[20.0,3.441],[30.0,1.128],[40.0,0.5685],[50.0,0.3681],[60.0,0.2778],[80.0,0.2018],[100.0,0.1704],[150.0,0.1378],[200.0,0.1223],[300.0,0.1042],[400.0,0.09276],[500.0,0.08445],[600.0,0.07802],[800.0,0.06841],[1000.0,0.06146],[1250.0,0.05496],[1500.0,0.05006],[2000.0,0.04324],[3000.0,0.03541],[4000.0,0.03106],[5000.0,0.02836],[6000.0,0.02655],[8000.0,0.02437],[10000.0,0.02318],[15000.0,0.02195],[20000.0,0.02168]]
    else:
        rho=8.90
        tdelta=5.4*1e-6 * lambda1**2 * rho * 28/58.693
        muoverrhos=[[1.0,9855.0],[1.00404,9753.0],[1.0081,9654.0],[1.0081,10990.0],[1.5,4234.0],[2.0,2049.0],[3.0,709.4],[4.0,328.2],[5.0,179.3],[6.0,109.0],[8.0,49.52],[8.3328,44.28],[8.3328,329.4],[10.0,209.0],[15.0,70.81],[20.0,32.2],[30.0,10.34],[40.0,4.6],[50.0,2.474],[60.0,1.512],[80.0,0.7306],[100.0,0.444],[150.0,0.2208],[200.0,0.1582],[300.0,0.1154],[400.0,0.09765],[500.0,0.08698],[600.0,0.07944],[800.0,0.06891],[1000.0,0.0616],[1250.0,0.05494],[1500.0,0.05015],[2000.0,0.04387],[3000.0,0.03745],[4000.0,0.03444],[5000.0,0.03289],[6000.0,0.0321],[8000.0,0.03164],[10000.0,0.03185],[15000.0,0.0332],[20000.0,0.03476]]

    muoverrhos = numpy.array(muoverrhos)
    l= muoverrhos.shape[0]


    if plot_cross_sections:
        XRL_MU = numpy.zeros(l)
        XRL_MU_E = numpy.zeros(l)

        for i in range(l):
            XRL_MU[i] = rho * xraylib.CS_Total(xraylib.SymbolToAtomicNumber(MATERIAL), muoverrhos[i,0])
            XRL_MU_E[i] = rho * xraylib.CS_Energy(xraylib.SymbolToAtomicNumber(MATERIAL), muoverrhos[i,0])
            # print("     >>",muoverrhos[i,0],rho*muoverrhos[i,1], XRL_MU[i], XRL_MU_E[i], XRL_MU_E[i] / XRL_MU[i])

        from srxraylib.plot.gol import plot
        plot(muoverrhos[:,0], rho*muoverrhos[:,1],
             muoverrhos[:,0], XRL_MU,
             muoverrhos[:, 0], XRL_MU_E,
             xlog=1,ylog=1,title=MATERIAL,xrange=[1,1e3],xtitle="Photon energy [keV]",ytitle="mu [cm-1]",
             legend=["mu xro","mu xraylib","energy loss mu"])

    for i in range(l):
        if muoverrhos[i, 0] <= e:
            if i == 0:
                mu = rho * muoverrhos[0, 1]
            elif i == (l - 1):
                mu = rho * muoverrhos[l - 1, 1]
            else:
                mu = rho * (
                        muoverrhos[i, 1] * (e - muoverrhos[i - 1, 0]) / (muoverrhos[i, 0] - muoverrhos[i - 1, 0]) - \
                        muoverrhos[i - 1, 1] * (e - muoverrhos[i, 0]) / (muoverrhos[i, 0] - muoverrhos[i - 1, 0])
                )




    xrl_mu = rho * xraylib.CS_Total(xraylib.SymbolToAtomicNumber(MATERIAL), e)
    xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", e, rho)).real

    if use_xraylib:
        mu = xrl_mu
        tdelta = 2 * xrl_delta

    # //compute desired focal length and real number of lenses
    dist=dist*1e+3 # //dist in mm

    if presf: # //compute focal length
        f=desf
    else:
        f=dist*desf/(dist+desf)

    if force_n_lenses:
        N = int(force_n_lenses)
        exf = numpy.sqrt(R * F / tdelta) * 1 / (numpy.sin(N * numpy.sqrt(tdelta * F / R)))

    else:
        N = numpy.sqrt(R/(tdelta*F)) * numpy.arcsin( numpy.sqrt(R*F/tdelta) / f)

        if numpy.isnan(N):
            raise Exception("The focus lies within the CRL.")

        # //compute integer number of lenses and true focal length
        N=numpy.floor(N)
        print("N: %f" % (N) )



        exfa = numpy.sqrt(R*F/tdelta) * 1/(numpy.sin(N    *numpy.sqrt(tdelta*F/R))) # //in mm
        exfb = numpy.sqrt(R*F/tdelta) * 1/(numpy.sin((N+1)*numpy.sqrt(tdelta*F/R))) # //in mm

        if (numpy.abs(exfa-f) > numpy.abs(exfb-f) and (N+1)<=numpy.sqrt(R/(tdelta*F))*numpy.pi/2):
            exf = exfb
            N = N + 1
        else:
            exf = exfa

    print("N: %f" % (N))

    if exf >= dist:
        raise Exception("The source lies within the focal distance.")



    lcrl= N * F  # length of crl
    ol = N * numpy.sqrt(F * tdelta / R) # \omega*L
    feff = exf * numpy.cos(ol) # effective focal length in mm
    pp = exf - feff - lcrl / 2 # position of principal plane in mm

    idist = exf * dist / (dist - exf) #image distance in mm
    Neff = (N / 2 + 0.25 * numpy.sin(2*ol) / numpy.sin(ol / N)) # effective number of lenses
    gapp = 2 * numpy.sqrt(R * 1e+3 * (th * 1e+3 - d)) #geometric aperture in um


    # compute effective aperture

    ap = 0.1 * mu * Neff * numpy.power(gapp / 2e+3, 2) /(2 * R)
    twod = D2 # document.getElementById("D2").checked;
    if twod:
        Deff = gapp * numpy.sqrt((1 - numpy.exp(-ap)) / ap) #effective aperture in um for 2D lenses
    else:
        Deff = gapp * 0.5 * numpy.sqrt(numpy.pi/ap) * erf(numpy.sqrt(ap)) # effective aperture in um for 1D lenses



    corr = 2 * numpy.sqrt( 2 * numpy.log(2))
    res = corr * lambda1 * 1e+2 * exf / (numpy.pi * Deff) # resolution limit in nm
    fspotw = corr * idist * 1e+6 * numpy.sqrt(numpy.power(w/(corr*dist*1e+3),2) +
                    numpy.power(lambda1*1e-4/(numpy.pi * Deff),2)) # focal spot size in nm
    if twod:
        fspoth = corr * idist * 1e+6 * numpy.sqrt(numpy.power(h/(corr*dist*1e+3),2) +
                    numpy.power(lambda1*1e-4/(numpy.pi * Deff),2)) # focal spot size (height) in nm
    else:
        fspoth = 0


    # compute transmission and gain
    aD = 0.1 * mu * Neff * numpy.power(Deff / 2e+3,2) / (2 * R)
    if twod:
        transm = numpy.exp(-mu * N * d * 1e-4) * (1 - numpy.exp(-2 * aD)) / (2*aD) # transmission for 2D lenses (ill. of eff. ap.)
        transmgapp = numpy.exp(-mu * N * d * 1e-4) * (1 - numpy.exp(-2 * ap)) / (2*ap) # transmission for 2D lenses (ill. of geom. ap.)
    else:
        transm = numpy.exp(-mu * N * d * 1e-4) * numpy.sqrt(numpy.pi) * 0.5 * erf(numpy.sqrt(2*aD)) /\
                numpy.sqrt(2 * aD)# transmission for 1D lenses (ill. of eff. ap.)
        transmgapp = numpy.exp(-mu * N * d * 1e-4) * numpy.sqrt(numpy.pi) * 0.5 * erf(numpy.sqrt(2*ap)) / \
                numpy.sqrt(2 * ap) # transmission for 1D lenses (ill. of geom. ap.)


    if twod:
        gain = transmgapp * numpy.power(1e+3 * gapp,2) / (fspotw * fspoth) # gain for 2D lenses
    else:
        gain = transmgapp * 1e+3 * gapp / fspotw # gain for 1D lenses

    eidist = idist - pp - N * F / 2


    tkt = {}

    tkt["e"] = e
    tkt["presf"] = presf
    tkt["desf"] = desf
    tkt["dist"] = dist
    tkt["R"] = R
    tkt["w"] = w
    tkt["h"] = h
    tkt["F"] = F
    tkt["d"] = d
    tkt["D2"] = D2
    tkt["th"] = th
    tkt["MATERIAL"] = MATERIAL
    tkt["use_xraylib"] = use_xraylib
    tkt["tdelta"] = tdelta
    tkt["mu"] = mu
    tkt["N"] = N
    tkt["lcrl"] = lcrl
    tkt["gapp"] = gapp
    tkt["Deff"] =  Deff
    tkt["pp"] = pp
    tkt["feff"] = feff
    tkt["idist"] = idist
    tkt["eidist"] = eidist
    tkt["fspotw"] = fspotw
    tkt["fspoth"] = fspoth
    tkt["res"] =  res
    tkt["transmgapp"] = transmgapp
    tkt["gain"] =  gain


    return tkt

def rxoptics_display(tkt):
    presf       = tkt["presf"]
    desf        = tkt["desf"]
    dist        = tkt["dist"]
    R           = tkt["R"]
    w           = tkt["w"]
    h           = tkt["h"]
    F           = tkt["F"]
    d           = tkt["d"]
    D2          = tkt["D2"]
    th          = tkt["th"]
    MATERIAL    = tkt["MATERIAL"]
    use_xraylib = tkt["use_xraylib"]
    e = tkt["e"]

    tdelta = tkt["tdelta"]
    mu = tkt["mu"]
    N = tkt["N"]
    lcrl = tkt["lcrl"]
    gapp = tkt["gapp"]
    Deff = tkt["Deff"]
    pp = tkt["pp"]
    feff = tkt["feff"]
    idist = tkt["idist"]
    eidist = tkt["eidist"]
    fspotw = tkt["fspotw"]
    fspoth = tkt["fspoth"]
    res = tkt["res"]
    transmgapp = tkt["transmgapp"]
    gain = tkt["gain"]

    print("===== Inputs: =====")
    print("MATERIAL: %s (Lens material is: %s)" % (MATERIAL, MATERIAL))
    print("R: %g (Radius of curvature R at the apex of the individual lenses)" % (R))

    print("D2: %d (Lens type: %dD)" % (D2, D2+1))

    print("d: %g μm (Average web thickness d, that is, the distance of the apices of each lens' paraboloids. The most common values are 30 μm and 50 μm.)" % d)
    print("th: %g mm (Thickness D of the individual lenses. The most common value is 1 mm)" % th)
    print("F: %g mm (Thickness F of the individual lenses' frames. The most common value is 2 mm)" % F)

    print("presf: %d (Please choose whether you wish to specify the focal length (0) of the CRL or the distance of the CRL to the secondary source (1))" % presf)
    if presf == 0:
        print("desf: %g mm (Focal length value)" % desf)
    else:
        print("desf: %g mm (Distance of the CRL to the secondary source value)" % desf)

    print("e: %g keV (photon energy of the x-ray beam)" % (e))

    print("dist: %g m (Distance of the x-ray source to the CRL (more precisely, to its secondary principal plane))" % (dist))

    if D2:
        print("w: %g um (width (FWHM) of the x-ray source)" % (w))
        print("h: %g um (height (FWHM) of the x-ray source)" % (h))
    else:
        print("w: %g um (lateral extent (FWHM) of the x-ray source in the focused direction)" % (w))

    print("==========\n\n")


    print("\n====== Results ========")

    print("Refractive decrement: %g" % (tdelta/2))
    print("Attenuation coefficient: %g cm-1" % (mu))
    print("Number of lenses: %d" % (N))
    print("Length of CRL: %g mm" % (lcrl))
    print("Geometric aperture 2R0 of the individual lenses: %g μm " % (gapp))
    print("Effective aperture: %g μm" % (Deff))
    # print("Principal plane (The distance of the CRL's secondary principal plane to its center. The secondary principal plane lies to the left of the center if the focus lies to the right)): %g mm Focal length: %g mm" % (pp,feff))
    print(
        "Principal plane (The distance of the CRL's secondary principal plane to its center. The secondary principal plane lies to the left of the center if the focus lies to the right)): %g mm" % (pp))
    print(
        "Focal length: %g mm" % (feff))

    print("Image distance: %g mm " % (idist))
    print("Effective image distance: %g mm" % (eidist))
    if D2 == 0:
        print("lateral extent (FWHM) of the secondary source (geometric plus diffraction contribution) in the focused direction: %g nm" % (fspotw) )
    else:
        print("Focal spot size: %g nm × %g nm" % (fspotw,fspoth))

    print(
        "Lateral resolution: %g nm" % (res))
    print("Transmission: %g Gain: %g" % (transmgapp,gain))

if __name__ == "__main__":

    # tkt = rxoptics_compute()
    tkt = rxoptics_compute(
        MATERIAL="Be",  # *
        R=0.2,        # * radius of curvature R at the apex of the individual lenses in the CRL [mm]
        D2=1,         # * lens type: 0=1D, 1=2D
        d=30,         # * web thickness [um]
        th=1.0,       # * thickness D of the individual lenses [mm]
        F=2.0,        # * thickness of the lens frame [mm]
        presf=0,      # * 0=focal length, 1=image distance,
        desf=28200.0, # * desired focal lengt - focal length (presf=0) in mm or image distance (presf=1) in mm
        e=10.0,  # * photon energy in keV
        dist=65.0,  # * distance source-lens [m]
        w=15.3,  # * width (FWHM) of the x-ray source [um]
        h=15.3,  # * height (FWHM) of the x-ray source [um]
        use_xraylib=0,  # internal(0) or xraylib(1) optics library
        force_n_lenses=None,
        )

    rxoptics_display(tkt)