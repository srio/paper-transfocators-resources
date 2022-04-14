from srxraylib.plot.gol import plot, plot_image, plot_show
import numpy
import h5py
from tools import fwhm
import matplotlib.pyplot as plt
from barc4plots.barc4plots import Image2Plot, ESRF_colors_2D

# see https://stackoverflow.com/questions/33159134/matplotlib-y-axis-label-with-multiple-colors
def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='x',anchorpad=0,bbox_to_anchor=None,**kw ):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw))
                    for text,color in zip(list_of_strings,list_of_colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        if bbox_to_anchor is None:
            bbox_to_anchor = (0.3, -0.15)  # (0.2, -0.09),
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=bbox_to_anchor,
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',rotation=90,**kw))
                     for text,color in zip(list_of_strings[::-1],list_of_colors) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        if bbox_to_anchor is None:
            bbox_to_anchor=(-0.15, 0.375) #(-0.10, 0.2)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=bbox_to_anchor,
                                          bbox_transform=ax.transAxes, borderpad=0.)

        ax.add_artist(anchored_ybox)



if __name__ == "__main__":
    directory = "./"
    # directory = "/scisoft/users/srio/COMSYL-SLURM/doubleslit"
    filename = "%s/wofry1d_h_doubleslit.dat" % directory
    cf = numpy.loadtxt(filename)
    APERTURES = cf[:,0]
    CF = cf[:,1]
    DoC = cf[:,2]


    # plot(APERTURES*1e6,CF, xtitle="aperture(outer) [um]", ytitle="CF")
    # plot(APERTURES*1e6,DoC, xtitle="aperture(outer) [um]", ytitle="DoC fitted")



    f = h5py.File("%s/wofry1d_h_doubleslit.h5" % directory, 'r')
    SD = f["images/Intensity/image_data"][()].T
    APERTURES = f["images/Intensity/axis_x"][()]
    x = f["images/Intensity/axis_y"][()]
    f.close()






    ii = 9

    ######################################################################################################################

    # SD[ii, :] = 0
    plt.rcParams.update({'font.size': 30})
    g = plot_image((SD / SD.max()).T,
                   x * 1e3,
                   APERTURES*1e6 - 2.5,
                   xtitle="",ytitle="",title="",
                   aspect='auto',
                   yrange=[18,120], xrange=[-1.5,1.5], #xrange=[-.750,.750],
                   add_colorbar=False,figsize=(12,10), show=0,
                   cmap=ESRF_colors_2D(8))


    multicolor_ylabel(g[1], ("$s_A$ [$\mu$m]", ""),
                      ('k', 'b'),
                      axis='y', size=None, weight=None)


    multicolor_ylabel(g[1], ("$x$ [mm] ($z$=46 m)", ""),
                      ('b', 'b'),
                      axis='x', size=None, weight=None)




    filename = "doubleslit_scan.pdf"
    plt.savefig(filename)
    print("File written to disk: %s" % filename)
    plot_show()



    ######################################################################################################################
    plt.rcParams.update({'font.size': 20})
    g = plot(x*1e3, SD[ii, :] / SD.max(), xtitle="", ytitle="",
         xrange=[-1.500,1.500], title = "$s_A$ = %3.1f $\mu$m" % (APERTURES[ii]*1e6 - 2.5), show=0, figsize=(10, 8) )

    multicolor_ylabel(g[1], ("$x$ [mm]", "($z$=46 m)"),
                      ('k', 'k'),
                      axis='x', size=None, weight=None)

    multicolor_ylabel(g[1], ("$\mathcal{I}$","[arbitrary units]" ),
                      ('k', 'k'),
                      axis='y', size=None, weight=None, bbox_to_anchor=(-0.12, 0.25))


    filename = "doubleslit_profile.pdf"
    plt.savefig(filename)
    print("File written to disk: %s" % filename)
    plot_show()


    ######################################################################################################################
    plt.rcParams.update({'font.size': 18})

    DoC_wofry1 = numpy.loadtxt("DoCprofile_wofry1.dat")
    DoC_srwme = numpy.loadtxt("DoCprofile_srwme.dat")

    APERTURES2 = []
    DoC2 = []
    i = 0
    while i < APERTURES.size:
        APERTURES2.append(APERTURES[i])
        DoC2.append(DoC[i])
        i += 3
    APERTURES2 = numpy.array(APERTURES2)
    DoC2 = numpy.array(DoC2)

    plt.rcParams.update({'font.size': 28})
    g = plot(
        DoC_wofry1[:, 0], DoC_wofry1[:, 1],
        DoC_srwme[:,0], DoC_srwme[:,1],
        1e6 * APERTURES2 - 2.50, DoC2,
        title="", legend=[
            "WOFRY (CL: %3.1f $\mu$m) " %   (fwhm(DoC_srwme[:,0], DoC_srwme[:,1])),
            "SRW (CL: %3.1f $\mu$m) "   %   (fwhm(DoC_wofry1[:,0], DoC_wofry1[:,1])),
            "$\mathcal{V}(s_A)$",
            ],
        marker=[None, None, '.'], linestyle=[None,None,""],
        xtitle="", ytitle="", xrange=[0,150], yrange=[0, 1.1], show=0, figsize=(10,8))


    multicolor_ylabel(g[1], ("$x_2-x_1$", "or", "$s_A$ [$\mu$m]"),
                      ('b', 'k', 'g'),
                      axis='x', size=None, weight=None)

    multicolor_ylabel(g[1], ("|DoC|", "or", " $\mathcal{V}$ "),
                      ('g', 'k', 'b'),
                      axis='y', size=None, weight=None, bbox_to_anchor=(-0.15, 0.35))





    filename = "doubleslit_DoC.pdf"
    plt.savefig(filename)
    print("File written to disk: %s" % filename)
    plot_show()